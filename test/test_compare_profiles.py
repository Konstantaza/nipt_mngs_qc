"""
Юнит-тесты compare_profiles.py (comparison model v2, контракт profile v2).

Проверяют гарантии новой модели сравнения фракций:
  1. идентичные профили -> Jaccard 1, ноль таксонов only_*;
  2. равные абсолютные числа прочтений -> таксон НЕ enriched
     (регрессия бага v1.0: fold по долям внутри фракции давал ~50000x);
  3. явное обогащение only_A ловится, шум ниже порогов отсекается;
  4. несовместимые контракты -> ошибка (data integrity);
  5. пара dehost:all — синтез all = dehost + host;
  6. скрипт не содержит ссылок на колонки классификатора.

Фикстуры — контракт profile v2 (TSV + sidecar JSON), пишутся через
conftest.write_profile_v2.

Запуск: python3 -m pytest test/test_compare_profiles.py -v
"""

import json
import subprocess
import sys
from pathlib import Path

import pytest

from conftest import write_profile_v2

ROOT = Path(__file__).resolve().parents[1]
COMPARE = ROOT / "bin" / "compare_profiles.py"


def _taxa(*rows):
    """rows: (taxid, name, reads)."""
    return [(t, n, float(r)) for t, n, r in rows]


@pytest.fixture
def profiles(tmp_path):
    """Пишет пару профилей v2 и возвращает фабрику запуска compare."""
    def write(prof_d, prof_h):
        # prof_*: kwargs для write_profile_v2 (sample/fraction/role общие)
        pd_tsv, pd_json = write_profile_v2(
            tmp_path, "TEST", "dehost", "full", **prof_d)
        ph_tsv, ph_json = write_profile_v2(
            tmp_path, "TEST", "host", "full", **prof_h)
        return (pd_tsv, pd_json), (ph_tsv, ph_json)

    def run(pair_files, pair="dehost:host", expect_ok=True):
        (pd_tsv, pd_json), (ph_tsv, ph_json) = pair_files
        out_tsv = tmp_path / f"out_{pair.replace(':', '_')}.tsv"
        out_json = tmp_path / f"out_{pair.replace(':', '_')}.json"
        cmd = [
            sys.executable, str(COMPARE),
            "--dehost_tsv", str(pd_tsv), "--dehost_json", str(pd_json),
            "--host_tsv", str(ph_tsv), "--host_json", str(ph_json),
            "--sample", "TEST", "--role", "full", "--pair", pair,
            "--min_reads", "10", "--host_min_reads", "10",
            "--enrichment_fold", "10",
            "--out_tsv", str(out_tsv), "--out_json", str(out_json),
        ]
        res = subprocess.run(cmd, capture_output=True, text=True)
        if expect_ok:
            assert res.returncode == 0, f"compare упал: {res.stderr}"
            return json.loads(out_json.read_text()), out_tsv
        return res

    return write, run


# ---------------------------------------------------------------- test 1
def test_identical_profiles_jaccard_1(profiles):
    """Идентичные профили: Jaccard 1, 0 таксонов only_*, ничего не enriched."""
    write, run = profiles
    taxa = _taxa(("111", "Taxon one", 1000), ("222", "Taxon two", 500),
                 ("333", "Taxon three", 100))
    pair_files = write(dict(input_reads=100_000, taxa=taxa),
                       dict(input_reads=100_000, taxa=taxa))
    summary, _ = run(pair_files)
    assert summary["jaccard"] == 1.0
    assert summary["n_only_A"] == 0
    assert summary["n_only_B"] == 0
    assert summary["n_shared"] == 3
    assert summary["bray_curtis_dissimilarity"] == 0.0
    assert summary["spearman_rho"] == 1.0
    assert summary["n_enriched"] == 0


# ---------------------------------------------------------------- test 2
def test_equal_absolute_reads_not_enriched(profiles):
    """Регрессия бага v1.0: равные абсолютные прочтения -> НЕ enriched.

    Размеры фракций различаются в 100 раз: старая модель (fold по долям
    внутри профиля) давала бы fold=100 и статус enriched. Новая: таксон
    присутствует в host (>= host_min_reads) -> категория both -> не enriched.
    """
    write, run = profiles
    pair_files = write(
        dict(input_reads=100_000, taxa=_taxa(("111", "Taxon", 5000))),
        dict(input_reads=10_000_000, taxa=_taxa(("111", "Taxon", 5000))))
    summary, tsv = run(pair_files)
    assert summary["n_enriched"] == 0
    assert summary["n_shared"] == 1
    row = tsv.read_text().splitlines()[1].split("\t")
    assert row[6] == "both"          # category
    assert row[8] == "false"         # enriched
    assert float(row[7]) == 100.0    # RPM-fold сам по себе велик — гейтинг решает

    # Тот же случай при равных input_reads: fold ровно 1.0
    pair_files = write(
        dict(input_reads=100_000, taxa=_taxa(("111", "Taxon", 5000))),
        dict(input_reads=100_000, taxa=_taxa(("111", "Taxon", 5000))))
    summary, tsv = run(pair_files)
    row = tsv.read_text().splitlines()[1].split("\t")
    assert float(row[7]) == 1.0
    assert row[8] == "false"
    assert summary["n_enriched"] == 0


# ---------------------------------------------------------------- test 3
def test_clear_enrichment_and_noise_filtering(profiles):
    """only_A с нулём в host -> enriched (fold=null в JSON); шум отсечён."""
    write, run = profiles
    pair_files = write(
        dict(input_reads=100_000,
             taxa=_taxa(("111", "Pathogen", 1000), ("999", "Noise", 5))),
        dict(input_reads=10_000_000,
             taxa=_taxa(("777", "HostOnly", 200))))
    summary, tsv = run(pair_files)
    assert summary["n_enriched"] == 1
    assert summary["n_only_A"] == 1      # Noise (5 < min_reads) исключён
    assert summary["n_only_B"] == 1
    assert summary["n_shared"] == 0
    assert summary["jaccard"] == 0.0
    hit = summary["enriched_taxa"][0]
    assert hit["taxid"] == "111"
    assert hit["reads"] == 1000.0
    assert hit["rpm"] == 10_000.0
    assert hit["enrichment_fold"] is None      # inf -> null в JSON
    lines = tsv.read_text().splitlines()
    assert len(lines) == 3                     # заголовок + 2 таксона
    assert lines[1].split("\t")[7] == "Inf"    # inf -> "Inf" в TSV


# ---------------------------------------------------------------- test 4
def test_contract_mismatch_fails(profiles, tmp_path):
    """Несовпадение db_id / версии контракта -> осмысленная ошибка."""
    write, run = profiles
    pair_files = write(
        dict(input_reads=100_000, taxa=_taxa(("111", "T", 100)),
             db_id="db_full"),
        dict(input_reads=100_000, taxa=_taxa(("111", "T", 100)),
             db_id="db_OTHER"))
    res = run(pair_files, expect_ok=False)
    assert res.returncode != 0
    assert "db_id" in res.stderr

    # битая версия контракта в sidecar host-профиля
    pair_files = write(
        dict(input_reads=100_000, taxa=_taxa(("111", "T", 100))),
        dict(input_reads=100_000, taxa=_taxa(("111", "T", 100))))
    (_, _), (ph_tsv, ph_json) = pair_files
    bad = json.loads(ph_json.read_text())
    bad["contract_version"] = "profile.v0"
    ph_json.write_text(json.dumps(bad))
    res = run(pair_files, expect_ok=False)
    assert res.returncode != 0
    assert "contract_version" in res.stderr


# ---------------------------------------------------------------- test 5
def test_pair_dehost_all_synthesis(profiles):
    """Пара dehost:all — 'all' синтезируется как dehost + host."""
    write, run = profiles
    pair_files = write(
        dict(input_reads=100_000, taxa=_taxa(("111", "Taxon", 1000))),
        dict(input_reads=900_000,
             taxa=_taxa(("111", "Taxon", 9000), ("222", "Other", 500))))
    summary, tsv = run(pair_files, pair="dehost:all")
    assert summary["input_reads"] == {"dehost": 100_000, "all": 1_000_000}
    rows = {ln.split("\t")[0]: ln.split("\t")
            for ln in tsv.read_text().splitlines()[1:]}
    # t1: 1000 в dehost, 10000 в all -> RPM равны (10000), fold 1.0, both
    assert rows["111"][2] == "1000" and rows["111"][3] == "10000"
    assert float(rows["111"][7]) == 1.0
    assert rows["111"][6] == "both"
    # t2: есть только в host (=> в all) -> only_B
    assert rows["222"][6] == "only_B"
    assert summary["n_enriched"] == 0


# ---------------------------------------------------------------- test 6
def test_no_classifier_column_references():
    """Критерий 3: скрипт не знает колонок форматов классификатора."""
    src = COMPARE.read_text()
    for token in ("new_est_reads", "fraction_total_reads", "taxonomy_id",
                  "kraken_assigned_reads", "added_reads", "taxonomy_lvl"):
        assert token not in src, f"compare_profiles.py ссылается на {token}"
    assert "kraken" not in src.lower()
    assert "bracken" not in src.lower()
