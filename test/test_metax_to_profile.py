"""
Юнит-тесты metax_to_profile.py — конвертера реального выхода MetaX
(hzi-bifo/Metax) в контракт profile v2 (NIPT-16).

Фикстуры — РЕАЛЬНЫЙ выход MetaX 2.0.3 на toy-данных (не синтетика):
  test/fixtures/metax_profile_full.txt     — community-профиль
      (13 колонок) прогона `metax profile` по dehost-FASTA TOY_NORMAL
      (2942 пары) против toy-БД db_metax_full;
  test/fixtures/metax_profile_pathogen.txt — прогон с --host 9606
      против db_metax_patho (16 колонок, host-поля "unknown" без
      pathogen-таблицы).
Провенанс: test/toy/db_metax_*/ (индексы + reference.metax.fa),
команды — в docs/profile_contract.md и report_pipeline_theory.md.

Проверяют:
  1. выход конвертера проходит validate_profile.py на обеих фикстурах;
  2. маппинг колонок: reads=col4 (строкой), rank=col3, contract
     abundance = reads/input_reads, col6 -> x_abundance_pct (проценты
     MetaX НЕ идут в abundance), col9/13 -> x_cov_prob/x_chunk_prob;
  3. sidecar: classified = min(input, ceil(sum(reads))) = 2941 при
     input 2942; unclassified = input - classified;
  4. --host-вариант: 16 колонок -> x_host_names/x_host_taxids/x_diseases;
  5. "NA" в вероятностях переносится строкой (синтетическая мини-
     фикстура — в реальном toy-выходе все вероятности числовые);
  6. пустой profile -> status=empty, валидный профиль;
  7. битая строка (12 колонок) -> громкая ошибка (exit != 0).

Запуск: python3 -m pytest test/test_metax_to_profile.py -v
"""

import csv
import json
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
CONVERTER = ROOT / "bin" / "metax_to_profile.py"
VALIDATOR = ROOT / "bin" / "validate_profile.py"
FIX_FULL = ROOT / "test" / "fixtures" / "metax_profile_full.txt"
FIX_PATHO = ROOT / "test" / "fixtures" / "metax_profile_pathogen.txt"

N_INPUT = 2942  # пар в dehost-FASTA TOY_NORMAL (см. провенанс фикстур)


@pytest.fixture
def inputs(tmp_path):
    """FASTA с N_INPUT записями (r1/r2) + пустой se-плейсхолдер."""
    r1 = tmp_path / "dehost_1.fa"
    r2 = tmp_path / "dehost_2.fa"
    se = tmp_path / "dehost_se.fa"
    recs = "".join(f">read{i}\nACGT\n" for i in range(N_INPUT))
    r1.write_text(recs)
    r2.write_text(recs)
    se.write_text("")
    return {"r1": r1, "r2": r2, "se": se, "tmp": tmp_path}


def run_converter(inputs, profile, role="full", expect_ok=True):
    tsv = inputs["tmp"] / f"out.{role}.profile.tsv"
    js = inputs["tmp"] / f"out.{role}.profile.json"
    cmd = [
        sys.executable, str(CONVERTER),
        "--profile", str(profile),
        "--r1", str(inputs["r1"]), "--r2", str(inputs["r2"]),
        "--se", str(inputs["se"]), "--single_end", "false",
        "--sample", "TEST", "--fraction", "dehost", "--role", role,
        "--db_id", "db_metax_full", "--annotator", "metax",
        "--annotator_version", "2.0.3",
        "--out_tsv", str(tsv), "--out_json", str(js),
    ]
    res = subprocess.run(cmd, capture_output=True, text=True)
    if expect_ok:
        assert res.returncode == 0, f"конвертер упал: {res.stderr}"
    return res, tsv, js


def run_validator(tsv, js):
    res = subprocess.run(
        [sys.executable, str(VALIDATOR), "--tsv", str(tsv),
         "--json", str(js)], capture_output=True, text=True)
    assert res.returncode == 0, f"валидатор отклонил профиль: {res.stderr}"


def read_tsv(tsv):
    with open(tsv, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


# --- 1-3. community-фикстура: валидность, маппинг, счётчики ----------------
def test_real_fixture_passes_validator(inputs):
    _, tsv, js = run_converter(inputs, FIX_FULL)
    run_validator(tsv, js)


def test_column_mapping_and_abundance(inputs):
    _, tsv, js = run_converter(inputs, FIX_FULL)
    rows = read_tsv(tsv)
    assert len(rows) == 3
    # сортировка по reads desc: 2097 (1911) -> 284813 (520) -> 10798 (510)
    assert [r["taxid"] for r in rows] == ["2097", "284813", "10798"]
    top = rows[0]
    assert top["name"] == "Mycoplasma genitalium"
    assert top["rank"] == "species"
    assert top["reads"] == "1911.00000"  # col4 строкой
    # contract abundance = reads/input, а НЕ col6 (проценты MetaX)
    assert float(top["abundance"]) == pytest.approx(1911 / N_INPUT,
                                                    rel=1e-5)
    assert top["x_abundance_pct"] == "3.39595"  # col6 сохранена как x_*
    assert top["x_cov_prob"] == "1.00000"       # col9
    assert top["x_chunk_prob"] == "1.00000"     # col13
    # идентичность проштампована в каждой строке
    assert top["sample"] == "TEST" and top["annotator"] == "metax"
    assert top["annotator_version"] == "2.0.3"
    assert top["db_id"] == "db_metax_full"


def test_sidecar_counts(inputs):
    _, tsv, js = run_converter(inputs, FIX_FULL)
    sidecar = json.loads(js.read_text())
    assert sidecar["input_reads"] == N_INPUT
    # sum(reads) = 2941 -> classified = min(2942, 2941)
    assert sidecar["classified_reads"] == 2941
    assert sidecar["unclassified_reads"] == 1
    assert sidecar["status"] == "ok"


# --- 4. --host-фикстура: 16 колонок ----------------------------------------
def test_host_fixture_16_columns(inputs):
    _, tsv, js = run_converter(inputs, FIX_PATHO, role="patho")
    run_validator(tsv, js)
    rows = read_tsv(tsv)
    assert len(rows) == 1
    row = rows[0]
    assert row["taxid"] == "284813"
    assert row["x_host_names"] == "unknown"
    assert row["x_host_taxids"] == "unknown"
    assert row["x_diseases"] == "unknown"
    # col6 = 100.00000 (ренормировка на один таксон) — в abundance НЕ идёт
    assert row["x_abundance_pct"] == "100.00000"
    # abundance пишется компактно (%.6g) — допуск на округление записи
    assert float(row["abundance"]) == pytest.approx(520 / N_INPUT,
                                                    rel=1e-5)
    sidecar = json.loads(js.read_text())
    assert sidecar["classified_reads"] == 520
    assert sidecar["unclassified_reads"] == N_INPUT - 520


# --- 5. NA в вероятностях (синтетическая строка формата MetaX) --------------
def test_na_probabilities_passthrough(inputs, tmp_path):
    prof = tmp_path / "na.profile.txt"
    prof.write_text(
        "Some taxon\t12345\tspecies\t10.00000\t0.1\t100.00000\t0.5\t0.5\t"
        "NA\t1.00000\t1.00000\t1.00000\tNA\n")
    _, tsv, js = run_converter(inputs, prof)
    run_validator(tsv, js)
    row = read_tsv(tsv)[0]
    assert row["x_cov_prob"] == "NA"
    assert row["x_chunk_prob"] == "NA"


# --- 6. пустой profile -> status=empty --------------------------------------
def test_empty_profile(inputs, tmp_path):
    prof = tmp_path / "empty.profile.txt"
    prof.write_text("")
    _, tsv, js = run_converter(inputs, prof)
    run_validator(tsv, js)
    sidecar = json.loads(js.read_text())
    assert sidecar["status"] == "empty"
    assert sidecar["classified_reads"] == 0
    assert sidecar["unclassified_reads"] == N_INPUT
    assert read_tsv(tsv) == []  # только заголовок


# --- 7. битая строка -> громкая ошибка --------------------------------------
def test_malformed_row_fails(inputs, tmp_path):
    prof = tmp_path / "bad.profile.txt"
    prof.write_text("Some taxon\t12345\tspecies\t10.0\t0.1\t100.0\t0.5\t"
                    "0.5\t1.0\t1.0\t1.0\t1.0\n")  # 12 колонок
    res, _, _ = run_converter(inputs, prof, expect_ok=False)
    assert res.returncode != 0
    assert "колонок" in res.stderr


@pytest.mark.parametrize("reads", ["nan", "inf", "-inf"])
def test_nonfinite_reads_rejected(inputs, tmp_path, reads):
    prof = tmp_path / "nonfinite.txt"
    cols = FIX_FULL.read_text().splitlines()[0].split("\t")
    cols[3] = reads
    prof.write_text("\t".join(cols) + "\n")
    res, _, _ = run_converter(inputs, prof, expect_ok=False)
    assert res.returncode != 0
    assert "конечным" in res.stderr


def test_counts_above_input_rejected(inputs, tmp_path):
    prof = tmp_path / "scaled.txt"
    cols = FIX_FULL.read_text().splitlines()[0].split("\t")
    cols[3] = str(N_INPUT + 1)
    prof.write_text("\t".join(cols) + "\n")
    res, _, _ = run_converter(inputs, prof, expect_ok=False)
    assert res.returncode != 0
    assert "sum(reads)" in res.stderr
