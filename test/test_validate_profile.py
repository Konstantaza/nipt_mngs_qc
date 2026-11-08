"""
Юнит-тесты validate_profile.py — валидатора контракта profile v2.

Проверяют инварианты контракта (docs/profile_contract.md):
  1. валидная пара profile.tsv + profile.json проходит (exit 0, OK),
     включая сверку идентичности --expect-*;
  2. отсутствие обязательной колонки TSV -> INVALID;
  3. неизвестная колонка без префикса x_ -> INVALID; колонка x_* —
     игнорируется (проходит);
  4. input_reads != classified + unclassified -> INVALID;
  5. status=empty при непустом TSV -> INVALID; пустой профиль
     (только заголовок, status=empty) — валиден;
  6. несовпадение --expect-role -> INVALID;
  7. сумма reads > classified_reads -> INVALID.

Запуск: python3 -m pytest test/test_validate_profile.py -v
"""

import subprocess
import sys
from pathlib import Path

from conftest import write_profile_v2

ROOT = Path(__file__).resolve().parents[1]
VALIDATOR = ROOT / "bin" / "validate_profile.py"

TAXA = [("2097", "Mycoplasma genitalium", 500),
        ("284813", "Encephalitozoon cuniculi", 100)]


def run_validator(tsv, json_path, extra=None):
    cmd = [sys.executable, str(VALIDATOR),
           "--tsv", str(tsv), "--json", str(json_path)]
    if extra:
        cmd += extra
    return subprocess.run(cmd, capture_output=True, text=True)


# ---------------------------------------------------------------- test 1
def test_valid_profile_passes(tmp_path):
    """Валидная пара -> exit 0 / OK, включая --expect-sample/fraction/role."""
    tsv, js = write_profile_v2(tmp_path, "TEST", "dehost", "full",
                               input_reads=20_000, classified_reads=10_000,
                               unclassified_reads=10_000, taxa=TAXA)
    res = run_validator(tsv, js, ["--expect-sample", "TEST",
                                  "--expect-fraction", "dehost",
                                  "--expect-role", "full"])
    assert res.returncode == 0, f"валидный профиль отклонён: {res.stderr}"
    assert res.stdout.startswith("OK:")


# ---------------------------------------------------------------- test 2
def test_missing_tsv_column_rejected(tmp_path):
    """TSV без обязательной колонки (rank) -> INVALID."""
    tsv, js = write_profile_v2(tmp_path, "TEST", "dehost", "full",
                               input_reads=20_000, classified_reads=10_000,
                               unclassified_reads=10_000, taxa=TAXA)
    lines = tsv.read_text().splitlines()
    header = lines[0].split("\t")
    drop = header.index("rank")
    kept = ["\t".join(f for i, f in enumerate(line.split("\t"))
                      if i != drop) for line in lines]
    tsv.write_text("\n".join(kept) + "\n")
    res = run_validator(tsv, js)
    assert res.returncode == 1
    assert "INVALID" in res.stderr


# ---------------------------------------------------------------- test 3
def test_unknown_column_rejected_x_ignored(tmp_path):
    """Колонка без префикса x_ -> INVALID; x_* — игнорируется."""
    tsv, js = write_profile_v2(tmp_path, "TEST", "dehost", "full",
                               input_reads=20_000, classified_reads=10_000,
                               unclassified_reads=10_000, taxa=TAXA)
    lines = tsv.read_text().splitlines()

    # x_* колонка адаптера — легальна, потребители игнорируют
    x_lines = [lines[0] + "\tx_kraken_pct"] + \
              [line + "\t1.23" for line in lines[1:]]
    tsv.write_text("\n".join(x_lines) + "\n")
    res = run_validator(tsv, js)
    assert res.returncode == 0, f"x_* колонка отклонена: {res.stderr}"

    # неизвестная колонка без префикса — нарушение контракта
    bad_lines = [lines[0] + "\tbogus"] + [line + "\t1" for line in lines[1:]]
    tsv.write_text("\n".join(bad_lines) + "\n")
    res = run_validator(tsv, js)
    assert res.returncode == 1
    assert "INVALID" in res.stderr


# ---------------------------------------------------------------- test 4
def test_counter_invariant_rejected(tmp_path):
    """input_reads != classified + unclassified -> INVALID."""
    tsv, js = write_profile_v2(tmp_path, "TEST", "dehost", "full",
                               input_reads=1_000, classified_reads=600,
                               unclassified_reads=300,  # 600+300 != 1000
                               taxa=[("2097", "Mycoplasma genitalium", 100)])
    res = run_validator(tsv, js)
    assert res.returncode == 1
    assert "INVALID" in res.stderr


# ---------------------------------------------------------------- test 5
def test_empty_status_rules(tmp_path):
    """status=empty с строками -> INVALID; пустой профиль валиден."""
    # пустой профиль (только заголовок, status=empty) — легален
    tsv, js = write_profile_v2(tmp_path, "TEST", "dehost", "full",
                               input_reads=2, classified_reads=2,
                               unclassified_reads=0, taxa=[])
    res = run_validator(tsv, js)
    assert res.returncode == 0, f"пустой профиль отклонён: {res.stderr}"
    assert "status=empty" in res.stdout

    # status=empty, но TSV содержит строки — нарушение
    tsv, js = write_profile_v2(tmp_path, "BAD", "dehost", "full",
                               input_reads=1_000, classified_reads=1_000,
                               unclassified_reads=0, taxa=TAXA,
                               status="empty")
    res = run_validator(tsv, js)
    assert res.returncode == 1
    assert "INVALID" in res.stderr


# ---------------------------------------------------------------- test 6
def test_expect_role_mismatch_rejected(tmp_path):
    """--expect-role не совпал с sidecar -> INVALID."""
    tsv, js = write_profile_v2(tmp_path, "TEST", "dehost", "full",
                               input_reads=20_000, classified_reads=10_000,
                               unclassified_reads=10_000, taxa=TAXA)
    res = run_validator(tsv, js, ["--expect-role", "patho"])
    assert res.returncode == 1
    assert "INVALID" in res.stderr


# ---------------------------------------------------------------- test 7
def test_sum_reads_exceeds_classified_rejected(tmp_path):
    """Сумма reads по таксонам > classified_reads -> INVALID."""
    tsv, js = write_profile_v2(tmp_path, "TEST", "dehost", "full",
                               input_reads=1_000, classified_reads=100,
                               unclassified_reads=900,
                               taxa=[("2097", "Mycoplasma genitalium", 500)])
    res = run_validator(tsv, js)
    assert res.returncode == 1
    assert "INVALID" in res.stderr
