"""
Юнит-тесты kraken2_bracken_to_profile.py — конвертера адаптера
Kraken2+Bracken (kraken2 report + bracken TSV -> контракт profile v2).

Проверяют:
  1. выход конвертера проходит validate_profile.py (включая --expect-*),
     счётчики sidecar — из report (input = U + R);
  2. содержимое строк TSV: reads = new_est_reads, abundance = reads/input,
     rank S -> species, сортировка по reads desc, компактная запись чисел;
  3. пустой bracken (только заголовок) -> status=empty, валидный профиль;
  4. report без строки U -> unclassified_reads = 0;
  5. строка human (taxid 9606) проходит в TSV (её читает qc_gate);
  6. таксоны с new_est_reads = 0 отбрасываются.

Запуск: python3 -m pytest test/test_kraken2_bracken_to_profile.py -v
"""

import csv
import json
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
CONVERTER = ROOT / "bin" / "kraken2_bracken_to_profile.py"
VALIDATOR = ROOT / "bin" / "validate_profile.py"

# kraken2 report: 6 колонок (pct, reads, direct, rank code, taxid, имя).
# Счётчики берутся только из строк U и R; кладные строки (R1/D/S)
# игнорируются.
REPORT = (
    "  0.34\t10\t10\tU\t0\tunclassified\n"
    " 99.66\t2932\t0\tR\t1\troot\n"
    " 99.66\t2932\t0\tR1\t131567\t  cellular organisms\n"
    " 61.18\t1800\t1800\tS\t2097\t    Mycoplasma genitalium\n"
)

BRACKEN_HEADER = ("name\ttaxonomy_id\ttaxonomy_lvl\tkraken_assigned_reads\t"
                  "added_reads\tnew_est_reads\tfraction_total_reads")

BRACKEN_ROWS = (
    "Mycoplasma genitalium\t2097\tS\t1800\t0\t1800\t0.6118\n"
    "Encephalitozoon cuniculi\t284813\tS\t517\t0\t517\t0.1757\n"
    "Human parvovirus B19\t10798\tS\t510\t0\t510\t0.1734\n"
    "Homo sapiens\t9606\tS\t50\t0\t50\t0.0170\n"
    "Bacillus subtilis\t1423\tS\t10\t0\t0\t0.0000\n"  # new_est=0 -> отброс
)


@pytest.fixture
def inputs(tmp_path):
    report = tmp_path / "k2.report"
    report.write_text(REPORT)
    bracken = tmp_path / "bracken.tsv"
    bracken.write_text(BRACKEN_HEADER + "\n" + BRACKEN_ROWS)
    return {"report": report, "bracken": bracken, "tmp": tmp_path}


def run_converter(inputs, report=None, bracken=None):
    tsv = inputs["tmp"] / "out.profile.tsv"
    js = inputs["tmp"] / "out.profile.json"
    cmd = [
        sys.executable, str(CONVERTER),
        "--report", str(report or inputs["report"]),
        "--bracken", str(bracken or inputs["bracken"]),
        "--sample", "TEST", "--fraction", "dehost", "--role", "full",
        "--db_id", "db_full", "--annotator", "kraken2_bracken",
        "--annotator_version", "2.17.2+3.0.1",
        "--out_tsv", str(tsv), "--out_json", str(js),
    ]
    res = subprocess.run(cmd, capture_output=True, text=True)
    assert res.returncode == 0, f"конвертер упал: {res.stderr}"
    return tsv, js


def run_validator(tsv, js, extra=None):
    cmd = [sys.executable, str(VALIDATOR), "--tsv", str(tsv),
           "--json", str(js)] + (extra or [])
    return subprocess.run(cmd, capture_output=True, text=True)


def read_taxa(tsv):
    with open(tsv, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


# ---------------------------------------------------------------- test 1
def test_output_passes_validator(inputs):
    """Выход конвертера валиден; счётчики sidecar = U+R из report."""
    tsv, js = run_converter(inputs)
    res = run_validator(tsv, js, ["--expect-sample", "TEST",
                                  "--expect-fraction", "dehost",
                                  "--expect-role", "full"])
    assert res.returncode == 0, f"выход конвертера невалиден: {res.stderr}"
    sidecar = json.loads(js.read_text())
    assert sidecar["input_reads"] == 2942        # 10 U + 2932 R
    assert sidecar["classified_reads"] == 2932
    assert sidecar["unclassified_reads"] == 10
    assert sidecar["status"] == "ok"


# ---------------------------------------------------------------- test 2
def test_row_content_and_abundance(inputs):
    """Строки TSV: reads/abundance/rank/сортировка/формат чисел."""
    tsv, _ = run_converter(inputs)
    taxa = read_taxa(tsv)
    # new_est_reads=0 отброшен, сортировка по reads desc
    assert [t["taxid"] for t in taxa] == ["2097", "284813", "10798", "9606"]
    myco = taxa[0]
    assert myco["name"] == "Mycoplasma genitalium"
    assert myco["rank"] == "species"            # код S -> species
    assert myco["reads"] == "1800"              # компактная запись float
    assert float(myco["abundance"]) == pytest.approx(1800 / 2942)
    for t in taxa:                              # идентичность в каждой строке
        assert t["sample"] == "TEST"
        assert t["fraction"] == "dehost"
        assert t["annotator"] == "kraken2_bracken"
        assert t["annotator_version"] == "2.17.2+3.0.1"
        assert t["db_id"] == "db_full"


# ---------------------------------------------------------------- test 3
def test_empty_bracken_gives_empty_status(inputs):
    """Bracken только с заголовком -> status=empty, валидный профиль."""
    empty = inputs["tmp"] / "empty.tsv"
    empty.write_text(BRACKEN_HEADER + "\n")
    tsv, js = run_converter(inputs, bracken=empty)
    sidecar = json.loads(js.read_text())
    assert sidecar["status"] == "empty"
    assert read_taxa(tsv) == []                 # только заголовок
    res = run_validator(tsv, js)
    assert res.returncode == 0, f"пустой профиль невалиден: {res.stderr}"


# ---------------------------------------------------------------- test 4
def test_report_without_u_line(inputs):
    """Report без строки U -> unclassified_reads = 0."""
    no_u = inputs["tmp"] / "no_u.report"
    no_u.write_text("\n".join(REPORT.splitlines()[1:]) + "\n")
    tsv, js = run_converter(inputs, report=no_u)
    sidecar = json.loads(js.read_text())
    assert sidecar["unclassified_reads"] == 0
    assert sidecar["input_reads"] == sidecar["classified_reads"] == 2932
    res = run_validator(tsv, js)
    assert res.returncode == 0, f"профиль без U-строки невалиден: {res.stderr}"


# ---------------------------------------------------------------- test 5
def test_human_row_propagates(inputs):
    """Строка taxid 9606 (human) попадает в TSV — её читает qc_gate."""
    tsv, _ = run_converter(inputs)
    human = [t for t in read_taxa(tsv) if t["taxid"] == "9606"]
    assert len(human) == 1
    assert human[0]["name"] == "Homo sapiens"
    assert float(human[0]["reads"]) == 50.0
