"""
Юнит-тесты stub_to_profile.py — fixture-driven stub-адаптера интерфейса
ANNOTATE (NIPT-16).

stub синтезирует профиль из fixture-файла (taxid, name, rank, weight),
адресуемого через dbsheet (<path>/<db_id>.tsv) — тесты реестра
аннотаторов без внешней БД и без классификатора.

Проверяют:
  1. выход конвертера проходит validate_profile.py;
  2. семантика весов: reads_i = input*weight_i, abundance_i = weight_i,
     classified = round(input*sum_w), остаток -> unclassified;
  3. сортировка по reads desc и проштамповка идентичности;
  4. пустой fixture (только заголовок) -> status=empty;
  5. пустой вход (0 записей FASTA) -> status=empty даже при строках;
  6. ошибки fixture фатальны: нет файла, неверный заголовок, weight
     вне [0,1], сумма weight > 1.

Запуск: python3 -m pytest test/test_stub_to_profile.py -v
"""

import csv
import json
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
CONVERTER = ROOT / "bin" / "stub_to_profile.py"
VALIDATOR = ROOT / "bin" / "validate_profile.py"

HEADER = "taxid\tname\trank\tweight\n"
N_INPUT = 100


@pytest.fixture
def inputs(tmp_path):
    """FASTA с N_INPUT записями (r1) + пустые r2/se-плейсхолдеры."""
    r1 = tmp_path / "dehost_1.fa"
    r1.write_text("".join(f">read{i}\nACGT\n" for i in range(N_INPUT)))
    r2 = tmp_path / "dehost_2.fa"
    r2.write_text("")
    se = tmp_path / "dehost_se.fa"
    se.write_text("")
    return {"r1": r1, "r2": r2, "se": se, "tmp": tmp_path}


def write_fixture(tmp_path, body, name="fixture.tsv"):
    fx = tmp_path / name
    fx.write_text(HEADER + body)
    return fx


def run_converter(inputs, fixture, expect_ok=True):
    tsv = inputs["tmp"] / "out.profile.tsv"
    js = inputs["tmp"] / "out.profile.json"
    native = inputs["tmp"] / "out.stub.txt"
    cmd = [
        sys.executable, str(CONVERTER),
        "--fixture", str(fixture),
        "--r1", str(inputs["r1"]), "--r2", str(inputs["r2"]),
        "--se", str(inputs["se"]), "--single_end", "false",
        "--sample", "TEST", "--fraction", "dehost", "--role", "full",
        "--db_id", "stub_full", "--annotator", "stub",
        "--annotator_version", "1.0-stub",
        "--out_tsv", str(tsv), "--out_json", str(js),
        "--out_native", str(native),
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


# --- 1-3. веса, валидность, сортировка --------------------------------------
def test_weights_and_validator(inputs, tmp_path):
    fx = write_fixture(
        tmp_path,
        "2097\tMycoplasma genitalium\tspecies\t0.6\n"
        "284813\tEncephalitozoon cuniculi\tspecies\t0.3\n")
    _, tsv, js = run_converter(inputs, fx)
    run_validator(tsv, js)
    rows = read_tsv(tsv)
    # сортировка по reads desc
    assert [r["taxid"] for r in rows] == ["2097", "284813"]
    assert float(rows[0]["reads"]) == pytest.approx(60.0)
    assert float(rows[0]["abundance"]) == pytest.approx(0.6)
    assert float(rows[1]["reads"]) == pytest.approx(30.0)
    assert rows[0]["sample"] == "TEST" and rows[0]["annotator"] == "stub"
    sidecar = json.loads(js.read_text())
    assert sidecar["input_reads"] == N_INPUT
    assert sidecar["classified_reads"] == 90   # round(100 * 0.9)
    assert sidecar["unclassified_reads"] == 10  # остаток весов
    assert sidecar["status"] == "ok"


def test_full_weight_replicates_legacy_stub(inputs, tmp_path):
    """weight=1.0 по 2097 — в точности поведение stub'а v1.4.0."""
    fx = write_fixture(
        tmp_path, "2097\tMycoplasma genitalium\tspecies\t1.0\n")
    _, tsv, js = run_converter(inputs, fx)
    run_validator(tsv, js)
    row = read_tsv(tsv)[0]
    assert float(row["reads"]) == pytest.approx(N_INPUT)
    assert float(row["abundance"]) == pytest.approx(1.0)
    sidecar = json.loads(js.read_text())
    assert sidecar["classified_reads"] == N_INPUT
    assert sidecar["unclassified_reads"] == 0


# --- 4. пустой fixture -> status=empty --------------------------------------
def test_empty_fixture(inputs, tmp_path):
    fx = write_fixture(tmp_path, "")
    _, tsv, js = run_converter(inputs, fx)
    run_validator(tsv, js)
    sidecar = json.loads(js.read_text())
    assert sidecar["status"] == "empty"
    assert sidecar["classified_reads"] == 0
    assert sidecar["unclassified_reads"] == N_INPUT
    assert read_tsv(tsv) == []


# --- 5. пустой вход -> status=empty даже при строках fixture -----------------
def test_empty_input(inputs, tmp_path):
    inputs["r1"].write_text("")  # 0 записей
    fx = write_fixture(
        tmp_path, "2097\tMycoplasma genitalium\tspecies\t1.0\n")
    _, tsv, js = run_converter(inputs, fx)
    run_validator(tsv, js)
    sidecar = json.loads(js.read_text())
    assert sidecar["status"] == "empty"
    assert sidecar["input_reads"] == 0
    assert read_tsv(tsv) == []


# --- 6. ошибки fixture фатальны ----------------------------------------------
def test_missing_fixture_fails(inputs, tmp_path):
    res, _, _ = run_converter(inputs, tmp_path / "absent.tsv",
                              expect_ok=False)
    assert res.returncode != 0
    assert "fixture не найден" in res.stderr


def test_bad_header_fails(inputs, tmp_path):
    fx = tmp_path / "bad_header.tsv"
    fx.write_text("taxid\tname\trank\tfraction\n"
                  "2097\tMycoplasma genitalium\tspecies\t1.0\n")
    res, _, _ = run_converter(inputs, fx, expect_ok=False)
    assert res.returncode != 0
    assert "заголовок" in res.stderr


def test_weight_out_of_range_fails(inputs, tmp_path):
    fx = write_fixture(
        tmp_path, "2097\tMycoplasma genitalium\tspecies\t1.5\n")
    res, _, _ = run_converter(inputs, fx, expect_ok=False)
    assert res.returncode != 0
    assert "вне [0,1]" in res.stderr


def test_weight_sum_above_one_fails(inputs, tmp_path):
    fx = write_fixture(
        tmp_path,
        "2097\tMycoplasma genitalium\tspecies\t0.7\n"
        "284813\tEncephalitozoon cuniculi\tspecies\t0.5\n")
    res, _, _ = run_converter(inputs, fx, expect_ok=False)
    assert res.returncode != 0
    assert "сумма weight" in res.stderr
