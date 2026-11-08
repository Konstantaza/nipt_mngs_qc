"""
Юнит-тесты qc_gate.py (контракт profile v2, версия пайплайна 1.10.0).

Проверяют контрактные гарантии QC gate:
  1. отсутствие обязательного артефакта (пара profile.tsv/.json / diff JSON /
     fastp) -> FAIL/PIPELINE_ERROR, а не молчаливый PASS_NEGATIVE;
  2. утечка human в dehost не оценивается при нехватке не-хост прочтений;
  3. позитивность — только по панели мишеней; без панели — not_evaluated;
  4. QC JSON пишется всегда (образец не «пропадает») и несёт annotator/db_id.

Фикстуры профилей — контракт profile v2 (TSV + sidecar JSON), пишутся через
conftest.write_profile_v2. Метрики dehost-ветки qc_gate вычисляет
классификаторо-независимо: human — строка taxid 9606 в TSV,
classified_frac — счётчики sidecar.

Запуск: python3 -m pytest test/test_qc_gate.py -v
"""

import json
import subprocess
import sys
from pathlib import Path

import pytest

from conftest import (ANNOTATOR, ANNOTATOR_COMPOSED, ANNOTATOR_VERSION,
                      write_profile_v2)

ROOT = Path(__file__).resolve().parents[1]
QC_GATE = ROOT / "bin" / "qc_gate.py"

# Пороги под синтетические числа (как в test-профиле пайплайна)
THRESHOLDS = {
    "min_input_pairs": 10_000,
    "min_q30": 0.50,
    "min_trim_survival": 0.50,
    "host_frac_min": 0.80,
    "min_nonhost_reads": 500,
    "max_dehost_human_frac": 0.05,
    "min_classified_frac": 0.0001,
    "max_consistency_dev": 0.01,
}


def _write_counts(path, primary_reads, paired_reads=None):
    """counts JSON стадии BAM2FASTA (bam_counts.py, primary-only)."""
    paired_reads = primary_reads if paired_reads is None else paired_reads
    path.write_text(json.dumps({
        "sample_id": "TEST", "fraction": "dehost",
        "primary_reads": primary_reads, "paired_reads": paired_reads,
        "unpaired_reads": primary_reads - paired_reads,
        "single_end": paired_reads == 0,
        "fragments": paired_reads // 2 + (primary_reads - paired_reads),
    }))


def _write_fastp(path, before=2_100_000, after=2_000_000, q30=0.95):
    path.write_text(json.dumps({
        "summary": {
            "before_filtering": {"total_reads": before, "q30_rate": q30},
            "after_filtering": {"total_reads": after},
        }
    }))


def _diff(sample, role, enriched):
    """Сводка comparison model v2 (dehost:host) в схеме compare JSON v1.4.0."""
    return {
        "sample_id": sample, "role": role, "db_id": f"db_{role}",
        "annotator": ANNOTATOR, "annotator_version": ANNOTATOR_VERSION,
        "pair": "dehost:host",
        "input_reads": {"dehost": 20_000, "host": 1_980_000},
        "n_shared": 1, "n_only_A": 3, "n_only_B": 1,
        "jaccard": 0.2, "bray_curtis_dissimilarity": 0.9,
        "spearman_rho": 1.0,
        "n_taxa_considered": 10, "n_enriched": len(enriched),
        "n_background": 10 - len(enriched),
        "enriched_taxa": enriched,
    }


ENRICHED_FULL = [{"taxid": 2097, "name": "Mycoplasma genitalium",
                  "dehost_reads": 500, "enrichment_fold": 20.0}]
ENRICHED_PATHO = [{"taxid": 284813, "name": "Encephalitozoon cuniculi",
                   "dehost_reads": 100, "enrichment_fold": 15.0}]


@pytest.fixture
def inputs(tmp_path):
    """Полный валидный набор входов; тесты ломают ровно одну часть."""
    sample = "TEST"
    dirs = {}
    for d in ("counts", "profiles", "diffs"):
        dirs[d] = tmp_path / d
        dirs[d].mkdir()

    cnt_d = dirs["counts"] / "TEST.dehost.counts.json"
    cnt_h = dirs["counts"] / "TEST.host.counts.json"
    _write_counts(cnt_d, 20_000)         # 10k не-хост пар
    _write_counts(cnt_h, 1_980_000)      # итого 1M пар, host_frac 0.99

    fastp = tmp_path / "fastp.json"
    _write_fastp(fastp)

    # 4 пары контракта profile v2. dehost.full несёт строку human
    # (taxid 9606, 200 прочтений из 20000 -> human_frac 0.01) и
    # classified_frac 0.5 (10000/20000); остальные профили валидны
    # и человека не содержат.
    write_profile_v2(dirs["profiles"], sample, "dehost", "full",
                     input_reads=20_000, classified_reads=10_000,
                     unclassified_reads=10_000,
                     taxa=[("9606", "Homo sapiens", 200),
                           ("2097", "Mycoplasma genitalium", 500)])
    for frac, role in (("dehost", "patho"), ("host", "full"),
                       ("host", "patho")):
        write_profile_v2(dirs["profiles"], sample, frac, role,
                         input_reads=20_000, classified_reads=10_000,
                         unclassified_reads=10_000,
                         taxa=[("2097", "Mycoplasma genitalium", 100)])
    (dirs["diffs"] / f"{sample}.full.dehost_vs_host.compare.json").write_text(
        json.dumps(_diff(sample, "full", ENRICHED_FULL)))
    (dirs["diffs"] / f"{sample}.patho.dehost_vs_host.compare.json").write_text(
        json.dumps(_diff(sample, "patho", ENRICHED_PATHO)))

    return {"sample": sample, "cnt_d": cnt_d, "cnt_h": cnt_h, "fastp": fastp,
            "dirs": dirs, "tmp": tmp_path}


def run_gate(inputs, extra_args=None, mode="fastq"):
    out_json = inputs["tmp"] / "out_qc.json"
    out_mqc = inputs["tmp"] / "out_mqc.tsv"
    cmd = [
        sys.executable, str(QC_GATE),
        "--sample", inputs["sample"], "--mode", mode,
        "--fastp_json", str(inputs["fastp"]),
        "--counts_dehost", str(inputs["cnt_d"]),
        "--counts_host", str(inputs["cnt_h"]),
        "--profile_dir", str(inputs["dirs"]["profiles"]),
        "--diff_dir", str(inputs["dirs"]["diffs"]),
        "--out_json", str(out_json), "--out_mqc", str(out_mqc),
    ]
    for k, v in THRESHOLDS.items():
        cmd += [f"--{k}", str(v)]
    if extra_args:
        cmd += extra_args
    res = subprocess.run(cmd, capture_output=True, text=True)
    assert res.returncode == 0, f"qc_gate упал: {res.stderr}"
    assert out_json.exists(), "QC JSON не записан: образец «пропал» без статуса"
    return json.loads(out_json.read_text()), out_mqc


def crit(report, cid):
    return next(c for c in report["criteria"] if c["id"] == cid)


# ---------------------------------------------------------------- test 1
def test_missing_profile_fails(inputs):
    """Отсутствие sidecar JSON одного профиля -> FAIL/PIPELINE_ERROR."""
    (inputs["dirs"]["profiles"] / "TEST.dehost.patho.profile.json").unlink()
    report, _ = run_gate(inputs)
    assert report["status"] == "FAIL"
    assert "PIPELINE_ERROR" in report["fallback_reasons"]
    assert crit(report, "artifacts_profiles")["pass"] is False


# ---------------------------------------------------------------- test 2
def test_missing_patho_diff_fails(inputs):
    """Отсутствие patho diff -> FAIL/PIPELINE_ERROR, а НЕ PASS_NEGATIVE."""
    (inputs["dirs"]["diffs"]
        / "TEST.patho.dehost_vs_host.compare.json").unlink()
    panel = inputs["tmp"] / "panel.txt"
    panel.write_text("284813\n")
    report, _ = run_gate(inputs, ["--target_panel", str(panel)])
    assert report["status"] == "FAIL"
    assert report["status"] != "PASS_NEGATIVE"
    assert "PIPELINE_ERROR" in report["fallback_reasons"]
    assert crit(report, "artifacts_diffs")["pass"] is False


# ---------------------------------------------------------------- test 3
def test_leak_not_evaluated_on_tiny_nonhost(inputs):
    """human_frac=1.0 при 2 не-хост парах: утечка НЕ оценивается.

    Ожидание: FALLBACK только по ZERO_MICROBIOME, критерий
    dehost_human_frac пропущен (pass=None), DEHOST_FAILURE отсутствует.
    """
    _write_counts(inputs["cnt_d"], 4)             # 2 не-хост пары
    _write_counts(inputs["cnt_h"], 1_999_996)     # итого 1M пар
    # dehost.full: все 4 прочтения — human (taxid 9606) -> human_frac 1.0
    write_profile_v2(inputs["dirs"]["profiles"], "TEST", "dehost", "full",
                     input_reads=4, classified_reads=4, unclassified_reads=0,
                     taxa=[("9606", "Homo sapiens", 4)])

    report, _ = run_gate(inputs)
    assert report["status"] == "FALLBACK"
    assert report["fallback_reasons"] == ["ZERO_MICROBIOME"]
    assert "DEHOST_FAILURE" not in report["fallback_reasons"]
    assert crit(report, "dehost_human_frac")["pass"] is None


# ---------------------------------------------------------------- test 4
def test_malformed_fastp_fails(inputs):
    """Битый fastp JSON в режиме fastq -> FAIL/PIPELINE_ERROR (оба варианта)."""
    # 4a: невалидный JSON
    inputs["fastp"].write_text("{not a json")
    report, _ = run_gate(inputs)
    assert report["status"] == "FAIL"
    assert "PIPELINE_ERROR" in report["fallback_reasons"]
    assert crit(report, "artifacts_fastp")["pass"] is False

    # 4b: валидный JSON без обязательных ключей summary.*
    inputs["fastp"].write_text(json.dumps({"unexpected": {}}))
    report, _ = run_gate(inputs)
    assert report["status"] == "FAIL"
    assert "PIPELINE_ERROR" in report["fallback_reasons"]
    assert crit(report, "artifacts_fastp")["pass"] is False


# ---------------------------------------------------------------- test 5
def test_panel_positivity(inputs):
    """Позитивность — только по панели; без панели — not_evaluated."""
    panel_hit = inputs["tmp"] / "panel_hit.txt"
    panel_hit.write_text("284813\n")
    report, _ = run_gate(inputs, ["--target_panel", str(panel_hit)])
    assert report["status"] == "PASS_POSITIVE"
    assert report["positivity_source"] == "target_panel"
    assert report["target_panel"]["hits"] == ["Encephalitozoon cuniculi"]

    panel_miss = inputs["tmp"] / "panel_miss.txt"
    panel_miss.write_text("999999\n")
    report, _ = run_gate(inputs, ["--target_panel", str(panel_miss)])
    assert report["status"] == "PASS_NEGATIVE"
    assert report["positivity_source"] == "target_panel"
    assert report["target_panel"]["hits"] == []

    # панель не передана вовсе
    report, _ = run_gate(inputs)
    assert report["status"] == "PASS_NEGATIVE"
    assert report["positivity_source"] == "not_evaluated"
    assert report["target_panel"]["file"] is None


# ---------------------------------------------------------------- test 6
def test_happy_path_contract(inputs):
    """Полный валидный вход -> PASS_POSITIVE + annotator/db_id/comparison_info."""
    panel = inputs["tmp"] / "panel.txt"
    panel.write_text("284813\n")
    report, out_mqc = run_gate(inputs, ["--target_panel", str(panel)])

    assert report["status"] == "PASS_POSITIVE"
    assert report["pipeline_version"] == "1.10.0"
    # profile v2: annotator в QC JSON — "<annotator> <annotator_version>"
    assert report["annotator"] == ANNOTATOR_COMPOSED
    assert report["db_id"] == {"full": "db_full", "patho": "db_patho"}
    assert report["target_panel"]["n_targets"] == 1
    # метрики dehost-ветки вычислены из контракта (human = taxid 9606)
    assert report["metrics"]["dehost_human_reads"] == 200
    assert report["metrics"]["dehost_human_frac"] == 0.01
    assert report["metrics"]["dehost_classified_frac"] == 0.5
    # comparison_info — информационная сводка по обеим ролям БД
    assert set(report["comparison_info"]) == {"full", "patho"}
    assert report["comparison_info"]["patho"]["n_enriched"] == 1
    # все критерии пройдены
    assert all(c["pass"] in (True, None) for c in report["criteria"])
    # MultiQC custom content записан из того же контракта
    mqc = out_mqc.read_text()
    assert "FINAL_STATUS\tPASS_POSITIVE" in mqc


# ---------------------------------------------------------------- test 7
def test_consistency_with_half_counts(inputs):
    """v1.5.0: consistency = dehost+host+half против trimmed (критерий 3).

    half=1000 прочтений, fastp after = 2_001_000 -> dev=0 -> pass;
    информационные метрики half_mapped_reads/half_mapped_frac заполнены.
    """
    cnt_half = inputs["dirs"]["counts"] / "TEST.half.counts.json"
    _write_counts(cnt_half, 1_000)
    _write_fastp(inputs["fastp"], after=2_001_000)
    report, _ = run_gate(inputs, ["--counts_half", str(cnt_half)])
    assert crit(report, "consistency")["pass"] is True
    assert report["metrics"]["consistency_dev"] == 0.0
    assert report["metrics"]["half_mapped_reads"] == 1_000
    assert report["metrics"]["half_mapped_frac"] == round(1_000 / 2_001_000, 6)


# ---------------------------------------------------------------- test 8
def test_consistency_half_mismatch_fails(inputs):
    """half входит в сумму партиции: рассинхрон -> SPLIT_INCONSISTENT."""
    cnt_half = inputs["dirs"]["counts"] / "TEST.half.counts.json"
    _write_counts(cnt_half, 100_000)   # лишние 100k против fastp after
    report, _ = run_gate(inputs, ["--counts_half", str(cnt_half)])
    assert crit(report, "consistency")["pass"] is False
    assert report["status"] == "FALLBACK"
    assert "SPLIT_INCONSISTENT" in report["fallback_reasons"]
