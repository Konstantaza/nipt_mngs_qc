"""FASTQ SE schema failures and real mixed PE/SE integration routing."""
import csv
import gzip
import json
from pathlib import Path
import shutil
import subprocess
import pytest

ROOT = Path(__file__).resolve().parents[1]
TOOLS = ("nextflow", "fastp", "hisat2", "samtools", "kraken2", "bracken", "multiqc")
pytestmark = pytest.mark.skipif(not all(shutil.which(t) for t in TOOLS),
                                reason="Native Nextflow integration toolchain unavailable")


def nf(tmp_path, rows):
    sheet = tmp_path/"samples.csv"
    with open(sheet, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["sample", "fastq_1", "fastq_2", "single_end"])
        w.writerows(rows)
    out = tmp_path/"results"
    proc = subprocess.run(
        ["nextflow", "-log", str(tmp_path/"nextflow.log"), "run", str(ROOT/"main.nf"),
         "-profile", "test", "-ansi-log", "false", "-work-dir", str(tmp_path/"work"),
         "--mode", "fastq", "--input", str(sheet), "--outdir", str(out),
         "--publish_validation_artifacts", "true"],
        cwd=tmp_path, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
        timeout=180)
    return proc, out


@pytest.mark.parametrize("row,message", [
    (["BAD", "r1", "", "false"], "fastq_2"),
    (["BAD", "r1", "r2", "true"], "single_end=true"),
    (["BAD", "r1", "", "maybe"], "single_end"),
])
def test_fastq_schema_rejects_contradictions(tmp_path, row, message):
    p, out = nf(tmp_path, [row])
    assert p.returncode != 0
    assert message in p.stdout
    assert not list(out.glob("*/qc/*_qc_report.json"))


def test_real_mixed_fastq_routing(tmp_path):
    r1 = ROOT/"test/toy/TOY_NORMAL_R1.fastq.gz"
    r2 = ROOT/"test/toy/TOY_NORMAL_R2.fastq.gz"
    se = tmp_path/"se50.fastq"
    # 400 records are enough to validate routing, not biological acceptance.
    with gzip.open(r1, "rt") as inp, open(se, "w") as out:
        for i in range(400):
            h, seq, plus, q = [inp.readline().rstrip() for _ in range(4)]
            out.write(f"{h}\n{seq[:50]}\n+\n{q[:50]}\n")
    p, out = nf(tmp_path, [["ROUTE_PE", r1, r2, ""], ["ROUTE_SE", se, "", "true"]])
    assert p.returncode == 0, p.stdout
    for s, is_se in (("ROUTE_PE", False), ("ROUTE_SE", True)):
        qc = json.loads((out/s/f"qc/{s}_qc_report.json").read_text())
        assert qc["metrics"]["single_end"] is is_se
        fastp = json.loads((out/s/f"qc/fastp/{s}.fastp.json").read_text())
        total = 0
        for fr in ("host", "dehost", "half"):
            cj = json.loads((out/s/f"qc/counts/{s}.{fr}.counts.json").read_text())
            total += cj["primary_reads"]
            if is_se and fr == "half":
                assert cj["primary_reads"] == 0
        assert total == fastp["summary"]["after_filtering"]["total_reads"]
        assert (out/s/"validation/trimmed/trim_R2.fastq.gz").exists() is (not is_se)
