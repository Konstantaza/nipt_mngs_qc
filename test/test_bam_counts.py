"""
Юнит-тесты bam_counts.py (v1.3.0) — primary-only счётчики фракций BAM.

Ключевой контракт (критерий приёмки 3): счётчики НЕ зависят от secondary
(0x100) и supplementary (0x800) записей — в отличие от строки «in total»
samtools flagstat, которая их включает.

Фикстуры — крошечные BAM, собранные из SAM-текста через samtools view.
Запуск: python3 -m pytest test/test_bam_counts.py -v
"""

import json
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
BAM_COUNTS = ROOT / "bin" / "bam_counts.py"

pytestmark = pytest.mark.skipif(shutil.which("samtools") is None,
                                reason="samtools не найден на PATH")

HEADER = "@HD\tVN:1.6\tSO:unknown\n@SQ\tSN:ref1\tLN:1000000\n"
SEQ50 = "A" * 50
QUAL50 = "I" * 50


def _sam_record(name, flag, rname="*", pos=0, mapq=0, cigar="*"):
    return (f"{name}\t{flag}\t{rname}\t{pos}\t{mapq}\t{cigar}\t*\t0\t0\t"
            f"{SEQ50}\t{QUAL50}\n")


def _write_bam(path, records):
    """SAM-текст -> BAM через samtools view."""
    sam = HEADER + "".join(records)
    res = subprocess.run(["samtools", "view", "-b", "-o", str(path), "-"],
                         input=sam, capture_output=True, text=True)
    assert res.returncode == 0, f"samtools view failed: {res.stderr}"


def _run_counts(bam):
    out = bam.with_suffix(".counts.json")
    res = subprocess.run([sys.executable, str(BAM_COUNTS),
                          "--bam", str(bam), "--sample", "T",
                          "--fraction", "dehost", "--out", str(out)],
                         capture_output=True, text=True)
    assert res.returncode == 0, f"bam_counts упал: {res.stderr}"
    return json.loads(out.read_text())


# флаги
PE_R1_UNMAPPED = 0x1 | 0x4 | 0x8 | 0x40          # 77
PE_R2_UNMAPPED = 0x1 | 0x4 | 0x8 | 0x80          # 141
SE_UNMAPPED = 0x4
SECONDARY = 0x100 | 0x1 | 0x40                   # 321, картированная
SUPPLEMENTARY = 0x800 | 0x1 | 0x40               # 2113, картированная


@pytest.fixture
def pe_bam(tmp_path):
    """2 пары primary-прочтений."""
    bam = tmp_path / "pe.bam"
    _write_bam(bam, [
        _sam_record("p1", PE_R1_UNMAPPED), _sam_record("p1", PE_R2_UNMAPPED),
        _sam_record("p2", PE_R1_UNMAPPED), _sam_record("p2", PE_R2_UNMAPPED),
    ])
    return bam


def test_paired_counts(pe_bam):
    """PE BAM: primary=4, paired=4, fragments=2, single_end=False."""
    c = _run_counts(pe_bam)
    assert c["primary_reads"] == 4
    assert c["paired_reads"] == 4
    assert c["unpaired_reads"] == 0
    assert c["fragments"] == 2
    assert c["single_end"] is False


def test_secondary_supplementary_invariant(pe_bam):
    """Критерий 3: secondary/supplementary записи НЕ меняют счётчики.

    Контроль инварианта: полное число записей (как «in total» flagstat)
    при этом растёт — именно поэтому старый подсчёт был завышен.
    """
    before = _run_counts(pe_bam)
    bam2 = pe_bam.with_name("pe_sec_supp.bam")
    _write_bam(bam2, [
        _sam_record("p1", PE_R1_UNMAPPED), _sam_record("p1", PE_R2_UNMAPPED),
        _sam_record("p2", PE_R1_UNMAPPED), _sam_record("p2", PE_R2_UNMAPPED),
        _sam_record("p1", SECONDARY, "ref1", 100, 60, "50M"),
        _sam_record("p2", SUPPLEMENTARY, "ref1", 5000, 60, "50M"),
    ])
    after = _run_counts(bam2)
    assert after == before
    # sanity: записей в файле действительно стало больше
    total = subprocess.run(["samtools", "view", "-c", str(bam2)],
                           capture_output=True, text=True).stdout.strip()
    assert int(total) == 6


def test_single_end_counts(tmp_path):
    """SE BAM: paired=0 -> single_end=True, fragments = числу прочтений."""
    bam = tmp_path / "se.bam"
    _write_bam(bam, [_sam_record(f"r{i}", SE_UNMAPPED) for i in range(3)])
    c = _run_counts(bam)
    assert c["primary_reads"] == 3
    assert c["paired_reads"] == 0
    assert c["unpaired_reads"] == 3
    assert c["fragments"] == 3
    assert c["single_end"] is True


def test_empty_bam(tmp_path):
    """Пустой BAM (только заголовок): нули, single_end=True."""
    bam = tmp_path / "empty.bam"
    _write_bam(bam, [])
    c = _run_counts(bam)
    assert c["primary_reads"] == 0
    assert c["paired_reads"] == 0
    assert c["fragments"] == 0
    assert c["single_end"] is True
