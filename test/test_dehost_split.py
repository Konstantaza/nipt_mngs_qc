"""
Юнит-тесты dehost_split.sh (v1.5.0) — партиция primary-записей BAM по
SAM-флагам на фракции host / dehost / half.

Ключевой контракт (критерий приёмки 3): каждая primary-запись попадает
РОВНО в одну фракцию — host (оба мейта картированы), dehost (оба
некартированы), half (полукартированная пара: ровно один мейт).
Secondary (0x100) и supplementary (0x800) записи исключаются (-F 0x900).
В hocort_split.nf (<= v1.4.0) полукартированные пары молча терялись
между двумя прогонами hocort.

Фикстуры — крошечные BAM, собранные из SAM-текста через samtools view.
Запуск: python3 -m pytest test/test_dehost_split.py -v
"""

import shutil
import subprocess
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SPLIT = ROOT / "bin" / "dehost_split.sh"

pytestmark = pytest.mark.skipif(shutil.which("samtools") is None,
                                reason="samtools не найден на PATH")

HEADER = "@HD\tVN:1.6\tSO:unknown\n@SQ\tSN:ref1\tLN:1000000\n"
SEQ50 = "A" * 50
QUAL50 = "I" * 50


def _sam(name, flag, mapped=False):
    rname, pos, mapq, cigar = (("ref1", 100, 60, "50M") if mapped
                               else ("*", 0, 0, "*"))
    return (f"{name}\t{flag}\t{rname}\t{pos}\t{mapq}\t{cigar}\t*\t0\t0\t"
            f"{SEQ50}\t{QUAL50}\n")


# флаги: 0x1 paired, 0x4 unmap, 0x8 munmap, 0x40 first, 0x80 second
BOTH_MAPPED = (0x1 | 0x2 | 0x20 | 0x40, 0x1 | 0x2 | 0x10 | 0x80)   # 99/147
BOTH_UNMAPPED = (0x1 | 0x4 | 0x8 | 0x40, 0x1 | 0x4 | 0x8 | 0x80)   # 77/141
R1_ONLY = (0x1 | 0x8 | 0x40, 0x1 | 0x4 | 0x80)                     # 73/133
R2_ONLY = (0x1 | 0x4 | 0x40, 0x1 | 0x8 | 0x80)                     # 69/137
SECONDARY = BOTH_MAPPED[0] | 0x100                                 # 355
SUPPLEMENTARY = BOTH_MAPPED[1] | 0x800                             # 2195


def _write_bam(path, records):
    sam = HEADER + "".join(records)
    res = subprocess.run(["samtools", "view", "-b", "-o", str(path), "-"],
                         input=sam, capture_output=True, text=True)
    assert res.returncode == 0, f"samtools view failed: {res.stderr}"


def _count(bam, extra=()):
    res = subprocess.run(["samtools", "view", "-c", *extra, str(bam)],
                         capture_output=True, text=True)
    assert res.returncode == 0
    return int(res.stdout.strip())


def _run_split(bam, out_prefix):
    res = subprocess.run(["bash", str(SPLIT), str(bam), str(out_prefix), "1"],
                         capture_output=True, text=True)
    assert res.returncode == 0, f"dehost_split упал: {res.stderr}"
    return {f: Path(f"{out_prefix}.{f}.bam") for f in ("host", "dehost", "half")}


@pytest.fixture
def mixed_bam(tmp_path):
    """4 пары: оба-mapped, оба-unmapped, R1-only, R2-only + sec/supp."""
    bam = tmp_path / "mixed.bam"
    _write_bam(bam, [
        _sam("p1", BOTH_MAPPED[0], mapped=True),
        _sam("p1", BOTH_MAPPED[1], mapped=True),
        _sam("p2", BOTH_UNMAPPED[0]), _sam("p2", BOTH_UNMAPPED[1]),
        _sam("p3", R1_ONLY[0], mapped=True), _sam("p3", R1_ONLY[1]),
        _sam("p4", R2_ONLY[0]), _sam("p4", R2_ONLY[1], mapped=True),
        _sam("p1", SECONDARY, mapped=True),
        _sam("p1", SUPPLEMENTARY, mapped=True),
    ])
    return bam


def test_categories_and_partition(mixed_bam, tmp_path):
    """Категории корректны и образуют точную партицию primary-записей."""
    outs = _run_split(mixed_bam, tmp_path / "split")
    n_host = _count(outs["host"])
    n_dehost = _count(outs["dehost"])
    n_half = _count(outs["half"])
    assert n_host == 2        # p1: оба мейта картированы
    assert n_dehost == 2      # p2: оба некартированы
    assert n_half == 4        # p3 (R1-only) + p4 (R2-only)
    # партиция: сумма фракций == числу primary-записей входа
    n_primary = _count(mixed_bam, extra=["-F", "0x900"])
    assert n_primary == 8     # 10 записей минус secondary и supplementary
    assert n_host + n_dehost + n_half == n_primary


def test_secondary_supplementary_excluded(mixed_bam, tmp_path):
    """Secondary/supplementary не попадают ни в одну фракцию."""
    outs = _run_split(mixed_bam, tmp_path / "split")
    for bam in outs.values():
        assert _count(bam, extra=["-f", "0x100"]) == 0
        assert _count(bam, extra=["-f", "0x800"]) == 0


def test_empty_categories_valid(tmp_path):
    """Только картированные пары: dehost/half — валидные пустые BAM."""
    bam = tmp_path / "mapped_only.bam"
    _write_bam(bam, [
        _sam("p1", BOTH_MAPPED[0], mapped=True),
        _sam("p1", BOTH_MAPPED[1], mapped=True),
    ])
    outs = _run_split(bam, tmp_path / "split")
    assert _count(outs["host"]) == 2
    for f in ("dehost", "half"):
        assert outs[f].exists() and outs[f].stat().st_size > 0
        assert _count(outs[f]) == 0
        # заголовок BAM читается (файл не битый)
        res = subprocess.run(["samtools", "view", "-H", str(outs[f])],
                             capture_output=True, text=True)
        assert res.returncode == 0 and "@HD" in res.stdout


def test_unpaired_records_partitioned(tmp_path):
    """SE-записи тоже покрыты партицией: mapped->host, unmapped->dehost."""
    bam = tmp_path / "se.bam"
    _write_bam(bam, [
        _sam("r1", 0x0, mapped=True),    # SE mapped -> host
        _sam("r2", 0x4),                 # SE unmapped -> dehost
    ])
    outs = _run_split(bam, tmp_path / "split")
    assert _count(outs["host"]) == 1
    assert _count(outs["dehost"]) == 1
    assert _count(outs["half"]) == 0
