#!/usr/bin/env python3
"""
bam_counts.py — первичные счётчики прочтений фракции BAM (v1.3.0).

Заменяет чтение «in total» из samtools flagstat: та строка включает
secondary (0x100) и supplementary (0x800) записи и завышает счётчики QC.
Здесь считаются ТОЛЬКО primary-записи:

  primary_reads   = samtools view -c -F 0x900  bam
  paired_reads    = samtools view -c -F 0x900 -f 0x1  bam
  unpaired_reads  = primary_reads - paired_reads
  single_end      = (paired_reads == 0)   # автоопределение per fraction
  fragments       = paired_reads // 2 + unpaired_reads

«fragments» — число секвенированных фрагментов: для парных данных это пары,
для single-end — сами прочтения. QC gate использует fragments там, где
раньше брал total/2.

Выход: <sample>.<fraction>.counts.json — контракт стадии BAM2FASTA для
QC gate (наряду с profile.json и compare JSON).
"""

import argparse
import json
import shutil
import subprocess
import sys

REQUIRED_KEYS = ("sample_id", "fraction", "primary_reads", "paired_reads",
                 "unpaired_reads", "single_end", "fragments")


def view_count(bam, extra_flags):
    """samtools view -c <flags> bam -> int."""
    cmd = ["samtools", "view", "-c", *extra_flags, bam]
    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.returncode != 0:
        raise RuntimeError(f"samtools view failed: {res.stderr.strip()}")
    return int(res.stdout.strip())


def count_bam(bam):
    """-> dict со счётчиками primary-only."""
    primary = view_count(bam, ["-F", "0x900"])
    paired = view_count(bam, ["-F", "0x900", "-f", "0x1"])
    unpaired = primary - paired
    return {
        "primary_reads": primary,
        "paired_reads": paired,
        "unpaired_reads": unpaired,
        "single_end": paired == 0,
        "fragments": paired // 2 + unpaired,
    }


def main():
    ap = argparse.ArgumentParser(
        description="Primary-only счётчики прочтений фракции BAM")
    ap.add_argument("--bam", required=True)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--fraction", required=True,
                    choices=["dehost", "host", "half", "all"])
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    if shutil.which("samtools") is None:
        sys.exit("bam_counts.py: samtools не найден на PATH")

    counts = count_bam(args.bam)
    record = {"sample_id": args.sample, "fraction": args.fraction, **counts}
    with open(args.out, "w") as fh:
        json.dump(record, fh, indent=2, ensure_ascii=False)
    print(f"[bam_counts] {args.sample}/{args.fraction}: "
          f"primary={counts['primary_reads']} paired={counts['paired_reads']} "
          f"single_end={counts['single_end']} fragments={counts['fragments']}")


if __name__ == "__main__":
    main()
