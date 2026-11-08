#!/usr/bin/env python3
"""
sum_counts.py — слияние counts JSON фракций в один (ветка all, v1.5.0).

Фракция all («все прочтения без dehost») синтезируется конкатенацией
FASTA конституентов (MAKE_ALL), без нового прогона выравнивателя. Этот
скрипт строит для неё counts JSON той же схемы, что выдаёт
bam_counts.py:

  primary_reads / paired_reads / unpaired_reads — суммы по конституентам
  fragments  — пересчитывается из сумм (paired // 2 + unpaired), чтобы
               инвариант схемы выполнялся точно
  single_end — AND по конституентам (all single-end <=> все single-end)

Конституенты: fastq-режим — dehost + host + half (партиция trimmed);
bam-режим — dehost + host (партиция NIPD принимается полной).

Использование:
  sum_counts.py --sample S --fraction all --out S.all.counts.json \
      S.dehost.counts.json S.host.counts.json [S.half.counts.json]
"""

import argparse
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from bam_counts import REQUIRED_KEYS  # noqa: E402

NUMERIC_KEYS = ("primary_reads", "paired_reads", "unpaired_reads")


def load_counts(path):
    """counts JSON -> dict; ValueError при битом/неполном файле."""
    with open(path) as fh:
        data = json.load(fh)
    if not isinstance(data, dict) or \
            any(k not in data for k in REQUIRED_KEYS):
        raise ValueError(f"{path}: нет обязательных ключей {REQUIRED_KEYS}")
    return data


def main():
    ap = argparse.ArgumentParser(
        description="Сумма counts JSON фракций (синтез фракции all)")
    ap.add_argument("--sample", required=True)
    ap.add_argument("--fraction", required=True,
                    choices=["dehost", "host", "half", "all"])
    ap.add_argument("--out", required=True)
    ap.add_argument("inputs", nargs="+",
                    help="counts JSON конституентов (>= 1)")
    args = ap.parse_args()

    parts = [load_counts(p) for p in args.inputs]
    primary = sum(p["primary_reads"] for p in parts)
    paired = sum(p["paired_reads"] for p in parts)
    unpaired = sum(p["unpaired_reads"] for p in parts)
    record = {
        "sample_id": args.sample,
        "fraction": args.fraction,
        "primary_reads": primary,
        "paired_reads": paired,
        "unpaired_reads": unpaired,
        "single_end": all(p["single_end"] for p in parts),
        "fragments": paired // 2 + unpaired,
    }
    with open(args.out, "w") as fh:
        json.dump(record, fh, indent=2, ensure_ascii=False)
    print(f"[sum_counts] {args.sample}/{args.fraction}: "
          f"{len(parts)} конституентов -> primary={primary} "
          f"paired={paired} single_end={record['single_end']} "
          f"fragments={record['fragments']}")


if __name__ == "__main__":
    main()
