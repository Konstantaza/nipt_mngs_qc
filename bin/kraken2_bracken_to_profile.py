#!/usr/bin/env python3
"""
kraken2_bracken_to_profile.py — конвертер адаптера Kraken2+Bracken.

Kraken2 report + Bracken TSV -> контракт profile v2
(profile.tsv + sidecar profile.json, см. docs/profile_contract.md).
ЕДИНСТВЕННОЕ место адаптера, знающее форматы Kraken2/Bracken
(парсинг перенесён из profile_contract.py, v1.3).

Формат kraken2 report (TSV, 6 колонок):
  1) pct — доля прочтений клады, %;  2) reads — кладные прочтения;
  3) direct;  4) rank code (U/R/D/.../S);  5) taxid;  6) имя.
Счётчики sidecar: input = U + R, classified = R, unclassified = U
(строки U может не быть, если unclassified = 0).

Формат bracken TSV (уровень S, заголовок обязателен):
  name, taxonomy_id, taxonomy_lvl, kraken_assigned_reads,
  added_reads, new_est_reads, fraction_total_reads
Per-taxon сигнал — переоценённые прочтения new_est_reads (float).
Пустой bracken (только заголовок — edge case «нулевой микробиом»)
=> пустой профиль: TSV только с заголовком, status="empty".

abundance = reads / input_reads (0 при input_reads = 0).
rank: коды Kraken2 мапятся в человекочитаемые (S -> species, ...).
"""

import argparse
import csv
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from profile_io import CONTRACT_VERSION, REQUIRED_TSV_COLUMNS  # noqa: E402

RANK_MAP = {"D": "domain", "K": "kingdom", "P": "phylum", "C": "class",
            "O": "order", "F": "family", "G": "genus", "S": "species"}


def parse_kraken_report(path):
    """kraken2 report -> счётчики sidecar."""
    total = classified = unclassified = 0
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 6:
                continue
            try:
                reads = int(f[1])
            except ValueError:
                continue
            code = f[3].strip()
            if code == "U":
                unclassified += reads
                total += reads
            elif code == "R":
                classified += reads
                total += reads
    return {"input_reads": total,
            "classified_reads": classified,
            "unclassified_reads": unclassified}


def parse_bracken_taxa(path):
    """Bracken TSV -> [(taxid, name, rank, reads)] с reads > 0."""
    taxa = []
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        header = reader.fieldnames or []
        for col in ("taxonomy_id", "name", "taxonomy_lvl", "new_est_reads"):
            if col not in header:
                sys.exit(f"Bracken-файл {path}: нет обязательной колонки "
                         f"{col!r}")
        for row in reader:
            try:
                reads = float(row["new_est_reads"])
            except (TypeError, ValueError):
                sys.exit(f"Ошибка парсинга Bracken-файла {path}: "
                         f"new_est_reads={row.get('new_est_reads')!r}")
            if reads > 0:
                taxa.append((row["taxonomy_id"].strip(),
                             row["name"].strip(),
                             RANK_MAP.get(row["taxonomy_lvl"].strip(),
                                          row["taxonomy_lvl"].strip()
                                          .lower()),
                             reads))
    taxa.sort(key=lambda t: t[3], reverse=True)
    return taxa


def fmt_num(x):
    """Компактная запись числа: 3794.0 -> '3794', 0.5 -> '0.5'."""
    return f"{x:.6f}".rstrip("0").rstrip(".") if x % 1 else f"{x:.0f}"


def main():
    ap = argparse.ArgumentParser(
        description="Kraken2 report + Bracken TSV -> profile v2 "
                    "(TSV + sidecar JSON)")
    ap.add_argument("--report", required=True, help="kraken2 report (TSV)")
    ap.add_argument("--bracken", required=True,
                    help="bracken TSV (источник per-taxon reads)")
    ap.add_argument("--sample", required=True)
    ap.add_argument("--fraction", required=True,
                    choices=["dehost", "host", "half", "all"])
    ap.add_argument("--role", required=True,
                    help="роль БД из dbsheet (full|patho)")
    ap.add_argument("--db_id", required=True, help="идентификатор БД из dbsheet")
    ap.add_argument("--annotator", required=True,
                    help="идентификатор адаптера (kraken2_bracken)")
    ap.add_argument("--annotator_version", required=True,
                    help="версии, напр. '2.17.2+3.0.1'")
    ap.add_argument("--out_tsv", required=True)
    ap.add_argument("--out_json", required=True)
    args = ap.parse_args()

    if not os.path.exists(args.report):
        sys.exit(f"kraken2 report не найден: {args.report}")
    if not os.path.exists(args.bracken):
        sys.exit(f"bracken TSV не найден: {args.bracken}")

    m = parse_kraken_report(args.report)
    taxa = parse_bracken_taxa(args.bracken)
    input_reads = m["input_reads"]
    status = "ok" if taxa else "empty"

    with open(args.out_tsv, "w", newline="") as out:
        w = csv.writer(out, delimiter="\t", lineterminator="\n")
        w.writerow(REQUIRED_TSV_COLUMNS)
        for taxid, name, rank, reads in taxa:
            abundance = reads / input_reads if input_reads else 0.0
            w.writerow([args.sample, args.fraction, args.annotator,
                        args.annotator_version, args.db_id, taxid, name,
                        rank, fmt_num(reads), fmt_num(abundance)])

    sidecar = {
        "contract_version": CONTRACT_VERSION,
        "sample_id": args.sample,
        "fraction": args.fraction,
        "role": args.role,
        "db_id": args.db_id,
        "annotator": args.annotator,
        "annotator_version": args.annotator_version,
        **m,
        "status": status,
    }
    with open(args.out_json, "w") as out:
        json.dump(sidecar, out, indent=2, ensure_ascii=False)

    print(f"[kraken2_bracken_to_profile] {args.sample}/{args.fraction}/"
          f"{args.role}: status={status}, taxa={len(taxa)}, "
          f"input_reads={input_reads}")


if __name__ == "__main__":
    main()
