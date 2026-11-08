#!/usr/bin/env python3
"""
metax_to_profile.py — конвертер реального выхода MetaX в profile v2.

NIPT-16: второй адаптер интерфейса ANNOTATE поверх классификатора
MetaX (hzi-bifo/Metax). Процесс ANNOTATE_METAX запускает
`metax profile`, а этот скрипт переводит его таблицу в контракт
profile v2 (docs/profile_contract.md): profile.tsv + sidecar
profile.json. native-отчётом служит сам сырой выход MetaX (его
копирует процесс; --out_native здесь нет).

Формат входа (проверен эмпирически на MetaX 2.0.3 и по документации
закреплённого релиза v2.0.3): tab-delimited, БЕЗ заголовка,
13 колонок (16 с --host):
  1  name            5  depth             9   cov_prob
  2  taxid           6  abundance_pct     10  fixed_chunk
  3  rank            7  breadth           11  flex_chunk
  4  reads           8  expected_breadth  12  expected_flex_chunk
                     (%.5f)               13  chunk_prob
  с --host: 14 host_names, 15 host_taxids, 16 diseases
  ("unknown", если таблица патоген-хозяин не задана).
Вероятности (9, 13): "NA" либо число — переносим строкой как есть.

ВАЖНО про колонку 6 (отклонение от буквы постановки, согласовано
с контрактом): MetaX РЕНОРМИРУЕТ abundance на сумму 100 (проценты
по прошедшим фильтры таксонам) — это НЕ доля от input_reads, и прямой
перенос нарушил бы контракт (abundance in [0,1], доля от input).
Поэтому:
  * contract abundance = reads / input_reads (вычисляем);
  * колонка 6 сохраняется как x_abundance_pct (строкой, для
    диагностики; потребители контракта её игнорируют).

classified_reads = min(input_reads, ceil(sum(reads))) — семантика
«прочтений, вошедших в отчёт» (post-filter). MetaX консервирует массу
EM, поэтому на ПОЛНОМ индексе sum(reads) <= input_reads; дробный
(fractional) индекс масштабирует counts и сломал бы инвариант —
адаптер требует полный индекс (metax_prepare.py проверяет метадату);
конвертер отвергает counts > input ещё до VALIDATE_PROFILE.

Пустой profile (0 строк — все отфильтровано или пустой вход; процесс
на пустом входе MetaX не вызывает, а создаёт пустой файл) =>
status="empty", TSV только с заголовком.

Совместимость: python >= 3.9; в процессе явно закреплён Python 3.14.7.
"""

import argparse
import csv
import json
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from profile_io import CONTRACT_VERSION, REQUIRED_TSV_COLUMNS  # noqa: E402

N_COLS_COMMUNITY = 13
N_COLS_HOST = 16
X_COLS = ("x_abundance_pct", "x_cov_prob", "x_chunk_prob")
X_HOST_COLS = ("x_host_names", "x_host_taxids", "x_diseases")


def fail(msg):
    sys.exit(f"[metax_to_profile] ОШИБКА: {msg}")


def count_fasta_records(path):
    """Число записей FASTA (строки '>')."""
    n = 0
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                n += 1
    return n


def parse_metax_profile(path):
    """Сырой profile.txt MetaX -> (rows, ncols).

    rows: [(name, taxid, rank, reads_str, reads_float,
            x_abundance_pct, x_cov_prob, x_chunk_prob,
            x_host_names|None, x_host_taxids|None, x_diseases|None)]
    """
    if not os.path.exists(path):
        fail(f"выход MetaX не найден: {path} "
             f"(процесс обязан передать metax_run[.pathogen].profile.txt)")
    rows = []
    ncols_seen = None
    with open(path, newline="") as fh:
        for n, raw in enumerate(fh, start=1):
            line = raw.rstrip("\n")
            if not line.strip():
                continue
            cols = line.split("\t")
            if len(cols) not in (N_COLS_COMMUNITY, N_COLS_HOST):
                fail(f"{path}:{n}: ожидается {N_COLS_COMMUNITY} колонок "
                     f"({N_COLS_HOST} с --host), получено {len(cols)}")
            if ncols_seen is None:
                ncols_seen = len(cols)
            elif len(cols) != ncols_seen:
                fail(f"{path}:{n}: смешанное число колонок "
                     f"({len(cols)} после {ncols_seen})")
            name, taxid, rank = cols[0].strip(), cols[1].strip(), \
                cols[2].strip()
            if not name or not taxid or not rank:
                fail(f"{path}:{n}: пустые name/taxid/rank недопустимы")
            reads_str = cols[3].strip()
            try:
                reads = float(reads_str)
            except ValueError:
                fail(f"{path}:{n}: reads не число: {reads_str!r}")
            if not math.isfinite(reads) or reads < 0:
                fail(f"{path}:{n}: reads должно быть конечным >= 0 ({reads})")
            host = cols[13:16] if len(cols) == N_COLS_HOST else \
                (None, None, None)
            rows.append((name, taxid, rank, reads_str, reads,
                         cols[5].strip(), cols[8].strip(), cols[12].strip(),
                         host[0], host[1], host[2]))
    return rows, (ncols_seen or N_COLS_COMMUNITY)


def main():
    ap = argparse.ArgumentParser(
        description="MetaX profile.txt -> profile v2 (TSV + sidecar)")
    ap.add_argument("--profile", required=True,
                    help="сырой выход MetaX: metax_run.profile.txt либо "
                         "metax_run.pathogen.profile.txt (правило суффикса "
                         ".pathogen при --host кодирует процесс)")
    ap.add_argument("--r1", required=True)
    ap.add_argument("--r2", required=True)
    ap.add_argument("--se", required=True)
    ap.add_argument("--single_end", required=True,
                    help="'true' — читаем se; иначе пары = записи r1")
    ap.add_argument("--sample", required=True)
    ap.add_argument("--fraction", required=True,
                    choices=["dehost", "host", "half", "all"])
    ap.add_argument("--role", required=True)
    ap.add_argument("--db_id", required=True)
    ap.add_argument("--annotator", required=True)
    ap.add_argument("--annotator_version", required=True)
    ap.add_argument("--out_tsv", required=True)
    ap.add_argument("--out_json", required=True)
    args = ap.parse_args()

    if args.single_end == "true":
        input_reads = count_fasta_records(args.se)
    else:
        input_reads = count_fasta_records(args.r1)

    rows, ncols = parse_metax_profile(args.profile)
    with_host = ncols == N_COLS_HOST

    status = "ok" if rows else "empty"
    if status == "ok" and input_reads == 0:
        fail(f"{args.profile}: {len(rows)} строк при input_reads=0 — "
             f"несогласованность входа и профиля")

    total_reads = sum(r[4] for r in rows)
    # Small tolerance for rounding the native %.5f per-taxon counts.
    if total_reads > input_reads + len(rows) * 0.000005 + 1e-8:
        fail(f"{args.profile}: sum(reads)={total_reads:g} > "
             f"input_reads={input_reads}; проверьте полный индекс и вход")
    classified = min(input_reads, int(math.ceil(total_reads))) \
        if status == "ok" else 0
    unclassified = input_reads - classified

    header = list(REQUIRED_TSV_COLUMNS) + list(X_COLS)
    if with_host:
        header += X_HOST_COLS

    # детерминированный порядок: по reads убыв., затем по taxid
    rows.sort(key=lambda r: (-r[4], r[1]))

    with open(args.out_tsv, "w", newline="") as out:
        w = csv.writer(out, delimiter="\t", lineterminator="\n")
        w.writerow(header)
        for (name, taxid, rank, reads_str, reads,
             x_ab, x_cov, x_chunk, h_names, h_taxids, h_dis) in rows:
            row = [args.sample, args.fraction, args.annotator,
                   args.annotator_version, args.db_id,
                   taxid, name, rank, reads_str,
                   "%.6g" % (reads / input_reads),
                   x_ab, x_cov, x_chunk]
            if with_host:
                row += [h_names, h_taxids, h_dis]
            w.writerow(row)

    sidecar = {
        "contract_version": CONTRACT_VERSION,
        "sample_id": args.sample,
        "fraction": args.fraction,
        "role": args.role,
        "db_id": args.db_id,
        "annotator": args.annotator,
        "annotator_version": args.annotator_version,
        "input_reads": input_reads,
        "classified_reads": classified,
        "unclassified_reads": unclassified,
        "status": status,
    }
    with open(args.out_json, "w") as out:
        json.dump(sidecar, out, indent=2, ensure_ascii=False)

    print(f"[metax_to_profile] {args.sample}/{args.fraction}/{args.role}: "
          f"status={status}, input_reads={input_reads}, "
          f"reported={total_reads:g}, rows={len(rows)}")


if __name__ == "__main__":
    main()
