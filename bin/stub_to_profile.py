#!/usr/bin/env python3
"""
stub_to_profile.py — fixture-driven stub-адаптер интерфейса ANNOTATE.

Назначение (NIPT-16): тестовый аннотатор без внешней БД и без
классификатора. Доказывает расширяемость реестра аннотаторов: добавление
аннотатора = процесс + конвертер + строка в dbsheet, без правок main.nf
и потребителей (compare, qc_gate).

Источник «классификации» — FIXTURE-файл, адресуемый через dbsheet:
процесс передаёт --fixture <path_из_dbsheet>/<db_id>.tsv. Формат fixture
(TSV с заголовком):

    taxid<TAB>name<TAB>rank<TAB>weight
    2097<TAB>Mycoplasma genitalium<TAB>species<TAB>1.0

weight — доля входных прочтений, относимая к таксону (0..1); сумма
weight <= 1 (остаток уходит в unclassified). Так fixture может
воспроизводить любой детерминированный профиль, в т.ч. с неполной
классификацией.

Семантика контракта profile v2 (docs/profile_contract.md):
  * input_reads  — число записей входного FASTA (PE: записи r1 —
    пара = одно событие, как kraken2 --paired; SE: записи se);
  * reads_i      = input_reads * weight_i (формат %.6g);
  * abundance_i  = weight_i (доля от input_reads — по определению
    контракта);
  * classified   = round(input_reads * sum(weight));
  * unclassified = input_reads - classified;
  * status       = ok, если есть строки и input_reads > 0, иначе empty
    (пустой профиль — TSV только с заголовком).

Ошибки fixture (нет файла, неверный заголовок, weight вне [0,1],
сумма > 1, пустые обязательные поля) — фатальны: адаптер обязан
упасть громко, а не эмитить сомнительный профиль.

native-отчёт — текстовый дамп stub'а (MultiQC его не распознаёт и
молча игнорирует).
"""

import argparse
import csv
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from profile_io import CONTRACT_VERSION, REQUIRED_TSV_COLUMNS  # noqa: E402

FIXTURE_HEADER = ("taxid", "name", "rank", "weight")


def fail(msg):
    sys.exit(f"[stub_to_profile] ОШИБКА: {msg}")


def count_fasta_records(path):
    """Число записей FASTA (строки '>')."""
    n = 0
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                n += 1
    return n


def load_fixture(path):
    """fixture TSV -> [(taxid, name, rank, weight)]. Строгая валидация."""
    if not os.path.exists(path):
        fail(f"fixture не найден: {path} "
             f"(ожидается <path_из_dbsheet>/<db_id>.tsv)")
    rows = []
    with open(path, newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader, None)
        if header is None or tuple(h.strip() for h in header) != FIXTURE_HEADER:
            fail(f"{path}: заголовок должен быть "
                 f"{'|'.join(FIXTURE_HEADER)}, получено {header!r}")
        for n, row in enumerate(reader, start=2):
            if not row or all(not c.strip() for c in row):
                continue  # пустые строки допустимы
            if len(row) != 4:
                fail(f"{path}:{n}: ожидается 4 колонки, получено {len(row)}")
            taxid, name, rank, weight_s = (c.strip() for c in row)
            if not taxid or not name or not rank:
                fail(f"{path}:{n}: пустые taxid/name/rank недопустимы")
            try:
                weight = float(weight_s)
            except ValueError:
                fail(f"{path}:{n}: weight не число: {weight_s!r}")
            if not 0.0 <= weight <= 1.0:
                fail(f"{path}:{n}: weight вне [0,1]: {weight}")
            rows.append((taxid, name, rank, weight))
    total = sum(r[3] for r in rows)
    if total > 1.0 + 1e-9:
        fail(f"{path}: сумма weight ({total:g}) > 1 — доли входных "
             f"прочтений не могут суммироваться выше единицы")
    return rows


def main():
    ap = argparse.ArgumentParser(
        description="stub-аннотатор: fixture TSV + FASTA -> profile v2 "
                    "(TSV + sidecar + native)")
    ap.add_argument("--fixture", required=True,
                    help="TSV: taxid, name, rank, weight "
                         "(<path_из_dbsheet>/<db_id>.tsv)")
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
    ap.add_argument("--out_native", required=True)
    args = ap.parse_args()

    fixture = load_fixture(args.fixture)

    if args.single_end == "true":
        input_reads = count_fasta_records(args.se)
    else:
        input_reads = count_fasta_records(args.r1)

    weight_sum = sum(r[3] for r in fixture)
    classified = int(round(input_reads * weight_sum))
    unclassified = input_reads - classified
    status = "ok" if (fixture and input_reads > 0) else "empty"

    # детерминированный порядок: по reads убыв., затем по taxid
    rows = sorted(fixture, key=lambda r: (-r[3], r[0]))

    with open(args.out_tsv, "w", newline="") as out:
        w = csv.writer(out, delimiter="\t", lineterminator="\n")
        w.writerow(REQUIRED_TSV_COLUMNS)
        if status == "ok":
            for taxid, name, rank, weight in rows:
                w.writerow([args.sample, args.fraction, args.annotator,
                            args.annotator_version, args.db_id,
                            taxid, name, rank,
                            "%.6g" % (input_reads * weight),
                            "%.6g" % weight])

    sidecar = {
        "contract_version": CONTRACT_VERSION,
        "sample_id": args.sample,
        "fraction": args.fraction,
        "role": args.role,
        "db_id": args.db_id,
        "annotator": args.annotator,
        "annotator_version": args.annotator_version,
        "input_reads": input_reads,
        "classified_reads": classified if status == "ok" else 0,
        "unclassified_reads": unclassified if status == "ok" else input_reads,
        "status": status,
    }
    with open(args.out_json, "w") as out:
        json.dump(sidecar, out, indent=2, ensure_ascii=False)

    with open(args.out_native, "w") as out:
        out.write("# stub annotator report\n")
        out.write(f"sample: {args.sample}\n")
        out.write(f"fraction: {args.fraction}\n")
        out.write(f"role: {args.role}\n")
        out.write(f"db_id: {args.db_id}\n")
        out.write(f"fixture: {args.fixture}\n")
        out.write(f"input_reads: {input_reads}\n")
        out.write(f"classified_reads: {sidecar['classified_reads']}\n")
        for taxid, name, rank, weight in rows:
            out.write(f"taxon: {name} (taxid {taxid}, rank {rank}, "
                      f"weight {weight:g})\n")

    print(f"[stub_to_profile] {args.sample}/{args.fraction}/{args.role}: "
          f"status={status}, input_reads={input_reads}, "
          f"fixture_rows={len(fixture)}")


if __name__ == "__main__":
    main()
