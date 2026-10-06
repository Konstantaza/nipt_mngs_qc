#!/usr/bin/env python3
"""
compare_profiles.py — «comparison model» пайплайна nipt_mngs_qc.

Сравнивает таксономические профили DEHOST- и HOST-фракций одного образца
(выходы Bracken на уровне видов) и строит дифференциальный профиль.

Научная логика
--------------
HOST-фракция — это прочтения, которые выровнялись на геном человека. Если
какой-то «микробный» таксон находится и там, это почти наверняка
миссклассифицированная человеческая ДНК (консервативные/повторяющиеся
участки) или контаминация. То есть HOST-профиль играет роль встроенного
негативного контроля конкретного образца.

Эвристика обогащения: таксон считается достоверно присутствующим
(«enriched»), если
  1) в DEHOST-профиле у него >= min_reads прочтений (по оценке Bracken), И
  2) его доля в DEHOST превышает долю в HOST минимум в enrichment_fold раз
     (если в HOST таксона нет вовсе — обогащение считается бесконечным).

Выходы: TSV с полным дифференциальным профилем + компактный JSON-сводка
(используется QC gate).
"""

import argparse
import json
import sys


def parse_bracken(path):
    """Читает выход Bracken (уровень S) -> dict taxid -> record."""
    taxa = {}
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        # Ожидаемые колонки Bracken:
        # name, taxonomy_id, taxonomy_lvl, kraken_assigned_reads,
        # added_reads, new_est_reads, fraction_total_reads
        idx = {name: i for i, name in enumerate(header)}
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            f = line.split("\t")
            try:
                taxid = f[idx["taxonomy_id"]]
                taxa[taxid] = {
                    "name": f[idx["name"]],
                    "new_est_reads": float(f[idx["new_est_reads"]]),
                    "fraction": float(f[idx["fraction_total_reads"]]),
                }
            except (KeyError, ValueError, IndexError) as e:
                sys.exit(f"Ошибка парсинга Bracken-файла {path}: {e}\nСтрока: {line}")
    return taxa


def main():
    ap = argparse.ArgumentParser(description="Differential profile: dehost vs host (Bracken)")
    ap.add_argument("--dehost", required=True, help="Bracken-файл DEHOST-фракции")
    ap.add_argument("--host", required=True, help="Bracken-файл HOST-фракции")
    ap.add_argument("--sample", required=True)
    ap.add_argument("--db", required=True, choices=["full", "patho"],
                    help="Тип профилирования (полная БД или панель патогенов)")
    ap.add_argument("--min_reads", type=float, default=10,
                    help="Минимум прочтений в DEHOST для рассмотрения таксона")
    ap.add_argument("--enrichment_fold", type=float, default=10,
                    help="Минимальное обогащение dehost/host для статуса enriched")
    ap.add_argument("--out_tsv", required=True)
    ap.add_argument("--out_json", required=True)
    args = ap.parse_args()

    dehost = parse_bracken(args.dehost)
    host = parse_bracken(args.host)

    rows = []
    n_enriched = 0
    for taxid, d in sorted(dehost.items(),
                           key=lambda kv: kv[1]["new_est_reads"], reverse=True):
        if d["new_est_reads"] < args.min_reads:
            continue  # таксоны с единичными прочтениями не рассматриваем
        h = host.get(taxid, {"new_est_reads": 0.0, "fraction": 0.0})
        if h["fraction"] > 0:
            fold = d["fraction"] / h["fraction"]
        else:
            fold = float("inf")  # в host-фракции таксона нет вовсе
        status = "enriched" if fold >= args.enrichment_fold else "background"
        if status == "enriched":
            n_enriched += 1
        rows.append({
            "taxid": taxid,
            "name": d["name"],
            "dehost_reads": d["new_est_reads"],
            "host_reads": h["new_est_reads"],
            "dehost_frac": d["fraction"],
            "host_frac": h["fraction"],
            "enrichment_fold": fold,
            "status": status,
        })

    with open(args.out_tsv, "w") as out:
        out.write("taxid\tname\tdehost_reads\thost_reads\tdehost_frac\t"
                  "host_frac\tenrichment_fold\tstatus\n")
        for r in rows:
            fold = "Inf" if r["enrichment_fold"] == float("inf") \
                else f'{r["enrichment_fold"]:.2f}'
            out.write(f'{r["taxid"]}\t{r["name"]}\t{r["dehost_reads"]:.0f}\t'
                      f'{r["host_reads"]:.0f}\t{r["dehost_frac"]:.6f}\t'
                      f'{r["host_frac"]:.6f}\t{fold}\t{r["status"]}\n')

    summary = {
        "sample_id": args.sample,
        "db": args.db,
        "n_taxa_considered": len(rows),
        "n_enriched": n_enriched,
        "n_background": len(rows) - n_enriched,
        "enriched_taxa": [
            {"taxid": r["taxid"], "name": r["name"],
             "dehost_reads": r["dehost_reads"],
             "enrichment_fold": (None if r["enrichment_fold"] == float("inf")
                                 else round(r["enrichment_fold"], 2))}
            for r in rows if r["status"] == "enriched"
        ],
    }
    with open(args.out_json, "w") as out:
        json.dump(summary, out, indent=2, ensure_ascii=False)


if __name__ == "__main__":
    main()
