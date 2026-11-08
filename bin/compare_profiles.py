#!/usr/bin/env python3
"""
compare_profiles.py — сравнение профилей фракций (comparison model v2).

Отвечает на вопрос: «находим ли мы одно и то же в двух фракциях одного
образца?» (для dehost:host ожидаем различия). Работает ТОЛЬКО на контракте
profile v2 (profile.tsv + sidecar profile.json, docs/profile_contract.md):
форматов классификатора здесь нет и быть не должно — они инкапсулированы
в конвертерах адаптеров ANNOTATE, а чтение контракта — в profile_io.py.

Почему не доли внутри профиля (баг v1.0)
----------------------------------------
v1.0 считала обогащение как отношение ДОЛЕЙ таксона внутри каждого
профиля. Доли нормированы на размер своей фракции, а фракции различаются
на порядки (dehost ~1–5% образца). Поэтому при РАВНЫХ абсолютных числах
прочтений в dehost и host (типичная картина утечки/миссклассификации
человеческих повторов) fold раздувался до десятков тысяч и таксон
ошибочно получал статус enriched.

Модель v2
---------
* Нормировка — RPM на input_reads фракции (прочтения на входе
  классификации этой фракции), а не доля внутри профиля.
* Присутствие таксона — по абсолютным прочтениям: в A при
  reads >= min_reads, в B при reads >= host_min_reads (отдельный порог
  на абсолютный сигнал в host-фракции).
* Категории: both / only_A / only_B.
* Флаг enriched (строгий): present_A И НЕ present_B И
  (rpm_B == 0 ИЛИ rpm_A/rpm_B >= enrichment_fold). Общий таксон (both)
  никогда не enriched — равные абсолютные числа прочтений обогащения
  не дают.
* Сводка сходства: n_shared, n_only_A, n_only_B, Jaccard,
  Bray–Curtis (dissimilarity по RPM, 0 = идентичны), Spearman по RPM
  общих таксонов.

Пары фракций (--pair A:B): dehost:host (дефолт) или dehost:all;
«all» синтезируется как dehost+host (сплит — партиция прочтений образца).

Выходы:
  <...>.compare.tsv  — per-taxon таблица
  <...>.compare.json — сводка (читает QC gate)
"""

import argparse
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from profile_io import CONTRACT_VERSION, ProfileError, read_profile  # noqa: E402


def load_profile(tsv_path, json_path):
    """Пара profile.tsv + sidecar -> dict; строгая валидация контракта."""
    try:
        return read_profile(tsv_path, json_path)
    except ProfileError as e:
        sys.exit(str(e))


def synthesize_all(prof_a, prof_b):
    """Фракция 'all' = сумма двух фракций (сплит — партиция прочтений)."""
    taxa = {}
    for p in (prof_a, prof_b):
        for t in p["taxa"]:
            rec = taxa.setdefault(t["taxid"], {"taxid": t["taxid"],
                                               "name": t["name"],
                                               "reads": 0.0})
            rec["reads"] += t["reads"]
    return {
        "fraction": "all",
        "input_reads": prof_a["input_reads"] + prof_b["input_reads"],
        "taxa": sorted(taxa.values(), key=lambda t: t["reads"],
                       reverse=True),
    }


def ranks(xs):
    """Средние ранги (связи делят ранг)."""
    order = sorted(range(len(xs)), key=lambda i: xs[i])
    r = [0.0] * len(xs)
    i = 0
    while i < len(xs):
        j = i
        while j + 1 < len(xs) and xs[order[j + 1]] == xs[order[i]]:
            j += 1
        avg = (i + j) / 2.0 + 1.0
        for k in range(i, j + 1):
            r[order[k]] = avg
        i = j + 1
    return r


def pearson(xs, ys):
    n = len(xs)
    mx = sum(xs) / n
    my = sum(ys) / n
    sxx = sum((x - mx) ** 2 for x in xs)
    syy = sum((y - my) ** 2 for y in ys)
    if sxx == 0 or syy == 0:
        return None
    sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    return sxy / (sxx * syy) ** 0.5


def spearman(xs, ys):
    """Rho Спирмена; None при < 2 точках или нулевой дисперсии рангов."""
    if len(xs) < 2:
        return None
    return pearson(ranks(xs), ranks(ys))


def main():
    ap = argparse.ArgumentParser(
        description="Fraction profile comparison on profile v2 contracts")
    ap.add_argument("--dehost_tsv", required=True)
    ap.add_argument("--dehost_json", required=True,
                    help="sidecar JSON dehost-профиля")
    ap.add_argument("--host_tsv", required=True)
    ap.add_argument("--host_json", required=True,
                    help="sidecar JSON host-профиля")
    ap.add_argument("--sample", required=True)
    ap.add_argument("--role", required=True,
                    help="роль БД из dbsheet (full|patho)")
    ap.add_argument("--pair", default="dehost:host",
                    help="пара фракций A:B; 'all' = синтез dehost+host")
    ap.add_argument("--min_reads", type=float, default=10,
                    help="присутствие в A: reads >= min_reads")
    ap.add_argument("--host_min_reads", type=float, default=10,
                    help="присутствие в B: reads >= host_min_reads "
                         "(отдельный абсолютный порог для host-фракции)")
    ap.add_argument("--enrichment_fold", type=float, default=10,
                    help="минимальное RPM-обогащение A/B для флага enriched")
    ap.add_argument("--out_tsv", required=True)
    ap.add_argument("--out_json", required=True)
    args = ap.parse_args()

    frac_a, _, frac_b = args.pair.partition(":")
    if not frac_a or not frac_b:
        sys.exit(f"--pair должен быть вида A:B, получено: {args.pair!r}")

    prof_d = load_profile(args.dehost_tsv, args.dehost_json)
    prof_h = load_profile(args.host_tsv, args.host_json)
    if prof_d["fraction"] != "dehost" or prof_h["fraction"] != "host":
        sys.exit("Перепутаны фракции профилей: "
                 f"{prof_d['fraction']!r} и {prof_h['fraction']!r}")
    for key in ("sample_id", "role", "db_id", "annotator",
                "annotator_version"):
        if prof_d[key] != prof_h[key]:
            sys.exit(f"Несовпадение {key} между профилями: "
                     f"{prof_d[key]!r} != {prof_h[key]!r}")
    if prof_d["sample_id"] != args.sample or prof_d["role"] != args.role:
        sys.exit("Профили не соответствуют --sample/--role")

    profiles = {"dehost": prof_d, "host": prof_h}
    if frac_a == "all" or frac_b == "all":
        profiles["all"] = synthesize_all(prof_d, prof_h)
    if frac_a not in profiles or frac_b not in profiles:
        sys.exit(f"Пара {args.pair!r}: неизвестная фракция "
                 f"(доступны: dehost, host, all)")
    pa, pb = profiles[frac_a], profiles[frac_b]

    input_a = pa["input_reads"]
    input_b = pb["input_reads"]
    taxa_a = {t["taxid"]: t for t in pa["taxa"]}
    taxa_b = {t["taxid"]: t for t in pb["taxa"]}

    def rpm(reads, input_reads):
        return reads / input_reads * 1e6 if input_reads else 0.0

    rows = []
    for taxid in sorted(set(taxa_a) | set(taxa_b),
                        key=lambda t: -(taxa_a.get(t, {}).get("reads", 0.0)
                                        + taxa_b.get(t, {}).get("reads", 0.0))):
        ta = taxa_a.get(taxid)
        tb = taxa_b.get(taxid)
        reads_a = ta["reads"] if ta else 0.0
        reads_b = tb["reads"] if tb else 0.0
        present_a = reads_a >= args.min_reads
        present_b = reads_b >= args.host_min_reads
        if not (present_a or present_b):
            continue  # ниже порога в обеих фракциях — шум
        rpm_a = rpm(reads_a, input_a)
        rpm_b = rpm(reads_b, input_b)
        category = ("both" if present_b else "only_A") if present_a \
            else "only_B"
        fold = (rpm_a / rpm_b) if rpm_b > 0 else float("inf")
        # Строгий enriched: общий таксон (both) никогда не enriched —
        # равные абсолютные числа прочтений обогащения не дают (баг v1.0)
        enriched = (present_a and not present_b
                    and (rpm_b == 0 or fold >= args.enrichment_fold))
        rows.append({
            "taxid": taxid,
            "name": (ta or tb)["name"],
            "reads_a": reads_a, "reads_b": reads_b,
            "rpm_a": rpm_a, "rpm_b": rpm_b,
            "category": category,
            "enrichment_fold": fold,
            "enriched": enriched,
        })

    n_shared = sum(1 for r in rows if r["category"] == "both")
    n_only_a = sum(1 for r in rows if r["category"] == "only_A")
    n_only_b = sum(1 for r in rows if r["category"] == "only_B")
    union = n_shared + n_only_a + n_only_b
    jaccard = (n_shared / union) if union else None
    denom = sum(r["rpm_a"] + r["rpm_b"] for r in rows)
    bray_curtis = (sum(abs(r["rpm_a"] - r["rpm_b"]) for r in rows) / denom) \
        if denom else None
    shared = [r for r in rows if r["category"] == "both"]
    rho = spearman([r["rpm_a"] for r in shared],
                   [r["rpm_b"] for r in shared])
    enriched_rows = [r for r in rows if r["enriched"]]

    # ---------- per-taxon TSV ----------
    with open(args.out_tsv, "w") as out:
        out.write(f"taxid\tname\treads_{frac_a}\treads_{frac_b}\t"
                  f"rpm_{frac_a}\trpm_{frac_b}\tcategory\t"
                  "enrichment_fold\tenriched\n")
        for r in rows:
            fold = "Inf" if r["enrichment_fold"] == float("inf") \
                else f'{r["enrichment_fold"]:.2f}'
            out.write(f'{r["taxid"]}\t{r["name"]}\t'
                      f'{r["reads_a"]:.0f}\t{r["reads_b"]:.0f}\t'
                      f'{r["rpm_a"]:.2f}\t{r["rpm_b"]:.2f}\t'
                      f'{r["category"]}\t{fold}\t'
                      f'{"true" if r["enriched"] else "false"}\n')

    # ---------- сводка JSON ----------
    summary = {
        "sample_id": args.sample,
        "role": args.role,
        "db_id": prof_d["db_id"],
        "annotator": prof_d["annotator"],
        "annotator_version": prof_d["annotator_version"],
        "pair": args.pair,
        "input_reads": {frac_a: input_a, frac_b: input_b},
        "thresholds": {"min_reads": args.min_reads,
                       "host_min_reads": args.host_min_reads,
                       "enrichment_fold": args.enrichment_fold},
        "n_shared": n_shared,
        "n_only_A": n_only_a,
        "n_only_B": n_only_b,
        "jaccard": (round(jaccard, 4) if jaccard is not None else None),
        "bray_curtis_dissimilarity": (round(bray_curtis, 6)
                                      if bray_curtis is not None else None),
        "spearman_rho": (round(rho, 4) if rho is not None else None),
        "n_taxa_considered": union,
        "n_enriched": len(enriched_rows),
        "n_background": union - len(enriched_rows),
        "enriched_taxa": [
            {"taxid": r["taxid"], "name": r["name"],
             "reads": r["reads_a"], "rpm": round(r["rpm_a"], 2),
             "enrichment_fold": (None if r["enrichment_fold"] == float("inf")
                                 else round(r["enrichment_fold"], 2))}
            for r in enriched_rows
        ],
    }
    with open(args.out_json, "w") as out:
        json.dump(summary, out, indent=2, ensure_ascii=False)

    print(f"[compare] {args.sample}/{args.role} {args.pair}: "
          f"shared={n_shared} only_A={n_only_a} only_B={n_only_b} "
          f"enriched={len(enriched_rows)}")


if __name__ == "__main__":
    main()
