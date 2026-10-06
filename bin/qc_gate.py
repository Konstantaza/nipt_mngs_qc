#!/usr/bin/env python3
"""
qc_gate.py — QC gate пайплайна nipt_mngs_qc.

Собирает метрики всех стадий одного образца, проверяет формальный чеклист
QC/адекватности метагеномной части и выносит вердикт:

  PASS_POSITIVE — QC пройден, в patho-профиле есть обогащённые таксоны
  PASS_NEGATIVE — QC пройден, патогены не детектированы (валидный «ноль»)
  FALLBACK      — образец не прошёл QC: результат НЕ отдаём, уходит
                  по отдельному пути обработки (коды причин в отчёте)
  FAIL          — входные данные неадекватны (анализ невозможен в принципе)

Принцип edge-case: «нулевой микробиом» / провал анализа — это НЕ пустой
отчёт, а явный статус FALLBACK с машиночитаемой причиной.

Выходы:
  <sample>_qc_report.json — контракт для внешней системы отчётов
  <sample>_qc_mqc.tsv     — custom content для MultiQC
"""

import argparse
import glob
import json
import os
import sys

HUMAN_TAXID = "9606"


def parse_flagstat(path):
    """samtools flagstat -> dict с ключевыми числами (paired-end)."""
    total = mapped = 0
    with open(path) as fh:
        for line in fh:
            if "in total" in line:
                total = int(line.split("+")[0].strip())
            elif ") mapped (" in line and "primary" not in line:
                mapped = int(line.split("+")[0].strip())
    return {"total_reads": total, "mapped_reads": mapped}


def parse_kraken_report(path):
    """Kraken2 report -> (доля unclassified, доля human, всего прочтений)."""
    unclass_frac = human_frac = 0.0
    total = 0
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 6:
                continue
            pct = float(f[0])
            reads = int(f[1])
            code = f[3]
            taxid = f[4].strip()
            if code == "U":
                unclass_frac = pct / 100.0
                total += reads
            elif code == "R":
                total += reads
            if taxid == HUMAN_TAXID:
                human_frac = pct / 100.0
    return unclass_frac, human_frac, total


def find_file(pattern, where="."):
    hits = sorted(glob.glob(os.path.join(where, pattern)))
    return hits[0] if hits else None


def main():
    ap = argparse.ArgumentParser(description="Per-sample QC gate for NIPT metagenomics")
    ap.add_argument("--sample", required=True)
    ap.add_argument("--fastp_json", default=None,
                    help="fastp JSON (или пустой stub для entry=bam)")
    ap.add_argument("--flagstat_dehost", required=True)
    ap.add_argument("--flagstat_host", required=True)
    ap.add_argument("--kraken_dir", required=True,
                    help="Директория с 4 kraken2 report (имена *.<dehost|host>.<full|patho>.kraken2.report)")
    ap.add_argument("--diff_dir", required=True,
                    help="Директория с 2 JSON comparison model (*.<full|patho>.diff.json)")
    # --- пороги (все приходят из nextflow.config) ---
    ap.add_argument("--min_input_pairs", type=float, default=5_000_000)
    ap.add_argument("--min_q30", type=float, default=0.80)
    ap.add_argument("--min_trim_survival", type=float, default=0.90)
    ap.add_argument("--host_frac_min", type=float, default=0.95)
    ap.add_argument("--min_nonhost_reads", type=float, default=50_000)
    ap.add_argument("--max_dehost_human_frac", type=float, default=0.05)
    ap.add_argument("--min_classified_frac", type=float, default=0.0001)
    ap.add_argument("--max_consistency_dev", type=float, default=0.01)
    ap.add_argument("--out_json", required=True)
    ap.add_argument("--out_mqc", required=True)
    args = ap.parse_args()

    criteria = []   # (id, value, threshold_str, passed, severity, reason_code)
    metrics = {}

    def add(cid, value, thr, passed, severity, reason=None, desc=""):
        criteria.append({"id": cid, "value": value, "threshold": thr,
                         "pass": passed, "severity_on_fail": severity,
                         "reason_code": reason if not passed else None,
                         "description": desc})

    # ---------- 1. Метрики входа (flagstat) ----------
    fs_d = parse_flagstat(args.flagstat_dehost)
    fs_h = parse_flagstat(args.flagstat_host)
    dehost_reads = fs_d["total_reads"]
    host_reads = fs_h["total_reads"]
    total_reads = dehost_reads + host_reads
    input_pairs = total_reads / 2.0
    host_frac = host_reads / total_reads if total_reads else 0.0
    nonhost_reads = dehost_reads / 2.0
    metrics.update({"input_pairs": input_pairs, "dehost_reads": dehost_reads,
                    "host_reads": host_reads, "host_fraction": round(host_frac, 6),
                    "nonhost_pairs": nonhost_reads})

    add("input_pairs", input_pairs, f">= {args.min_input_pairs:.0f}",
        input_pairs >= args.min_input_pairs, "FAIL", "LOW_INPUT",
        "Достаточность входа: число пар прочтений")

    add("host_fraction", round(host_frac, 4), f">= {args.host_frac_min}",
        host_frac >= args.host_frac_min, "FALLBACK", "HOST_FRACTION_LOW",
        "Доля хост-прочтений в ожидаемом диапазоне (ниже — аномалия wet-lab)")

    add("nonhost_pairs", nonhost_reads, f">= {args.min_nonhost_reads:.0f}",
        nonhost_reads >= args.min_nonhost_reads, "FALLBACK", "ZERO_MICROBIOME",
        "Достаточность не-хост сигнала (edge case: нулевой микробиом)")

    # ---------- 2. fastp (только entry=fastq) ----------
    fastp_ok = False
    if args.fastp_json and os.path.exists(args.fastp_json):
        try:
            with open(args.fastp_json) as fh:
                fj = json.load(fh)
            before = fj["summary"]["before_filtering"]["total_reads"]
            after = fj["summary"]["after_filtering"]["total_reads"]
            q30 = fj["summary"]["before_filtering"]["q30_rate"]
            survival = after / before if before else 0.0
            fastp_ok = True
            metrics.update({"q30_rate": round(q30, 4),
                            "trim_survival": round(survival, 4),
                            "fastp_total_before": before})
            add("q30_rate", round(q30, 4), f">= {args.min_q30}",
                q30 >= args.min_q30, "FALLBACK", "LOW_QUALITY",
                "Качество сырых прочтений (Q30)")
            add("trim_survival", round(survival, 4), f">= {args.min_trim_survival}",
                survival >= args.min_trim_survival, "FALLBACK", "LOW_QUALITY",
                "Доля пар, переживших тримминг")
            # консистентность: dehost+host против выхода fastp
            dev = abs(total_reads - after) / after if after else 1.0
            metrics["consistency_dev"] = round(dev, 6)
            add("consistency", round(dev, 6), f"<= {args.max_consistency_dev}",
                dev <= args.max_consistency_dev, "FALLBACK", "SPLIT_INCONSISTENT",
                "Консистентность разбиения dehost+host против trimmed")
        except (KeyError, json.JSONDecodeError):
            pass
    if not fastp_ok:
        for cid in ("q30_rate", "trim_survival", "consistency"):
            criteria.append({"id": cid, "value": None, "threshold": None,
                             "pass": None, "severity_on_fail": None,
                             "reason_code": None,
                             "description": "Не оценивается (entry=bam, fastp upstream)"})

    # ---------- 3. Kraken2-санити ----------
    kr_d_full = find_file(f"*dehost.full*.report", args.kraken_dir) or \
        find_file("*dehost*full*", args.kraken_dir)
    if kr_d_full:
        uncl, human, _ = parse_kraken_report(kr_d_full)
        classified = 1.0 - uncl
        metrics.update({"dehost_classified_frac": round(classified, 6),
                        "dehost_human_frac": round(human, 6)})
        add("dehost_human_frac", round(human, 4), f"<= {args.max_dehost_human_frac}",
            human <= args.max_dehost_human_frac, "FALLBACK", "DEHOST_FAILURE",
            "Утечка человеческой ДНК в dehost-фракцию (сбой dehost)")
        add("classified_frac", round(classified, 6), f">= {args.min_classified_frac}",
            classified >= args.min_classified_frac, "FALLBACK", "NO_CLASSIFIED_SIGNAL",
            "Санити: классифицируемая доля dehost-прочтений > 0")
    else:
        add("kraken_report", None, "present", False, "FAIL", "PIPELINE_ERROR",
            "Kraken2 report для dehost/full не найден")

    # ---------- 4. Детекция патогенов (comparison model) ----------
    diff_patho = find_file("*.patho.diff.json", args.diff_dir)
    n_patho_enriched = 0
    if diff_patho:
        with open(diff_patho) as fh:
            dp = json.load(fh)
        n_patho_enriched = dp.get("n_enriched", 0)
        metrics["patho_enriched_taxa"] = n_patho_enriched
        metrics["patho_enriched_names"] = [t["name"] for t in dp.get("enriched_taxa", [])]

    # ---------- 5. Вердикт ----------
    failed = [c for c in criteria if c["pass"] is False]
    hard_fail = [c for c in failed if c["severity_on_fail"] == "FAIL"]
    soft_fail = [c for c in failed if c["severity_on_fail"] == "FALLBACK"]
    if hard_fail:
        status = "FAIL"
    elif soft_fail:
        status = "FALLBACK"
    else:
        status = "PASS_POSITIVE" if n_patho_enriched > 0 else "PASS_NEGATIVE"
    reasons = sorted({c["reason_code"] for c in failed if c["reason_code"]})

    report = {
        "sample_id": args.sample,
        "pipeline": "nipt_mngs_qc",
        "pipeline_version": "1.0.0",
        "status": status,
        "fallback_reasons": reasons,
        "criteria": criteria,
        "metrics": metrics,
    }
    with open(args.out_json, "w") as out:
        json.dump(report, out, indent=2, ensure_ascii=False)

    # ---------- 6. MultiQC custom content ----------
    with open(args.out_mqc, "w") as out:
        out.write("# id: 'nipt_mngs_qc'\n")
        out.write("# section_name: 'Metagenome QC gate'\n")
        out.write("# description: 'Per-sample QC criteria for the NIPT metagenomic adjunct'\n")
        out.write("# plot_type: 'table'\n")
        out.write("Sample\tCriterion\tValue\tThreshold\tPass\n")
        for c in criteria:
            val = "NA" if c["value"] is None else c["value"]
            thr = "NA" if c["threshold"] is None else c["threshold"]
            p = "NA" if c["pass"] is None else ("PASS" if c["pass"] else "FAIL")
            out.write(f"{args.sample}\t{c['id']}\t{val}\t{thr}\t{p}\n")
        out.write(f"{args.sample}\tFINAL_STATUS\t{status}\t-\t"
                  f"{'PASS' if status.startswith('PASS') else 'FAIL'}\n")

    print(f"[qc_gate] {args.sample}: {status}"
          + (f" ({', '.join(reasons)})" if reasons else ""))


if __name__ == "__main__":
    main()
