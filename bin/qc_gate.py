#!/usr/bin/env python3
"""
qc_gate.py — QC gate пайплайна nipt_mngs_qc (v1.9.0).

v1.9.0: второй адаптер аннотатора (NIPT-16) — реестр ANNOTATE пополнен
fixture-driven stub (тесты без БД) и реальным MetaX 2.0.3 (процесс +
конвертер bin/metax_to_profile.py); сменяемость доказана mixed-dbsheet
(kraken2 + stub в одном прогоне) без правок main.nf/compare/qc;
сам gate без изменений.

v1.8.0: воспроизводимость окружения (NIPT-7) — точные пины всех
инструментов в environment.yml, per-process директивы conda/container
в модулях, исправлен ретрай Kraken2 по памяти (удвоение через
MemoryUnit вместо умножения строки), process.resourceLimits,
manifest.nextflowVersion; сам gate без изменений.

v1.7.0: сквозной smoke-harness на toy-данных (test/run_smoke.sh +
test/check_smoke.py), самодостаточный профиль test (input/dbsheet/
host_index/панель внутри профиля) и docker-профиль с Dockerfile;
сам gate без изменений.

v1.6.0: fail-fast валидация входов перенесена в начало main.nf
(обязательные параметры по режиму, колонки и уникальность sample в
samplesheet/dbsheet, glob индекса *.ht2); сам gate без изменений.

v1.5.0: опциональный вход --counts_half (полукартированные пары DEHOST,
mode=fastq); критерий consistency проверяет точную партицию
dehost+host+half против trimmed; информационные метрики
half_mapped_reads / half_mapped_frac.

Собирает метрики одного образца и выносит вердикт:

  PASS_POSITIVE — QC пройден, обогащённые таксоны пересекли панель мишеней
  PASS_NEGATIVE — QC пройден, позитивности нет (см. positivity_source)
  FALLBACK      — мягкий критерий не пройден: результат НЕ отдаём,
                  образец уходит по fallback-ветке (коды причин в отчёте)
  FAIL          — вход неадекватен или отсутствует обязательный артефакт
                  (PIPELINE_ERROR): анализ невозможен/неполон

Принципы этой версии:
  * Форматов отчётов классификатора здесь НЕТ: метрики веток читаются из
    контракта profile v2 (profile.tsv + sidecar profile.json стадии
    ANNOTATE, docs/profile_contract.md; чтение — через profile_io.py),
    сводки сравнения фракций — из compare JSON comparison model v2.
    Размеры фракций — из counts JSON стадии BAM2FASTA (bam_counts.py:
    primary-only, -F 0x900; secondary и supplementary записи НЕ
    считаются, в отличие от строки «in total» flagstat). Парсится также
    fastp JSON.
  * Отсутствие обязательного артефакта — FAIL/PIPELINE_ERROR, а не молчаливый
    PASS_NEGATIVE и не пропажа образца.
  * Утечка human в dehost оценивается только при достаточном числе не-хост
    прочтений (nonhost_pairs >= min_nonhost_reads): на единицах прочтений
    доля human бессмысленна.
  * Позитивность — только по панели мишеней (--target_panel). Без панели
    позитивность не оценивается: PASS_NEGATIVE + positivity_source =
    "not_evaluated". Сводка сравнения dehost:host — информационная
    (comparison_info), на вердикт не влияет.

Выходы:
  <sample>_qc_report.json — контракт для внешней системы отчётов
  <sample>_qc_mqc.tsv     — custom content для MultiQC
"""

import argparse
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from profile_io import try_read_profile  # noqa: E402

EXPECTED_FRACTIONS = ("dehost", "host")
EXPECTED_ROLES = ("full", "patho")

HUMAN_TAXID = "9606"


COUNTS_REQUIRED_KEYS = ("primary_reads", "paired_reads", "unpaired_reads",
                        "single_end", "fragments")


def load_counts(path):
    """counts JSON (bam_counts.py) -> dict или None, если отсутствует/битый.

    Пустой/битый файл или JSON без обязательных ключей — это отсутствие
    обязательного артефакта (PIPELINE_ERROR), а не «0 прочтений».
    """
    data, err = load_json(path)
    if err:
        return None
    if not isinstance(data, dict) or \
            any(k not in data for k in COUNTS_REQUIRED_KEYS):
        return None
    return data


def load_json(path):
    """-> (data, error_str). error_str is None при успехе."""
    if not path or not os.path.exists(path) or os.path.getsize(path) == 0:
        return None, "файл отсутствует или пуст"
    try:
        with open(path) as fh:
            return json.load(fh), None
    except json.JSONDecodeError as e:
        return None, f"битый JSON: {e}"


def read_panel(path):
    """Файл панели мишеней -> (taxids:set, names:set) или None, если не задан.

    Одна запись на строку; числовая строка = taxid, иначе имя вида.
    Пустой/отсутствующий файл (stub) = панель не передана.
    """
    if not path or not os.path.exists(path) or os.path.getsize(path) == 0:
        return None
    taxids, names = set(), set()
    with open(path) as fh:
        for line in fh:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            (taxids if s.isdigit() else names).add(s)
    return (taxids, names) if (taxids or names) else None


def main():
    ap = argparse.ArgumentParser(description="Per-sample QC gate for NIPT metagenomics")
    ap.add_argument("--sample", required=True)
    ap.add_argument("--mode", required=True, choices=["fastq", "bam"])
    ap.add_argument("--fastp_json", default=None,
                    help="fastp JSON; в режиме bam — stub {'stub': true}")
    ap.add_argument("--counts_dehost", required=True,
                    help="counts JSON dehost-фракции (bam_counts.py)")
    ap.add_argument("--counts_host", required=True,
                    help="counts JSON host-фракции (bam_counts.py)")
    ap.add_argument("--counts_half", default=None,
                    help="counts JSON half-фракции (полукартированные пары, "
                         "mode=fastq); в режиме bam — stub/отсутствует -> 0")
    ap.add_argument("--profile_dir", required=True,
                    help="Директория с 4 парами контракта profile v2 "
                         "(<sample>.<fraction>.<role>.profile.tsv/.json)")
    ap.add_argument("--diff_dir", required=True,
                    help="Директория с 2 сводками comparison model "
                         "(<sample>.<role>.dehost_vs_host.compare.json)")
    ap.add_argument("--target_panel", default=None,
                    help="Файл панели мишеней (taxid/имя на строку) или пустой stub")
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

    criteria = []   # dict: id, value, threshold, pass, severity_on_fail, reason_code
    metrics = {}

    def add(cid, value, thr, passed, severity, reason=None, desc=""):
        criteria.append({"id": cid, "value": value, "threshold": thr,
                         "pass": passed, "severity_on_fail": severity,
                         "reason_code": reason if passed is False else None,
                         "description": desc})

    def skip(cid, desc):
        criteria.append({"id": cid, "value": None, "threshold": None,
                         "pass": None, "severity_on_fail": None,
                         "reason_code": None, "description": desc})

    # ---------- 0. Обязательные артефакты ----------
    # контракт profile v2: 4 пары profile.tsv + sidecar profile.json
    profiles = {}
    missing_profiles = []
    for frac in EXPECTED_FRACTIONS:
        for role in EXPECTED_ROLES:
            base = os.path.join(args.profile_dir,
                                f"{args.sample}.{frac}.{role}.profile")
            data, err = try_read_profile(base + ".tsv", base + ".json")
            if err:
                missing_profiles.append(f"{frac}.{role} ({err})")
            else:
                profiles[(frac, role)] = data
    if missing_profiles:
        add("artifacts_profiles", None, "4x (profile.tsv + profile.json)",
            False, "FAIL", "PIPELINE_ERROR",
            "Отсутствуют/битые контракты профиля: " + "; ".join(missing_profiles))

    # compare JSON x2 (сводки comparison model обязательной пары dehost:host)
    diffs = {}
    missing_diffs = []
    for db in EXPECTED_ROLES:
        p = os.path.join(
            args.diff_dir,
            f"{args.sample}.{db}.dehost_vs_host.compare.json")
        data, err = load_json(p)
        if err:
            missing_diffs.append(f"{db} ({err})")
        else:
            diffs[db] = data
    if missing_diffs:
        add("artifacts_diffs", None, "2 diff JSON", False, "FAIL",
            "PIPELINE_ERROR",
            "Отсутствуют/битые сводки comparison model: " + "; ".join(missing_diffs))

    # ---------- 1. Метрики входа (counts JSON, primary-only) ----------
    cnt_d = load_counts(args.counts_dehost)
    cnt_h = load_counts(args.counts_host)
    # half-фракция (полукартированные пары) — опциональна: в режиме bam
    # её нет (партиция NIPD из двух фракций) -> stub -> 0, поведение
    # как в v1.4.0. В режиме fastq DEHOST выдаёт counts всегда (пустая
    # фракция — валидный JSON с нулями).
    cnt_half = load_counts(args.counts_half) if args.counts_half else None
    half_reads = cnt_half["primary_reads"] if cnt_half else 0
    if cnt_d is None or cnt_h is None:
        add("artifacts_counts", None, "2 counts JSON", False, "FAIL",
            "PIPELINE_ERROR",
            "Отсутствует/не парсится counts JSON dehost или host")
        skip("input_pairs", "Не вычисляется: нет counts JSON")
        skip("host_fraction", "Не вычисляется: нет counts JSON")
        skip("nonhost_pairs", "Не вычисляется: нет counts JSON")
        input_pairs = host_frac = nonhost_pairs = None
    else:
        # primary_reads — без secondary/supplementary (в отличие от «in
        # total» flagstat); fragments — пары для PE, прочтения для SE.
        dehost_reads = cnt_d["primary_reads"]
        host_reads = cnt_h["primary_reads"]
        total_reads = dehost_reads + host_reads
        input_pairs = cnt_d["fragments"] + cnt_h["fragments"]
        host_frac = host_reads / total_reads if total_reads else 0.0
        nonhost_pairs = cnt_d["fragments"]
        metrics.update({"input_pairs": input_pairs, "dehost_reads": dehost_reads,
                        "host_reads": host_reads,
                        "host_fraction": round(host_frac, 6),
                        "nonhost_pairs": nonhost_pairs,
                        "single_end": bool(cnt_d["single_end"]
                                           and cnt_h["single_end"])})
        # Полукартированные пары — информационные метрики (на вердикт не
        # влияют; влияют только на критерий consistency в режиме fastq)
        partition_reads = total_reads + half_reads
        metrics["half_mapped_reads"] = half_reads
        metrics["half_mapped_frac"] = (round(half_reads / partition_reads, 6)
                                       if partition_reads else 0.0)

        add("input_pairs", input_pairs, f">= {args.min_input_pairs:.0f}",
            input_pairs >= args.min_input_pairs, "FAIL", "LOW_INPUT",
            "Достаточность входа: число пар прочтений")

        add("host_fraction", round(host_frac, 4), f">= {args.host_frac_min}",
            host_frac >= args.host_frac_min, "FALLBACK", "HOST_FRACTION_LOW",
            "Доля хост-прочтений в ожидаемом диапазоне (ниже — аномалия wet-lab)")

        add("nonhost_pairs", nonhost_pairs, f">= {args.min_nonhost_reads:.0f}",
            nonhost_pairs >= args.min_nonhost_reads, "FALLBACK", "ZERO_MICROBIOME",
            "Достаточность не-хост сигнала (edge case: нулевой микробиом)")

    # ---------- 2. fastp ----------
    fastp_data, fastp_err = (None, "не передан")
    if args.fastp_json:
        fastp_data, fastp_err = load_json(args.fastp_json)
    is_stub = isinstance(fastp_data, dict) and fastp_data.get("stub") is True

    if is_stub:
        # штатный режим bam: fastp выполнен upstream в NIPD
        if args.mode == "fastq":
            add("artifacts_fastp", None, "fastp JSON", False, "FAIL",
                "PIPELINE_ERROR",
                "Режим fastq, но вместо отчёта fastp получен stub")
        for cid in ("q30_rate", "trim_survival", "consistency"):
            skip(cid, "Не оценивается (mode=bam, fastp upstream)")
    elif fastp_err:
        add("artifacts_fastp", None, "fastp JSON", False, "FAIL",
            "PIPELINE_ERROR", f"fastp JSON недоступен: {fastp_err}")
        for cid in ("q30_rate", "trim_survival", "consistency"):
            skip(cid, "Не вычисляется: нет fastp JSON")
    else:
        try:
            before = fastp_data["summary"]["before_filtering"]["total_reads"]
            after = fastp_data["summary"]["after_filtering"]["total_reads"]
            q30 = fastp_data["summary"]["before_filtering"]["q30_rate"]
        except (KeyError, TypeError):
            add("artifacts_fastp", None, "fastp JSON", False, "FAIL",
                "PIPELINE_ERROR",
                "fastp JSON без обязательных ключей summary.*")
            for cid in ("q30_rate", "trim_survival", "consistency"):
                skip(cid, "Не вычисляется: битый fastp JSON")
        else:
            survival = after / before if before else 0.0
            metrics.update({"q30_rate": round(q30, 4),
                            "trim_survival": round(survival, 4),
                            "fastp_total_before": before})
            add("q30_rate", round(q30, 4), f">= {args.min_q30}",
                q30 >= args.min_q30, "FALLBACK", "LOW_QUALITY",
                "Качество сырых прочтений (Q30)")
            add("trim_survival", round(survival, 4), f">= {args.min_trim_survival}",
                survival >= args.min_trim_survival, "FALLBACK", "LOW_QUALITY",
                "Доля пар, переживших тримминг")
            if cnt_d is not None and cnt_h is not None:
                # v1.5.0: DEHOST — точная партиция trimmed-прочтений,
                # поэтому в сумму входит и half-фракция (критерий 3:
                # dehost + host + half == trimmed)
                total_reads = (cnt_d["primary_reads"] + cnt_h["primary_reads"]
                               + half_reads)
                dev = abs(total_reads - after) / after if after else 1.0
                metrics["consistency_dev"] = round(dev, 6)
                add("consistency", round(dev, 6), f"<= {args.max_consistency_dev}",
                    dev <= args.max_consistency_dev, "FALLBACK",
                    "SPLIT_INCONSISTENT",
                    "Консистентность разбиения dehost+host+half против trimmed")
            else:
                skip("consistency", "Не вычисляется: нет counts JSON")

    # ---------- 3. Метрики dehost-ветки из контракта профиля ----------
    prof_dehost_full = profiles.get(("dehost", "full"))
    if prof_dehost_full is None:
        skip("dehost_human_frac",
             "Не вычисляется: нет профиля dehost.full")
        skip("classified_frac",
             "Не вычисляется: нет профиля dehost.full")
    else:
        # profile v2: метрики вычисляются классификаторо-независимо.
        # human — строка taxid 9606 в TSV (0, если БД не содержит human
        # или таксон ниже порога отчётности адаптера); classified_frac —
        # из счётчиков sidecar.
        inp = prof_dehost_full["input_reads"]
        human_reads = sum(t["reads"] for t in prof_dehost_full["taxa"]
                          if t["taxid"] == HUMAN_TAXID)
        human = human_reads / inp if inp else 0.0
        classified = (prof_dehost_full["classified_reads"] / inp
                      if inp else 0.0)
        metrics.update({"dehost_classified_frac": round(classified, 6),
                        "dehost_human_frac": round(human, 6),
                        "dehost_human_reads": human_reads})
        # Утечка human оценивается только при достаточном не-хост сигнале:
        # на единицах прочтений доля human бессмысленна (см. ZERO_MICROBIOME)
        if nonhost_pairs is None:
            skip("dehost_human_frac", "Не вычисляется: нет counts JSON")
        elif nonhost_pairs < args.min_nonhost_reads:
            skip("dehost_human_frac",
                 f"Не оценивается: не-хост пар < {args.min_nonhost_reads:.0f} "
                 "(недостаточно для оценки утечки)")
        else:
            add("dehost_human_frac", round(human, 4),
                f"<= {args.max_dehost_human_frac}",
                human <= args.max_dehost_human_frac, "FALLBACK", "DEHOST_FAILURE",
                "Утечка человеческой ДНК в dehost-фракцию (сбой dehost)")
        add("classified_frac", round(classified, 6),
            f">= {args.min_classified_frac}",
            classified >= args.min_classified_frac, "FALLBACK",
            "NO_CLASSIFIED_SIGNAL",
            "Санити: классифицируемая доля dehost-прочтений > 0")

    # ---------- 4. Позитивность по панели мишеней ----------
    panel = read_panel(args.target_panel)
    diff_patho = diffs.get("patho")
    enriched_patho = (diff_patho or {}).get("enriched_taxa", [])
    panel_hits = []
    if panel is not None and diff_patho is not None:
        taxids, names = panel
        panel_hits = [t for t in enriched_patho
                      if str(t.get("taxid")) in taxids or t.get("name") in names]
        positivity_source = "target_panel"
        panel_info = {"file": os.path.basename(args.target_panel),
                      "n_targets": len(taxids) + len(names),
                      "hits": [t["name"] for t in panel_hits],
                      "positivity_source": positivity_source}
    elif panel is None:
        positivity_source = "not_evaluated"
        panel_info = {"file": None, "n_targets": 0, "hits": [],
                      "positivity_source": positivity_source}
    else:
        # панель задана, но patho diff отсутствует — FAIL уже зафиксирован
        positivity_source = "target_panel"
        panel_info = {"file": os.path.basename(args.target_panel),
                      "n_targets": len(panel[0]) + len(panel[1]), "hits": [],
                      "positivity_source": positivity_source}
    metrics["patho_enriched_taxa"] = (diff_patho or {}).get("n_enriched", None)
    metrics["patho_panel_hits"] = len(panel_hits)

    # ---------- 5. Вердикт ----------
    failed = [c for c in criteria if c["pass"] is False]
    hard_fail = [c for c in failed if c["severity_on_fail"] == "FAIL"]
    soft_fail = [c for c in failed if c["severity_on_fail"] == "FALLBACK"]
    if hard_fail:
        status = "FAIL"
    elif soft_fail:
        status = "FALLBACK"
    else:
        status = "PASS_POSITIVE" if panel_hits else "PASS_NEGATIVE"
    reasons = sorted({c["reason_code"] for c in failed if c["reason_code"]})

    # ---------- 6. Контракт QC JSON ----------
    annotator = None
    db_id = {}
    for (frac, role), p in profiles.items():
        if annotator is None and p.get("annotator"):
            # profile v2: annotator + annotator_version — раздельные поля
            annotator = " ".join(x for x in (p.get("annotator"),
                                             p.get("annotator_version"))
                                 if x)
        db_id.setdefault(role, p.get("db_id"))

    report = {
        "sample_id": args.sample,
        "pipeline": "nipt_mngs_qc",
        "pipeline_version": "1.10.0",
        "status": status,
        "fallback_reasons": reasons,
        "positivity_source": positivity_source,
        "annotator": annotator,
        "db_id": db_id,
        "target_panel": panel_info,
        "criteria": criteria,
        "metrics": metrics,
        # Информационная сводка comparison model v2; на вердикт НЕ влияет
        "comparison_info": {
            db: {"pair": d.get("pair"),
                 "n_shared": d.get("n_shared"),
                 "n_only_A": d.get("n_only_A"),
                 "n_only_B": d.get("n_only_B"),
                 "jaccard": d.get("jaccard"),
                 "bray_curtis_dissimilarity": d.get(
                     "bray_curtis_dissimilarity"),
                 "spearman_rho": d.get("spearman_rho"),
                 "n_taxa_considered": d.get("n_taxa_considered"),
                 "n_enriched": d.get("n_enriched"),
                 "n_background": d.get("n_background"),
                 "enriched_taxa": d.get("enriched_taxa", [])}
            for db, d in diffs.items()
        },
    }
    with open(args.out_json, "w") as out:
        json.dump(report, out, indent=2, ensure_ascii=False)

    # ---------- 7. MultiQC custom content (из того же контракта) ----------
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
