#!/usr/bin/env python3
"""
check_smoke.py — проверка результатов smoke-прогонов nipt_mngs_qc (NIPT-13).

Использование: check_smoke.py RUN RESULTS_DIR
  RUN         — bam | fastq | stub | mixed (какой прогон проверять);
  RESULTS_DIR — каталог --outdir прогона.

Проверки:
  * набор образцов в результатах совпадает с ожидаемым (никто не потерян
    и не лишний);
  * статус и fallback_reasons каждого образца — точное совпадение;
  * сводки сравнения dehost:host (compare JSON, роль full): у нормальных
    образцов kraken2-прогонов наборы таксонов фракций РАЗЛИЧАЮТСЯ
    (n_only_A + n_only_B > 0, Jaccard < 1), у TOY_LEAK — совпадают
    (Jaccard == 1.0: формальная детекция сбоя dehost);
  * обязательные артефакты на месте и непустые: профили dehost/host x
    роли БД (tsv+json), counts JSON, compare tsv/json, qc_mqc.tsv,
    multiqc_report.html; в fastq-прогоне — ещё fastp json/html и
    hisat2 summary.

Прогоны (NIPT-16): stub — оба аннотатора stub (fixture-driven);
mixed — role full kraken2_bracken + role patho stub (доказательство
сменяемости аннотаторов без правок main.nf/compare/qc). В mixed
позитивность гаснет (PASS_NEGATIVE): панельный хит идёт по роли patho
(qc_gate: панель ∩ enriched_taxa patho-диффа), а stub даёт 2097 вне
панели; при этом compare и детекция сбоя dehost (роль full, kraken2)
ведут себя в точности как в прогоне bam.

Любое несовпадение — понятное сообщение в stderr и exit 1 (критерий
приёмки 2). При успехе — одна сводная строка и exit 0.
"""
import json
import sys
from pathlib import Path

# Ожидаемые вердикты (baseline v1.5.0+, см. report_pipeline_theory.md 6.2)
EXPECT = {
    'bam': {
        'TOY_NORMAL': {'status': 'PASS_POSITIVE', 'reasons': []},
        'TOY_SORTED': {'status': 'PASS_POSITIVE', 'reasons': []},
        'TOY_SE':     {'status': 'PASS_POSITIVE', 'reasons': []},
        'TOY_LEAK':   {'status': 'FALLBACK',
                       'reasons': ['DEHOST_FAILURE', 'HOST_FRACTION_LOW']},
    },
    'fastq': {
        'TOY_NORMAL': {'status': 'PASS_POSITIVE', 'reasons': []},
        'TOY_ZERO':   {'status': 'FALLBACK', 'reasons': ['ZERO_MICROBIOME']},
    },
    'stub': {
        'TOY_NORMAL': {'status': 'PASS_NEGATIVE', 'reasons': []},
        'TOY_SORTED': {'status': 'PASS_NEGATIVE', 'reasons': []},
        'TOY_SE':     {'status': 'PASS_NEGATIVE', 'reasons': []},
        'TOY_LEAK':   {'status': 'FALLBACK', 'reasons': ['HOST_FRACTION_LOW']},
    },
    # mixed: role full = kraken2 (db_full), role patho = stub.
    # Позитивность определяется РОЛЬЮ patho (панель ∩ enriched_taxa
    # patho-диффа, см. qc_gate.py): stub в роли patho даёт 2097 вне
    # панели => PASS_NEGATIVE, хотя роль full (kraken2) ведёт себя
    # в точности как в прогоне bam (compare и детекция сбоя dehost
    # у TOY_LEAK — по роли full и совпадают с bam). Это и есть
    # доказательство NIPT-16: замена одной строки dbsheet меняет
    # только соответствующую ветку аннотации.
    'mixed': {
        'TOY_NORMAL': {'status': 'PASS_NEGATIVE', 'reasons': []},
        'TOY_SORTED': {'status': 'PASS_NEGATIVE', 'reasons': []},
        'TOY_SE':     {'status': 'PASS_NEGATIVE', 'reasons': []},
        'TOY_LEAK':   {'status': 'FALLBACK',
                       'reasons': ['DEHOST_FAILURE', 'HOST_FRACTION_LOW']},
    },
}

ROLES = ('full', 'patho')
FRACTIONS = {'bam': ('dehost', 'host'), 'stub': ('dehost', 'host'),
             'mixed': ('dehost', 'host'), 'fastq': ('dehost', 'host', 'half')}
# Сравнение наборов таксонов dehost:host проверяется только там, где
# роль full обслуживает kraken2: stub по построению даёт одинаковый
# профиль фракций, а в mixed роль full — kraken2 (как в bam).
COMPARE_DIFFERS = {'bam': {'TOY_NORMAL', 'TOY_SORTED', 'TOY_SE'},
                   'fastq': {'TOY_NORMAL', 'TOY_ZERO'},
                   'stub': set(),
                   'mixed': {'TOY_NORMAL', 'TOY_SORTED', 'TOY_SE'}}
COMPARE_SAME = {'bam': {'TOY_LEAK'}, 'fastq': set(), 'stub': set(),
                'mixed': {'TOY_LEAK'}}


def fail(errors, msg):
    errors.append(msg)


def check_file(errors, path):
    if not path.is_file():
        fail(errors, f"нет файла: {path}")
    elif path.stat().st_size == 0:
        fail(errors, f"пустой файл: {path}")


def main():
    if len(sys.argv) != 3 or sys.argv[1] not in EXPECT:
        sys.exit(f"usage: check_smoke.py {{{'|'.join(EXPECT)}}} RESULTS_DIR")
    run, outdir = sys.argv[1], Path(sys.argv[2])
    errors = []
    expect = EXPECT[run]

    # --- набор образцов: ни потерянных, ни лишних ---
    found = {p.name.split('_qc_report')[0]
             for p in outdir.glob('*/qc/*_qc_report.json')}
    if found != set(expect):
        fail(errors, f"набор образцов: ожидалось {sorted(expect)}, "
                     f"найдено {sorted(found)}")

    for sample, exp in expect.items():
        qc = outdir / sample / 'qc' / f'{sample}_qc_report.json'
        if not qc.is_file():
            continue  # уже зафиксировано выше
        rep = json.loads(qc.read_text())

        # --- вердикт и причины ---
        if rep.get('status') != exp['status']:
            fail(errors, f"{sample}: status {rep.get('status')!r} != "
                         f"ожидаемого {exp['status']!r}")
        if sorted(rep.get('fallback_reasons', [])) != sorted(exp['reasons']):
            fail(errors, f"{sample}: fallback_reasons "
                         f"{rep.get('fallback_reasons')} != ожидаемых "
                         f"{exp['reasons']}")
        if not rep.get('pipeline_version'):
            fail(errors, f"{sample}: пустой pipeline_version в QC JSON")

        # --- сводка сравнения dehost:host (роль full) ---
        cmp_json = (outdir / sample / 'profiles' /
                    f'{sample}.full.dehost_vs_host.compare.json')
        if sample in COMPARE_DIFFERS[run] | COMPARE_SAME[run]:
            if not cmp_json.is_file():
                fail(errors, f"{sample}: нет {cmp_json.name}")
            else:
                c = json.loads(cmp_json.read_text())
                jac, only = c.get('jaccard'), c.get('n_only_A', 0) + \
                    c.get('n_only_B', 0)
                if sample in COMPARE_DIFFERS[run]:
                    if not (jac is not None and jac < 1.0 and only > 0):
                        fail(errors, f"{sample}: наборы таксонов dehost/host "
                                     f"должны различаться (Jaccard={jac}, "
                                     f"only_A+only_B={only})")
                elif not (jac == 1.0 and only == 0):
                    fail(errors, f"{sample} (сбой dehost): наборы таксонов "
                                 f"должны совпадать (Jaccard={jac}, "
                                 f"only_A+only_B={only})")

        # --- обязательные артефакты ---
        for fr in FRACTIONS[run]:
            if fr != 'half':  # half не аннотируется — профилей нет
                for role in ROLES:
                    for ext in ('tsv', 'json'):
                        check_file(errors, outdir / sample / 'profiles' /
                                   f'{sample}.{fr}.{role}.profile.{ext}')
            check_file(errors, outdir / sample / 'qc' / 'counts' /
                       f'{sample}.{fr}.counts.json')
        for role in ROLES:
            for ext in ('tsv', 'json'):
                check_file(errors, outdir / sample / 'profiles' /
                           f'{sample}.{role}.dehost_vs_host.compare.{ext}')
        check_file(errors, outdir / sample / 'qc' / f'{sample}_qc_mqc.tsv')
        if run == 'fastq':
            check_file(errors, outdir / sample / 'qc' / 'fastp' /
                       f'{sample}.fastp.json')
            check_file(errors, outdir / sample / 'qc' / 'dehost' /
                       f'{sample}.hisat2.summary.txt')

    check_file(errors, outdir / 'multiqc' / 'multiqc_report.html')

    if errors:
        for e in errors:
            print(f"FAIL: {e}", file=sys.stderr)
        sys.exit(f"check_smoke [{run}]: {len(errors)} несовпадений")
    print(f"check_smoke [{run}]: OK — {len(expect)} образцов, вердикты и "
          f"артефакты соответствуют ожиданиям")


if __name__ == '__main__':
    main()
