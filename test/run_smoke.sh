#!/usr/bin/env bash
# =============================================================================
# run_smoke.sh — сквозной smoke-harness nipt_mngs_qc (NIPT-13)
#
# Одна команда проверяет пайплайн целиком:
#   bash test/run_smoke.sh
#
# Что делает:
#   1. pytest — все юнит-тесты bin/ (критерий приёмки 3);
#   2. nextflow — четыре прогона на toy-данных (входы заданы внутри
#      самодостаточного профиля test, v1.7.0):
#        bam   — mode=bam,   dbsheet kraken2_bracken (4 образца);
#        fastq — mode=fastq, dbsheet kraken2_bracken (2 образца);
#        stub  — mode=bam,   dbsheet со stub-аннотатором (4 образца);
#        mixed — mode=bam,   dbsheet kraken2_bracken (role full) +
#                stub (role patho) — доказательство сменяемости
#                аннотаторов без правок main.nf/compare/qc (NIPT-16);
#   3. test/check_smoke.py — сверка вердиктов QC JSON, сводок сравнения
#      dehost:host и обязательных артефактов с ожиданиями; любое
#      несовпадение — exit 1 (критерий приёмки 2).
#
# Окружение (критерий приёмки 1 — чистая машина с docker):
#   * есть docker  — pytest и три мульти-инструментальных процесса
#     (DEHOST, BAM2FASTA, ANNOTATE_KRAKEN2_BRACKEN) идут в контейнере
#     nipt_mngs_qc:1.10.0 (при отсутствии образ собирается из Dockerfile;
#     разовая операция, дальше — кэш), а одно-инструментные процессы
#     тянут точные biocontainers-образы с quay.io (NIPT-7, per-process
#     контейнеры); nextflow нужен на хосте (при отсутствии — бутстрап
#     ./nextflow через get.nextflow.io, нужна java);
#   * нет docker   — fallback: инструменты с PATH (conda-окружение из
#     environment.yml).
#
# Аннотатор metax в harness НЕ входит: его conda-канал (zldeng) не
# входит в стандартную матрицу окружений. Ручной прогон на toy-БД:
#   nextflow run main.nf -profile test --mode bam \
#     --dbsheet test/toy/dbsheet_metax_real.csv --outdir results_metax
# (требует zldeng::metax=2.0.3 и Python
# на PATH; см. docs/profile_contract.md).
#
# Результаты: test/results_smoke/{bam,fastq,stub,mixed}/ + логи рядом.
# Ожидаемое время: ~10 минут на ноутбуке (без учёта разовой сборки образа).
# =============================================================================
set -euo pipefail

cd "$(dirname "$0")/.."
REPO="$PWD"
IMAGE="nipt_mngs_qc:1.10.0"
SMOKE_DIR="test/results_smoke"

log() { echo "[run_smoke] $*"; }

# --- 1. runtime: docker или локальные инструменты -------------------------
if command -v docker >/dev/null 2>&1 && docker info >/dev/null 2>&1; then
    RUNTIME=docker
    if ! docker image inspect "$IMAGE" >/dev/null 2>&1; then
        log "образ $IMAGE не найден — собираю из Dockerfile (разово)"
        docker build -t "$IMAGE" .
    fi
    PROFILES="test,docker"
    run_py() {
        docker run --rm -u "$(id -u):$(id -g)" \
            -v "$REPO:/repo" -w /repo "$IMAGE" python3 "$@"
    }
else
    RUNTIME=local
    PROFILES="test"
    for t in fastp hisat2 samtools kraken2 bracken multiqc python3; do
        command -v "$t" >/dev/null 2>&1 || {
            echo "ERROR: '$t' не найден на PATH. Активируйте conda-окружение" \
                 "из environment.yml или запустите на машине с docker." >&2
            exit 1
        }
    done
    run_py() { python3 "$@"; }
fi
log "runtime: $RUNTIME (профили: $PROFILES)"

# --- 2. nextflow ------------------------------------------------------------
if command -v nextflow >/dev/null 2>&1; then
    NF=nextflow
elif [ -x "$REPO/nextflow" ]; then
    NF="$REPO/nextflow"
else
    command -v java >/dev/null 2>&1 || {
        echo "ERROR: nextflow не найден; для автобутстрапа нужна java." \
             "Установите nextflow: https://nextflow.io" >&2
        exit 1
    }
    log "nextflow не найден — бутстрап в $REPO/nextflow"
    curl -fsSL https://get.nextflow.io | bash
    NF="$REPO/nextflow"
fi
log "nextflow: $("$NF" -version 2>&1 | grep -m1 version || true)"

# --- 3. чистый старт --------------------------------------------------------
rm -rf "$REPO/$SMOKE_DIR"
mkdir -p "$REPO/$SMOKE_DIR"

# --- 4. юнит-тесты ----------------------------------------------------------
log "pytest (юнит-тесты bin/)"
run_py -m pytest test/ -q

# --- 5. сквозные прогоны ----------------------------------------------------
run_nf() {  # run_nf <имя> <extra args nextflow...>
    local name="$1"; shift
    log "nextflow run: $name"
    if ! "$NF" -log "$REPO/$SMOKE_DIR/$name.nextflow.log" run "$REPO/main.nf" \
            -profile "$PROFILES" \
            -work-dir "$REPO/$SMOKE_DIR/work" \
            --outdir "$REPO/$SMOKE_DIR/$name" "$@"; then
        echo "ERROR: прогон '$name' упал. Хвост лога:" >&2
        tail -n 30 "$REPO/$SMOKE_DIR/$name.nextflow.log" >&2
        exit 1
    fi
}

run_nf bam   --mode bam
run_nf fastq --mode fastq
run_nf stub  --mode bam --dbsheet "$REPO/test/toy/dbsheet_stub.csv"
run_nf mixed --mode bam --dbsheet "$REPO/test/toy/dbsheet_mixed.csv"

# --- 6. сверка с ожиданиями -------------------------------------------------
for name in bam fastq stub mixed; do
    run_py test/check_smoke.py "$name" "$SMOKE_DIR/$name"
done

log "ВСЕ ПРОВЕРКИ ПРОШЛИ (runtime: $RUNTIME, время: ${SECONDS} c)"
