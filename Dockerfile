# =============================================================================
# nipt_mngs_qc — all-in-one образ инструментов (v1.10.0)
#
# Собирается из environment.yml (те же точные пины, что и у разработческого
# conda-окружения). Репозиторий в образ НЕ входит: Nextflow сам монтирует
# main.nf / modules / bin в контейнеры задач, а входные файлы стейджатся
# каналами. NIPT-7: образ используется тремя мульти-инструментальными
# процессами (DEHOST, BAM2FASTA, ANNOTATE_KRAKEN2_BRACKEN) через
# params.pipeline_image — готовых mulled-образов с точными комбинациями
# версий на quay.io нет; одно-инструментные процессы тянут точные
# biocontainers-образы напрямую. Сборка и запуск (это делает harness
# test/run_smoke.sh на чистой машине):
#
#   docker build -t nipt_mngs_qc:1.10.0 .
#   nextflow run main.nf -profile test,docker
#
# NIPT-16: опциональный ANNOTATE_METAX тоже использует pipeline_image,
# но стандартный образ не содержит MetaX. План отдельного образа
# 1.9.0-metax с zldeng::metax=2.0.3 — docs/profile_contract.md §8.4.
# Внешний MA для MetaX 2 не нужен.
#
# MPLCONFIGDIR / XDG_CACHE_HOME + chmod 777: с docker.userEmulation задачи
# исполняются под UID пользователя хоста, HOME внутри контейнера может быть
# недоступен на запись — matplotlib/multiqc складывают кэш в /tmp.
# =============================================================================
FROM mambaorg/micromamba:1.5.10

COPY environment.yml /tmp/environment.yml
RUN micromamba create -y -n nipt_mngs_qc -f /tmp/environment.yml \
    && micromamba clean -a -y

ENV PATH=/opt/conda/envs/nipt_mngs_qc/bin:$PATH \
    MPLCONFIGDIR=/tmp/.mplconfig \
    XDG_CACHE_HOME=/tmp/.cache
RUN mkdir -p /tmp/.mplconfig /tmp/.cache \
    && chmod 777 /tmp/.mplconfig /tmp/.cache

# Сброс entrypoint базового образа micromamba: Nextflow вызывает
# `bash -c <task>` напрямую, обёртка _entrypoint.sh не нужна.
ENTRYPOINT []
CMD ["/bin/bash"]
