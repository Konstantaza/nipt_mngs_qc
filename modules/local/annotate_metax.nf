// =============================================================================
// ANNOTATE_METAX — адаптер классификатора MetaX (hzi-bifo/Metax) к
// интерфейсу ANNOTATE (NIPT-16)
//
// Реальный классификатор: `metax profile` (встроенный маппер + EM)
// по входному FASTA фракции; выход конвертируется в контракт profile v2
// скриптом bin/metax_to_profile.py. native-отчёт — сырой profile MetaX.
//
// CLI проверен эмпирически на MetaX 2.0.3 (conda-канал zldeng):
//   metax profile -d <db_dir> --dmp-dir <db_dir>/taxonomy \
//       -i <mate1> -I <mate2> -o metax_run -t <cpus> --no-progress
//   * вход FASTA/FASTQ (gzip тоже); PE задаётся парой -i/-I;
//   * role == 'patho' => добавляется --host 9606 (контекст NIPT —
//     хозяин человек), и MetaX пишет metax_run.pathogen.profile.txt
//     (16 колонок) вместо metax_run.profile.txt (13 колонок);
//   * на ПУСТОМ входе MetaX падает ("file is empty") — процесс
//     коротким замыканием создаёт пустой профиль, конвертер эмитит
//     status=empty.
//
// ТРЕБОВАНИЯ К БД (контракт каталога из колонки path dbsheet):
//   <path>/metax_db.mxi, <path>/metax_db.tsv.gz   (metax index)
//   <path>/taxonomy/{nodes,names,merged,delnodes}.dmp
// Индекс обязан быть ПОЛНЫМ (не fractional): дробный индекс
// масштабирует counts и ломает инварианты контракта (см.
// bin/metax_to_profile.py). План сборки БД и контейнера —
// docs/profile_contract.md.
//
// Окружение: metax есть только в личном conda-канале автора (zldeng),
// biocontainers-образа нет => container = общий образ пайплайна;
// для docker-профиля metax нужно добавить в environment.yml
// перед сборкой образа (точные строки — docs/profile_contract.md).
// =============================================================================
process ANNOTATE_METAX {
    tag { "${meta.id}/${meta.fraction}/${meta.role}" }
    label 'process_medium'

    // MetaX 2 — самодостаточный бинарник, maCMD больше не нужен.
    // Python нужен нашим конвертеру и подготовке SE, а не самому MetaX.
    conda 'zldeng::metax=2.0.3 conda-forge::python=3.14.7'
    container params.pipeline_image

    publishDir "${params.outdir}", mode: 'copy',
        saveAs: { fn -> "${meta.id}/profiles/${fn}" }

    input:
    // meta: [id, fraction, single_end, db_id, role, annotator]
    tuple val(meta), path(r1), path(r2), path(se), path(db)

    output:
    tuple val(meta), path("${meta.id}.${meta.fraction}.${meta.role}.profile.tsv"),
          path("${meta.id}.${meta.fraction}.${meta.role}.profile.json"),
          path("${meta.id}.${meta.fraction}.${meta.role}.metax.txt"),
          emit: annotated

    script:
    def prefix = "${meta.id}.${meta.fraction}.${meta.role}"
    def host_arg = meta.role == 'patho' ? '--host 9606' : ''
    def profile_out = meta.role == 'patho' ?
        'metax_run.pathogen.profile.txt' : 'metax_run.profile.txt'
    """
    metax_prepare.py --db ${db} --single_end ${meta.single_end} \
        --se ${se} --out_se metax_single_end.fasta

    # число входных прочтений (PE: записи r1 — пара = одно событие;
    # SE: записи se)
    if [ "${meta.single_end}" = "true" ]; then
        N_INPUT=\$(grep -c '^>' ${se} || true)
    else
        N_INPUT=\$(grep -c '^>' ${r1} || true)
    fi

    if [ "\$N_INPUT" -eq 0 ]; then
        # MetaX падает на пустом входе — эмитим пустой профиль напрямую
        : > ${profile_out}
    elif [ "${meta.single_end}" = "true" ]; then
        metax profile -d ${db} --dmp-dir ${db}/taxonomy \
            -i metax_single_end.fasta -o metax_run -t ${task.cpus} \
            --no-progress ${host_arg}
    else
        metax profile -d ${db} --dmp-dir ${db}/taxonomy \
            -i ${r1} -I ${r2} -o metax_run -t ${task.cpus} \
            --no-progress ${host_arg}
    fi

    METAX_VERSION=\$(metax --version | awk '{print \$2}')

    metax_to_profile.py \
        --profile ${profile_out} \
        --r1 ${r1} --r2 ${r2} --se ${se} \
        --single_end ${meta.single_end} \
        --sample ${meta.id} \
        --fraction ${meta.fraction} \
        --role ${meta.role} \
        --db_id ${meta.db_id} \
        --annotator ${meta.annotator} \
        --annotator_version "\$METAX_VERSION" \
        --out_tsv ${prefix}.profile.tsv \
        --out_json ${prefix}.profile.json

    # native-отчёт — сырой выход MetaX (MultiQC его не распознаёт)
    cp ${profile_out} ${prefix}.metax.txt
    """
}
