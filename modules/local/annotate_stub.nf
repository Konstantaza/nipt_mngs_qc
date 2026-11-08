// =============================================================================
// ANNOTATE_STUB — fixture-driven stub-адаптер интерфейса ANNOTATE (NIPT-16)
//
// Тестовый аннотатор без внешней БД и без классификатора: профиль
// синтезируется из fixture-файла <path_из_dbsheet>/<db_id>.tsv
// (taxid, name, rank, weight; см. bin/stub_to_profile.py). Назначение —
// доказательство расширяемости реестра: добавление аннотатора = этот
// процесс + конвертер + строка в dbsheet (+ одна ветка в реестре
// subworkflows/local/annotate.nf). main.nf и потребители (compare,
// qc_gate) не меняются.
//
// Выход — контракт profile v2; native — текстовый отчёт stub'а
// (MultiQC его не распознаёт и игнорирует).
// =============================================================================
process ANNOTATE_STUB {
    tag { "${meta.id}/${meta.fraction}/${meta.role}" }
    label 'process_low'

    // только python3 (stub_to_profile.py, stdlib)
    conda 'conda-forge::python=3.14.7'
    container 'quay.io/biocontainers/python:3.14'

    publishDir "${params.outdir}", mode: 'copy',
        saveAs: { fn -> "${meta.id}/profiles/${fn}" }

    input:
    // meta: [id, fraction, single_end, db_id, role, annotator]
    // db — каталог из колонки path dbsheet; процесс читает из него
    // fixture <db>/<db_id>.tsv.
    tuple val(meta), path(r1), path(r2), path(se), path(db)

    output:
    tuple val(meta), path("${meta.id}.${meta.fraction}.${meta.role}.profile.tsv"),
          path("${meta.id}.${meta.fraction}.${meta.role}.profile.json"),
          path("${meta.id}.${meta.fraction}.${meta.role}.stub.txt"),
          emit: annotated

    script:
    def prefix = "${meta.id}.${meta.fraction}.${meta.role}"
    """
    stub_to_profile.py \
        --fixture ${db}/${meta.db_id}.tsv \
        --r1 ${r1} --r2 ${r2} --se ${se} \
        --single_end ${meta.single_end} \
        --sample ${meta.id} \
        --fraction ${meta.fraction} \
        --role ${meta.role} \
        --db_id ${meta.db_id} \
        --annotator ${meta.annotator} \
        --annotator_version "1.0-stub" \
        --out_tsv ${prefix}.profile.tsv \
        --out_json ${prefix}.profile.json \
        --out_native ${prefix}.stub.txt
    """
}
