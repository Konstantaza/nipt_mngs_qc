// =============================================================================
// VALIDATE_PROFILE — framework-валидация контракта profile v2
//
// Прогоняет bin/validate_profile.py на каждый выход адаптера ANNOTATE:
// нарушение контракта = падение задачи, а не молчаливый битый профиль
// downstream. Адаптерам не нужно помнить про валидацию — её делает
// интерфейс. Файлы уже опубликованы адаптером, здесь — pass-through.
// =============================================================================
process VALIDATE_PROFILE {
    tag { "${meta.id}/${meta.fraction}/${meta.role}" }
    label 'process_low'

    // NIPT-7 (v1.8.0): только python3 (validate_profile.py, stdlib)
    conda 'conda-forge::python=3.14.7'
    container 'quay.io/biocontainers/python:3.14'

    input:
    tuple val(meta), path(tsv), path(json), path(nat)

    output:
    tuple val(meta), path(tsv), path(json), path(nat), emit: profiles

    script:
    """
    validate_profile.py \
        --tsv ${tsv} \
        --json ${json} \
        --expect-sample ${meta.id} \
        --expect-fraction ${meta.fraction} \
        --expect-role ${meta.role}
    """
}
