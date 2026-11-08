// =============================================================================
// COMPARE — comparison model v2: сравнение профилей фракций на контракте
// profile v2 (profile.tsv + sidecar JSON), НЕ на сырых выходах
// классификатора.
//
// Запускается на каждое (образец, роль БД, пара из params.compare_pairs).
// По умолчанию пара одна: dehost:host. Логика — в bin/compare_profiles.py.
// =============================================================================
process COMPARE {
    tag { "${sample_id}/${role}/${pair}" }
    label 'process_low'

    // NIPT-7 (v1.8.0): только python3 (compare_profiles.py, stdlib)
    conda 'conda-forge::python=3.14.7'
    container 'quay.io/biocontainers/python:3.14'

    publishDir "${params.outdir}", mode: 'copy',
        saveAs: { fn -> "${sample_id}/profiles/${fn}" }

    input:
    tuple val(sample_id), val(role),
          path(dehost_tsv), path(dehost_json),
          path(host_tsv), path(host_json), val(pair)

    output:
    tuple val(sample_id), val(role), val(pair),
          path("${sample_id}.${role}.${pair.replace(':', '_vs_')}.compare.tsv"),
          path("${sample_id}.${role}.${pair.replace(':', '_vs_')}.compare.json"),
          emit: diff

    script:
    def slug = pair.replace(':', '_vs_')
    def prefix = "${sample_id}.${role}.${slug}"
    """
    compare_profiles.py \
        --dehost_tsv ${dehost_tsv} \
        --dehost_json ${dehost_json} \
        --host_tsv ${host_tsv} \
        --host_json ${host_json} \
        --sample ${sample_id} \
        --role ${role} \
        --pair ${pair} \
        --min_reads ${params.compare_min_reads} \
        --host_min_reads ${params.compare_host_min_reads} \
        --enrichment_fold ${params.enrichment_fold} \
        --out_tsv ${prefix}.compare.tsv \
        --out_json ${prefix}.compare.json
    """
}
