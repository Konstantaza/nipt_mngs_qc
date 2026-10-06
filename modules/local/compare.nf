// =============================================================================
// COMPARE — comparison model: дифференциальный профиль dehost vs host
//
// Запускается 2 раза на образец: для full-БД и для patho-БД (EukPathDB).
// Логика эвристики обогащения — в bin/compare_profiles.py.
// =============================================================================
process COMPARE {
    tag { "${sample_id}/${dbname}" }
    label 'process_low'

    publishDir "${params.outdir}", mode: 'copy',
        saveAs: { fn -> "${sample_id}/profiles/${fn}" }

    input:
    tuple val(sample_id), val(dbname), path(bracken_dehost), path(bracken_host)

    output:
    tuple val(sample_id), val(dbname),
          path("${sample_id}.${dbname}.diff.tsv"),
          path("${sample_id}.${dbname}.diff.json"), emit: diff

    script:
    def prefix = "${sample_id}.${dbname}"
    """
    compare_profiles.py \
        --dehost ${bracken_dehost} \
        --host ${bracken_host} \
        --sample ${sample_id} \
        --db ${dbname} \
        --min_reads ${params.bracken_min_reads} \
        --enrichment_fold ${params.enrichment_fold} \
        --out_tsv ${prefix}.diff.tsv \
        --out_json ${prefix}.diff.json
    """
}
