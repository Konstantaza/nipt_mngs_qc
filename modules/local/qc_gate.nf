// =============================================================================
// QC_GATE — формальные критерии QC/адекватности на уровне одного образца
//
// Собирает метрики всех стадий (fastp, flagstat, kraken2, comparison model),
// проверяет чеклист и выносит вердикт PASS_POSITIVE / PASS_NEGATIVE /
// FALLBACK / FAIL. Главный артефакт — <sample>_qc_report.json для внешней
// системы отчётов + custom content для MultiQC.
// =============================================================================
process QC_GATE {
    tag { sample_id }
    label 'process_low'

    publishDir "${params.outdir}", mode: 'copy',
        saveAs: { fn -> "${sample_id}/qc/${fn}" }

    input:
    tuple val(sample_id), path(flagstat_dehost), path(flagstat_host),
          path(fastp_json), path(kraken_reports), path(diff_jsons)

    output:
    tuple val(sample_id), path("${sample_id}_qc_report.json"), emit: report
    path "${sample_id}_qc_mqc.tsv", emit: mqc

    script:
    """
    mkdir -p kraken_in diff_in
    # Имена файлов детерминированы, порядок в канале не гарантирован —
    # раскладываем по директориям и ищем glob-ом внутри qc_gate.py
    cp ${kraken_reports} kraken_in/
    cp ${diff_jsons} diff_in/

    qc_gate.py \
        --sample ${sample_id} \
        --fastp_json ${fastp_json} \
        --flagstat_dehost ${flagstat_dehost} \
        --flagstat_host ${flagstat_host} \
        --kraken_dir kraken_in \
        --diff_dir diff_in \
        --min_input_pairs ${params.min_input_pairs} \
        --min_q30 ${params.min_q30} \
        --min_trim_survival ${params.min_trim_survival} \
        --host_frac_min ${params.host_frac_min} \
        --min_nonhost_reads ${params.min_nonhost_reads} \
        --max_dehost_human_frac ${params.max_dehost_human_frac} \
        --min_classified_frac ${params.min_classified_frac} \
        --max_consistency_dev ${params.max_consistency_dev} \
        --out_json ${sample_id}_qc_report.json \
        --out_mqc ${sample_id}_qc_mqc.tsv
    """
}
