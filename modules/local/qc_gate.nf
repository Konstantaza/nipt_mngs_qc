// =============================================================================
// QC_GATE — формальные критерии QC/адекватности на уровне одного образца
//
// v1.4: работает только по контрактам — profile.tsv + sidecar JSON
// (стадия ANNOTATE, profile v2), diff JSON (COMPARE), counts JSON
// (BAM2FASTA/bam_counts.py, primary-only), fastp JSON. Форматы
// классификаторов здесь не парсятся. Отсутствие обязательного артефакта —
// FAIL/PIPELINE_ERROR, а не молчаливый PASS_NEGATIVE.
// v1.5: + counts_half (полукартированные пары DEHOST; в bam-режиме —
// stub -> 0); consistency = dehost+host+half против trimmed.
//
// Главный артефакт — <sample>_qc_report.json для внешней системы отчётов
// + custom content для MultiQC.
// =============================================================================
process QC_GATE {
    tag { sample_id }
    label 'process_low'

    // NIPT-7 (v1.8.0): только python3 (qc_gate.py, stdlib) + coreutils (mkdir/cp)
    conda 'conda-forge::python=3.14.7'
    container 'quay.io/biocontainers/python:3.14'

    publishDir "${params.outdir}", mode: 'copy',
        saveAs: { fn -> "${sample_id}/qc/${fn}" }

    input:
    tuple val(sample_id), path(counts_dehost), path(counts_host),
          path(counts_half), path(fastp_json), path(profile_tsvs),
          path(profile_jsons), path(diff_jsons)
    path target_panel   // файл панели мишеней или пустой stub assets/no_panel.txt

    output:
    tuple val(sample_id), path("${sample_id}_qc_report.json"), emit: report
    path "${sample_id}_qc_mqc.tsv", emit: mqc

    script:
    """
    mkdir -p profiles_in diff_in
    # Списки могут быть пустыми (образец без артефактов после сбоя upstream):
    # это штатная ситуация — qc_gate.py вынесет FAIL/PIPELINE_ERROR.
    if [ -n "${profile_tsvs}" ]; then cp ${profile_tsvs} profiles_in/; fi
    if [ -n "${profile_jsons}" ]; then cp ${profile_jsons} profiles_in/; fi
    if [ -n "${diff_jsons}" ]; then cp ${diff_jsons} diff_in/; fi

    qc_gate.py \
        --sample ${sample_id} \
        --mode ${params.mode} \
        --fastp_json ${fastp_json} \
        --counts_dehost ${counts_dehost} \
        --counts_host ${counts_host} \
        --counts_half ${counts_half} \
        --profile_dir profiles_in \
        --diff_dir diff_in \
        --target_panel ${target_panel} \
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
