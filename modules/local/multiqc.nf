// =============================================================================
// MULTIQC — агрегация QC-артефактов в единый HTML-отчёт
//
// v1.1: работает ТОЛЬКО по застейдженным входам (custom content из контракта
// QC gate, kraken reports, flagstat, fastp) — params.outdir не сканируется.
// MultiQC из коробки понимает fastp, samtools flagstat, kraken2; вердикты
// QC gate подключаются как custom content (*_qc_mqc.tsv).
// =============================================================================
process MULTIQC {
    label 'process_low'

    // NIPT-7 (v1.8.0): точный пин; container — для -profile docker
    conda 'bioconda::multiqc=1.35'
    container 'quay.io/biocontainers/multiqc:1.35--pyhdfd78af_2'

    publishDir "${params.outdir}/multiqc", mode: 'copy'

    input:
    path qc_files   // все QC-артефакты всех образцов (collect):
                    // имена уникальны за счёт префикса sample_id

    output:
    path 'multiqc_report.html'
    path 'multiqc_report_data', emit: data   // MultiQC называет data-папку по имени отчёта

    script:
    // Сканируем только свою work-директорию с застейдженными входами.
    """
    multiqc . \
        --config ${projectDir}/multiqc_config.yaml \
        --filename multiqc_report.html \
        --force
    """
}
