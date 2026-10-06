// =============================================================================
// MULTIQC — агрегация всех QC-артефактов в единый HTML-отчёт
//
// MultiQC из коробки понимает fastp, samtools flagstat, kraken2, bracken;
// вердикты QC gate подключаются как custom content (*_qc_mqc.tsv).
// =============================================================================
process MULTIQC {
    label 'process_low'

    publishDir "${params.outdir}/multiqc", mode: 'copy'

    input:
    path qc_files   // не читается напрямую: нужен только для порядка зависимостей
                    // (MULTIQC ждёт завершения всех QC_GATE / PROFILE)
    val outdir_abs  // params.outdir, заранее resolved в абсолютный путь
                    // (относительный путь в work-директории задачи не существует)

    output:
    path 'multiqc_report.html'
    path 'multiqc_report_data', emit: data   // MultiQC называет data-папку по имени отчёта

    script:
    // Сканируем опубликованные результаты, а не стейджинг: так MultiQC видит
    // fastp / flagstat / kraken / custom content всех образцов.
    // Собственную выходную подпапку multiqc/ игнорируем (важно при -resume).
    """
    multiqc ${outdir_abs} \
        --config ${projectDir}/multiqc_config.yaml \
        --filename multiqc_report.html \
        --ignore '*/multiqc/*' \
        --force
    """
}
