// =============================================================================
// FASTP — тримминг адаптеров и фильтрация по качеству (только mode=fastq)
//
// Теория: сырые прочтения DNBSeq содержат остатки адаптеров и низкокачественные
// концы. Адаптерные хвосты могут ложно картироваться на консервативные участки
// микробных геномов, поэтому тримминг — обязательная первая стадия.
// fastp отдаёт JSON/HTML-отчёты, которые читает и QC gate, и MultiQC.
// =============================================================================
process FASTP {
    tag { sample_id }
    label 'process_low'

    // Публикуем только компактные QC-отчёты; триммированные FASTQ — промежуточные.
    // saveAs вычисляется в контексте задачи — так делают per-sample подпапки.
    publishDir "${params.outdir}", mode: 'copy', pattern: '*.{json,html}',
        saveAs: { fn -> "${sample_id}/qc/fastp/${fn}" }

    input:
    tuple val(sample_id), path(r1), path(r2)

    output:
    tuple val(sample_id), path('trim_R1.fastq.gz'), path('trim_R2.fastq.gz'), emit: reads
    tuple val(sample_id), path("${sample_id}.fastp.json"), emit: json
    path "${sample_id}.fastp.html", emit: html

    script:
    // Имя файла содержит sample_id: иначе MultiQC схлопывает одноимённые
    // fastp.json разных образцов в один сэмпл.
    """
    fastp \
        --in1 ${r1} --in2 ${r2} \
        --out1 trim_R1.fastq.gz --out2 trim_R2.fastq.gz \
        --length_required ${params.fastp_min_len} \
        --thread ${task.cpus} \
        --json ${sample_id}.fastp.json --html ${sample_id}.fastp.html
    """
}
