// =============================================================================
// FASTP — тримминг адаптеров и фильтрация по качеству (только mode=fastq)
//
// Теория: сырые прочтения DNBSeq содержат остатки адаптеров и низкокачественные
// концы. Адаптерные хвосты могут ложно картироваться на консервативные участки
// микробных геномов, поэтому тримминг — обязательная первая стадия.
// fastp отдаёт JSON/HTML-отчёты, которые читает и QC gate, и MultiQC.
// =============================================================================
process FASTP {
    tag { meta.id }
    label 'process_low'

    // NIPT-7 (v1.8.0): точные пины окружения; container — для -profile docker
    conda 'bioconda::fastp=1.3.7'
    container 'quay.io/biocontainers/fastp:1.3.7'

    // Публикуем только компактные QC-отчёты; триммированные FASTQ — промежуточные.
    // saveAs вычисляется в контексте задачи — так делают per-sample подпапки.
    publishDir "${params.outdir}", mode: 'copy', pattern: '*.{json,html}',
        saveAs: { fn -> "${meta.id}/qc/fastp/${fn}" }
    publishDir "${params.outdir}", mode: 'copy', pattern: 'trim_R*.fastq.gz',
        enabled: params.publish_validation_artifacts,
        saveAs: { fn -> "${meta.id}/validation/trimmed/${fn}" }

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path('trim_R*.fastq.gz'), emit: reads
    tuple val(meta.id), path("${meta.id}.fastp.json"), emit: json
    path "${meta.id}.fastp.html", emit: html

    script:
    // Имя файла содержит sample_id: иначе MultiQC схлопывает одноимённые
    // fastp.json разных образцов в один сэмпл.
    def inputs = reads instanceof List ? reads : [reads]
    def mates = meta.single_end ? '' : "--in2 ${inputs[1]} --out2 trim_R2.fastq.gz"
    """
    fastp \
        --in1 ${inputs[0]} --out1 trim_R1.fastq.gz ${mates} \
        --length_required ${params.fastp_min_len} \
        --thread ${task.cpus} \
        --json ${meta.id}.fastp.json --html ${meta.id}.fastp.html
    """
}
