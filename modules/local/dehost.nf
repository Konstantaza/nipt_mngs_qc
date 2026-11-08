// =============================================================================
// DEHOST — отделение прочтений человека (mode=fastq), v1.5.0
//
// Теория: 95–99% прочтений НИПТ — ДНК человека (материнская + фетальная
// вкДНК). Один проход HISAT2 на геном человека с --no-spliced-alignment
// (вход — ДНК, сплайсированные выравнивания неуместны; hocort гонял
// HISAT2 со сплайсингом по умолчанию и делал это дважды). Разбиение по
// SAM-флагам — bin/dehost_split.sh: host (оба мейта картированы),
// dehost (оба некартированы), half (полукартированные — ровно один
// мейт; раньше молча терялись между двумя прогонами hocort). Три
// фракции — точная партиция trimmed-прочтений (критерий приёмки 3).
//
// --summary-file: статистика выравнивания HISAT2 публикуется и идёт в
// MultiQC (у hocort лог выравнивателя был скрыт в work/).
//
// Выходные BAM — промежуточные (не публикуются, storage hygiene);
// дальше каждая фракция идёт в per-fraction BAM2FASTA.
// =============================================================================
process DEHOST {
    tag { meta.id }
    label 'process_medium'

    // NIPT-7 (v1.8.0): hisat2+samtools в одном окружении. Готового mulled-образа
    // с этой точной комбинацией версий на quay.io нет (проверено по хешам
    // mulled-v2) — используем all-in-one образ пайплайна (params.pipeline_image,
    // собирается из Dockerfile в корне репозитория).
    conda 'bioconda::hisat2=2.2.3 bioconda::samtools=1.24'
    container params.pipeline_image

    publishDir "${params.outdir}", mode: 'copy', pattern: '*.hisat2.summary.txt',
        saveAs: { fn -> "${meta.id}/qc/dehost/${fn}" }
    publishDir "${params.outdir}", mode: 'copy', pattern: '{host,dehost,half}.bam',
        enabled: params.publish_validation_artifacts,
        saveAs: { fn -> "${meta.id}/validation/partitions/${fn}" }

    input:
    tuple val(meta), path(reads)
    path index_files   // все файлы HISAT2-индекса человека (стейджатся в cwd)

    output:
    tuple val(meta.id), path('dehost.bam'), path('host.bam'),
          path('half.bam'), emit: bams
    tuple val(meta.id), path("${meta.id}.hisat2.summary.txt"), emit: summary

    script:
    def idx = file(params.host_index).name   // префикс индекса (без пути)
    def inputs = reads instanceof List ? reads : [reads]
    def read_args = meta.single_end ? "-U ${inputs[0]}" : "-1 ${inputs[0]} -2 ${inputs[1]}"
    """
    hisat2 -x ${idx} ${read_args} \\
        --no-spliced-alignment \\
        --summary-file ${meta.id}.hisat2.summary.txt \\
        -p ${task.cpus} \\
        | samtools view -@ ${task.cpus} -b -o aligned.bam -

    dehost_split.sh aligned.bam split ${task.cpus}
    mv split.dehost.bam dehost.bam
    mv split.host.bam   host.bam
    mv split.half.bam   half.bam
    """
}
