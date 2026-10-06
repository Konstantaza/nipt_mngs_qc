// =============================================================================
// BAM2FASTA — конвертация DEHOST/HOST BAM в FASTA + сбор статистики (flagstat)
//
// Теория: Kraken2 классифицирует последовательности, а не выравнивания,
// поэтому BAM превращаем в парные FASTA (samtools fasta -1/-2).
// Заодно считаем samtools flagstat по обоим BAM — это базовые числа QC gate:
// размер фракций, доля хоста, консистентность разбиения.
// =============================================================================
process BAM2FASTA {
    tag { sample_id }
    label 'process_low'

    publishDir "${params.outdir}", mode: 'copy', pattern: '*.flagstat.txt',
        saveAs: { fn -> "${sample_id}/qc/flagstat/${fn}" }

    input:
    // stageAs с фиксированными именами: защита от коллизии, когда dehost_bam
    // и host_bam имеют одинаковое имя файла (или это один и тот же файл).
    tuple val(sample_id), path(dehost_bam, stageAs: 'dehost_input.bam'),
          path(host_bam, stageAs: 'host_input.bam')

    output:
    tuple val(sample_id), path('dehost_1.fasta'), path('dehost_2.fasta'),
          path('host_1.fasta'), path('host_2.fasta'), emit: seqs
    tuple val(sample_id), path("${sample_id}.dehost.flagstat.txt"),
          path("${sample_id}.host.flagstat.txt"), emit: stats

    script:
    // Имена flagstat содержат sample_id: защита от коллизий при стейджинге
    // и от слияния сэмплов в MultiQC.
    """
    samtools fasta -@ ${task.cpus} -1 dehost_1.fasta -2 dehost_2.fasta -s /dev/null ${dehost_bam}
    samtools fasta -@ ${task.cpus} -1 host_1.fasta   -2 host_2.fasta   -s /dev/null ${host_bam}
    samtools flagstat -@ ${task.cpus} ${dehost_bam} > ${sample_id}.dehost.flagstat.txt
    samtools flagstat -@ ${task.cpus} ${host_bam}   > ${sample_id}.host.flagstat.txt
    """
}
