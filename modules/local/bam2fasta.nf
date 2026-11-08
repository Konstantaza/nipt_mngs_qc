// =============================================================================
// BAM2FASTA — конвертация BAM одной фракции в FASTA + счётчики (v1.5.0)
//
// Теория: Kraken2 классифицирует последовательности, а не выравнивания,
// поэтому BAM превращаем в FASTA (samtools fasta). Процесс per-fraction:
// вход (meta[id, fraction], bam) — одна фракция одного образца; общее
// ядро для обеих точек входа (mode=bam: dehost/host от NIPD; mode=fastq:
// dehost/host/half от DEHOST).
//
// Логика скрипта без изменений от v1.3.0/v1.4.0 (та же семантика чисел):
//   1) samtools collate -u -O перед fasta: координатно-сортированный BAM
//      держит мейтов НЕ соседними; без collate часть пар молча теряется
//      в «singletons» или склеивается позиционно в чужие пары.
//   2) fasta с -0/-s: single-end прочтения и singletons не улетают молча
//      в stdout; объединённый <frac>_se.fasta — вход аннотации для SE.
//
// Счётчики фракций — bam_counts.py (primary-only, -F 0x900): публикуются
// в <sample>/qc/counts/*.counts.json и читаются QC gate. flagstat
// сохраняется только как вход MultiQC.
// =============================================================================
process BAM2FASTA {
    tag { "${meta.id}/${meta.fraction}" }
    label 'process_low'

    // NIPT-7 (v1.8.0): samtools + python3 (bam_counts.py). Точного mulled-образа
    // с этой комбинацией версий нет — all-in-one образ пайплайна.
    conda 'bioconda::samtools=1.24 conda-forge::python=3.14.7'
    container params.pipeline_image

    publishDir "${params.outdir}", mode: 'copy', pattern: '*.flagstat.txt',
        saveAs: { fn -> "${meta.id}/qc/flagstat/${fn}" }
    publishDir "${params.outdir}", mode: 'copy', pattern: '*.counts.json',
        saveAs: { fn -> "${meta.id}/qc/counts/${fn}" }

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("${meta.fraction}_1.fasta"),
          path("${meta.fraction}_2.fasta"),
          path("${meta.fraction}_se.fasta"), emit: seqs
    tuple val(meta), path("${meta.id}.${meta.fraction}.counts.json"),
          emit: counts
    tuple val(meta), path("${meta.id}.${meta.fraction}.flagstat.txt"),
          emit: stats

    script:
    // Имена counts/flagstat содержат id образца: защита от коллизий при
    // стейджинге и от слияния сэмплов в MultiQC.
    // -F 0x900 у fasta — это дефолт samtools, оставляем явно: secondary и
    // supplementary записи не должны попадать в аннотацию.
    """
    samtools collate -u -O -@ ${task.cpus} ${bam} | \\
        samtools fasta -@ ${task.cpus} -F 0x900 \\
            -1 ${meta.fraction}_1.fasta -2 ${meta.fraction}_2.fasta \\
            -0 ${meta.fraction}_0.fasta -s ${meta.fraction}_s.fasta -
    cat ${meta.fraction}_0.fasta ${meta.fraction}_s.fasta \\
        > ${meta.fraction}_se.fasta

    bam_counts.py --bam ${bam} --sample ${meta.id} \\
        --fraction ${meta.fraction} \\
        --out ${meta.id}.${meta.fraction}.counts.json

    samtools flagstat -@ ${task.cpus} ${bam} \\
        > ${meta.id}.${meta.fraction}.flagstat.txt
    """
}
