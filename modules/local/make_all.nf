// =============================================================================
// MAKE_ALL — синтез фракции all (все прочтения без dehost), v1.5.0
//
// all = конкатенация конституентов: fastq-режим — dehost+host+half
// (= trimmed, точная партиция DEHOST); bam-режим — dehost+host
// (партиция NIPD принимается полной — политика half у NIPD недоступна).
// Новый прогон выравнивателя НЕ нужен: FASTA конституентов уже построены
// per-fraction BAM2FASTA.
//
// Синхронность мейтов: cat по алфавитному glob (*_1.fasta, *_2.fasta,
// *_se.fasta) — порядок конституентов (dehost < half < host) одинаков
// во всех трёх glob, а внутри каждого файла мейты стоят на одинаковых
// позициях => пары не перемешиваются.
//
// counts JSON — сумма счётчиков конституентов (bin/sum_counts.py;
// single_end = AND). Процесс включается в граф только при 'all' в
// params.fractions (критерий приёмки 4: ветка включается одним
// параметром, без правок кода).
// =============================================================================
process MAKE_ALL {
    tag { sample_id }
    label 'process_low'

    // NIPT-7 (v1.8.0): только python3 (sum_counts.py) + coreutils (cat)
    conda 'conda-forge::python=3.14.7'
    container 'quay.io/biocontainers/python:3.14'

    publishDir "${params.outdir}", mode: 'copy', pattern: '*.counts.json',
        saveAs: { fn -> "${sample_id}/qc/counts/${fn}" }

    input:
    // Списки файлов конституентов (имена не пересекаются: префикс —
    // фракция), сгруппированные по образцу через groupKey в main.nf
    tuple val(sample_id), path(r1_files), path(r2_files), path(se_files),
          path(counts_files)

    output:
    tuple val(sample_id), path('all_1.fasta'), path('all_2.fasta'),
          path('all_se.fasta'), emit: seqs
    tuple val(sample_id), path("${sample_id}.all.counts.json"), emit: counts

    script:
    """
    cat *_1.fasta  > all_1.fasta
    cat *_2.fasta  > all_2.fasta
    cat *_se.fasta > all_se.fasta

    sum_counts.py --sample ${sample_id} --fraction all \\
        --out ${sample_id}.all.counts.json ${counts_files}
    """
}
