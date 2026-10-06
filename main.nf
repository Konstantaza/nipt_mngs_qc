#!/usr/bin/env nextflow
// =============================================================================
// nipt_mngs_qc — метагеномный анализ поверх данных НИПТ с per-sample QC gate
//
// Схема (по ТЗ задачи 3):
//   RAW fastq -> fastp -> TRIM fastq -> hocort -> DEHOST bam + HOST bam
//     -> samtools bam2fa -> DEHOST/HOST fasta
//     -> kraken2+bracken x {full DB, patho DB (EukPathDB)}
//     -> comparison model -> Differential profiles -> results
//
// Две точки входа (params.mode):
//   fastq — полный путь от сырых прочтений (fastp + hocort внутри)
//   bam   — старт от готовых DEHOST/HOST bam NIPD-системы (основной режим)
//
// Запуск:
//   nextflow run main.nf --mode bam   --input samplesheet_bam.csv \
//        --db_full /path/kraken2_std --db_patho /path/eukpathdb
//   nextflow run main.nf --mode fastq --input samplesheet_fastq.csv \
//        --host_index /path/hisat2/human --db_full ... --db_patho ...
// =============================================================================

nextflow.enable.dsl = 2

include { FASTP }        from './modules/local/fastp'
include { HOCORT_SPLIT } from './modules/local/hocort_split'
include { BAM2FASTA }    from './modules/local/bam2fasta'
include { PROFILE }      from './modules/local/profile'
include { COMPARE }      from './modules/local/compare'
include { QC_GATE }      from './modules/local/qc_gate'
include { MULTIQC }      from './modules/local/multiqc'

workflow {

    // ------------------------------------------------------------------
    // 1. Вход: samplesheet CSV -> канал образцов
    // ------------------------------------------------------------------
    if (params.mode == 'fastq') {
        // Относительные пути в samplesheet резолвим от директории самого
        // samplesheet (а не от launchDir): так тестовые данные работают
        // из любого каталога запуска.
        def sheet_dir = file(params.input, checkIfExists: true).parent
        ch_input = Channel.fromPath(params.input, checkIfExists: true)
            .splitCsv(header: true)
            .map { row -> [row.sample,
                           file("${sheet_dir}/${row.fastq_1}", checkIfExists: true),
                           file("${sheet_dir}/${row.fastq_2}", checkIfExists: true)] }

        FASTP(ch_input)

        // HISAT2-индекс человека: стейджим все файлы префикса (idx.*.ht2)
        ch_index = Channel.fromPath("${params.host_index}.*").collect()
        HOCORT_SPLIT(FASTP.out.reads, ch_index)

        ch_bams  = HOCORT_SPLIT.out.bams          // (sample, dehost.bam, host.bam)
        ch_fastp = FASTP.out.json                 // (sample, fastp.json)
    }
    else if (params.mode == 'bam') {
        def sheet_dir = file(params.input, checkIfExists: true).parent
        ch_bams = Channel.fromPath(params.input, checkIfExists: true)
            .splitCsv(header: true)
            .map { row -> [row.sample,
                           file("${sheet_dir}/${row.dehost_bam}", checkIfExists: true),
                           file("${sheet_dir}/${row.host_bam}", checkIfExists: true)] }

        // fastp выполнен upstream в NIPD — подставляем пустой stub,
        // QC gate пометит fastp-критерии как «не оцениваются»
        ch_fastp = ch_bams.map { s, d, h ->
            [s, file("${projectDir}/assets/no_fastp.json")]
        }
    }
    else {
        error "params.mode должен быть 'fastq' или 'bam', получено: '${params.mode}'"
    }

    // ------------------------------------------------------------------
    // 2. BAM -> FASTA + flagstat (общее ядро для обеих точек входа)
    // ------------------------------------------------------------------
    BAM2FASTA(ch_bams)

    // ------------------------------------------------------------------
    // 3. Профилирование: fork на {dehost, host} x {full, patho} = 4 запуска
    // ------------------------------------------------------------------
    ch_frac = BAM2FASTA.out.seqs.flatMap { s, d1, d2, h1, h2 ->
        [[s, 'dehost', d1, d2], [s, 'host', h1, h2]]
    }
    ch_dbs = Channel.of(['full', params.db_full], ['patho', params.db_patho])

    ch_prof_in = ch_frac.combine(ch_dbs).map { s, frac, r1, r2, dbname, dbpath ->
        [s, frac, dbname, r1, r2, file(dbpath, checkIfExists: true)]
    }
    PROFILE(ch_prof_in)

    // ------------------------------------------------------------------
    // 4. Comparison model: join dehost+host профилей по (sample, db)
    // ------------------------------------------------------------------
    ch_br = PROFILE.out.profiles.map { s, frac, dbname, krep, brk ->
        ["${s}|||${dbname}", s, frac, dbname, brk]
    }
    ch_dehost = ch_br.filter { it[2] == 'dehost' }
        .map { key, s, f, db, brk -> [key, s, db, brk] }
    ch_host = ch_br.filter { it[2] == 'host' }
        .map { key, s, f, db, brk -> [key, brk] }

    COMPARE(ch_dehost.join(ch_host)
        .map { key, s, db, brd, brh -> [s, db, brd, brh] })

    // ------------------------------------------------------------------
    // 5. QC gate: собрать все метрики образца в один кортеж
    // ------------------------------------------------------------------
    ch_kreports = PROFILE.out.profiles
        .map { s, frac, dbname, krep, brk -> [s, krep] }
        .groupTuple()                                   // (sample, [4 report])
    ch_diffs = COMPARE.out.diff
        .map { s, db, tsv, js -> [s, js] }
        .groupTuple()                                   // (sample, [2 json])

    ch_qc = BAM2FASTA.out.stats                         // (s, fd, fh)
        .join(ch_fastp)                                 // + fastp.json
        .join(ch_kreports)                              // + [kraken reports]
        .join(ch_diffs)                                 // + [diff jsons]

    QC_GATE(ch_qc)

    // ------------------------------------------------------------------
    // 6. MultiQC по опубликованным QC-артефактам
    // ------------------------------------------------------------------
    ch_mqc_deps = QC_GATE.out.mqc
        .mix(PROFILE.out.profiles.map { s, f, d, krep, brk -> krep })
        .collect()
    // params.outdir может быть относительным — MULTIQC работает в своей
    // work-директории, поэтому передаём абсолютный путь (resolve от launchDir)
    MULTIQC(ch_mqc_deps, file(params.outdir).toAbsolutePath().toString())
}
