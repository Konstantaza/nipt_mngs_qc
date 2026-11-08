#!/usr/bin/env nextflow
// =============================================================================
// nipt_mngs_qc — метагеномный анализ поверх данных НИПТ с per-sample QC gate
//
// Схема (v1.5.0, по GNP-3):
//   mode=fastq: RAW fastq -> fastp -> DEHOST (один проход HISAT2
//     --no-spliced-alignment; samtools view по флагам ->
//     dehost | host | half bam — точная партиция trimmed-прочтений)
//   mode=bam:   готовые DEHOST/HOST bam от NIPD (как в GNP-3)
//     -> BAM2FASTA per-fraction (collate|fasta + counts + flagstat)
//     -> [если 'all' в params.fractions] MAKE_ALL (конкатенация FASTA
//        конституентов + сумма counts; без нового выравнивания)
//     -> ANNOTATE — ОДИН И ТОТ ЖЕ subworkflow на всех фракциях из
//        params.fractions x строки dbsheet -> контракт profile v2
//     -> comparison model -> Differential profiles -> QC gate -> MultiQC
//
// Две точки входа (params.mode):
//   fastq — полный путь от сырых прочтений (fastp + DEHOST внутри)
//   bam   — старт от готовых DEHOST/HOST bam NIPD-системы (основной режим)
//
// Запуск:
//   nextflow run main.nf --mode bam   --input samplesheet_bam.csv \
//        --dbsheet dbsheet.csv
//   nextflow run main.nf --mode fastq --input samplesheet_fastq.csv \
//        --host_index /path/hisat2/human --dbsheet dbsheet.csv
//   # ветка all (все прочтения без dehost) включается одним параметром:
//   ... --fractions dehost,host,all
//
// dbsheet CSV (заменяет --db_full/--db_patho, v1.4.0):
//   db_id,annotator,path,role
//   kraken2_std,kraken2_bracken,/path/kraken2_std,full
//   eukpathdb,kraken2_bracken,/path/eukpathdb,patho
// Контракт профиля и интерфейс ANNOTATE — docs/profile_contract.md.
//
// v1.6.0: fail-fast валидация входов — обязательные параметры по режиму
//   (input, dbsheet, host_index для mode=fastq), обязательные колонки,
//   непустые значения и уникальность sample проверяются eager, ДО старта
//   задач; пути из CSV: абсолютные и URI (s3://, https://, ...) — как есть,
//   относительные — от каталога самого CSV; HISAT2-индекс стейджится
//   glob'ом <prefix>.*.ht2 с checkIfExists (посторонние файлы рядом с
//   префиксом не захватываются).
// =============================================================================

nextflow.enable.dsl = 2

include { FASTP }     from './modules/local/fastp'
include { DEHOST }    from './modules/local/dehost'
include { BAM2FASTA } from './modules/local/bam2fasta'
include { MAKE_ALL }  from './modules/local/make_all'
include { ANNOTATE }  from './subworkflows/local/annotate'
include { COMPARE }   from './modules/local/compare'
include { QC_GATE }   from './modules/local/qc_gate'
include { MULTIQC }   from './modules/local/multiqc'

// =============================================================================
// resolveCsvPath (v1.6.0) — резолвинг пути из samplesheet/dbsheet CSV:
//   абсолютный (/...) и URI (s3://, https://, file://, ...) — как есть;
//   относительный — от каталога самого CSV (а не от launchDir): запуск
//   работает из любого каталога. checkIfExists: отсутствующий файл —
//   понятная ошибка вместо молчаливого пустого канала.
// =============================================================================
def resolveCsvPath(base_dir, String p) {
    def is_uri = p ==~ /^[A-Za-z][A-Za-z0-9+.-]*:\/\/.*/
    return (p.startsWith('/') || is_uri) ?
        file(p, checkIfExists: true) :
        file("${base_dir}/${p}", checkIfExists: true)
}

workflow {

    // ------------------------------------------------------------------
    // 0. Обязательные параметры по режиму — fail fast ДО построения
    //    каналов и старта задач (v1.6.0). Без этого: file(null) с
    //    невнятной ошибкой; пустой канал индекса -> граф молча не
    //    запускается и пайплайн «успешно» завершается без результатов.
    // ------------------------------------------------------------------
    if (!(params.mode in ['fastq', 'bam'])) {
        error "params.mode должен быть 'fastq' или 'bam', получено: '${params.mode}'"
    }
    if (!params.input) {
        error "не задан обязательный параметр --input (samplesheet CSV)"
    }
    if (!params.dbsheet) {
        error "не задан обязательный параметр --dbsheet (CSV: db_id,annotator,path,role)"
    }
    if (params.mode == 'fastq' && !params.host_index) {
        error "mode=fastq требует --host_index (префикс HISAT2-индекса человека для стадии DEHOST)"
    }

    // ------------------------------------------------------------------
    // 0b. Список фракций: params.fractions (default 'dehost,host').
    //    dehost и host обязательны (их профили/сводки читает QC gate);
    //    all — опциональная ветка «все прочтения без dehost». half
    //    (полукартированные) — QC-категория: считается и проверяется в
    //    consistency, но в аннотацию не идёт, поэтому в списке
    //    фракций её нет.
    // ------------------------------------------------------------------
    def fractions = params.fractions.split(',').collect { it.trim() }
    def allowed_fractions = ['dehost', 'host', 'all']
    if (!fractions || fractions.any { !(it in allowed_fractions) }) {
        error "params.fractions: допустимы только ${allowed_fractions}, " +
              "получено: '${params.fractions}'"
    }
    if (fractions.unique().size() != fractions.size()) {
        error "params.fractions: дубликаты в '${params.fractions}'"
    }
    if (!('dehost' in fractions && 'host' in fractions)) {
        error "params.fractions: 'dehost' и 'host' обязательны " +
              "(их профили читает QC gate), получено: '${params.fractions}'"
    }

    // ------------------------------------------------------------------
    // 1. Вход: samplesheet CSV — eager-парсинг и валидация (v1.6.0):
    //    обязательные колонки режима, непустые значения, уникальность
    //    sample — всё до старта задач. Path.splitCsv — тот же парсер,
    //    что у канального splitCsv: валидация и downstream видят одни и
    //    те же строки. Пути резолвятся через resolveCsvPath (см. выше).
    // ------------------------------------------------------------------
    def sheet_file = file(params.input, checkIfExists: true)
    def sheet_dir  = sheet_file.parent
    def req_cols = (params.mode == 'fastq') ?
        ['sample', 'fastq_1', 'fastq_2'] : ['sample', 'dehost_bam', 'host_bam']
    def sheet_rows = sheet_file.splitCsv(header: true)
    if (!sheet_rows) {
        error "samplesheet ${params.input}: нет строк данных"
    }
    def header_cols = sheet_rows[0].keySet() as List
    def missing_cols = req_cols - header_cols
    if (missing_cols) {
        error "samplesheet ${params.input} (mode=${params.mode}): нет " +
              "обязательных колонок ${missing_cols}; заголовок: ${header_cols}"
    }
    def sample_ids = sheet_rows.collect { r -> r.sample?.trim() }
    if (sample_ids.any { !it }) {
        error "samplesheet ${params.input}: пустое значение в колонке 'sample'"
    }
    def dup_ids = sample_ids.countBy { it }.findAll { k, v -> v > 1 }.keySet()
    if (dup_ids) {
        error "samplesheet ${params.input}: дубликаты sample: ${dup_ids.sort()}"
    }
    sheet_rows.eachWithIndex { r, i ->
        def value = (r.get('single_end') ?: '').trim().toLowerCase()
        if (params.mode == 'fastq' && !(value in ['', 'true', 'false'])) {
            error "samplesheet ${params.input}: строка ${i + 2}: single_end должен быть true или false"
        }
        def se = params.mode == 'fastq' && value == 'true'
        def cells = (req_cols - 'sample') - (se ? ['fastq_2'] : [])
        cells.each { c ->
            if (!r[c]?.trim()) {
                error "samplesheet ${params.input}: строка ${i + 2}: " +
                      "пустое значение в колонке '${c}'"
            }
        }
        if (se && r.fastq_2?.trim()) {
            error "samplesheet ${params.input}: строка ${i + 2}: single_end=true требует пустой fastq_2"
        }
        r.sample = r.sample.trim()
        if (params.mode == 'fastq') r.single_end = se
    }

    ch_samples = Channel.from(sheet_rows).map { row -> [row.sample] }

    if (params.mode == 'fastq') {
        ch_input = Channel.from(sheet_rows)
            .map { row ->
                def reads = [resolveCsvPath(sheet_dir, row.fastq_1.trim())]
                if (!row.single_end) reads << resolveCsvPath(sheet_dir, row.fastq_2.trim())
                [[id: row.sample, single_end: row.single_end], reads]
            }

        FASTP(ch_input)

        // HISAT2-индекс человека: стейджим ТОЛЬКО файлы префикса *.ht2
        // (посторонние файлы рядом с префиксом не захватываются). Eager-
        // проверка glob'а — до старта задач; checkIfExists — страховка
        // на уровне канала. Пустой glob = понятная ошибка, а не
        // молчаливый пустой канал (поведение «null.*» до v1.6.0).
        def hp = file(params.host_index)
        def hglob = "${hp.name}.*.ht2{,l}"
        def hmatches = []
        try {
            ['ht2', 'ht2l'].each { ext ->
                java.nio.file.Files.newDirectoryStream(hp.parent, "${hp.name}.*.${ext}")
                    .withCloseable { ds -> hmatches.addAll(ds.collect { it }) }
            }
        } catch (java.nio.file.NoSuchFileException e) {
            error "--host_index '${params.host_index}': каталог '${hp.parent}' не существует"
        }
        if (!hmatches) {
            error "--host_index '${params.host_index}': не найдено файлов " +
                  "'${hglob}' (HISAT2-индекс человека)"
        }
        ch_index = Channel.fromPath("${params.host_index}.*.ht2{,l}", checkIfExists: true).collect()
        DEHOST(FASTP.out.reads, ch_index)

        // Три фракции — точная партиция trimmed-прочтений (критерий 3)
        ch_bams = DEHOST.out.bams.flatMap { s, d, h, hf ->
            [[[id: s, fraction: 'dehost'], d],
             [[id: s, fraction: 'host'],   h],
             [[id: s, fraction: 'half'],   hf]]
        }
        ch_fastp = FASTP.out.json                 // (sample, fastp.json)
    }
    else if (params.mode == 'bam') {
        ch_bams = Channel.from(sheet_rows)
            .flatMap { row ->
                def s = row.sample
                [[[id: s, fraction: 'dehost'],
                  resolveCsvPath(sheet_dir, row.dehost_bam)],
                 [[id: s, fraction: 'host'],
                  resolveCsvPath(sheet_dir, row.host_bam)]]
            }

        // fastp выполнен upstream в NIPD — подставляем пустой stub,
        // QC gate пометит fastp-критерии как «не оцениваются»
        ch_fastp = Channel.from(sheet_rows)
            .map { row -> [row.sample,
                           file("${projectDir}/assets/no_fastp.json")] }
    }
    else {
        error "params.mode должен быть 'fastq' или 'bam', получено: '${params.mode}'"
    }

    // ------------------------------------------------------------------
    // 2. BAM -> FASTA per-fraction (collate восстанавливает пары)
    //    + counts JSON (общее ядро для обеих точек входа)
    // ------------------------------------------------------------------
    BAM2FASTA(ch_bams)
    // seqs:   (meta[id,fraction], <f>_1.fasta, <f>_2.fasta, <f>_se.fasta)
    // counts: (meta[id,fraction], <id>.<f>.counts.json)
    // stats:  (meta[id,fraction], <id>.<f>.flagstat.txt)

    // ------------------------------------------------------------------
    // 2b. Ветка all (только при 'all' в params.fractions): конкатенация
    //     FASTA конституентов + сумма counts. Конституенты: fastq —
    //     dehost+host+half (партиция trimmed), bam — dehost+host
    //     (партиция NIPD). groupKey с явным size — groupTuple не
    //     зависает на неполных группах.
    // ------------------------------------------------------------------
    if ('all' in fractions) {
        def n_constituents = (params.mode == 'fastq') ? 3 : 2
        ch_all_seqs = BAM2FASTA.out.seqs
            .map { meta, r1, r2, se ->
                [groupKey(meta.id, n_constituents), r1, r2, se]
            }
            .groupTuple()
        ch_all_counts = BAM2FASTA.out.counts
            .map { meta, cj -> [groupKey(meta.id, n_constituents), cj] }
            .groupTuple()
        MAKE_ALL(ch_all_seqs.join(ch_all_counts, failOnMismatch: true))
    }

    // ------------------------------------------------------------------
    // 3. Профилирование: subworkflow ANNOTATE (интерфейс v1.4.0) — один
    //    и тот же для всех фракций (критерий 1). reads: [meta(id,
    //    fraction, single_end), r1, r2, se]; single_end фракции читаем
    //    из её counts JSON (автоопределение BAM2FASTA по флагу PAIRED).
    //    join seqs~counts по кортежному ключу meta с failOnMismatch:
    //    обе стороны обязаны существовать (внутренний join).
    // ------------------------------------------------------------------
    ch_reads = BAM2FASTA.out.seqs
        .join(BAM2FASTA.out.counts, failOnMismatch: true)
        .filter { meta, r1, r2, se, cj -> meta.fraction in fractions }
        .map { meta, r1, r2, se, cj ->
            def j = new groovy.json.JsonSlurper().parseText(cj.text)
            [meta + [single_end: j.single_end], r1, r2, se]
        }
    if ('all' in fractions) {
        ch_reads = ch_reads.mix(
            MAKE_ALL.out.seqs
                .join(MAKE_ALL.out.counts, failOnMismatch: true)
                .map { s, r1, r2, se, cj ->
                    def j = new groovy.json.JsonSlurper().parseText(cj.text)
                    [[id: s, fraction: 'all', single_end: j.single_end],
                     r1, r2, se]
                })
    }

    // dbsheet CSV — eager-парсинг и валидация (v1.6.0): обязательные
    // колонки, непустые значения, уникальность ролей, существование
    // каталогов БД — всё до старта задач. Пути резолвятся тем же
    // resolveCsvPath: относительные — от каталога dbsheet, абсолютные
    // и URI — как есть.
    def dbsheet_file = file(params.dbsheet, checkIfExists: true)
    def dbsheet_dir  = dbsheet_file.parent
    def dbs_req  = ['db_id', 'annotator', 'path', 'role']
    def dbs_rows = dbsheet_file.splitCsv(header: true)
    if (!dbs_rows) {
        error "dbsheet ${params.dbsheet}: нет строк данных"
    }
    def missing_dbs = dbs_req - (dbs_rows[0].keySet() as List)
    if (missing_dbs) {
        error "dbsheet ${params.dbsheet}: нет обязательных колонок " +
              "${missing_dbs}; заголовок: ${dbs_rows[0].keySet() as List}"
    }
    dbs_rows = dbs_rows.collect { row ->
        dbs_req.each { col ->
            if (!row[col]?.trim()) {
                error "dbsheet ${params.dbsheet}: пустое значение в колонке '${col}'"
            }
        }
        [db_id: row.db_id, annotator: row.annotator, role: row.role,
         path: resolveCsvPath(dbsheet_dir, row.path)]
    }
    def roles = dbs_rows.collect { it.role }
    if (roles.unique().size() != roles.size()) {
        error "dbsheet ${params.dbsheet}: роли должны быть уникальны, получено: ${roles}"
    }
    ch_dbs = Channel.from(dbs_rows)
    // Число строк dbsheet — размер для groupKey ниже; value-канал
    // допускает несколько подписчиков (combine в ch_contracts и ch_diffs)
    ch_n_dbs = Channel.value(dbs_rows.size())

    ANNOTATE(ch_reads, ch_dbs)
    // ANNOTATE.out.profiles: (meta[id, fraction, single_end, db_id, role,
    //                          annotator], profile.tsv, profile.json, native)

    // ------------------------------------------------------------------
    // 4. Comparison model v2: сравнение профилей фракций на контракте
    //    profile v2. COMPARE питается от profile.tsv + sidecar JSON (НЕ
    //    от сырых выходов классификатора). Пары — из params.compare_pairs;
    //    dehost:host обязательна (её сводки читает QC gate). Join по
    //    кортежному ключу [id, role] с failOnMismatch (строковые ключи
    //    v1.4.0 хрупки; обе фракции обязаны существовать для роли).
    // ------------------------------------------------------------------
    def compare_pairs = params.compare_pairs.split(',').collect { it.trim() }
    if (!('dehost:host' in compare_pairs)) {
        error "params.compare_pairs должен содержать 'dehost:host' " +
              "(её сводки читает QC gate), получено: '${params.compare_pairs}'"
    }

    ch_ctr = ANNOTATE.out.profiles.map { meta, tsv, js, nat ->
        [[meta.id, meta.role], meta.id, meta.fraction, meta.role, tsv, js]
    }
    ch_ctr_d = ch_ctr.filter { it[2] == 'dehost' }
        .map { key, s, f, role, tsv, js -> [key, s, role, tsv, js] }
    ch_ctr_h = ch_ctr.filter { it[2] == 'host' }
        .map { key, s, f, role, tsv, js -> [key, tsv, js] }

    COMPARE(ch_ctr_d.join(ch_ctr_h, failOnMismatch: true)
        .map { key, s, role, dtsv, djs, htsv, hjs ->
            [s, role, dtsv, djs, htsv, hjs]
        }
        .combine(Channel.from(compare_pairs)))

    // ------------------------------------------------------------------
    // 5. QC gate: КАЖДЫЙ образец samplesheet получает QC JSON.
    //    Сборка идёт ОТ авторитетного канала образцов левыми join
    //    (remainder: true): отсутствующие артефакты -> stub/пустые списки,
    //    qc_gate.py вынесет FAIL/PIPELINE_ERROR — образец не пропадает молча.
    //    groupKey с явным size (n_fractions x n_dbs / n_dbs) вместо
    //    groupTuple без size: группа завершается по счётчику, а не по
    //    закрытию канала.
    // ------------------------------------------------------------------
    ch_contracts = ANNOTATE.out.profiles
        .combine(ch_n_dbs)
        .map { meta, tsv, js, nat, ndbs ->
            [groupKey(meta.id, fractions.size() * ndbs), tsv, js]
        }
        .groupTuple()                       // (sample, [tsv], [json])
    // В QC gate идут только сводки обязательной пары dehost:host
    ch_diffs = COMPARE.out.diff
        .filter { s, role, pair, tsv, js -> pair == 'dehost:host' }
        .combine(ch_n_dbs)
        .map { s, role, pair, tsv, js, ndbs -> [groupKey(s, ndbs), js] }
        .groupTuple()                       // (sample, [diff json])

    // counts per фракция -> per-sample каналы (half в bam-режиме
    // отсутствует -> remainder даст null -> stub -> half=0)
    ch_counts = BAM2FASTA.out.counts
        .map { meta, cj -> [meta.id, meta.fraction, cj] }
    ch_cd    = ch_counts.filter { it[1] == 'dehost' }
        .map { s, f, cj -> [s, cj] }
    ch_ch    = ch_counts.filter { it[1] == 'host' }
        .map { s, f, cj -> [s, cj] }
    ch_chalf = ch_counts.filter { it[1] == 'half' }
        .map { s, f, cj -> [s, cj] }

    ch_qc = ch_samples
        .join(ch_cd, remainder: true)           // (s, cd | null)
        .join(ch_ch, remainder: true)           // + ch | null
        .join(ch_chalf, remainder: true)        // + chalf | null
        .join(ch_fastp, remainder: true)        // + fastp json | null
        .join(ch_contracts, remainder: true)    // + [tsv],[json] | null
        .join(ch_diffs, remainder: true)        // + [diff json] | null
        .map { s, cd, ch, chalf, fj, tsvs, jsons, diffs ->
            [ s,
              cd    ?: file("${projectDir}/assets/missing.txt"),
              ch    ?: file("${projectDir}/assets/missing.txt"),
              chalf ?: file("${projectDir}/assets/missing.txt"),
              fj    ?: file("${projectDir}/assets/no_fastp.json"),
              tsvs ?: [],
              jsons ?: [],
              diffs ?: [] ]
        }

    // Панель мишеней необязательна; пустой stub = «не задана»
    def panel_file = params.target_panel ?
        file(params.target_panel, checkIfExists: true) :
        file("${projectDir}/assets/no_panel.txt")
    QC_GATE(ch_qc, Channel.value(panel_file))

    // ------------------------------------------------------------------
    // 6. MultiQC по застейдженным входам (params.outdir НЕ читается).
    //    native адаптеров миксуется целиком: MultiQC распознаёт kraken
    //    report, неродные форматы (напр. metax stub) молча игнорирует.
    //    В fastq-режиме добавляется summary HISAT2 стадии DEHOST.
    // ------------------------------------------------------------------
    ch_mqc_files = QC_GATE.out.mqc
        .mix(ANNOTATE.out.profiles.map { meta, tsv, js, nat -> nat })
        .mix(BAM2FASTA.out.stats.map { meta, f -> f })
    if (params.mode == 'fastq') {
        ch_mqc_files = ch_mqc_files
            .mix(FASTP.out.json.map { s, j -> j })
            .mix(FASTP.out.html)
            .mix(DEHOST.out.summary.map { s, f -> f })
    }
    MULTIQC(ch_mqc_files.collect())
}
