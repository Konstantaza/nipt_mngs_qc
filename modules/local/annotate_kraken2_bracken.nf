// =============================================================================
// ANNOTATE_KRAKEN2_BRACKEN — адаптер Kraken2+Bracken интерфейса ANNOTATE
//
// Теория:
//   Kraken2 разбивает каждое прочтение на k-меры (k=35) и ищет точные
//   совпадения в БД; таксон прочтения = наименьший общий предок (LCA)
//   всех совпавших k-меров. --confidence задаёт минимальную долю k-меров,
//   согласных с классификацией (снижает ложноположительные вызовы).
//   Bracken затем байесовски переоценивает видовые изобилия: перераспределяет
//   прочтения, которые Kraken2 смог отнести только к роду/семейству,
//   пропорционально уникальным k-мерам видов.
//
// Адаптер = этот процесс + bin/kraken2_bracken_to_profile.py.
// Выход — контракт profile v2 (profile.tsv + sidecar profile.json,
// docs/profile_contract.md); native — kraken2 report (в MultiQC).
// Публикуем только компактные report/TSV; сырой per-read вывод Kraken2
// не сохраняем (storage hygiene).
// =============================================================================
process ANNOTATE_KRAKEN2_BRACKEN {
    tag { "${meta.id}/${meta.fraction}/${meta.role}" }
    label 'process_kraken'

    // NIPT-7 (v1.8.0): kraken2+bracken+python3. Готового mulled-образа с этой
    // точной комбинацией версий на quay.io нет — all-in-one образ пайплайна.
    conda 'bioconda::kraken2=2.17.2 bioconda::bracken=3.1 conda-forge::python=3.14.7'
    container params.pipeline_image

    publishDir "${params.outdir}", mode: 'copy',
        saveAs: { fn -> fn.endsWith('.report') ?
            "${meta.id}/qc/kraken/${fn}" : "${meta.id}/profiles/${fn}" }

    input:
    // meta: [id, fraction, single_end, db_id, role, annotator]
    // single_end=true -> kraken2 на se (БЕЗ --paired); иначе --paired r1 r2.
    // se — объединённый файл непарных прочтений (-0 + -s из BAM2FASTA).
    tuple val(meta), path(r1), path(r2), path(se), path(db)

    output:
    // Контракт profile v2 + native (kraken2 report для MultiQC)
    tuple val(meta), path("${meta.id}.${meta.fraction}.${meta.role}.profile.tsv"),
          path("${meta.id}.${meta.fraction}.${meta.role}.profile.json"),
          path("${meta.id}.${meta.fraction}.${meta.role}.kraken2.report"),
          emit: annotated

    script:
    def prefix = "${meta.id}.${meta.fraction}.${meta.role}"
    """
    # Single-end фракция (автоопределение BAM2FASTA: нет primary-прочтений
    # с флагом PAIRED) аннотируется без --paired — иначе SE-прочтения
    # молча терялись бы (edge case «ложный ZERO_MICROBIOME»).
    if [ "${meta.single_end}" = "true" ]; then
        kraken2 --db ${db} \
            --threads ${task.cpus} \
            --confidence ${params.kraken_confidence} \
            --report ${prefix}.kraken2.report \
            --output /dev/null \
            ${se}
    else
        kraken2 --db ${db} \
            --threads ${task.cpus} \
            --confidence ${params.kraken_confidence} \
            --paired \
            --report ${prefix}.kraken2.report \
            --output /dev/null \
            ${r1} ${r2}
    fi

    # Edge case «нулевой микробиом»: bracken падает (exit 1, "no reads found"),
    # если в report НЕТ ни одного таксона уровня species (-l S) с числом
    # кладных прочтений >= порога (-t). Это шире, чем «пустой report»:
    # report может содержать прочтения, но все ниже порога Bracken.
    # Воспроизводим точное условие Bracken: строки с кодом ранга 'S'
    # и clade-прочтениями (колонка 2 report) >= threshold. Если таких нет — пишем пустой
    # bracken-файл с заголовком и идём дальше: конвертер выпишет
    # status="empty", а вердикт вынесет QC gate (FALLBACK/ZERO_MICROBIOME).
    N_CLASS=\$(awk -F'\\t' -v thr=${params.bracken_min_reads} '\$4 == "S" && \$2+0 >= thr {s += 1} END {print s+0}' ${prefix}.kraken2.report)
    if [ "\$N_CLASS" -gt 0 ]; then
        bracken -d ${db} \
            -i ${prefix}.kraken2.report \
            -o ${prefix}.bracken.tsv \
            -l S -t ${params.bracken_min_reads} -r ${params.bracken_read_len}
    else
        printf 'name\\ttaxonomy_id\\ttaxonomy_lvl\\tkraken_assigned_reads\\tadded_reads\\tnew_est_reads\\tfraction_total_reads\\n' \
            > ${prefix}.bracken.tsv
    fi

    # Конвертация в контракт profile v2; версии аннотаторов — в рантайме
    KRAKEN_VER=\$(kraken2 --version 2>&1 | grep -oE '[0-9]+\\.[0-9.]+' | head -1)
    BRACKEN_VER=\$(bracken -v 2>&1 | grep -oE '[0-9]+\\.[0-9.]+' | head -1)
    kraken2_bracken_to_profile.py \
        --report ${prefix}.kraken2.report \
        --bracken ${prefix}.bracken.tsv \
        --sample ${meta.id} \
        --fraction ${meta.fraction} \
        --role ${meta.role} \
        --db_id ${meta.db_id} \
        --annotator ${meta.annotator} \
        --annotator_version "\${KRAKEN_VER:-unknown}+\${BRACKEN_VER:-unknown}" \
        --out_tsv ${prefix}.profile.tsv \
        --out_json ${prefix}.profile.json
    """
}
