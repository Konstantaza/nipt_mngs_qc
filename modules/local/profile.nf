// =============================================================================
// PROFILE — таксономическое профилирование Kraken2 + Bracken
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
// Процесс запускается 4 раза на образец: {dehost, host} x {full, patho} БД.
// Публикуем только компактные report/TSV; сырой per-read вывод Kraken2
// не сохраняем (storage hygiene).
// =============================================================================
process PROFILE {
    tag { "${sample_id}/${fraction}/${dbname}" }
    label 'process_kraken'

    publishDir "${params.outdir}", mode: 'copy',
        saveAs: { fn -> fn.endsWith('.report') ?
            "${sample_id}/qc/kraken/${fn}" : "${sample_id}/profiles/${fn}" }

    input:
    tuple val(sample_id), val(fraction), val(dbname),
          path(r1), path(r2), path(db)

    output:
    tuple val(sample_id), val(fraction), val(dbname),
          path("${sample_id}.${fraction}.${dbname}.kraken2.report"),
          path("${sample_id}.${fraction}.${dbname}.bracken.tsv"), emit: profiles

    script:
    def prefix = "${sample_id}.${fraction}.${dbname}"
    """
    kraken2 --db ${db} \
        --threads ${task.cpus} \
        --confidence ${params.kraken_confidence} \
        --paired \
        --report ${prefix}.kraken2.report \
        --output /dev/null \
        ${r1} ${r2}

    # Edge case «нулевой микробиом»: bracken падает (exit 1, "no reads found"),
    # если в report НЕТ ни одного таксона уровня species (-l S) с числом
    # кладных прочтений >= порога (-t). Это шире, чем «пустой report»:
    # report может содержать прочтения, но все ниже порога Bracken.
    # Воспроизводим точное условие Bracken: строки с кодом ранга 'S'
    # и clade-прочтениями (колонка 2 report) >= threshold. Если таких нет — пишем пустой
    # bracken-файл с заголовком и идём дальше: вердикт вынесет QC gate
    # (FALLBACK / ZERO_MICROBIOME), а не падение стадии.
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
    """
}
