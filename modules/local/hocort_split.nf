// =============================================================================
// HOCORT_SPLIT — разделение прочтений на DEHOST и HOST фракции (mode=fastq)
//
// Теория: 95–99% прочтений НИПТ — ДНК человека (материнская + фетальная вкДНК).
// HoCoRT картирует прочтения на геном человека (здесь — HISAT2-индекс) и
// раздаёт их по SAM-флагам:
//   -f true  -> unmapped пары  (samtools view -f 13: paired+unmapped+mate-unmapped)
//   -f false -> mapped пары    (samtools view -F 12 -f 1)
// Мы прогоняем hocort дважды (dehost / host) и конвертируем FASTQ -> BAM
// (samtools import), чтобы интерфейс стадии совпадал с контрактом
// NIPD-системы (DEHOST/HOST bam).
//
// Замечание по эффективности: два прогона = двойное картирование. В проде,
// где вход — готовые BAM от NIPD, эта стадия не выполняется (mode=bam).
// Однопроходная альтернатива (hisat2 + samtools view -f 13 / -F 12)
// описана в отчёте.
// =============================================================================
process HOCORT_SPLIT {
    tag { sample_id }
    label 'process_medium'

    input:
    tuple val(sample_id), path(r1), path(r2)
    path index_files   // все файлы HISAT2-индекса человека (стейджатся в cwd)

    output:
    tuple val(sample_id), path('dehost.bam'), path('host.bam'), emit: bams

    script:
    def idx = file(params.host_index).name   // префикс индекса (без пути)
    """
    # DEHOST: некартированные пары (микробная фракция)
    hocort map hisat2 \
        -x ${idx} \
        -i ${r1} ${r2} \
        -o dehost_R1.fastq dehost_R2.fastq \
        -t ${task.cpus} -f true

    # HOST: картированные пары (человеческая фракция)
    hocort map hisat2 \
        -x ${idx} \
        -i ${r1} ${r2} \
        -o host_R1.fastq host_R2.fastq \
        -t ${task.cpus} -f false

    # FASTQ -> BAM (контракт интерфейса NIPD: на выходе стадии — BAM)
    samtools import -@ ${task.cpus} -1 dehost_R1.fastq -2 dehost_R2.fastq -o dehost.bam
    samtools import -@ ${task.cpus} -1 host_R1.fastq   -2 host_R2.fastq   -o host.bam
    """
}
