# nipt_mngs_qc

Метагеномный анализ поверх данных НИПТ с формальным per-sample QC gate.
Nextflow DSL2.

**Документация:**
- `report_pipeline_theory.md` — теория инструментов, устройство пайплайна,
  как читать/расширять (образовательный отчёт);
- `QC_checklist.md` — формальные критерии QC/адекватности, пороги,
  fallback-правила, контракт JSON для внешней системы отчётов.

## Установка

```bash
conda env create -f environment.yml
conda activate nipt_mngs
# требуются также: nextflow, samtools на PATH
```

## Запуск

> **Важно:** команды в этом разделе — **шаблоны для реальных данных**.
> `samplesheet_bam.csv` / `samplesheet_fastq.csv` и пути `/path/to/...`
> нужно заменить на свои файлы (в корне проекта этих файлов нет).
> Готовая к запуску «из коробки» команда — в разделе
> [Smoke-test на toy-данных](#smoke-test-на-toy-данных) ниже.

### Основной режим: от DEHOST/HOST BAM NIPD-системы

Samplesheet (`samplesheet_bam.csv`):

```csv
sample,dehost_bam,host_bam
SAMPLE1,/path/S1.dehost.bam,/path/S1.host.bam
```

```bash
nextflow run main.nf --mode bam \
    --input samplesheet_bam.csv \
    --db_full /path/to/kraken2_full_db \
    --db_patho /path/to/eukpathdb \
    --outdir results -resume
```

### Полный путь от сырых прочтений

Samplesheet (`samplesheet_fastq.csv`):

```csv
sample,fastq_1,fastq_2
SAMPLE1,/path/S1_R1.fastq.gz,/path/S1_R2.fastq.gz
```

```bash
nextflow run main.nf --mode fastq \
    --input samplesheet_fastq.csv \
    --host_index /path/to/hisat2/human_prefix \
    --db_full /path/to/kraken2_full_db \
    --db_patho /path/to/eukpathdb \
    --outdir results -resume
```

Обе БД должны быть собраны Kraken2 + `bracken-build` с длиной прочтения,
равной `params.bracken_read_len` (по умолчанию 150).

## Выходы (`results/`)

```
<sample>/
├── qc/
│   ├── <sample>_qc_report.json      # вердикт + критерии (контракт для внешней системы)
│   ├── <sample>_qc_mqc.tsv          # custom content для MultiQC
│   ├── fastp/                       # fastp JSON/HTML (только mode=fastq)
│   └── kraken/                      # 4 kraken2 report
├── profiles/                        # bracken TSV ×4 + diff TSV/JSON ×2
multiqc/multiqc_report.html          # сводный QC-отчёт
```

Статусы: `PASS_POSITIVE` / `PASS_NEGATIVE` / `FALLBACK` (образец на
отдельную обработку, причины в `fallback_reasons`) / `FAIL`.

## Smoke-test на toy-данных

```bash
# toy-данные уже в test/toy/ (регенерация: test/make_toy_data.py + test/build_toy_db.sh)
# Пути в samplesheet относительные (от директории samplesheet) — команда
# работает из любого каталога запуска и на любой машине.
nextflow run main.nf --mode fastq --input test/toy/samplesheet_fastq.csv \
    --host_index $PWD/test/toy/host_index/human \
    --db_full $PWD/test/toy/db_full --db_patho $PWD/test/toy/db_patho \
    --outdir test/results_fastq -profile test

nextflow run main.nf --mode bam --input test/toy/samplesheet_bam.csv \
    --host_index $PWD/test/toy/host_index/human \
    --db_full $PWD/test/toy/db_full --db_patho $PWD/test/toy/db_patho \
    --outdir test/results_bam -profile test
```

### Запуск на ноутбуке (WSL2, ~4 GB RAM)

Профиль `test` ужимает все задачи до 2 CPU / 2 GB — smoke-test проходит
на машине с 4 GB RAM. Дополнительные советы:

- **Не запускайте из `/mnt/c/...`** (диск Windows в WSL2): файловая
  система 9P очень медленная и ломает блокировки. Скопируйте пайплайн в
  домашний каталог WSL (`~/nipt_mngs_qc`) и запускайте оттуда.
- `NXF_SYNTAX_PARSER=v1` не нужен (это обходной флаг старых версий
  Nextflow; на 26.04.6 работает парсер по умолчанию).
- Если kraken-задача убита по памяти (exit 137), профиль `test` сам
  перезапустит её один раз с повышенным лимитом.

Проверенные результаты (см. `report_pipeline_theory.md`, часть 6):

| Режим | Образец | Вердикт |
|---|---|---|
| fastq | TOY_NORMAL (97% human + 3 микроба) | PASS_POSITIVE |
| fastq | TOY_ZERO (100% human) | FALLBACK: ZERO_MICROBIOME, DEHOST_FAILURE |
| bam | TOY_NORMAL | PASS_POSITIVE |
| bam | TOY_LEAK (dehost = host BAM) | FALLBACK: DEHOST_FAILURE, HOST_FRACTION_LOW |

## Калибровка порогов

Все пороги QC — параметры (`nextflow.config`, секция `params`),
переопределяются из командной строки, например:
`--min_nonhost_reads 100000 --host_frac_min 0.90`. Дефолты — стартовые
эвристики под low-pass НИПТ; калибровка на реальных когортах описана в
`QC_checklist.md`.
