# nipt_mngs_qc

Метагеномный анализ поверх данных НИПТ с формальным per-sample QC gate.
Nextflow DSL2. Текущая версия: **v1.10.0**.

**Новый in-silico тест Samovar 99,5% HOST:** `tests/insilico/toy/README.md`
содержит точные серверные команды для соседних репозиториев, PE126/SE50,
фиксированный seed и автоматическую приёмку. Фактический локальный итог:
**FAIL строгой приёмки**, а не «всё прошло»: оба workflow завершились,
HOST и микробный состав в допуске, но human остаётся в DEHOST;
SE получил `FALLBACK/DEHOST_FAILURE`. См. сохранённый отчёт теста.

**Документация:**
- `report_pipeline_theory.md` — теория инструментов, устройство пайплайна,
  как читать/расширять (образовательный отчёт);
- `QC_checklist.md` — формальные критерии QC/адекватности, пороги,
  fallback-правила, контракт JSON для внешней системы отчётов;
- `docs/profile_contract.md` — контракт таксономического профиля
  (profile v2) и интерфейс subworkflow ANNOTATE: формат TSV + sidecar
  JSON, правила пустого профиля, как добавить новый аннотатор.

## Установка

```bash
conda env create -f environment.yml
conda activate nipt_mngs_qc
# требуется также nextflow >= 26.04.6 (задано в manifest.nextflowVersion)
```

С v1.8.0 все инструменты в `environment.yml` закреплены точной версией
(fastp 1.3.7, kraken2 2.17.2, bracken 3.1, samtools 1.24, hisat2 2.2.3,
multiqc 1.35, python 3.14.7, pytest 9.1.1; с v1.10.0 также seqkit 2.9.0) — окружение воспроизводимо
(NIPT-7). hocort не используется и не ставится.

### Профили окружения (v1.8.0)

У каждого процесса — свои директивы `conda` и `container` с точными
версиями (см. `modules/local/*.nf`):

- `-profile conda` — per-process conda-окружения (Nextflow собирает по
  одному env на уникальную директиву);
- `-profile docker` — per-process контейнеры: одно-инструментные
  процессы тянут точные образы с quay.io/biocontainers, а три
  мульти-инструментальных (DEHOST, BAM2FASTA, ANNOTATE_KRAKEN2_BRACKEN)
  используют all-in-one образ `nipt_mngs_qc:1.10.0`
  (`params.pipeline_image`) — готовых mulled-образов с точными
  комбинациями версий на quay.io нет. Перед запуском с `docker` образ
  нужно собрать: `docker build -t nipt_mngs_qc:1.10.0 .`
  (`test/run_smoke.sh` делает это автоматически);
- без профиля (`standard`) — инструменты с PATH, как раньше.

**MetaX — опционально (v1.9.0, NIPT-16).** Реестр содержит три
адаптера: `kraken2_bracken`, fixture-driven `stub` (тесты без БД) и
настоящий `metax` **2.0.3**. Последнему нужны отдельное окружение
`zldeng::metax=2.0.3` + Python и индекс **MetaX 2**, а не старый
`metax_db.json`. У MetaX 2 встроенный маппер, внешний MA не нужен.
Per-process conda-директива задаёт оба пакета; стандартный Docker-образ
MetaX не содержит. Требования, сборка БД/варианта контейнера, точный
формат и происхождение fixtures: `docs/profile_contract.md`, §8.

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
    --dbsheet dbsheet.csv \
    --target_panel /path/to/target_panel.txt \   # опционально
    --outdir results -resume
```

Базы аннотаторов задаются **dbsheet** (CSV с заголовком
`db_id,annotator,path,role`; с v1.4.0 заменяет `--db_full/--db_patho`):

```csv
db_id,annotator,path,role
kraken2_full_db,kraken2_bracken,/path/to/kraken2_full_db,full
eukpathdb,kraken2_bracken,/path/to/eukpathdb,patho
```

`role` (`full`/`patho`) — ключ downstream-логики (роли уникальны);
`annotator` — идентификатор адаптера из реестра subworkflow ANNOTATE
(сейчас `kraken2_bracken`, `stub`, `metax`);
`path` — абсолютный, URI или относительный от директории dbsheet
(см. «Валидация входов и пути в CSV»). Смена
классификатора = правка dbsheet, без изменений кода потребителей
(см. `docs/profile_contract.md`).

`--target_panel` — файл мишеней (taxid или имя вида, по одному на строку).
Позитивность (`PASS_POSITIVE`) определяется только пересечением
обогащённых таксонов patho-БД с этой панелью. **Без панели позитивность не
оценивается**: QC-пройденные образцы получают `PASS_NEGATIVE` с
`positivity_source: "not_evaluated"` в QC JSON.

### Полный путь от сырых прочтений

Samplesheet (`samplesheet_fastq.csv`):

```csv
sample,fastq_1,fastq_2
SAMPLE1,/path/S1_R1.fastq.gz,/path/S1_R2.fastq.gz
```

С v1.10.0 `mode=fastq` поддерживает настоящий SE и смешанный samplesheet:

```csv
sample,fastq_1,fastq_2,single_end
PE_SAMPLE,/path/PE_R1.fastq.gz,/path/PE_R2.fastq.gz,false
SE_SAMPLE,/path/SE50.fastq.gz,,true
```

`single_end` необязателен: отсутствует/пустой/false → прежний PE с
обязательным R2; true → SE с пустым `fastq_2`. Другие значения и true
с непустым R2 — ошибка до старта задач. Колонка `fastq_2` остаётся
в заголовке даже для SE. FASTP использует один вход/выход, HISAT2 — `-U`;
mapped SE → HOST, unmapped SE → DEHOST, half остаётся пустым.
Вход BAM и интерфейс ANNOTATE не изменены.

```bash
nextflow run main.nf --mode fastq \
    --input samplesheet_fastq.csv \
    --host_index /path/to/hisat2/human_prefix \
    --dbsheet dbsheet.csv \
    --outdir results -resume
```

Dehost (с v1.5.0) — отдельная стадия DEHOST: **один** проход HISAT2 с
`--no-spliced-alignment` (вход — ДНК) + разбиение primary-записей по
SAM-флагам (`bin/dehost_split.sh`) на три фракции: `host` (оба мейта
картированы), `dehost` (оба некартированы), `half` (полукартированные —
ровно один мейт; учитываются явно в QC, но не аннотируются). Три фракции
— точная партиция trimmed-прочтений: QC-критерий `consistency` проверяет
`dehost + host + half == trimmed`. hocort удалён (раньше картировал
вход дважды и молча терял полукартированные пары).

Для адаптера `kraken2_bracken` обе БД должны быть собраны Kraken2 +
`bracken-build` с длиной прочтения, равной `params.bracken_read_len`
(по умолчанию 150).

### Валидация входов и пути в CSV (v1.6.0)

Пайплайн падает **до старта задач** с понятным сообщением, если:

- не задан обязательный параметр режима: `--input` и `--dbsheet` всегда,
  `--host_index` при `--mode fastq`;
- в samplesheet нет обязательных колонок режима (`sample,fastq_1,fastq_2`
  или `sample,dehost_bam,host_bam`), есть пустые значения в них, дубли
  `sample` или нет строк данных;
- в dbsheet нет колонок `db_id,annotator,path,role`, есть пустые
  значения, дубли `role` или нет строк данных;
- по префиксу `--host_index` не найдено файлов `<prefix>.*.ht2`.

С v1.7.0 профиль `test` самодостаточен (`input`, `dbsheet`,
`host_index`, `target_panel` заданы внутри профиля), поэтому fail-fast
на отсутствие параметра срабатывает только вне профиля `test`:
`-profile test` без `--input` — штатный toy-прогон, а не ошибка.

Пути внутри samplesheet и dbsheet резолвятся так: **абсолютные**
(`/data/S1.bam`) и **URI** (`s3://...`, `https://...`, `file://...`) —
как есть; **относительные** — от каталога самого CSV (а не от каталога
запуска), поэтому запуск работает из любой директории. Отсутствующий
файл — ошибка `checkIfExists` с явным путём.

HISAT2-индекс стейджится glob'ом `<prefix>.*.ht2`: посторонние файлы
рядом с префиксом (`human.fa`, `human.ht2.bak`, `humanity.1.ht2` и т.п.)
в задачу не попадают.

### Фракции аннотации (`--fractions`, v1.5.0)

Один и тот же subworkflow ANNOTATE запускается на всех фракциях из
`params.fractions` (по умолчанию `'dehost,host'`). Опциональная ветка
`all` — все прочтения без dehost — включается одним параметром, без
правок кода:

```bash
nextflow run main.nf --mode bam --input samplesheet_bam.csv \
    --dbsheet dbsheet.csv --fractions dehost,host,all \
    --outdir results -resume
```

`all` синтезируется конкатенацией FASTA конституентов (fastq:
dehost+host+half = trimmed; bam: dehost+host) — новый прогон
выравнивателя не нужен. Профили `all` публикуются рядом с
`dehost`/`host` (`<sample>.all.<role>.profile.tsv/.json`) и проходят ту
же VALIDATE_PROFILE; QC gate их игнорирует (вердикт не меняется).

## Выходы (`results/`)

```
<sample>/
├── qc/
│   ├── <sample>_qc_report.json      # вердикт + критерии (контракт для внешней системы)
│   ├── <sample>_qc_mqc.tsv          # custom content для MultiQC
│   ├── fastp/                       # fastp JSON/HTML (только mode=fastq)
│   ├── dehost/                      # HISAT2 summary DEHOST (только mode=fastq, v1.5)
│   ├── flagstat/                    # samtools flagstat по фракциям (вход MultiQC)
│   ├── counts/                      # counts JSON по фракциям (dehost, host, +half в
│   │                                #   mode=fastq, +all при --fractions ...all; v1.5)
│   └── kraken/                      # kraken2 report ×(фракции × БД) (native адаптера)
├── profiles/                        # контракт profile v2: profile.tsv + sidecar
│                                    #   profile.json — по паре на фракцию × роль БД
│                                    #   (dehost|host [|all] × full|patho)
│                                    # + compare TSV/JSON ×2 (dehost:host по каждой роли БД)
multiqc/multiqc_report.html          # сводный QC-отчёт (только из staged-входов)
```

Статусы: `PASS_POSITIVE` / `PASS_NEGATIVE` / `FALLBACK` (образец на
отдельную обработку, причины в `fallback_reasons`) / `FAIL`.
Отсутствие обязательного артефакта (пара profile.tsv/.json, compare JSON,
fastp в mode=fastq, counts JSON в mode=bam) — `FAIL` с кодом
`PIPELINE_ERROR`: образец никогда не «пропадает» без статуса и не
получает молчаливый `PASS_NEGATIVE`.

## Smoke-test на toy-данных

### Одна команда (v1.7.0)

```bash
bash test/run_smoke.sh
```

Сквозной harness: прогоняет pytest (**60 юнит-тестов**), затем четыре
end-to-end прогона — mode=bam и mode=fastq с Kraken2-dbsheet,
mode=bam со stub и смешанный **Kraken2(full)+stub(patho)** — и сверяет результаты с
ожиданиями (`test/check_smoke.py`): вердикты QC (TOY_NORMAL —
PASS_POSITIVE; TOY_ZERO — FALLBACK/ZERO_MICROBIOME; TOY_LEAK —
FALLBACK/DEHOST_FAILURE), различие наборов таксонов dehost:host у
нормальных образцов и их совпадение у TOY_LEAK (формальная детекция
сбоя dehost), наличие обязательных артефактов. Любое несовпадение —
понятное сообщение в stderr и exit 1; при успехе — exit 0.
В pure-stub и mixed-прогонах нормальные образцы имеют **PASS_NEGATIVE**:
позитивность определяется enriched-таксонами роли patho, где stub
выдаёт taxid 2097 вне панели. В mixed full-сравнение и детекция
DEHOST_FAILURE у TOY_LEAK сохраняются (full остаётся Kraken2).
Реальный MetaX проверен **отдельно** и не входит в стандартную
матрицу harness (не требует личного conda-канала для обычного smoke).

Профиль `test` самодостаточен: `input`, `dbsheet`, `host_index` и
`target_panel` заданы внутри профиля, пути в CSV относительные —
ничего, кроме `--mode`, передавать не нужно.

Два рантайма (автовыбор): если есть docker — собирается образ
`nipt_mngs_qc:1.10.0` (Dockerfile в корне, зависимости из закреплённого
`environment.yml`) и прогоны идут с `-profile test,docker`: три
мульти-инструментальных процесса — в этом образе, одно-инструментальные
тянут точные biocontainers-образы с quay.io (per-process контейнеры,
v1.8.0); иначе — локально, если на PATH есть fastp, hisat2, samtools,
kraken2, bracken, multiqc (например, из conda-окружения
`environment.yml`). Если nextflow не найден, скрипт скачивает его сам
(нужна java). Локальный прогон на toy-данных — около минуты; на чистой
машине с docker добавляется одноразовая сборка образа и скачивание
biocontainers-образов. Docker-рантайм здесь не проверялся, поэтому
время первого Docker-запуска не гарантируется.

### Ручной запуск отдельных режимов

Те же прогоны по отдельности, с явными путями (эквивалентны дефолтам
профиля `test`):

```bash
# toy-данные уже в test/toy/ (регенерация: test/make_toy_data.py + test/build_toy_db.sh;
# toy-BAM для mode=bam, включая координатно-сортированный TOY_SORTED и
# single-end TOY_SE: test/make_test_bams.py)
# Пути в samplesheet относительные (от директории samplesheet; с v1.6.0
# также абсолютные и URI) — команда работает из любого каталога запуска
# и на любой машине.
nextflow run main.nf --mode fastq --input test/toy/samplesheet_fastq.csv \
    --host_index $PWD/test/toy/host_index/human \
    --dbsheet test/toy/dbsheet.csv \
    --target_panel $PWD/test/toy/target_panel.txt \
    --outdir test/results_fastq -profile test

nextflow run main.nf --mode bam --input test/toy/samplesheet_bam.csv \
    --dbsheet test/toy/dbsheet.csv \
    --target_panel $PWD/test/toy/target_panel.txt \
    --outdir test/results_bam -profile test

# Прогон с fixture-driven stub (проверка расширяемости ANNOTATE:
# смена классификатора = другая строка dbsheet, код не меняется)
nextflow run main.nf --mode bam --input test/toy/samplesheet_bam.csv \
    --dbsheet test/toy/dbsheet_stub.csv \
    --target_panel $PWD/test/toy/target_panel.txt \
    --outdir test/results_stub -profile test

# Два аннотатора в одном прогоне: Kraken2 full + stub patho
nextflow run main.nf -profile test --mode bam \
    --dbsheet test/toy/dbsheet_mixed.csv --outdir test/results_mixed

# Настоящий MetaX 2.0.3 (отдельный запуск; бинарник metax на PATH,
# либо добавить профиль conda: -profile test,conda)
nextflow run main.nf -profile test --mode bam \
    --dbsheet test/toy/dbsheet_metax_real.csv --outdir test/results_metax_real

# Ветка all (все прочтения без dehost) — одним параметром (v1.5.0)
nextflow run main.nf --mode fastq --input test/toy/samplesheet_fastq.csv \
    --host_index $PWD/test/toy/host_index/human \
    --dbsheet test/toy/dbsheet.csv --fractions dehost,host,all \
    --target_panel $PWD/test/toy/target_panel.txt \
    --outdir test/results_fastq_all -profile test
```

Юнит-тесты (60 шт.: QC gate — контрактные гарантии, гейтинг утечки
human, логика панели, consistency с half-фракцией; comparison model v2;
счётчики BAM, включая инвариантность к secondary/supplementary записям;
партиция DEHOST по SAM-флагам; слияние counts JSON ветки all;
валидатор контракта profile v2; конвертер адаптера kraken2_bracken,
включая пустой профиль; fixture-driven stub; MetaX на реальных
13/16-колоночных fixtures, NA/empty/ошибках counts; запрет fractional
индексов и уникализация SE-имён). Входит первым шагом в `test/run_smoke.sh`:

```bash
python3 -m pytest test/ -v
```

### Запуск на ноутбуке (WSL2, ~4 GB RAM)

Профиль `test` ужимает все задачи до 2 CPU / 2 GB — smoke-test проходит
на машине с 4 GB RAM. Дополнительные советы:

- **Не запускайте из `/mnt/c/...`** (диск Windows в WSL2): файловая
  система 9P очень медленная и ломает блокировки. Скопируйте пайплайн в
  домашний каталог WSL (`~/nipt_mngs_qc`) и запускайте оттуда.
- `NXF_SYNTAX_PARSER=v1` не нужен (это обходной флаг старых версий
  Nextflow; на 26.04.6 работает парсер по умолчанию).
- Если kraken-задача убита по памяти (exit 137), она автоматически
  перезапускается один раз с удвоенным лимитом памяти (v1.8.0; в
  профиле `test`: 2 GB → 4 GB). Потолок ресурсов на задачу ограничен
  `process.resourceLimits`.

Проверенные результаты (см. `report_pipeline_theory.md`, часть 6):

| Режим | Образец | Вердикт |
|---|---|---|
| fastq | TOY_NORMAL (97% human + 3 микроба) | PASS_POSITIVE (хит панели: *E. cuniculi*) |
| fastq | TOY_ZERO (100% human) | FALLBACK: ZERO_MICROBIOME |
| bam | TOY_NORMAL | PASS_POSITIVE (хит панели: *E. cuniculi*) |
| bam | TOY_LEAK (dehost = host BAM) | FALLBACK: DEHOST_FAILURE, HOST_FRACTION_LOW |
| bam | TOY_SORTED (координатно-сортированный, +70k secondary) | PASS_POSITIVE (профиль и счётчики ≡ TOY_NORMAL) |
| bam | TOY_SE (single-end) | PASS_POSITIVE (непустой профиль, `single_end=true`) |
| bam + dbsheet_stub | NORMAL/SORTED/SE; LEAK | PASS_NEGATIVE; FALLBACK (HOST_FRACTION_LOW) — fixture-driven stub |
| bam + dbsheet_mixed | NORMAL/SORTED/SE; LEAK | PASS_NEGATIVE; FALLBACK (DEHOST_FAILURE, HOST_FRACTION_LOW) — два аннотатора в одном workflow |
| bam + dbsheet_metax_real | NORMAL/SORTED/SE; LEAK | PASS_POSITIVE; FALLBACK (DEHOST_FAILURE, HOST_FRACTION_LOW) — настоящий MetaX 2.0.3 |

> С v1.1 утечка human в dehost (`DEHOST_FAILURE`) оценивается только при
> достаточном не-хост сигнале (`nonhost_pairs >= min_nonhost_reads`):
> на единицах прочтений доля human бессмысленна. Поэтому TOY_ZERO
> (0 не-хост пар) больше не получает `DEHOST_FAILURE` — только
> `ZERO_MICROBIOME`. TOY_LEAK (92k не-хост пар) по-прежнему ловится.

## Samovar toy: 99,5% HOST, PE126 и SE50

Подробный runbook: `tests/insilico/toy/README.md`. После разовой подготовки
соседнего Samovar и его отдельного Python 3.11:

```bash
bash tests/insilico/toy/run.sh all \
  --samovar-dir ../samovar \
  --outdir insilico_runs/toy_995
```

600000 пар, seed=42, три равные микробные компоненты; настоящий Samovar
generate/prepare/exec, синхронный seqkit shuffle, SE из R1[:50].
По умолчанию host — штатный укороченный toy FASTA Samovar; для полного
человеческого референса нужны явные `--host-fasta` и подходящие ресурсы
(подробности в runbook). `--hisat2-index` проверяется по последовательностям.
Разные распределения Bracken строятся для 126 и 50 nt.

Скрипт пишет `acceptance.tsv/json`; **exit 0 только при всех PASS**.
Exit 1 может означать корректно выявленное несоответствие биологическому
критерию, а не падение Nextflow. Проверенный запуск: HOST PE 99.41113%,
SE 99.40630%; все микробные доли в допуске ±20% относительно 1/3,
но full-профиль DEHOST содержит human: 20 PE-пар / 224 SE-прочтения.
PE QC PASS_POSITIVE, SE QC FALLBACK/DEHOST_FAILURE; оба MultiQC сформированы.
Полный человеческий референс и Docker здесь не проверялись.

`toy.config` меняет только тестовые пороги глубины/HOST; production
defaults сохранены. `--publish_validation_artifacts true` публикует
trimmed FASTQ и BAM партиций для GT-аудита, по умолчанию выключен.
QC-ключ `input_pairs` исторический: для SE он измеряется в прочтениях;
классификатор считает PE-пары, тогда как `primary_reads` — отдельные мейты.
HOST в приёмке = H/(H+D+half), в QC = H/(H+D).

Регрессии:
`python -m pytest test/ tests/insilico/toy/test_toy.py -q` — 87 тестов;
`bash test/run_smoke.sh` — прежние четыре сквозных прогона.

## Калибровка порогов

Все пороги QC — параметры (`nextflow.config`, секция `params`),
переопределяются из командной строки, например:
`--min_nonhost_reads 100000 --host_frac_min 0.90`. Дефолты — стартовые
эвристики под low-pass НИПТ; калибровка на реальных когортах описана в
`QC_checklist.md`.

Параметры comparison model v2: `--compare_pairs 'dehost:host'` (можно
несколько через запятую, напр. `'dehost:host,dehost:all'`),
`--compare_min_reads` (присутствие в dehost), `--compare_host_min_reads`
(отдельный абсолютный порог в host), `--enrichment_fold` (RPM-обогащение
для флага enriched).
