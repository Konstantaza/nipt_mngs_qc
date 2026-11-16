# MetaX с нуля на существующем Samovar toy

## Что делаем

Устанавливаем **MetaX 2.0.3 отдельно** и собираем две маленькие БД из уже
имеющихся референсов. Перезапускаем `nipt_mngs_qc` v1.10.0 на исходных
PE126/SE50 без новой генерации, изменения dehost или подгонки QC-порогов.

**Это технический тест интеграции, не анализ реальных образцов и не
production-БД.** Переход к реальным данным пока не выполнен.
В MetaX 2 требуются новые индексы; Kraken2-БД не конвертируются простым
переименованием. Используем полные индексы **без `-f`**. [1, 2]

### Ваш сервер и входы

В присланной диагностике: 125 GiB RAM всего / 110 GiB доступно, 32 CPU,
около 25 T свободного места на сетевом томе; Python 3.14.7,
Nextflow 26.10.5, Java 21. MetaX на PATH отсутствовал.
Это снимок ресурсов, не их эксклюзивное резервирование.

Нужные файлы найдены:

| Назначение | Файл относительно репозитория |
|---|---|
| Human FASTA | `insilico_runs/toy_995/refs/9606.fna` |
| PhiX FASTA | `../samovar/data/test_genomes/meta/2886930.fna` |
| Yeast FASTA | `../samovar/data/test_genomes/meta/4932.fna` |
| E. coli FASTA | `../samovar/data/test_genomes/meta/562.fna` |
| PE | `insilico_runs/toy_995/reads/toy_PE_R1.fastq`, `toy_PE_R2.fastq` |
| SE | `insilico_runs/toy_995/reads/toy_SE50.fastq` |
| Manifests | `insilico_runs/toy_995/generation.json`, `build.json` |
| Samplesheets | `insilico_runs/toy_995/samplesheet_pe.csv`, `samplesheet_se.csv` |
| Taxonomy | `tests/insilico/toy/taxonomy/nodes.dmp`, `names.dmp` |

Четыре FASTA занимают суммарно **11 562 082 bytes**. Размер файлов с
заголовками и переносами строк — не число нуклеотидов.

## 1. Разместить комплект

Скачайте `metax_samovar_server_v1.tar.gz` из ответа и перенесите его **в
корень серверного репозитория** удобным вам способом. Далее все команды
выполняйте на сервере в одной сессии. Git commit/push не требуется.

```bash
cd /mnt/tank/scratch/kyamschikov/nipt_mngs_qc || exit 1
conda activate nipt_mngs_qc

# Не извлекать поверх ранее распакованного комплекта.
test ! -e metax_samovar_server_v1 &&
tar -xzf metax_samovar_server_v1.tar.gz
```

После распаковки появится каталог `metax_samovar_server_v1`:
`install_metax.sh`, `server.py`, `test_server.py` и эта инструкция.
Core pipeline не меняется. Все следующие команды выполняйте по очереди;
при `STOP` или ненулевом exit code не переходите к следующему этапу.

## 2. Установить MetaX отдельно

```bash
bash metax_samovar_server_v1/install_metax.sh
```

Будет создано отдельное окружение:
`/mnt/tank/scratch/kyamschikov/metax_env_2_0_3`.
Скрипт устанавливает `metax=2.0.3`, `python=3.14.7` из `zldeng` /
`conda-forge`, записывает `install.log` и `metax_environment.explicit.txt`.
Существующий каталог окружения не перезаписывается: если он уже есть,
проверяется наличие именно MetaX 2.0.3.

**Не выполняйте `conda activate` для нового MetaX-окружения.**
Основной Python и остальные инструменты остаются из `nipt_mngs_qc`.
При pipeline-run на PATH добавляется только symlink бинарника MetaX.

Если установка не удалась, пришлите текст ошибки и
`metax_samovar_server_v1/install.log`, если он появился. Не устанавливайте
другую версию MetaX и не меняйте основное окружение для обхода ошибки.

## 3. Проверить входы и подготовить FASTA/config

```bash
python metax_samovar_server_v1/server.py prepare
```

Скрипт проверит до запуска workflow:

- наличие основных инструментов и MetaX 2.0.3;
- manifest v1.10.0 и присутствие закреплённого MetaX-адаптера;
- hashes исходных PE/SE FASTQ и human FASTA по `generation.json`;
- исходные samplesheets и режим PE/SE, без изменения их содержимого;
- реальный HISAT2 prefix из `build.json`, все восемь частей и их hashes;
- совпадение восстановленных `hisat2-inspect` последовательностей с human
  FASTA — не только наличие файлов;
- toy taxonomy: 43 узла, корректные ancestry и species labels;
- исходную target panel и безопасные для shell пути/имена.

Новый рабочий каталог:
`/mnt/tank/scratch/kyamschikov/nipt_mngs_qc/insilico_runs/metax_samovar_v1`.
Если он уже существует, `prepare` остановится, а не перезапишет его.
При ошибке не удаляйте каталог вслепую — сначала сохраните диагностику.

Новые FASTA содержат те же последовательности с заголовками
`genome_id|contig_id|s:species_taxid`. Все contigs одного исходного FASTA
имеют общий genome ID, отдельные contig IDs и сохранённые исходные
заголовки / sequence hashes в `reference_manifest.json`. Это соответствует
формату MetaX 2. [2]

| БД | Состав |
|---|---|
| `db_full` | Human 9606 + phiX 10847 + yeast 4932 + E. coli 562 |
| `db_patho` | Только phiX, yeast и E. coli |
| `db_pilot` | Только исходный E. coli FASTA, для измерения индексации |

Для phiX используется явное соответствие source taxid `2886930`
(`no rank`) → species taxid `10847` из исходного toy taxdump.
Таксономические записи не выдумываются. `merged.dmp` и `delnodes.dmp`
пустые для этого неизменённого тестового набора, не замена полному
production taxdump.

**Граница контроля происхождения:** старый `generation.json` содержит
hashes reads и human, но не исходных трёх microbial FASTA. Их текущие
hashes фиксируются сейчас; это само по себе не доказывает, что файлы
не менялись после старой генерации. Поэтому этот этап проверяет прежде
всего интеграцию, а не даёт новый независимый benchmark точности.

### QC и ресурсы

Подключается ваш исходный `tests/insilico/toy/toy.config`:

- minimum input: **100 000** пар / SE reads;
- minimum nonhost: **1 000**;
- minimum host fraction: **0,993**;
- minimum Q30: **0,80**;
- minimum trim survival: **0,90**;
- maximum human fraction в dehost: **0,05**;
- minimum classified fraction: **0,0001**;
- `compare_min_reads=10`, `compare_host_min_reads=10`.

Остальные QC-параметры наследуются от v1.10.0 без изменения, включая
`max_consistency_dev=0.01`. Сохраняются исходная target panel, PE126/SE50,
`publish_validation_artifacts=true` и `kraken_mem='2 GB'` из toy settings.
Параметр Bracken read length сохранён для соответствия исходной команде,
но **Bracken в MetaX-ветке не выполняется**.

Переопределены только инфраструктурные настройки: **4 потока, 8 GB на
задачу, очередь local executor = 1**, отдельные work/outdirs; явно задан
`resourceLimits`, чтобы toy `kraken_mem` не урезал MetaX. Это запросы /
потолки Nextflow, не cgroup-ограничение потребления на Linux.
Conda/Docker/Singularity execution отключены: используем инструменты PATH.
Профиль `test` **не используется**.

Развёрнутые значения сохраняются в `logs/nextflow_resolved.config`.

## 4. Собрать две БД

```bash
python metax_samovar_server_v1/server.py index
```

Порядок:

1. Индекс E. coli pilot с записью wall time и peak RSS.
2. Оценка масштаба full-индекса; при прогнозе выше 5 GiB — остановка,
   без автоматического изменения ресурсов.
3. Последовательная сборка `db_full`, затем `db_patho`, **без `-f`**.
4. Копирование taxonomy, проверка `#fractional false` в gzip metadata,
   `metax info` для каждого индекса.
5. Сохранение `build_metax.json`, hashes и измерений.

Оценка pilot — грубый screening с запасом 1,5×, не гарантия линейного
масштабирования. Она относится к **индексации**, а не к последующей
памяти профилирования большого набора reads.

Ожидаются `BUILT: full and patho`, `dbsheet_metax.csv`,
`pilot_estimate.json`, `build_metax.json` и логи. Повторная сборка поверх
этих каталогов запрещена. Kraken2 `toy_995/db_full`, `toy_995/db_patho`
и результаты прежнего запуска остаются нетронутыми.

Первое профилирование создаёт рядом с индексом
`metax_db.mxi.metamcache`; новые БД должны оставаться доступными на запись.
`metax info` не прогревает этот cache и не является проверкой полноценного
профилирования. [1, 2]

## 5. Запустить PE, затем SE

```bash
python metax_samovar_server_v1/server.py run pe
```

После завершения PE:

```bash
python metax_samovar_server_v1/server.py run se
```

Не запускайте обе команды параллельно для первого теста.
Каждая повторно проверяет hashes входов, core pipeline и построенных БД.
Используются **исходные** `toy_995/samplesheet_pe.csv` /
`samplesheet_se.csv`; фактический `--host_index` берётся из `build.json`.
Нет новых reads, GT-фильтрации, исправления dehost или замены входов.

Новые outputs:

```text
insilico_runs/metax_samovar_v1/
  prepare.json
  reference_manifest.json
  dbsheet_metax.csv
  execution.config
  pilot_estimate.json
  build_metax.json
  db_pilot/
  db_full/
  db_patho/
  logs/
    pe.nextflow.log
    pipeline_pe.log
    pe.trace.tsv
    se.nextflow.log
    pipeline_se.log
    se.trace.tsv
  results_pe/
  results_se/
  work_pe/
  work_se/
  run_pe.json
  run_se.json
```

`run_pe.json` / `run_se.json` появляются только после exit 0 и проверки
четырёх profile JSON: host/dehost × full/patho. Их смысл — **workflow
выполнился**, не «QC прошёл» и не «профиль точен».

В каждом `results_pe` / `results_se` будут:

- `<sample>/profiles/`: native `*.metax.txt`, contract `*.profile.tsv`,
  `*.profile.json`, `*.compare.tsv`, `*.compare.json`;
- `<sample>/qc/`: `*_qc_report.json`, counts JSON, flagstat и QC-таблицы;
- `multiqc/multiqc_report.html` и data directory;
- дополнительные FASTQ/BAM validation artifacts исходного toy-режима.

Точные имена sample берутся из ваших samplesheets, не подменяются.

### Как читать результат

1. **Исполнение:** exit code, `run_*.json`, отсутствие failed tasks.
2. **QC:** содержимое `*_qc_report.json`, включая FALLBACK и причины.
3. **Профили:** таксоны, количества, human в dehost, compare full/patho.
4. **Сравнение с Kraken2:** только после разбора единиц и семантики
   reported counts, не как автоматическая read-level confusion matrix.

MetaX-adapter использует post-filter reported counts; это не тот же
показатель, что raw mapping rate или Kraken root counts. Native abundance
MetaX и контрактная доля `reads/input` также различаются. Старый
Kraken-specific `toy.py evaluate` здесь **не запускается**.

Название роли `patho` означает ветку workflow, **не доказанную патогенность**
трёх тестовых организмов. Адаптер добавляет `--host 9606`, но не передаёт
pathogen table; host/disease поля без неё могут быть `unknown`.
Host-mode также отличается поведением fallback от full-role.
Клинических выводов из этих меток делать нельзя. [3]

Не ожидайте заранее PASS, особенно для SE50: замена аннотатора не
исправляет выполненный ранее этап удаления human. Не меняйте пороги
после просмотра результата ради зелёного отчёта.

## 6. Что прислать мне

После PE/SE выполните:

```bash
python metax_samovar_server_v1/server.py collect
```

Скрипт напечатает **`SEND THIS FILE:`**, точное имя и SHA256 архива
`metax_samovar_results_<UTC timestamp>_<suffix>.tar.gz` в новом рабочем каталоге.
Загрузите этот **один архив** в чат.

Внутри: preparation/build manifests, versions, conda explicit list,
pilot/index measurements, resolved config, dbsheet, Nextflow logs/trace,
native/contract profiles, compare, counts, QC JSON, MultiQC и небольшие
task command/error logs; `collection.json` перечисляет выбранные файлы.
FASTQ, FASTA, BAM, индексы, cache и весь work в архив не включаются.

Краткий контроль без ручного поиска имён sample:

```bash
ls -lh insilico_runs/metax_samovar_v1/run_{pe,se}.json
find insilico_runs/metax_samovar_v1/results_pe \
     insilico_runs/metax_samovar_v1/results_se \
     -name '*_qc_report.json' -print -exec cat {} \;
```

Скопируйте также последние строки терминала PE/SE. Не присылайте только
скриншот MultiQC: он не заменяет native profiles, QC JSON и логи.

## Если что-то упало

Не перезапускайте `prepare` / `index` поверх существующих каталогов
и не удаляйте work. Если новый рабочий каталог уже появился:

```bash
python metax_samovar_server_v1/server.py collect
```

Пришлите созданный архив и текст `STOP`. Для pipeline-ошибки посмотрите:

```bash
tail -n 100 insilico_runs/metax_samovar_v1/logs/pipeline_pe.log
tail -n 100 insilico_runs/metax_samovar_v1/logs/pipeline_se.log
```

Лог отсутствующего / ещё не начатого запуска может не существовать.
Если `prepare` остановился до создания каталога, пришлите полный вывод
терминала; сборщику пока нечего архивировать.
При ошибке индексации пришлите также соответствующий `logs/index_*.log`
через общий архив. Частичная БД не достраивается автоматически.

Только после разбора причины и её устранения допустим явный resume
конкретного неуспешного pipeline-run:

```bash
python metax_samovar_server_v1/server.py run pe --resume
# Или, если неуспешным был SE:
python metax_samovar_server_v1/server.py run se --resume
```

Успешный run защищён marker-файлом от случайной перезаписи.
PE/SE используют отдельные launch/work каталоги, чтобы не смешивать
Nextflow history при resume. Если входы или pipeline изменились,
hash check остановит run — это не случай для обхода проверки.

## Что уже проверено при подготовке комплекта

**Локально, не на вашем сервере:**

- shell syntax installer и Python syntax;
- 14 software-тестов: FASTA integrity, taxonomy, full-index rejection,
  безопасные пути, команды PE/SE, запрет overwrite, measurement и
  исключение больших входов из диагностического архива, PATH isolation
  и success/failure/resume launcher;
- preflight на существующих Samovar-файлах, включая hashes и совпадение
  human FASTA / HISAT2 sequences;
- чтение merged Nextflow config: исходные QC settings, 4 CPU / 8 GB,
  очередь 1;
- реальные сборки MetaX 2.0.3: pilot, full, patho;
- native PE/SE × full/patho на **первых 4096 исходных записей** каждого
  FASTQ по порядку, без отбора по GT; все четыре профиля прошли штатные
  конвертер и валидатор pipeline;
- сбор и просмотр состава диагностического архива.

Измерения индексации этой локальной копии, не прогноз для сервера:

| Индекс | Wall time, s | Peak RSS, MiB |
|---|---:|---:|
| Pilot E. coli | 0,244 | 34,16 |
| Full | 0,472 | 60,82 |
| Patho | 0,366 | 51,75 |

В full FASTA **11 418 456 нуклеотидов**. Время зависит от CPU, файловой
системы и cache; на сервере будут свои записи `*.resources.json`.
Smoke на 4096 reads — проверка программной совместимости, не оценка
чувствительности, специфичности, QC или стабильности на 600 000 парах.

**Не выполнены здесь:** новая conda-установка на сервере, полный
Nextflow E2E этим комплектом на всех 600 000 PE / 600 000 SE и обработка
реальных образцов. Не следует считать локальные smoke результаты вашими
серверными PE/SE результатами.

## Следующий этап — реальные данные

После разбора присланного архива отдельно определим:

1. Какие реальные FASTQ или обе фракции BAM доступны и разрешены к анализу.
2. Их происхождение, human reference, read layout, контроли и metadata.
3. Состав production MetaX-БД и согласованный taxonomy snapshot.
4. Pilot, ресурсные требования и отдельные outdirs с production QC defaults.

Малую toy-БД не переносить в клиническую интерпретацию. Публичные данные
или большие production downloads без отдельного решения не подставляются.

## References

1. GitHub - hzi-bifo/Metax at v2.0.3 · GitHub. [cited 2026 Nov 9]. Available from: https://github.com/hzi-bifo/Metax/tree/v2.0.3
2. MANUAL.md. [cited 2026 Nov 9]. Available from: https://github.com/hzi-bifo/Metax/blob/v2.0.3/MANUAL.md
3. https://raw.githubusercontent.com/hzi-bifo/Metax/v2.0.3/src/profiler/community.rs. [cited 2026 Nov 9]. Available from: https://raw.githubusercontent.com/hzi-bifo/Metax/v2.0.3/src/profiler/community.rs
