# Контракт таксономического профиля (profile v2)

Единый классификаторо-независимый формат профиля фракции образца и
интерфейс стадии аннотации пайплайна nipt_mngs_qc (v1.4.0; набор значений
`fraction` расширен значением `all` в v1.5.0, сам контракт неизменен).

**Зачем.** До v1.4.0 контракт профиля (profile.v1, один JSON) был спроектирован
вокруг семантики Kraken2/Bracken. profile v2 отделяет контракт от
классификатора: профиль — это пара файлов **TSV + sidecar JSON**, а
Kraken2+Bracken — первый адаптер интерфейса ANNOTATE. Замена или добавление
классификатора не требует правок main.nf и потребителей (comparison model,
QC gate): они читают только этот контракт.

## 1. Файлы профиля

Один прогон аннотации (образец × фракция × БД) производит два файла:

```
<sample>.<fraction>.<role>.profile.tsv     # per-taxon таблица
<sample>.<fraction>.<role>.profile.json    # sidecar со счётчиками и статусом
```

- `sample` — идентификатор образца (samplesheet);
- `fraction` — `dehost` | `host` | `all` (`all` — с v1.5, при
  `--fractions dehost,host,all`; фракция `half` существует только как
  QC-счётчик и никогда не аннотируется, поэтому в профилях не встречается);
- `role` — роль БД из dbsheet (`full` | `patho`; см. ниже). Идентификатор
  самой БД (`db_id`) живёт внутри контракта, а не в имени файла.

Оба файла обязательны. Отсутствие или битость любого из них —
`FAIL`/`PIPELINE_ERROR` в QC gate (критерий «наличие обязательных
артефактов»).

## 2. profile.tsv — per-taxon таблица

TSV с заголовком, одна строка = один таксон. Обязательные колонки
(порядок фиксирован для читаемости, но потребители обязаны читать
**по имени колонки**):

| Колонка | Тип | Смысл |
|---|---|---|
| `sample` | str | идентификатор образца (дублирует sidecar) |
| `fraction` | str | `dehost` \| `host` \| `all` (`all` — с v1.5) |
| `annotator` | str | идентификатор аннотатора из dbsheet (напр. `kraken2_bracken`) |
| `annotator_version` | str | версия(и) инструмента, снятые в рантайме (напр. `2.17.2+3.0.1`) |
| `db_id` | str | идентификатор БД из dbsheet |
| `taxid` | str | NCBI taxid таксона |
| `name` | str | научное имя таксона |
| `rank` | str | человекочитаемый ранг lowercase: `species`, `genus`, `family`, ... (адаптер Kraken2 мапит коды: S→species, G→genus и т.д.) |
| `reads` | float ≥ 0 | **оценка числа прочтений** таксона. Допускается дробное значение: Bracken байесовски переоценивает прочтения и выдаёт float. Классификаторы без переоценки пишут целые. Единицы — прочтения (для SE-фракций) или пары-фрагменты (для PE: Kraken2 `--paired` считает пару одним событием) |
| `abundance` | float, 0..1 | **доля от input_reads фракции** (= reads / input_reads; семантика `fraction_total_reads` Bracken). Это НЕ относительная доля внутри классифицированных |

**Дополнительные колонки `x_*`.** Адаптер может добавлять любые колонки
с префиксом `x_` (напр. `x_kraken_assigned_reads`). Потребители обязаны
их игнорировать; валидатор их не проверяет. Колонки без префикса `x_`,
не входящие в обязательный список, запрещены (валидатор падает) — это
защита от опечаток в обязательных именах.

**Сортировка**: по убыванию `reads` (рекомендация, не проверяется).

## 3. profile.json — sidecar

```json
{
  "contract_version": "profile.v2",
  "sample_id": "TOY_NORMAL",
  "fraction": "dehost",
  "role": "full",
  "db_id": "db_full",
  "annotator": "kraken2_bracken",
  "annotator_version": "2.17.2+3.0.1",
  "input_reads": 2942,
  "classified_reads": 2900,
  "unclassified_reads": 42,
  "status": "ok"
}
```

Обязательные поля:

| Поле | Тип | Смысл |
|---|---|---|
| `contract_version` | str | ровно `"profile.v2"` |
| `input_reads` | int ≥ 0 | прочтения на входе классификации фракции (база RPM-нормировки comparison model) |
| `classified_reads` | int ≥ 0 | прочтения, получившие таксон |
| `unclassified_reads` | int ≥ 0 | прочтения без таксона |
| `status` | str | `ok` \| `empty` (см. §4) |
| `sample_id`, `fraction`, `role`, `db_id`, `annotator`, `annotator_version` | str | поля идентичности — дублируют колонки TSV **намеренно**: у пустого профиля в TSV нет ни одной строки, и без sidecar файл терял бы самоописание |

Инвариант: `input_reads = classified_reads + unclassified_reads`
(проверяется валидатором).

## 4. Пустой профиль

Пустой профиль — **валидное** состояние (edge case «нулевой микробиом»:
прочтения на входе есть, но ни один таксон не прошёл пороги классификатора,
либо прочтений нет вовсе). Правила:

- `status = "empty"` ⇔ TSV содержит только заголовок (0 строк данных);
- `status = "ok"` ⇔ TSV содержит ≥ 1 строку данных;
- при `empty` допустимы и `classified_reads = 0`, и `classified_reads > 0`
  (классификатор что-то отнёс, но ниже порога отчётности — случай Bracken);
- вердикт по пустому профилю выносит QC gate (`ZERO_MICROBIOME` /
  `NO_CLASSIFIED_SIGNAL`), а не падение стадии аннотации.

## 5. Интерфейс ANNOTATE и dbsheet

Subworkflow `subworkflows/local/annotate.nf`:

```
ANNOTATE(reads, dbsheet) -> (meta, profile.tsv, profile.json, native)
```

- `reads` — канал `[meta, r1, r2, se]`, где `meta = [id, fraction,
  single_end]`; для PE-фракций читаются `r1/r2`, для SE (`single_end=true`)
  — `se`;
- `dbsheet` — CSV с колонками **`db_id, annotator, path, role`**
  (параметр `--dbsheet`, заменяет `--db_full/--db_patho`):
  - `db_id` — идентификатор БД (попадает в контракт);
  - `annotator` — идентификатор адаптера (реестр в annotate.nf;
    сейчас: `kraken2_bracken`, `stub`, `metax`);
  - `path` — путь к БД (относительные пути резолвятся от директории
    dbsheet); для `stub` это каталог тестовых fixture, а не индекс БД;
  - `role` — `full` | `patho`; ключ downstream-логики (QC gate требует
    обе роли, comparison model работает per role). Роли уникальны в
    пределах dbsheet;
- выходной `meta` дополняется `db_id, role, annotator`;
- `native` — сырой вывод классификатора (для Kraken2 — report; для
  MetaX — настоящий profile.txt; для stub — текстовый тестовый отчёт).
  Весь `native` миксуется в MultiQC: неродные
  форматы MultiQC молча игнорирует.

Каждый выход адаптера проходит framework-процесс VALIDATE_PROFILE
(`bin/validate_profile.py`): нарушение контракта = падение задачи,
а не молчаливый битый профиль downstream.

**Добавление аннотатора** = новый процесс-адаптер + конвертер в bin/ +
строка в dbsheet + одна ветка в реестре annotate.nf. main.nf и
потребители (compare_profiles.py, qc_gate.py) не меняются.

## 6. Потребители контракта

- `bin/compare_profiles.py` (comparison model v2): читает taxid/name/reads
  из TSV и `input_reads` из sidecar; RPM = reads / input_reads × 1e6.
- `bin/qc_gate.py`: `classified_frac = classified_reads / input_reads`
  (dehost.full); `human_frac` = reads строки `taxid = 9606` в TSV
  dehost.full / input_reads (классификаторо-независимо; если БД не
  содержит human или таксон ниже порога отчётности — human_frac = 0).
- Оба читают колонки по имени через общий reader `bin/profile_io.py` и
  игнорируют `x_*`.

## 7. Отличия от profile.v1 (миграция v1.3.0 → v1.4.0)

| Было (v1) | Стало (v2) |
|---|---|
| один JSON с embedded `taxa` | пара TSV (taxa) + sidecar JSON (счётчики) |
| `db` (full\|patho) + `db_id` | `role` (full\|patho) + `db_id` из dbsheet |
| `annotator` — свободная строка версий | раздельные `annotator` + `annotator_version` |
| `human_reads/human_frac` — clade-pct репорта Kraken2 | вычисляются потребителем из строки TSV taxid 9606 |
| `classified_frac` в контракте | вычисляется из счётчиков sidecar |
| явного статуса пустого профиля не было | `status: ok\|empty` + правила §4 |
| параметры `--db_full/--db_patho` | `--dbsheet` (db_id, annotator, path, role) |

## 8. Второй адаптер: stub и настоящий MetaX (NIPT-16, v1.9.0)

### 8.1 Stub без классификатора и внешней БД

`ANNOTATE_STUB` + `bin/stub_to_profile.py` читают
`<path_из_dbsheet>/<db_id>.tsv` — имя fixture выбирается по `db_id`,
никаких захардкоженных таксонов в конвертере нет:

```tsv
taxid	name	rank	weight
2097	Mycoplasma genitalium	species	1.0
```

`weight` — тестовая доля входных записей FASTA (PE: число записей r1,
SE: число записей se). `reads = input_reads × weight`,
`abundance = weight`, `classified_reads = round(input_reads × sum(weight))`;
остаток — unclassified. Веса в [0,1], сумма ≤ 1; отсутствие/битость
fixture — явная ошибка. Пустой fixture или вход — `status=empty`.
Это исключительно тестовый аннотатор, не биологическая классификация.

Примеры: `test/toy/dbsheet_stub.csv` (обе роли stub),
`test/toy/dbsheet_mixed.csv` (full=Kraken2, patho=stub).
Mixed-прогон — критерий сменяемости: один workflow использует два
аннотатора, потребители читают тот же контракт. Позитивность в нём
**PASS_NEGATIVE**, поскольку панель пересекается с enriched-таксонами
роли **patho**, а не full; TOY_LEAK сохраняет оба fallback по full/counts.

### 8.2 Версия MetaX и отличия CLI от исходной постановки

Адаптер закреплён на **MetaX 2.0.3** (`hzi-bifo/Metax`, tag `v2.0.3`).
Версия реально установлена и прогнана на toy-БД, не эмулируется.
Документация релиза: [README](https://github.com/hzi-bifo/Metax/blob/v2.0.3/README.md)
и [MANUAL](https://github.com/hzi-bifo/Metax/blob/v2.0.3/MANUAL.md).

**Не смешивать с MetaX 0.9:** упомянутые в постановке `metax_db.json`,
заголовок `genome|txid|species_txid|sequence|size` и PE-аргумент `-p`
относятся к старому интерфейсу. В проверенной 2.0.3 индекс —
`metax_db.mxi` + `metax_db.tsv.gz`, PE — `-i R1 -I R2`, маппер
встроен (внешний MA/`maCMD` не требуется). Индекс 0.9 нужно пересобрать;
совместимость со старым JSON-индексом этот адаптер **не заявляет**.

Контракт каталога `path` в dbsheet:

```text
db/
  metax_db.mxi
  metax_db.tsv.gz
  taxonomy/
    nodes.dmp
    names.dmp
    merged.dmp
    delnodes.dmp
```

`metax_prepare.py` до запуска проверяет `#fractional false` в сжатой
метадате: fractional/неопознанный индекс отклоняется (масштабированные
counts несовместимы с reads/input). Для SE скрипт назначает уникальный
ID **каждой** записи FASTA, не меняя последовательности: у BAM-derived
SE две записи могут иметь одинаковый QNAME, и MetaX иначе объединяет
их. PE-вход не переименовывается, pairing задаётся `-i/-I`.

Процесс запускает:

```bash
metax profile --db db --dmp-dir db/taxonomy \
    -i fraction_1.fasta -I fraction_2.fasta -o metax_run -t 2 --no-progress
# SE: только -i metax_single_end.fasta
# role=patho: дополнительно --host 9606
```

Обычный выход: `metax_run.profile.txt`; с `--host 9606`:
`metax_run.pathogen.profile.txt`. Процесс выбирает правильный файл,
конвертирует его и публикует копию как `<sample>.<fraction>.<role>.metax.txt`.
Пустой FASTA сам MetaX отвергает: процесс не запускает классификацию
при 0 записей и создаёт пустой native-файл, затем валидный empty-профиль.
Без pathogen-таблицы host-поля могут быть `unknown`: это **не**
подтверждение ассоциации патоген–человек. Здесь роль patho определяется
выбором отдельной БД, как и у Kraken2.

### 8.3 Конвертация и провенанс реальных fixtures

`bin/metax_to_profile.py --profile ...` читает TSV **без заголовка**,
13 колонок (16 с host-полями):

| Колонка MetaX (1-based) | Поле profile v2 |
|---|---|
| 1 / 2 / 3 | name / taxid / rank |
| 4 | reads (дробное значение сохраняется) |
| 6 | **x_abundance_pct**, не contract abundance |
| 9 / 13 | x_cov_prob / x_chunk_prob (`NA` переносится без изменения) |
| 14 / 15 / 16, если есть | x_host_names / x_host_taxids / x_diseases |

**Осознанное отличие от буквального маппинга col6→abundance:**
col6 — относительное изобилие MetaX в процентах (на toy-файле сумма ≈100),
а контракт требует **reads/input_reads в [0,1]**. Даже деление col6
на 100 не даст долю входных прочтений. Поэтому исходная величина
сохранена без потерь как `x_abundance_pct`, а `abundance` вычисляется.
Пример: B19 имеет 510 прочтений из 2942 (abundance≈0.173351), но
col6=93.88893. Нельзя подменять одно другим.

Счётчики MetaX здесь имеют **post-filter** семантику:
`classified_reads = min(input_reads, ceil(sum(reported reads)))`;
unclassified — остаток (включает отфильтрованные таксоны).
Это не сырой mapping-rate и не эквивалент Kraken2 root-count.
Сумма reads > input (кроме допуска округления нативных %.5f) и
нечисловые/бесконечные reads отвергаются.

Fixtures `test/fixtures/metax_profile_full.txt` и
`metax_profile_pathogen.txt` — **реальные** выходы MetaX 2.0.3
по `TOY_NORMAL.dehost.bam` (2942 пары), не выдуманные строки.
Получены через `samtools collate` → `samtools fasta -1/-2`,
затем команды выше с toy-БД full/patho. Финальный сквозной прогон
`test/results_metax_real/` воспроизвёл оба файла байт-в-байт.
В full: myco=1911, E. cuniculi=520, B19=510; всего classified=2941,
unclassified=1. SE после уникализации имён: input=5884, classified=5882.
Дополнительные тестовые строки `NA` и ошибочных counts явно
синтетические и используются только для проверки парсера.

### 8.4 Сборка БД, окружения и контейнера

Отдельное окружение (стандартный Kraken2-стек не меняется):

```bash
conda create -n metax2 -c zldeng -c conda-forge metax=2.0.3 python=3.14.7
conda activate metax2
```

В процессе задано
`conda 'zldeng::metax=2.0.3 conda-forge::python=3.14.7'`.
Для полного workflow либо используйте `-profile test,conda`, либо
добавьте бинарник metax на PATH основного окружения пайплайна.
Приведённые выше проверки v1.9 выполнены локально с инструментами
на PATH; новая per-process conda-директива MetaX задана, но отдельный
end-to-end прогон именно с ней не выполнялся.
Conda-пакеты MetaX 2 доступны для Linux x86-64/ARM64 и macOS ARM64;
в этом проекте проверен **Linux x86-64**. Исходная сборка:
Rust ≥1.89 + C compiler, tag v2.0.3, `cargo build --release`.

План сборки production-БД (референсы и taxonomy должны быть версионированы):

```text
>genome_id|sequence_id|s:species_taxid
ACGT...
```

```bash
metax index reference.metax.fa -o db_new -t 2   # полный индекс, без -f
mkdir -p db_new/taxonomy
# Скопировать nodes/names/merged/delnodes.dmp из закреплённого NCBI taxdump.
```

В toy-БД референсы уже подготовлены:
`test/toy/db_metax_{full,patho}/reference.metax.fa`; индексы собраны
MetaX 2.0.3, taxonomy — минимальная тестовая (merged/delnodes пусты).
Для пересборки используйте каждый reference.metax.fa с новым `-o`,
затем скопируйте соответствующий taxonomy/. Первый запуск создаёт
`metax_db.mxi.metamcache` рядом с индексом: обеспечьте права записи
при подготовке БД, а для read-only production-монтажа подготовьте
кэш заранее. Runtime-кэши в поставку не включаются.

**Docker:** стандартный `environment.yml` намеренно не тянет личный
канал автора. Для MetaX создайте отдельный вариант environment.yml,
добавив `zldeng` в channels и `zldeng::metax=2.0.3` в dependencies
(Python уже закреплён). Соберите тем же Dockerfile отдельный образ:
`docker build -t nipt_mngs_qc:1.9.0-metax .` и передайте
`--pipeline_image nipt_mngs_qc:1.9.0-metax -profile test,docker`.
Это **план**, Docker в среде проверки отсутствует; сборка этого образа
и его end-to-end запуск не проверены. MA/maCMD добавлять не нужно.
