# QC-чеклист метагеномной части НИПТ (задача 3)

Формальные критерии адекватности на уровне **одного образца**. Реализованы в
`bin/qc_gate.py`; все пороги — параметры `nextflow.config` (калибруются без
правки кода). Вердикты и причины машиночитаемы (JSON) и видны в MultiQC.

## In-silico Samovar 99,5% HOST (v1.10.0)

Источники/команды/проверки — `tests/insilico/toy/README.md`.
Этот тест **не прошёл строгую приёмку**, хотя Nextflow PE/SE успешно
завершился. Не путать exit status workflow с QC образца.

| Метрика/критерий | PE126 | SE50 |
|---|---:|---:|
| HOST = H/(H+D+half), требование 99.3–99.7% | 99.41113%, PASS | 99.40630%, PASS |
| Три микробных species; относительные доли ±20% | PASS | PASS |
| Human в DEHOST full, при min_reads=10 | 20 пар, FAIL | 224 reads, FAIL |
| QC verdict | PASS_POSITIVE | FALLBACK/DEHOST_FAILURE |
| Human classified / DEHOST input | 0.6658% | 6.3709% |
| dehost:host + MultiQC | сформированы | сформированы |

Тестовый `toy.config`: min_input_pairs=100000, min_nonhost_reads=1000,
host_frac_min=0.993; остальные пороги качества остаются production:
Q30≥0.80, survival≥0.90, human leakage≤0.05.
Production defaults не менялись. Три вида в роли patho — software-test
panel, не клиническая панель патогенов. Изменять пороги ради PASS нельзя.

Приёмка дополнительно сверяет ground truth taxid в ID и точное разбиение
trimmed IDs по HOST/DEHOST/half. Наличие human в DEHOST не приписывается
ошибке классификатора автоматически: по GT там 40 человеческих PE-мейтов
и 516 SE-прочтений. В closed toy DB не проверяется специфичность к
широкому набору отсутствующих организмов.

С v1.10.0 SE поддержан и для FASTQ: `single_end=true`, пустой fastq_2,
HISAT2 `-U`, пустая half-фракция. Исторические QC-ключи `input_pairs`
и `nonhost_pairs` означают фрагменты: PE-пары / SE-прочтения.
Для консистентности используются **read records** в обоих режимах.
`publish_validation_artifacts` включает аудитные FASTQ/BAM, иначе
они остаются промежуточными.

## Статусы образца

| Статус | Смысл | Что отдаём наружу |
|---|---|---|
| `PASS_POSITIVE` | QC пройден, обогащённые таксоны patho-БД пересекли **панель мишеней** | Полный отчёт + differential profiles |
| `PASS_NEGATIVE` | QC пройден, позитивности нет (см. `positivity_source`) | Полный отчёт с валидным «нулём» |
| `FALLBACK` | Образец не прошёл QC | **НЕ отдаём результат**; JSON с кодами причин уходит по отдельному пути |
| `FAIL` | Вход неадекватен или отсутствует обязательный артефакт (анализ невозможен/неполон) | JSON с причиной; образец на пересеквенирование/отбраковку |

Ключевой принцип edge-case: **нулевой микробиом / провал анализа — это не
пустой отчёт, а явный статус FALLBACK с причиной**. «Нет гарантии результата»
контролируется на уровне образца, а не батча.

Второй принцип (v1.1): **отсутствие обязательного артефакта — это
`FAIL`/`PIPELINE_ERROR`, а не молчаливый `PASS_NEGATIVE` и не пропажа
образца**. Каждый образец samplesheet получает QC JSON со статусом.

## Чеклист критериев

| # | Критерий | Метрика (источник) | Порог по умолчанию | Severity | Код причины |
|---|---|---|---|---|---|
| 1 | Достаточность входа | фрагменты = fragments_dehost + fragments_host (counts JSON) | ≥ 5 000 000 | FAIL | `LOW_INPUT` |
| 2 | Доля хоста в норме | primary_host/(primary_dehost+primary_host) (counts JSON) | ≥ 0.95 | FALLBACK | `HOST_FRACTION_LOW` |
| 3 | Достаточность не-хост сигнала | fragments dehost (counts JSON) | ≥ 50 000 | FALLBACK | `ZERO_MICROBIOME` |
| 4 | Качество сырых прочтений | Q30 (fastp json) | ≥ 0.80 | FALLBACK | `LOW_QUALITY` |
| 5 | Выживаемость тримминга | after/before (fastp json) | ≥ 0.90 | FALLBACK | `LOW_QUALITY` |
| 6 | Консистентность разбиения | \|dehost+host+half − trimmed\|/trimmed (v1.5: точная партиция DEHOST) | ≤ 0.01 | FALLBACK | `SPLIT_INCONSISTENT` |
| 7 | Утечка human в dehost | доля прочтений taxid 9606 в **profile.tsv** dehost.full (контракт profile v2) | ≤ 0.05 | FALLBACK | `DEHOST_FAILURE` |
| 8 | Санити классификации | `classified_reads / input_reads` из sidecar **profile.json** dehost.full | > 0.0001 | FALLBACK | `NO_CLASSIFIED_SIGNAL` |
| 9 | Наличие обязательных артефактов | 4 пары profile.tsv + profile.json, compare JSON ×2, counts JSON ×2 (mode=bam) / ×3 (+half, mode=fastq), fastp JSON (mode=fastq) | — | FAIL | `PIPELINE_ERROR` |

Критерии 4–6 оцениваются только при `mode=fastq`; при `mode=bam` fastp
выполнен upstream в NIPD и критерии помечаются `pass=null` («не оценивается»).

**Гейтинг критерия 7 (v1.1):** утечка human оценивается только при
достаточном не-хост сигнале (`nonhost_pairs ≥ min_nonhost_reads`). На
единицах прочтений доля human бессмысленна (любой единичный хит даёт
«100%»), поэтому при нехватке сигнала критерий помечается `pass=null`,
а образец и так падает по критерию 3 (`ZERO_MICROBIOME`).

**Источник метрик веток (v1.4.0, контракт profile v2):** все метрики
классификации QC gate получает из контракта profile v2 — пары
`profile.tsv` (per-taxon строки) + sidecar `profile.json` (счётчики
и статус), производимой subworkflow ANNOTATE (`docs/profile_contract.md`).
Чтение — только через общий reader `bin/profile_io.py`; каждый выход
адаптера валидируется процессом VALIDATE_PROFILE
(`bin/validate_profile.py`). QC gate форматов классификатора не знает
и файлы не глобит: все входы — явные пути. Метрики
классификаторо-независимы: human — строка `taxid 9606` в TSV (0, если
БД не содержит human или таксон ниже порога отчётности адаптера),
`classified_frac` — счётчики sidecar.

**Счётчики из counts JSON (v1.3):** в режиме bam счётчики критериев
1–3 и 6 берутся из per-fraction counts JSON
(`<sample>.<fraction>.counts.json`, `bin/bam_counts.py`):
`primary_reads` = `samtools view -c -F 0x900` (без secondary/
supplementary — строка «in total» flagstat их включала и завышала
счётчики); `fragments = paired_reads//2 + unpaired_reads`. Для
single-end BAM (автоопределение: `paired_reads == 0` per fraction)
«пары» = прочтения (фрагменты), в QC JSON добавляется метрика
`single_end`. Flagstat по-прежнему публикуется (`qc/flagstat/`), но
используется только MultiQC.

**Партиция DEHOST и half-фракция (v1.5.0):** в режиме fastq стадия
DEHOST (один проход HISAT2 `--no-spliced-alignment` + разбиение по
SAM-флагам, `bin/dehost_split.sh`) раскладывает primary-записи в три
фракции: `host` (оба мейта картированы), `dehost` (оба некартированы),
`half` (полукартированные — ровно один мейт). Полукартированные пары
**учитываются явно**: counts JSON фракции `half`, информационные
метрики `half_mapped_reads` / `half_mapped_frac` в QC JSON (на вердикт
не влияют) и критерий 6, который теперь проверяет точную партицию
`dehost + host + half == trimmed` (в v1.4 hocort молча терял
полукартированные пары между двумя прогонами). В аннотацию `half` не
идёт (это QC-категория, а не биологическая фракция). В режиме bam
half-фракции нет (партиция NIPD из двух фракций): counts_half = stub →
half = 0, критерий 6 вырождается в форму v1.4 (dehost+host против
trimmed не проверяется — fastp upstream).

**Позитивность (v1.1):** определяется только пересечением обогащённых
таксонов patho-БД с **панелью мишеней** (`--target_panel`, taxid или имя
вида на строку):
- панель задана и есть пересечение → `PASS_POSITIVE`
  (`positivity_source: "target_panel"`);
- панель задана, пересечения нет → `PASS_NEGATIVE`
  (`positivity_source: "target_panel"`);
- панель не задана → позитивность **не оценивается**: `PASS_NEGATIVE`
  с `positivity_source: "not_evaluated"` (строгое поведение: без панели
  нельзя отличить «нет мишени» от «не смотрели»).

Сводка сравнения dehost:host (`comparison_info` в QC JSON) —
**информационная**: на вердикт не влияет, нужна внешней системе отчётов
для контекста.

Критерии отчётности таксонов (не блокирующие, а формирующие содержимое;
comparison model v2, см. `report_pipeline_theory.md` 2.6):
- таксон присутствует в dehost при `reads ≥ compare_min_reads` (10),
  в host — при `reads ≥ compare_host_min_reads` (10, отдельный
  абсолютный порог);
- **enriched (строго)**: присутствует в dehost И отсутствует в host И
  RPM-обогащение ≥ `enrichment_fold` (10×) или ноль прочтений в host.
  Общий таксон (both) никогда не enriched — равные абсолютные числа
  прочтений обогащения не дают (баг v1.0 исправлен);
- вид попадает в отчёт из списка `enriched_taxa` patho-сравнения;
- если обогащённых таксонов в patho-профиле нет — это `PASS_NEGATIVE`
  (валидный отрицательный результат), а не ошибка.

## Обоснование дефолтов

- **≥ 5 млн пар**: low-pass НИПТ обычно даёт 10–13 млн прочтений; 5 млн пар —
  нижняя граница, при которой не-хост остаток (~1–5%) ещё даёт ≥ 50 тыс. пар.
- **≥ 50 тыс. не-хост пар**: при типичной доле не-хоста 1–5% от 5 млн пар это
  50–250 тыс.; ниже — статистика по видам вырождается (LOD Kim et al. 2026:
  сотни копий для CMV при обычной глубине).
- **host_frac ≥ 0.95**: в НИПТ 95–99%+ прочтений — человек; меньшая доля
  означает аномалию библиотеки/экстракции, а не «богатый микробиом».
- **human в dehost ≤ 5%**: при корректном dehost человеческий остаток единичен;
  > 5% — признак сбоя картирования/индекса хоста.
- **enrichment_fold = 10**: host-фракция как встроенный негативный контроль;
  10× — консервативный порог против миссклассификации человеческих повторов.

## Матрица fallback

| Код | Интерпретация | Путь обработки |
|---|---|---|
| `LOW_INPUT` | Мало прочтений | Пересеквенирование с большим покрытием |
| `LOW_QUALITY` | Низкое Q30 / агрессивный тримминг | Пересеквенирование; проверка библиотеки |
| `ZERO_MICROBIOME` | Не-хост сигнал ниже LOD | Не отдавать результат; при клиническом подозрении — выделенный тест с обогащением (ssDNA/сайз-селекция) |
| `HOST_FRACTION_LOW` | Аномально мало хоста | Ручной разбор wet-lab; проверка образца |
| `DEHOST_FAILURE` | Human просачивается в dehost | Проверка индекса/референса хоста; повтор dehost |
| `SPLIT_INCONSISTENT` | dehost+host+half ≠ trimmed | Расследование стадии DEHOST/сплита |
| `NO_CLASSIFIED_SIGNAL` | Kraken2 ничего не классифицировал | Проверка БД; ручной разбор |
| `PIPELINE_ERROR` | Отсутствует/битый обязательный артефакт (пара profile.tsv/.json, compare JSON, counts JSON, fastp) | Перезапуск; эскалация биоинформатику |

## JSON-контракт для внешней системы отчётов

`<sample>_qc_report.json` (v1.10.0):
```json
{
  "sample_id": "SAMPLE_001",
  "pipeline": "nipt_mngs_qc",
  "pipeline_version": "1.10.0",
  "status": "PASS_POSITIVE | PASS_NEGATIVE | FALLBACK | FAIL",
  "fallback_reasons": ["ZERO_MICROBIOME"],
  "positivity_source": "target_panel | not_evaluated",
  "annotator": "kraken2_bracken 2.17.2+3.0.1",
  "db_id": {"full": "kraken2_full_db", "patho": "eukpathdb"},
  "target_panel": {"file": "target_panel.txt", "n_targets": 12,
                   "hits": ["Encephalitozoon cuniculi"],
                   "positivity_source": "target_panel"},
  "criteria": [
    {"id": "input_pairs", "value": 6231451.0, "threshold": ">= 5000000",
     "pass": true, "severity_on_fail": "FAIL", "reason_code": null,
     "description": "Достаточность входа: число пар прочтений"}
  ],
  "metrics": {"input_pairs": 0, "host_fraction": 0, "nonhost_pairs": 0,
              "q30_rate": 0, "dehost_human_frac": 0, "patho_enriched_taxa": 0,
              "patho_panel_hits": 0, "single_end": false,
              "half_mapped_reads": 0, "half_mapped_frac": 0},
  "comparison_info": {
    "full":  {"pair": "dehost:host", "n_shared": 0, "n_only_A": 3,
              "n_only_B": 1, "jaccard": 0.0,
              "bray_curtis_dissimilarity": 1.0, "spearman_rho": null,
              "n_taxa_considered": 4, "n_enriched": 3, "n_background": 1,
              "enriched_taxa": []},
    "patho": {"pair": "dehost:host", "...": "..."}
  }
}
```

Новые поля v1.1: `positivity_source`, `annotator` (инструмент+версия из
контракта профиля), `db_id` (идентификаторы БД), `target_panel`,
`comparison_info` (информационная сводка dehost:host, на вердикт не
влияет). В v1.2 `comparison_info` несёт метрики сходства профилей
comparison model v2 (`n_shared`, `n_only_A/B`, `jaccard`,
`bray_curtis_dissimilarity`, `spearman_rho`). В v1.3 счётчики критериев
1–3 и 6 переведены с flagstat на per-fraction counts JSON (primary-only,
`samtools view -c -F 0x900`); в `metrics` добавлено `single_end`
(автоопределение SE по флагам BAM; для SE-образцов «pairs»-метрики — в
прочтениях/фрагментах). В v1.4.0 профили — контракт profile v2
(`docs/profile_contract.md`): `annotator` собирается из раздельных полей
`annotator` + `annotator_version` sidecar (например,
`"kraken2_bracken 2.17.2+3.0.1"`), `db_id` ключируется ролью из dbsheet
(`full`/`patho`); `dehost_human_frac` вычисляется по строке taxid 9606
в TSV (классификаторо-независимо), а не по clade-проценту отчёта
Kraken2 — на toy-данных это даёт 0.0 вместо 0.0003 у TOY_NORMAL
(единичные human-прочтения ниже порога отчётности Bracken; на вердикт
не влияет — порог критерия 0.05). В v1.5.0 критерий 6 проверяет точную
партицию DEHOST `dehost+host+half == trimmed`; в `metrics` добавлены
`half_mapped_reads` / `half_mapped_frac` (информационные, на вердикт не
влияют; в режиме bam — 0). В v1.6.0 изменилась только валидация входов
в main.nf (fail-fast обязательных параметров и колонок samplesheet/
dbsheet, резолвинг абсолютных/URI/относительных путей, glob индекса
`*.ht2`) — набор критериев, пороги и контракт QC JSON не менялись.
В v1.7.0 добавлены сквозной smoke-harness на toy-данных
(`test/run_smoke.sh` + `test/check_smoke.py`: pytest, три end-to-end
прогона — bam/fastq с Kraken2 и bam со stub-аннотатором metax — и
сверка вердиктов, сравнений dehost:host и обязательных артефактов),
самодостаточный профиль `test` (input/dbsheet/host_index/панель внутри
профиля) и docker-профиль с Dockerfile — критерии, пороги и контракт
QC JSON снова не менялись.
В v1.8.0 (NIPT-7) изменилось только окружение и ресурсы: точные пины
всех инструментов в `environment.yml`, per-process директивы
`conda`/`container` в модулях, исправлен ретрай Kraken2 по памяти
(удвоение через `MemoryUnit` вместо умножения строки), добавлены
`process.resourceLimits` и `manifest.nextflowVersion` — критерии,
пороги и контракт QC JSON не менялись. NB: пакет bracken 3.1
самоидентифицируется строкой «3.0.1» (не бампнутая константа
upstream), поэтому в `annotator` остаётся `2.17.2+3.0.1` — это
ожидаемо, а не признак неправильной версии.
В v1.9.0 (NIPT-16) добавлены fixture-driven `stub` и настоящий
`metax` 2.0.3: один dbsheet может смешивать Kraken2 и stub.
Критерии, пороги и структура QC JSON **не менялись**. В mixed-прогоне
full=Kraken2 сохраняет контроль human, patho=stub даёт PASS_NEGATIVE
(панель проверяется только по enriched-таксонам patho). Настоящий
MetaX на toy дал PASS_POSITIVE у NORMAL/SORTED/SE и FALLBACK по
обоим прежним причинам у LEAK. MetaX classified_reads здесь —
ceil суммы **отчётных** counts после фильтрации, а не сырой mapping
rate: min_classified_frac требует отдельной калибровки при смене
аннотатора. В contract abundance идёт reads/input; проценты MetaX
сохранены как x_abundance_pct. Full-индекс обязателен; SE-имена
уникализируются адаптером. 60/60 тестов, 4 smoke-прогона и отдельный
реальный MetaX-прогон прошли. Форматы, fixtures и план DB/container:
`docs/profile_contract.md`, §8. Поле верхнего уровня `annotator`
в QC JSON остаётся одной строкой первого доступного профиля:
для mixed оно не является полным перечнем аннотаторов. Точный
annotator/version каждой роли хранится в sidecar профиля.
`pass=null` у критерия
означает «не оценивался» (режим bam, гейтинг утечки, отсутствие входа
для вычисления).

## Стыковка с MultiQC

- MultiQC получает **только staged-входы** канала (custom content QC gate,
  kraken2 reports, flagstat, fastp json/html и summary HISAT2 стадии
  DEHOST в mode=fastq); `params.outdir` не сканируется — в отчёт
  попадает ровно то, что произвёл текущий запуск.
- Стандартные модули: fastp (json), samtools flagstat, kraken2 (report),
  hisat2 (summary, mode=fastq).
- Вердикты QC gate — custom content `<sample>_qc_mqc.tsv` (таблица
  критериев + `FINAL_STATUS`), генерируется из того же QC JSON.
