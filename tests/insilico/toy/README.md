# Samovar toy: 99,5% HOST — PE126 и настоящий SE50

**Фактический статус v1.10.0: строгая приёмка FAIL.** Сам Samovar и оба
Nextflow-прогона успешно выполнены здесь, но обнаружена human-утечка
в DEHOST. Инструкции ниже воспроизводят измерение, а не обещают PASS.

## 1. Репозитории рядом

```text
work/
  samovar/
  nipt_mngs_qc/
```

Запускать команды ниже из `work/nipt_mngs_qc`. Нужны Linux, Java,
Nextflow ≥26.04.6, git, uv и conda/mamba. Не нужен C++-аннотатор Samovar,
Kaiju, MetaX, Docker или полный установщик Samovar для генерации.

Проверенная ревизия Samovar:
`d6bb8246f9e31c7890f4393aa92bc7e0aa53d593` (0.12.150).
Неполного набора папок недостаточно: нужен `src/samovar`.
Если соседний репозиторий ещё не получен:

```bash
git clone https://github.com/CTLab/samovar.git ../samovar
git -C ../samovar checkout d6bb8246f9e31c7890f4393aa92bc7e0aa53d593
```

Если там уже есть свои изменения, **не делайте checkout поверх них**:
используйте отдельный checkout/worktree или существующую совместимую версию.
Скрипт сохраняет фактическую ревизию и dirty state, проверяет интерфейсы
до генерации; другая ревизия не считается идентичным воспроизведением.

## 2. Разовая подготовка окружений

Пайплайн использует Python 3.14; Samovar требует отдельный Python 3.11.
Не понижайте Python существующего окружения пайплайна.

```bash
conda env create -f environment.yml
conda activate nipt_mngs_qc
# Если окружение уже есть: conda env update -f environment.yml
# Если uv отсутствует:
conda install -c conda-forge uv

bash tests/insilico/toy/setup_samovar_env.sh ../samovar
nextflow -version
seqkit version
```

Setup создаёт `../samovar/.venv-nipt-toy`, использует Python 3.11
(он должен быть доступен uv; при необходимости `uv python install 3.11`)
и `samovar_requirements.lock.txt`, затем ставит сам соседний Samovar.
Существующий venv не затирается. Скрипт проверен созданием второго
чистого окружения здесь. Можно передать иной путь venv вторым аргументом;
тогда используйте `--samovar-python /absolute/venv/bin/python`.
Важны ISS и Snakemake рядом с этим Python.

Точные pipeline-пины — `environment.yml`; Bracken conda package 3.1
в проверенном окружении сообщает CLI version 3.0.1. Это записано,
а не скрыто. Локально был Nextflow 26.10.4.

## 3. Воспроизвести штатный toy одной командой

```bash
bash tests/insilico/toy/run.sh all \
  --samovar-dir ../samovar \
  --outdir insilico_runs/toy_995
```

По умолчанию: seed=42, 600000 пар, HOST=0.995, n_samples=1,
ISS hiseq, генерация одним потоком, остальные шаги двумя потоками.
HOST — **укороченный штатный toy** `../samovar/data/test_genomes/host/9606.fna`.
Индекс dehost строится из того же FASTA.

Порядок:
1. Проверить CLI/версии/референсы; создать run-local `samovar_config.json`.
2. Настоящий `samovar generate` с явным `--host_fraction 0.995`;
   `prepare --startpoint setup_reads --endpoint setup_reads`; `exec`.
   Prepare нужен, поскольку один generate не создаёт основной скрипт exec.
3. Проверить фактические PE126, taxid в ID, количества и pairing.
4. `seqkit shuffle --rand-seed 42` для обоих mates; проверить совпадение
   **каждой** пары после перемешивания и сохранение набора IDs.
5. Создать SE из R1, обрезав sequence **и quality** до 50 nt.
6. Собрать full Kraken2 (human + три микроба) и software-test patho panel
   (только три микроба); Bracken distributions на 126 и 50 nt.
7. Отдельные `mode=fastq` PE/SE прогоны; GT-аудит trimmed и BAM-фракций.

Штатный toy-builder Samovar с мутациями/пропусками здесь не используется.
Taxid source 2886930 сохраняется в GT, species 10847 — в Bracken;
сопоставление записано в `taxid_map.tsv`, реальные ранги не изменяются.

`samplesheet_pe.csv`, `samplesheet_se.csv`, `dbsheet.csv` в этом каталоге —
переносимые примеры для стандартного outdir. Скрипт также генерирует
реальные CSV внутри выбранного outdir с точными абсолютными путями.

## 4. Запуск с вашим человеческим референсом и индексом на сервере

Окружения — как выше. Укажите **тот же FASTA**, из которого построен индекс:

```bash
bash tests/insilico/toy/run.sh all \
  --samovar-dir ../samovar \
  --host-fasta /absolute/path/to/human.fa \
  --hisat2-index /absolute/path/to/hisat2/human \
  --threads 4 \
  --task-memory '16 GB' \
  --kraken-memory '16 GB' \
  --outdir insilico_runs/human_995
```

Это шаблон: замените оба абсолютных пути; префикс индекса — без `.1.ht2`.
`.ht2` и `.ht2l` поддержаны. Содержимое индекса извлекается
`hisat2-inspect` и сравнивается с нормализованным FASTA по каждому контигу;
несовпадение останавливает build. Можно не передавать индекс: он будет
построен заново из заданного FASTA. Представление `refs/9606.fna`
нужно для корректной human-метки Samovar, последовательности те же.

**Полный геном здесь не запускался.** 16 GB — начальный запрос памяти,
не гарантия достаточности; согласуйте его с размером реально построенной
Kraken2 БД и доступной RAM. Человеческий FASTA включается также в full
Kraken2, поэтому время, диск и память существенно выше toy.
Для первого полного запуска предусмотрите десятки GB свободного места;
порог preflight 4 GiB — только минимальный toy-check, не sizing полного генома.
Скрипт не скачивает GRCh38 автоматически и не резервирует сервер.

## 5. Где увидеть результат

В `insilico_runs/toy_995/` (или выбранном outdir):

- `acceptance.tsv` — каждая проверка, наблюдаемое/ожидаемое, PASS/FAIL.
- `acceptance.json` — те же проверки и подробные GT/QC/profile-счётчики.
- `results_pe/SAMOVAR_PE/qc/SAMOVAR_PE_qc_report.json`.
- `results_se/SAMOVAR_SE/qc/SAMOVAR_SE_qc_report.json`.
- `results_pe/multiqc/multiqc_report.html` и аналогичный SE HTML.
- `results_{pe,se}/SAMOVAR_{PE,SE}/profiles/` — profile TSV/JSON для
  DEHOST/HOST × full/patho, dehost_vs_host comparison TSV/JSON.
- `generation.json`, `manifest.json`, `build.json`, `tool_versions.json`,
  `samovar_packages.txt`, `*.log`, `*.trace.tsv` — provenance и диагностика.
- `validation/trimmed` и `validation/partitions` внутри sample —
  данные для проверки GT (в обычном пайплайне публикация выключена).

```bash
column -s $'\t' -t insilico_runs/toy_995/acceptance.tsv
python -m json.tool insilico_runs/toy_995/acceptance.json
```

**Exit 0 означает все критерии PASS. Exit 1 требует просмотра отчёта/лога.**
Если Nextflow оба раза завершился, но checker вернул 1 — это может быть
настоящий биологический FAIL, как в проверенном запуске ниже.
При падении до проверки смотрите названный stage log; прежний PASS
не сохраняется как действительный результат прерванного rerun.

## 6. Фактический результат здесь

| Показатель | PE126 | SE50 |
|---|---:|---:|
| Исходный объём | 600000 пар | 600000 reads R1[:50] |
| GT HOST до QC | 99.5% | 99.5% |
| HOST / (HOST+DEHOST+half), read records | 99.41113% | 99.40630% |
| HOST / (HOST+DEHOST), QC-метрика | 99.4882% | 99.4063% |
| DEHOST, primary read records | 6008 | 3516 |
| GT микробных read records в DEHOST | 5968 | 3000 |
| phiX174 / yeast / E. coli, среди трёх microbes | 33.3445 / 33.2775 / 33.3780% | 34.0292 / 32.0807 / 33.8900% |
| Human в full DEHOST profile | 20 пар | 224 reads |
| QC | PASS_POSITIVE | FALLBACK / DEHOST_FAILURE |
| Comparison и MultiQC | есть | есть |

Строгая приёмка: **31/34 checks PASS, общий FAIL**.
PE нарушает запрет лишних species с reads≥10 из-за human.
SE нарушает тот же запрет и порог human leakage≤5%: 224/3516=6.3709%.
Это не исчезновение микробов и не доказанная мисклассификация:
по GT в DEHOST остались 40 human-мейтов PE и 516 human SE reads.
Пороги/профили/данные не подгонялись ради PASS.

Локальный all занял около 664 s. Максимальный child RSS, сообщённый
Linux resource accounting, около 1.15 GiB; это не сумма одновременно
работавших процессов. Рабочие файлы около 1.8 GB без окружений.
Два независимых pilot по 20000 пар дали побайтово идентичные raw,
shuffled PE и SE FASTQ; приёмка pilot намеренно не считается PASS
по глубине. Полный человеческий индекс и Docker не тестировались.

## 7. Стадии, повторный запуск, проверки

Те же аргументы можно использовать с `preflight`, `generate`, `build`,
`run`, `check` вместо `all`. `generate` и `build` проверяют hashes готовых
результатов и переиспользуют их. Изменение seed/depth/reference/версии
Samovar требует нового outdir. Nextflow `run` использует `-resume`;
не удаляйте `work/` и `.nextflow/` до завершения/проверки.
Незавершённые db_full/db_patho без build.json нужно перенести в сторону
перед новым build (обе папки), чтобы библиотека не добавлялась повторно.

После изменения кода rerun `run`, затем `check`: sentinel содержит hashes
кода/конфигурации/входов. Старый результат не считается проверкой новой версии.

```bash
python -m pytest test/ tests/insilico/toy/test_toy.py -q
bash test/run_smoke.sh
```

Проверено 87/87 тестов и четыре старых сквозных прогона. В том числе
mixed PE+SE samplesheet, пустой SE half, counts/FASTP record units,
некорректные SE-флаги в samplesheet, парный shuffle, GT IDs, extra taxa,
отсутствующий вид и граница относительных 20%.

Таксономия в `taxonomy/` включена для offline build:
43 реальных узла и оригинальные ENA JSON-ответы с URL.
`fetch_taxonomy.py` — отдельное явное обновление, не автоматический
шаг при запуске; после обновления нужен новый outdir и review mapping.

Closed toy не проверяет клиническую чувствительность, контаминацию
реагентов или специфичность к организмам, которых нет в маленькой БД.
