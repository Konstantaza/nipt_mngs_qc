#!/usr/bin/env python3
"""
profile_io.py — общий reader контракта профиля (profile v2).

Единственное место потребителей, знающее устройство пары
profile.tsv + profile.json (sidecar). Спецификация —
docs/profile_contract.md. Используется compare_profiles.py, qc_gate.py
и validate_profile.py.

Принципы:
  * TSV читается ПО ИМЕНАМ колонок; колонки x_* и их порядок
    игнорируются (адаптеры могут добавлять свои);
  * ошибки контракта — ProfileError с человекочитаемым сообщением;
  * try_read_profile — небросающая обёртка для QC gate: отсутствие
    артефакта там — критерий PIPELINE_ERROR, а не падение скрипта.

read_profile возвращает dict:
  contract_version, sample_id, fraction, role, db_id, annotator,
  annotator_version, input_reads, classified_reads, unclassified_reads,
  status, taxa — [{taxid, name, rank, reads, abundance}]
(форма taxa совместима с profile.v1: taxid/name/reads — потребители
мигрируют без смены логики).
"""

import csv
import json
import os

CONTRACT_VERSION = "profile.v2"

REQUIRED_TSV_COLUMNS = ("sample", "fraction", "annotator",
                        "annotator_version", "db_id", "taxid", "name",
                        "rank", "reads", "abundance")

SIDECAR_REQUIRED_KEYS = ("contract_version", "sample_id", "fraction",
                         "role", "db_id", "annotator",
                         "annotator_version", "input_reads",
                         "classified_reads", "unclassified_reads",
                         "status")

STATUSES = ("ok", "empty")

# Колонки TSV, дублирующие поля идентичности sidecar (проверяется
# согласованность для непустого профиля)
_IDENTITY_TSV_TO_SIDECAR = {"sample": "sample_id", "fraction": "fraction",
                            "annotator": "annotator",
                            "annotator_version": "annotator_version",
                            "db_id": "db_id"}


class ProfileError(Exception):
    """Нарушение контракта profile v2."""


def load_sidecar(path):
    """sidecar JSON -> dict; ProfileError при отсутствии/битости/типах."""
    if not path or not os.path.exists(path) or os.path.getsize(path) == 0:
        raise ProfileError(f"sidecar JSON отсутствует или пуст: {path}")
    try:
        with open(path) as fh:
            data = json.load(fh)
    except json.JSONDecodeError as e:
        raise ProfileError(f"sidecar JSON битый ({path}): {e}")
    if not isinstance(data, dict):
        raise ProfileError(f"sidecar JSON не объект: {path}")
    for key in SIDECAR_REQUIRED_KEYS:
        if key not in data:
            raise ProfileError(f"{path}: нет обязательного поля {key!r}")
    if data["contract_version"] != CONTRACT_VERSION:
        raise ProfileError(
            f"{path}: ожидается contract_version={CONTRACT_VERSION}, "
            f"получено {data['contract_version']!r}")
    for key in ("input_reads", "classified_reads", "unclassified_reads"):
        v = data[key]
        if not isinstance(v, int) or isinstance(v, bool) or v < 0:
            raise ProfileError(
                f"{path}: поле {key!r} должно быть int >= 0, "
                f"получено {v!r}")
    if data["input_reads"] != \
            data["classified_reads"] + data["unclassified_reads"]:
        raise ProfileError(
            f"{path}: input_reads ({data['input_reads']}) != "
            f"classified ({data['classified_reads']}) + "
            f"unclassified ({data['unclassified_reads']})")
    if data["status"] not in STATUSES:
        raise ProfileError(
            f"{path}: status должен быть одним из {STATUSES}, "
            f"получено {data['status']!r}")
    return data


def load_taxa(path):
    """profile.tsv -> (header, taxa). taxa: [{taxid, name, rank, reads,
    abundance, _identity}] — типизированные значения."""
    if not path or not os.path.exists(path) or os.path.getsize(path) == 0:
        raise ProfileError(f"profile.tsv отсутствует или пуст: {path}")
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        header = reader.fieldnames or []
        missing = [c for c in REQUIRED_TSV_COLUMNS if c not in header]
        if missing:
            raise ProfileError(
                f"{path}: нет обязательных колонок: {', '.join(missing)}")
        unknown = [c for c in header
                   if c not in REQUIRED_TSV_COLUMNS
                   and not c.startswith("x_")]
        if unknown:
            raise ProfileError(
                f"{path}: неизвестные колонки без префикса x_: "
                + ", ".join(unknown))
        taxa = []
        for n, row in enumerate(reader, start=2):  # строка 1 — заголовок
            try:
                reads = float(row["reads"])
                abundance = float(row["abundance"])
            except (TypeError, ValueError):
                raise ProfileError(
                    f"{path}:{n}: reads/abundance не числа: "
                    f"{row.get('reads')!r}/{row.get('abundance')!r}")
            if reads < 0:
                raise ProfileError(f"{path}:{n}: reads < 0 ({reads})")
            if not 0.0 <= abundance <= 1.0:
                raise ProfileError(
                    f"{path}:{n}: abundance вне [0,1] ({abundance})")
            for col in ("taxid", "name", "rank"):
                if not (row.get(col) or "").strip():
                    raise ProfileError(
                        f"{path}:{n}: пустая обязательная колонка {col!r}")
            taxa.append({
                "taxid": row["taxid"].strip(),
                "name": row["name"].strip(),
                "rank": row["rank"].strip(),
                "reads": reads,
                "abundance": abundance,
                "_identity": {col: (row.get(col) or "").strip()
                              for col in _IDENTITY_TSV_TO_SIDECAR},
            })
    return header, taxa


def read_profile(tsv_path, json_path):
    """Пара TSV+sidecar -> dict контракта. Полная проверка инвариантов."""
    sidecar = load_sidecar(json_path)
    _, taxa = load_taxa(tsv_path)

    # status <=> число строк (см. docs/profile_contract.md §4)
    if sidecar["status"] == "empty" and taxa:
        raise ProfileError(
            f"{tsv_path}: status=empty, но TSV содержит {len(taxa)} строк")
    if sidecar["status"] == "ok" and not taxa:
        raise ProfileError(
            f"{tsv_path}: status=ok, но TSV без строк данных "
            "(пустой профиль обязан иметь status=empty)")

    # согласованность идентичности TSV <-> sidecar
    for t in taxa:
        for col, key in _IDENTITY_TSV_TO_SIDECAR.items():
            if t["_identity"][col] != str(sidecar[key]):
                raise ProfileError(
                    f"{tsv_path}: колонка {col!r} ({t['_identity'][col]!r}) "
                    f"не совпадает с sidecar {key!r} ({sidecar[key]!r})")

    # санити: сумма reads не превышает classified_reads
    # (доли от input; переоценщики типа Bracken не превышают classified)
    total = sum(t["reads"] for t in taxa)
    if total > sidecar["classified_reads"] * 1.01 + 1:
        raise ProfileError(
            f"{tsv_path}: сумма reads ({total}) > classified_reads "
            f"({sidecar['classified_reads']})")

    for t in taxa:
        t.pop("_identity")
    return {**sidecar, "taxa": taxa}


def try_read_profile(tsv_path, json_path):
    """-> (profile|None, error|None). Для QC gate: отсутствие артефакта —
    критерий PIPELINE_ERROR, а не исключение."""
    try:
        return read_profile(tsv_path, json_path), None
    except ProfileError as e:
        return None, str(e)
