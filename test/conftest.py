"""
Общие helpers юнит-тестов: запись фикстур контракта profile v2
(profile.tsv + sidecar profile.json, docs/profile_contract.md).
"""

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "bin"))

from profile_io import CONTRACT_VERSION, REQUIRED_TSV_COLUMNS  # noqa: E402

ANNOTATOR = "kraken2_bracken"
ANNOTATOR_VERSION = "2.17.2+3.0.1"
# Строка annotator в QC JSON собирается как "annotator annotator_version"
ANNOTATOR_COMPOSED = f"{ANNOTATOR} {ANNOTATOR_VERSION}"


def write_profile_v2(directory, sample, fraction, role, input_reads, taxa,
                     db_id=None, annotator=ANNOTATOR,
                     annotator_version=ANNOTATOR_VERSION,
                     classified_reads=None, unclassified_reads=0,
                     status=None):
    """Пишет пару <sample>.<fraction>.<role>.profile.{tsv,json}.

    taxa: [(taxid, name, reads)] — rank пишется 'species'.
    По умолчанию: classified_reads = input_reads, status = ok при
    непустом taxa. Возвращает (tsv_path, json_path).
    """
    directory = Path(directory)
    db_id = db_id or f"db_{role}"
    if classified_reads is None:
        classified_reads = input_reads
    if status is None:
        status = "ok" if taxa else "empty"
    base = directory / f"{sample}.{fraction}.{role}.profile"
    tsv_path = Path(str(base) + ".tsv")
    json_path = Path(str(base) + ".json")
    with open(tsv_path, "w") as fh:
        fh.write("\t".join(REQUIRED_TSV_COLUMNS) + "\n")
        for taxid, name, reads in taxa:
            abundance = reads / input_reads if input_reads else 0.0
            fh.write("\t".join([sample, fraction, annotator,
                                annotator_version, db_id, str(taxid), name,
                                "species", f"{reads:g}",
                                f"{abundance:g}"]) + "\n")
    sidecar = {
        "contract_version": CONTRACT_VERSION,
        "sample_id": sample,
        "fraction": fraction,
        "role": role,
        "db_id": db_id,
        "annotator": annotator,
        "annotator_version": annotator_version,
        "input_reads": input_reads,
        "classified_reads": classified_reads,
        "unclassified_reads": unclassified_reads,
        "status": status,
    }
    json_path.write_text(json.dumps(sidecar))
    return tsv_path, json_path
