"""Adapter safeguards; test-only FASTA strings are not biological findings."""

import gzip

import pytest

from metax_prepare import require_full_index, unique_fasta_ids


def test_full_index_allowed(tmp_path):
    with gzip.open(tmp_path / "metax_db.tsv.gz", "wt") as out:
        out.write("#metax_db_version\t2\n#fractional\tfalse\n")
    require_full_index(tmp_path)


@pytest.mark.parametrize("metadata", ["#fractional\ttrue\n", "#version\t2\n"])
def test_fractional_or_unknown_index_rejected(tmp_path, metadata):
    with gzip.open(tmp_path / "metax_db.tsv.gz", "wt") as out:
        out.write(metadata)
    with pytest.raises(ValueError):
        require_full_index(tmp_path)


def test_duplicate_se_names_get_unique_ids_preserving_sequences(tmp_path):
    source, dest = tmp_path / "se.fa", tmp_path / "unique.fa"
    source.write_text(">same\nAC\nGT\n>same\nTGCA\n")
    assert unique_fasta_ids(source, dest) == 2
    assert dest.read_text() == ">metax_se_1\nAC\nGT\n>metax_se_2\nTGCA\n"
