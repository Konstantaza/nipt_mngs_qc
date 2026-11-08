"""Unit checks: fixtures here are deliberately tiny software-test records."""
import importlib.util
import json
from pathlib import Path
from types import SimpleNamespace
import pytest

SPEC = importlib.util.spec_from_file_location("toy_acceptance", Path(__file__).with_name("toy.py"))
toy = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(toy)


@pytest.mark.parametrize("header,expected", [
    ("@taxid:562_9/1", "562"),
    ("@taxid:4932_2/2 extra", "4932"),
    ("@taxid:2886930_100/1", "2886930"),
    ("@fragment|taxid:9606/1", "9606"),
    ("fragment|taxid:9606", "9606"),
])
def test_truth_id(header, expected):
    assert toy.taxid(header) == expected


@pytest.mark.parametrize("header", ["@x_99/1", "@taxid:999_0/1",
                                    "@taxid:562_0|taxid:9606/1"])
def test_bad_truth(header):
    with pytest.raises(ValueError):
        toy.taxid(header)


def make_fq(path, ids, n=50):
    path.write_text("".join(f"@{i}\n{'A'*n}\n+\n{'I'*n}\n" for i in ids))
    return path


def test_pair_shuffle_must_be_synchronised(tmp_path):
    r1 = make_fq(tmp_path/"r1.fq", ["taxid:562_1/1", "taxid:4932_2/1"])
    r2 = make_fq(tmp_path/"r2.fq", ["taxid:4932_2/2", "taxid:562_1/2"])
    with pytest.raises(ValueError, match="Mate mismatch"):
        toy.check_pairs(r1, r2)


def test_pair_count_mismatch(tmp_path):
    r1 = make_fq(tmp_path/"r1.fq", ["taxid:562_1/1", "taxid:4932_2/1"])
    r2 = make_fq(tmp_path/"r2.fq", ["taxid:562_1/2"])
    with pytest.raises(ValueError):
        toy.check_pairs(r1, r2)


def test_pairs_good(tmp_path):
    r1 = make_fq(tmp_path/"r1.fq", ["taxid:562_1/1", "taxid:4932_2/1"])
    r2 = make_fq(tmp_path/"r2.fq", ["taxid:562_1/2", "taxid:4932_2/2"])
    assert toy.check_pairs(r1, r2) == 2


def test_quality_length_mismatch(tmp_path):
    path = tmp_path/"bad.fq"
    path.write_text("@taxid:562_0/1\nAAAA\n+\nIII\n")
    with pytest.raises(ValueError):
        list(toy.fq(path))


def test_duplicate_id(tmp_path):
    path = make_fq(tmp_path/"bad.fq", ["taxid:562_1/1"]*2)
    with pytest.raises(ValueError, match="Duplicate"):
        toy.inspect_fastq(path)


def test_wrong_read_length(tmp_path):
    path = make_fq(tmp_path/"se.fq", ["taxid:562_1/1"], 50)
    assert toy.inspect_fastq(path, 50)[0]["reads"] == 1
    with pytest.raises(ValueError):
        toy.inspect_fastq(path, 126)


@pytest.mark.parametrize("counts,ok", [
    ({"10847": 1000, "4932": 1000, "562": 1000}, True),
    ({"10847": 1200, "4932": 1000, "562": 800}, True),
    ({"10847": 1500, "4932": 1000, "562": 500}, False),
    ({"10847": 1000, "4932": 1000}, False),
    ({"10847": 1000, "4932": 1000, "562": 1000, "9606": 10}, False),
    ({"10847": 1000, "4932": 1000, "562": 1000, "9606": 9}, True),
])
def test_species_acceptance(tmp_path, counts, ok):
    path = tmp_path/"profile.tsv"
    path.write_text("taxid\treads\n" + "".join(f"{t}\t{n}\n" for t, n in counts.items()))
    assert toy.profile_metrics(path)["pass"] == ok


def test_fasta_comparison_ignores_wrap_not_sequence(tmp_path):
    a, b = tmp_path/"a.fa", tmp_path/"b.fa"
    a.write_text(">chr1 description\nacgt\nacgt\n")
    b.write_text(">chr1\nACGTACGT\n")
    assert toy.fasta_hashes(a) == toy.fasta_hashes(b)
    b.write_text(">chr1\nACGTACGA\n")
    assert toy.fasta_hashes(a) != toy.fasta_hashes(b)


def test_taxonomy_preserves_phix_rank():
    nodes = {}
    for line in (toy.HERE/"taxonomy/nodes.dmp").read_text().splitlines():
        tid, parent, rank, *_ = [x.strip() for x in line.split("|")]
        nodes[tid] = (parent, rank)
    assert nodes["2886930"] == ("10847", "no rank")
    assert nodes["10847"][1] == "species"
    for t in toy.TAXIDS:
        seen = set()
        while t != "1":
            assert t not in seen
            seen.add(t)
            t = nodes[t][0]


def test_missing_run_never_passes(tmp_path):
    assert not toy.acceptance(SimpleNamespace(outdir=tmp_path, pairs=600000))
    result = json.loads((tmp_path/"acceptance.json").read_text())
    assert result["pass"] is False
    assert any(not x["pass"] for x in result["checks"])
