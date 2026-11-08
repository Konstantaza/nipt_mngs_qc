"""
Юнит-тесты sum_counts.py (v1.5.0) — слияние counts JSON конституентов
в counts JSON фракции all (схема bam_counts.py сохраняется).

Запуск: python3 -m pytest test/test_sum_counts.py -v
"""

import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SUM_COUNTS = ROOT / "bin" / "sum_counts.py"


def _counts(path, primary, paired, single_end):
    """counts JSON схемы bam_counts.py."""
    path.write_text(json.dumps({
        "sample_id": "T", "fraction": "x",
        "primary_reads": primary, "paired_reads": paired,
        "unpaired_reads": primary - paired,
        "single_end": single_end,
        "fragments": paired // 2 + (primary - paired),
    }))


def _run(tmp_path, inputs, expect_ok=True):
    out = tmp_path / "T.all.counts.json"
    res = subprocess.run([sys.executable, str(SUM_COUNTS),
                          "--sample", "T", "--fraction", "all",
                          "--out", str(out)]
                         + [str(p) for p in inputs],
                         capture_output=True, text=True)
    if expect_ok:
        assert res.returncode == 0, f"sum_counts упал: {res.stderr}"
        return json.loads(out.read_text())
    assert res.returncode != 0, "sum_counts должен был упасть"
    return None


def test_merge_two_pe(tmp_path):
    """dehost+host (PE): суммы счётчиков, single_end = AND (False)."""
    d = tmp_path / "T.dehost.counts.json"
    h = tmp_path / "T.host.counts.json"
    _counts(d, primary=100, paired=100, single_end=False)
    _counts(h, primary=200, paired=200, single_end=False)
    c = _run(tmp_path, [d, h])
    assert c["fraction"] == "all"
    assert c["primary_reads"] == 300
    assert c["paired_reads"] == 300
    assert c["unpaired_reads"] == 0
    assert c["fragments"] == 150
    assert c["single_end"] is False


def test_merge_three_se_with_zeros(tmp_path):
    """dehost+host+half (SE, half пустой): нули не ломают сумму,
    single_end = AND (True)."""
    d = tmp_path / "T.dehost.counts.json"
    h = tmp_path / "T.host.counts.json"
    hf = tmp_path / "T.half.counts.json"
    _counts(d, primary=10, paired=0, single_end=True)
    _counts(h, primary=20, paired=0, single_end=True)
    _counts(hf, primary=0, paired=0, single_end=True)
    c = _run(tmp_path, [d, h, hf])
    assert c["primary_reads"] == 30
    assert c["paired_reads"] == 0
    assert c["unpaired_reads"] == 30
    assert c["fragments"] == 30
    assert c["single_end"] is True


def test_broken_input_rejected(tmp_path):
    """JSON без обязательных ключей — ошибка, а не молчаливые нули."""
    bad = tmp_path / "bad.counts.json"
    bad.write_text(json.dumps({"primary_reads": 5}))
    _run(tmp_path, [bad], expect_ok=False)
