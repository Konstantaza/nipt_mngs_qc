"""Software fixtures only; no generated reads or fabricated scientific results."""
import argparse
import contextlib
import gzip
import io
import json
from pathlib import Path
import subprocess
import sys
import tarfile
import tempfile
import unittest
from unittest.mock import patch

import server as s


class KitTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(prefix="metax_kit_tests_")
        self.root = Path(self.tmp.name)

    def tearDown(self):
        self.tmp.cleanup()

    def test_unsafe_path(self):
        for name in ("a b", "a;touch_bad", "a$b", "a'quote"):
            with self.subTest(name=name), self.assertRaises(ValueError):
                s.safe_path(self.root / name)
        self.assertEqual(s.safe_path(self.root / "abc-123"), self.root / "abc-123")

    def test_fasta_invalid(self):
        for text in ("ACGT\n", ">\nACGT\n", ">a\n", ">a\nACGTX\n"):
            with self.subTest(text=text):
                p = self.root / "bad.fa"
                p.write_text(text)
                with self.assertRaises(ValueError):
                    list(s.records(p))

    def test_reheader_preserves_sequences(self):
        refs = {}
        for source in s.SOURCE_SPECIES:
            p = self.root / f"{source}.fa"
            # Tiny arbitrary nucleotide strings are software fixtures, not data.
            p.write_text(">original description\nAcgTN\n>a2\nTTGC\n")
            refs[source] = p
        rows = s.reheader(refs, self.root)
        self.assertEqual(len(rows), 8)
        self.assertEqual([seq for _, seq in s.records(self.root / "reference_full.fasta")],
                         ["AcgTN", "TTGC"] * 4)
        self.assertEqual(len(list(s.records(self.root / "reference_patho.fasta"))), 6)
        self.assertEqual(len(list(s.records(self.root / "reference_pilot.fasta"))), 2)
        self.assertEqual({r["species_taxid"] for r in rows}, {"9606", "10847", "4932", "562"})

    def test_duplicate_contigs_rejected(self):
        p = self.root / "bad.fa"
        p.write_text(">same first\nACGT\n>same second\nTGCA\n")
        with self.assertRaises(ValueError):
            s.reheader({"562": p}, self.root)

    def test_seq_map_ignores_case_not_identity(self):
        p, q = self.root / "a.fa", self.root / "b.fa"
        p.write_text(">same\nacgT\n")
        q.write_text(">same\nACGT\n")
        self.assertEqual(s.seq_map(p), s.seq_map(q))
        q.write_text(">different\nACGT\n")
        self.assertNotEqual(s.seq_map(p), s.seq_map(q))

    def test_fractional_rejected(self):
        (self.root / "metax_db.mxi").write_bytes(b"software-fixture")
        for value in ("true", "false", None):
            with self.subTest(value=value):
                with gzip.open(self.root / "metax_db.tsv.gz", "wt") as f:
                    f.write("#fractional\t" + value + "\n" if value else "#other\tvalue\n")
                if value == "false":
                    s.full_index(self.root)
                else:
                    with self.assertRaises(ValueError):
                        s.full_index(self.root)

    def test_real_taxonomy(self):
        if not hasattr(self.__class__, "taxonomy"):
            self.skipTest("No reference taxonomy supplied")
        result = s.taxonomy_check(self.__class__.taxonomy)
        self.assertEqual(result["nodes"], 43)
        self.assertEqual(result["phiX_source_to_species"], {"2886930": "10847"})

    def test_pipeline_command_preserves_toy(self):
        args = argparse.Namespace(kind="pe", repo=self.root / "repo", toy=self.root / "toy",
                                  out=self.root / "new", resume=False)
        data = {"hisat2_index": str(self.root / "human"), "target_panel": str(self.root / "panel")}
        cmd = s.pipeline_command(args, data)
        self.assertEqual(cmd[cmd.index("-profile") + 1], "standard")
        self.assertEqual(cmd[cmd.index("--input") + 1], str(self.root / "toy/samplesheet_pe.csv"))
        self.assertEqual(cmd[cmd.index("--bracken_read_len") + 1], "126")
        self.assertNotIn("-resume", cmd)
        args.kind, args.resume = "se", True
        cmd = s.pipeline_command(args, data)
        self.assertEqual(cmd[cmd.index("--bracken_read_len") + 1], "50")
        self.assertIn("-resume", cmd)

    def test_collect_excludes_large_inputs(self):
        for folder in ("results_pe/SAMPLE/profiles", "results_pe/SAMPLE/validation",
                       "work_pe/ab/abc", "db_full", "logs"):
            (self.root / folder).mkdir(parents=True)
        for name in ("results_pe/SAMPLE/profiles/result.profile.json",
                     "logs/pipeline_pe.log", "work_pe/ab/abc/.command.err"):
            (self.root / name).write_text("{}")
        for name in ("results_pe/SAMPLE/validation/host.bam", "reads.fastq",
                     "db_full/metax_db.mxi", "work_pe/ab/abc/reads.fasta"):
            (self.root / name).write_text("not-to-be-collected")
        with contextlib.redirect_stdout(io.StringIO()):
            s.collect(argparse.Namespace(out=self.root))
        archive, = self.root.glob("metax_samovar_results_*.tar.gz")
        with tarfile.open(archive) as t:
            names = t.getnames()
            self.assertEqual(len(names), 4)
            self.assertTrue(any(n.endswith(".command.err") for n in names))
            self.assertFalse(any(n.endswith((".bam", ".fastq", ".mxi", ".fasta")) for n in names))

    def test_measure_reports_failure_and_memory(self):
        metrics = self.root / "measure.json"
        p = subprocess.run([sys.executable, str(Path(s.__file__)), "_measure", str(metrics),
                            sys.executable, "-c", "raise SystemExit(7)"], capture_output=True)
        self.assertEqual(p.returncode, 7)
        data = json.loads(metrics.read_text())
        self.assertEqual(data["exit_code"], 7)
        self.assertGreater(data["peak_rss_kib_linux"], 0)

    def test_prepare_never_overwrites(self):
        with self.assertRaisesRegex(ValueError, "overwrite"):
            s.prepare(argparse.Namespace(out=self.root))

    def launcher_fixture(self):
        (self.root / "logs").mkdir()
        metax = self.root / "pinned-metax"
        metax.write_text("unit-test fixture, never executed")
        args = argparse.Namespace(kind="pe", repo=self.root / "repo", toy=self.root / "toy",
                                  out=self.root, resume=False)
        data = {"metax": str(metax), "hisat2_index": str(self.root / "human"),
                "target_panel": str(self.root / "panel")}
        return args, data

    def test_launcher_path_and_success_marker(self):
        args, data = self.launcher_fixture()
        def fake_nextflow(command, log, env, cwd, measured):
            self.assertEqual(env["PATH"].split(":")[0], str(self.root / "toolbin"))
            self.assertEqual((self.root / "toolbin/metax").resolve(), Path(data["metax"]))
            self.assertEqual(cwd, self.root / "launch_pe")
            self.assertTrue(measured)
            folder = self.root / "results_pe/unit_fixtures"
            folder.mkdir(parents=True)
            for n in range(4):
                s.write_json(folder / f"{n}.profile.json", {"annotator": "metax"})
        with patch.object(s, "state", return_value=data), patch.object(s, "check_tools"), \
                patch.object(s, "validate_db"), patch.object(s, "run_logged", side_effect=fake_nextflow), \
                contextlib.redirect_stdout(io.StringIO()):
            s.pipeline_run(args)
        self.assertEqual(s.read_json(self.root / "run_pe.json")["exit_code"], 0)
        self.assertEqual(s.read_json(self.root / "attempt_pe.json")["state"], "completed")

    def test_launcher_missing_profiles_is_failure(self):
        args, data = self.launcher_fixture()
        with patch.object(s, "state", return_value=data), patch.object(s, "check_tools"), \
                patch.object(s, "validate_db"), patch.object(s, "run_logged"), \
                self.assertRaisesRegex(ValueError, "Expected host/dehost"):
            s.pipeline_run(args)
        self.assertFalse((self.root / "run_pe.json").exists())
        self.assertEqual(s.read_json(self.root / "attempt_pe.json")["state"], "failed")

    def test_launcher_requires_explicit_resume(self):
        args, data = self.launcher_fixture()
        (self.root / "work_pe").mkdir()
        with patch.object(s, "state", return_value=data), patch.object(s, "check_tools"), \
                patch.object(s, "validate_db"), self.assertRaisesRegex(ValueError, "resume"):
            s.pipeline_run(args)


if __name__ == "__main__":
    if "--taxonomy" in sys.argv:
        i = sys.argv.index("--taxonomy")
        KitTests.taxonomy = Path(sys.argv[i + 1])
        del sys.argv[i:i + 2]
    unittest.main()
