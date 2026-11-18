#!/usr/bin/env python3
"""Portable SOFTWARE tests. Fixtures are NOT sequencing/biological results.

Install runner at tools/validation/run_metax_confirm_v1.py and this file at
test/test_run_metax_confirm_v1.py. Run: python -m unittest discover -s test
-p test_run_metax_confirm_v1.py -v
No aligner, MetaX, Nextflow workflow, installation or network is executed.
"""
import argparse
import contextlib
import gzip
import importlib.util
import io
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest import mock


HERE = Path(__file__).resolve()
CANDIDATES = (HERE.parent / "run_metax_confirm_v1.py",
              HERE.parents[1] / "tools/validation/run_metax_confirm_v1.py")
RUNNER = next((p for p in CANDIDATES if p.is_file()), None)
if RUNNER is None:
    raise RuntimeError("Install run_metax_confirm_v1.py in tools/validation/")
spec = importlib.util.spec_from_file_location("metax_confirm", RUNNER)
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


class RunnerTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(prefix="metax_software_test_")
        self.root = Path(self.tmp.name)
        self.repo = self.root / "repo"
        self.repo.mkdir()
        self.args = argparse.Namespace(
            repo=self.repo, toy=self.repo / "insilico_runs/toy_995",
            baseline=self.repo / "insilico_runs/metax_samovar_v1",
            out=self.repo / "insilico_runs/metax_confirm_v2")
        self.args.toy.mkdir(parents=True)
        self.args.baseline.mkdir(parents=True)

    def tearDown(self):
        self.tmp.cleanup()

    def put(self, path, text="software-test fixture\n"):
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(text)
        return path

    def fixture_manifests(self):
        a = self.args
        source = self.put(a.repo / "main.nf")
        index = a.toy / "host_index/human"
        inputs = {str(source): m.sha(source)}
        for i in range(1, 9):
            p = self.put(Path(f"{index}.{i}.ht2"))
            inputs[str(p)] = m.sha(p)
        for kind in ("pe", "se"):
            p = self.put(a.toy / f"samplesheet_{kind}.csv",
                         f"sample,fastq_1,fastq_2,single_end\nSAMOVAR_{kind.upper()},fixture,,true\n")
            inputs[str(p)] = m.sha(p)
        binary = self.put(self.root / "metax")
        config = self.put(a.baseline / "execution.config")
        sheet = self.put(a.baseline / "dbsheet_metax.csv",
                         "db_id,annotator,path,role\n" +
                         "".join(f"test_{r},metax,{a.baseline / ('db_'+r)},{r}\n"
                                 for r in ("full", "patho")))
        db_hashes = {}
        for role in ("full", "patho"):
            db = a.baseline / f"db_{role}"
            self.put(db / "metax_db.mxi")
            with gzip.open(db / "metax_db.tsv.gz", "wt") as handle:
                handle.write("#fractional\tfalse\n")
            for name in ("nodes.dmp", "names.dmp", "merged.dmp", "delnodes.dmp"):
                self.put(db / "taxonomy" / name, "" if name.startswith(("merged", "delnodes")) else "test\n")
            db_hashes.update({str(p): m.sha(p) for p in db.rglob("*") if p.is_file()})
        prep = {"prepare_completed": True, "human_index_sequences_verified": True,
                "repo": str(a.repo), "toy": str(a.toy), "out": str(a.baseline),
                "hisat2_index": str(index), "target_panel": "software-fixture",
                "metax": str(binary), "metax_sha256": m.sha(binary),
                "input_sha256": inputs,
                "prepared_sha256": {str(p): m.sha(p) for p in (sheet, config)}}
        build = {"index_mode": "full", "fractional": False, "files_sha256": db_hashes}
        manifests = {"prepare.json": prep, "build_metax.json": build}
        for name, data in manifests.items():
            m.write_json(a.baseline / name, data)
        return manifests, {n:m.canonical_sha(d) for n,d in manifests.items()}

    def test_canonical_manifest_independent_of_json_layout(self):
        self.assertEqual(m.canonical_sha({"b":2,"a":1}), m.canonical_sha({"a":1,"b":2}))

    def test_hash_mismatch_stops(self):
        p = self.put(self.root / "input")
        with self.assertRaisesRegex(ValueError, "SHA256 mismatch"):
            m.verify_hashes({p:"0"*64})

    def test_missing_hash_file_stops(self):
        with self.assertRaisesRegex(ValueError, "Missing pinned"):
            m.verify_hashes({self.root / "absent":"0"*64})

    def test_source_change_stops(self):
        p = self.put(self.repo / "main.nf")
        with mock.patch.dict(m.SOURCE_PINS, {"main.nf":"0"*64}, clear=True):
            with self.assertRaisesRegex(ValueError, "SHA256"):
                m.verify_sources(self.repo)

    def test_extra_runtime_file_stops_but_validation_tools_allowed(self):
        p = self.put(self.repo / "main.nf")
        for folder in ("bin", "modules", "subworkflows"):
            (self.repo / folder).mkdir()
        self.put(self.repo / "tools/validation/runner.py")
        with mock.patch.dict(m.SOURCE_PINS, {"main.nf":m.sha(p)}, clear=True):
            self.assertEqual(m.verify_sources(self.repo), {"main.nf":m.sha(p)})
            self.put(self.repo / "bin/unreviewed.py")
            with self.assertRaisesRegex(ValueError, "file set differs"):
                m.verify_sources(self.repo)

    def test_genuine_manifest_schema_fixture(self):
        _, pins = self.fixture_manifests()
        with mock.patch.dict(m.MANIFEST_PINS, pins, clear=True):
            self.assertEqual(set(m.load_manifests(self.args)), set(pins))

    def test_manifest_change_stops(self):
        self.fixture_manifests()
        with self.assertRaisesRegex(ValueError, "manifest mismatch"):
            m.load_manifests(self.args)

    def test_db_change_stops(self):
        _, pins = self.fixture_manifests()
        self.put(self.args.baseline / "db_full/metax_db.mxi", "changed\n")
        with mock.patch.dict(m.MANIFEST_PINS, pins, clear=True):
            with self.assertRaisesRegex(ValueError, "SHA256"):
                m.load_manifests(self.args)

    def test_added_cache_not_rejected(self):
        _, pins = self.fixture_manifests()
        self.put(self.args.baseline / "db_full/metax_db.mxi.metamcache")
        with mock.patch.dict(m.MANIFEST_PINS, pins, clear=True):
            m.load_manifests(self.args)

    def test_fractional_database_rejected(self):
        self.fixture_manifests()
        path = self.args.baseline / "db_full/metax_db.tsv.gz"
        with gzip.open(path, "wt") as handle:
            handle.write("#fractional\ttrue\n")
        with self.assertRaisesRegex(ValueError, "non-fractional"):
            m.full_index(path.parent)

    def test_missing_database_rejected(self):
        self.fixture_manifests()
        db = self.args.baseline / "db_patho"
        (db / "metax_db.mxi").unlink()
        with self.assertRaisesRegex(ValueError, "Missing/empty"):
            m.full_index(db)

    def test_wrong_index_mode_stops_even_with_fixture_manifest_pin(self):
        manifests, pins = self.fixture_manifests()
        manifests["build_metax.json"]["index_mode"] = "fractional"
        m.write_json(self.args.baseline / "build_metax.json", manifests["build_metax.json"])
        pins["build_metax.json"] = m.canonical_sha(manifests["build_metax.json"])
        with mock.patch.dict(m.MANIFEST_PINS, pins, clear=True):
            with self.assertRaisesRegex(ValueError, "index mode"):
                m.load_manifests(self.args)

    def test_version_parsers(self):
        cases = {"python":("Python 3.14.7","3.14.7"),
                 "nextflow":("N E X T F L O W\n version 26.10.5 build 1","26.10.5"),
                 "java":('openjdk version "21.0.9" 2025',"21"),
                 "samtools":("samtools 1.24\nUsing htslib 1.24","1.24"),
                 "metax":("metax 2.0.3\n","2.0.3")}
        for name, (text,want) in cases.items():
            self.assertEqual(m.version_value(name,text),want)
        with self.assertRaisesRegex(ValueError, "Cannot identify"):
            m.version_value("metax", "unknown")

    def test_changed_tool_version_stops(self):
        (self.args.out / "logs").mkdir(parents=True)
        self.put(self.args.baseline / "logs/version_python.log", "Python 3.14.7\n")
        with mock.patch.object(m.shutil,"which",return_value=sys.executable), \
             mock.patch.object(m,"capture",return_value="Python 3.14.8\n"):
            with self.assertRaisesRegex(ValueError, "version changed"):
                m.verify_tools(self.args, {})

    def test_argv_safety_and_independent_directories(self):
        prep = {"hisat2_index":"/valid/index","target_panel":"/valid/panel"}
        pe, se = [m.pipeline_command(self.args,prep,k) for k in ("pe","se")]
        for k, command in (("pe",pe),("se",se)):
            self.assertNotIn("-resume",command)
            self.assertNotIn("shell=True",command)
            self.assertEqual(command[command.index("-work-dir")+1],str(self.args.out / f"work_{k}"))
            self.assertEqual(command[command.index("--outdir")+1],str(self.args.out / f"results_{k}"))
            self.assertEqual(command[command.index("-profile")+1],"standard")
            self.assertEqual(command[command.index("--kraken_mem")+1],"2 GB")
            self.assertIn("-C",command)
        with self.assertRaises(ValueError):
            m.pipeline_command(self.args,prep,"pe;touch_exploit")
        with self.assertRaisesRegex(ValueError, "Unsafe"):
            m.safe_path(self.root / "semicolon;or space")

    def test_metax_symlink_does_not_replace_python(self):
        self.args.out.mkdir(parents=True)
        binary = self.put(self.root / "other_env/bin/metax")
        original = os.environ.get("PATH","")
        env = m.metax_environment(self.args, {"metax":str(binary)})
        self.assertEqual(env["PATH"].split(os.pathsep,1)[1],original)
        self.assertEqual((self.args.out / "toolbin/metax").resolve(),binary)
        self.assertFalse((self.args.out / "toolbin/python").exists())

    def test_existing_output_untouched(self):
        self.args.out.mkdir(parents=True)
        sentinel = self.put(self.args.out / "sentinel")
        before = sentinel.read_bytes()
        with self.assertRaisesRegex(ValueError,"Output exists"):
            m.execute(self.args)
        self.assertEqual(sentinel.read_bytes(),before)
        self.assertFalse((self.args.out / "provenance.json").exists())

    def test_failure_creates_diagnostics_without_running_pipeline(self):
        with mock.patch.object(m,"verify_sources",side_effect=ValueError("source changed")), \
             mock.patch.object(m,"run_kind") as launch, contextlib.redirect_stdout(io.StringIO()), \
             contextlib.redirect_stderr(io.StringIO()):
            self.assertEqual(m.execute(self.args),2)
            launch.assert_not_called()
        self.assertEqual(m.read_json(self.args.out / "provenance.json")["state"],"stopped")
        self.assertIn("source changed",(self.args.out / "review_minimal.txt").read_text())

    def test_collector_excludes_large_or_sensitive_artifacts(self):
        (self.args.out / "logs").mkdir(parents=True)
        self.put(self.args.out / "provenance.json",'{"software_test": true}\n')
        self.put(self.args.out / "logs/pipeline_pe.log","line\n"*40000)
        self.put(self.args.out / "results_pe/S/validation/partitions/dehost.bam","FORBIDDEN_BAM")
        self.put(self.args.out / "results_pe/S/validation/trimmed/R1.fastq","FORBIDDEN_FASTQ")
        self.put(self.args.out / "logs/install.log","FORBIDDEN_INSTALL")
        self.put(self.args.out / "work_pe/aa/bb/.command.run","FORBIDDEN_ENV")
        self.put(self.args.out / "work_pe/aa/bb/.command.err","meaningful error\n")
        text = m.collect(self.args).read_text()
        for token in ("FORBIDDEN_BAM","FORBIDDEN_FASTQ","FORBIDDEN_INSTALL","FORBIDDEN_ENV"):
            self.assertNotIn(token,text)
        self.assertIn("meaningful error",text)
        self.assertLess(len(text.encode()),m.MAX_BUNDLE+128*1024)
        self.assertIn("TRUNCATED",text)

    def test_command_failure_no_shell_and_diagnostics(self):
        (self.args.out / "logs").mkdir(parents=True)
        prep = {"hisat2_index":"/valid/index","target_panel":"/valid/panel"}
        with mock.patch.object(m.subprocess,"run",return_value=subprocess.CompletedProcess([],9)) as run, \
             contextlib.redirect_stdout(io.StringIO()):
            with self.assertRaisesRegex(ValueError,"exit=9"):
                m.run_kind(self.args,prep,"pe",{})
        self.assertNotIn("shell",run.call_args.kwargs)
        self.assertIsInstance(run.call_args.args[0],list)
        record = m.read_json(self.args.out / "run_pe.json")
        self.assertEqual(record["exit_code"],9)
        self.assertEqual(record["state"],"stopped")
        self.assertTrue((self.args.out / "launch_pe").is_dir())
        self.assertFalse((self.args.out / "launch_se").exists())

    def test_artifact_check_accepts_qc_fallback_without_claiming_accuracy(self):
        kind, sample = "se", "SAMOVAR_SE"
        for fraction in ("dehost","host"):
            for role in ("full","patho"):
                path = self.args.out / f"results_se/{sample}/profiles/{fraction}.{role}.profile.json"
                self.put(path,json.dumps({"annotator":"metax","annotator_version":"2.0.3",
                                         "sample_id":sample,"fraction":fraction,"role":role}))
        self.put(self.args.out / f"results_se/{sample}/qc/{sample}_qc_report.json",
                 json.dumps({"sample_id":sample,"status":"FALLBACK","fallback_reasons":["DEHOST_FAILURE"]}))
        self.put(self.args.out / "logs/se.trace.tsv","status\texit\n"+"COMPLETED\t0\n"*17)
        result = m.check_outputs(self.args,kind)
        self.assertEqual(result["qc_status"],"FALLBACK")
        self.assertIn("NOT asserted",result["biological_accuracy"])

    def test_check_only_never_runs_workflow(self):
        prep = {"metax":"/fixture/metax","hisat2_index":"/fixture/index","target_panel":"/fixture/panel"}
        with mock.patch.object(m,"verify_sources",return_value={}), \
             mock.patch.object(m,"load_manifests",return_value={"prepare.json":prep}), \
             mock.patch.object(m,"verify_tools",return_value={}), \
             mock.patch.object(m,"metax_environment",return_value={}), \
             mock.patch.object(m,"capture",return_value="software-test\n"), \
             mock.patch.object(m,"run_kind") as run, contextlib.redirect_stdout(io.StringIO()):
            self.assertEqual(m.execute(self.args,check_only=True),0)
            run.assert_not_called()
        self.assertEqual(m.read_json(self.args.out / "provenance.json")["state"],"preflight_passed")

    def test_execution_order_and_stop_after_pe_failure(self):
        prep = {"metax":"/fixture/metax","hisat2_index":"/fixture/index","target_panel":"/fixture/panel"}
        with mock.patch.object(m,"verify_sources",return_value={}), \
             mock.patch.object(m,"load_manifests",return_value={"prepare.json":prep}), \
             mock.patch.object(m,"verify_tools",return_value={}), \
             mock.patch.object(m,"metax_environment",return_value={}), \
             mock.patch.object(m,"capture",return_value="software-test\n"), \
             mock.patch.object(m,"run_kind",side_effect=ValueError("PE execution failed")) as run, \
             contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
            self.assertEqual(m.execute(self.args),2)
            self.assertEqual(run.call_count,1)
            self.assertEqual(run.call_args.args[2],"pe")
        self.assertIn("PE execution failed",(self.args.out / "review_minimal.txt").read_text())

    def test_collector_fastp_summary_without_full_curves(self):
        (self.args.out / "logs").mkdir(parents=True)
        self.put(self.args.out / "results_pe/S/qc/fastp/S.fastp.json",
                 json.dumps({"summary":{"before_filtering":{"total_reads":4}},
                             "filtering_result":{"passed_filter_reads":2},
                             "read1_before_filtering":{"quality_curves":"EXCLUDED_CURVE"}}))
        bundle=m.collect(self.args).read_text()
        self.assertIn("passed_filter_reads",bundle)
        self.assertNotIn("EXCLUDED_CURVE",bundle)


if __name__ == "__main__":
    unittest.main()
