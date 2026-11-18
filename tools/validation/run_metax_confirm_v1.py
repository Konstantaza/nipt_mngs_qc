#!/usr/bin/env python3
"""Confirm the reviewed nipt_mngs_qc 1.10.0 MetaX workflow, without changing it.

Python >=3.10, standard library only. Run from the repository in its ORIGINAL
server location, using the existing metax_samovar_v1 preparation and databases.
No install/build/resume/overwrite. Added validation tools/docs are allowed;
changed runtime sources or manifests require a separate review, not new pins.
Exit 0 = execution/artifact checks passed, NOT QC PASS or biological accuracy.
"""
import argparse
import csv
import gzip
import hashlib
import io
import json
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import sys
import time

VERSION = "1"
MANIFEST_PINS = {
    "prepare.json": "176fd4b344f473ff6d4709889cde570ef7cb6c454b256f45121dbfd4d4ef4cef",
    "build_metax.json": "6e1fec06fb63c3bb8f89819656a7ee3d4d46a9f13f901fc3c045804c83d16821",
}
SOURCE_PINS = {
    "main.nf": "db2f892ec8a0d562ad46b2e5e6414723692009f2f1afa8449ab1d12401c92749",
    "nextflow.config": "a8b89dcd17a729c181991e1ca64bf5cdfd8c1f4a0c224f90b400e3c8a12a6f4b",
    "multiqc_config.yaml": "0d4d0d295fe1909df2db18166d766c4fe52510347856b123d9ebe5631b25fa27",
    "modules/local/annotate_kraken2_bracken.nf": "96a38c33d0f1b7128f42a9f33b84fdf70be33986e6bc969a879c275be88747f3",
    "modules/local/annotate_metax.nf": "2f768f700dae47d680f68b07116e76bd82c73fc76d23d29a7eb20fb3bef11c20",
    "modules/local/annotate_stub.nf": "6bc2d88bdf17bb552956f6104275e0b5616399ece8e3783707d1c0f558127067",
    "modules/local/bam2fasta.nf": "e4ca6064e19427fa088f1971558b23e2721eebbf4f4cf7984f332b1c039d7941",
    "modules/local/compare.nf": "d22281c7b42beb33d3f39ffef27972a60de65080235078e10140cc90ba892db6",
    "modules/local/dehost.nf": "1f5c50ceb272d98939b46433f2786e9bc560f26c92c9b082189e2ad4f045c921",
    "modules/local/fastp.nf": "ec99010932c379cd00c633ca6b133357d012d2763f4f7d8e811328cb3642ca0f",
    "modules/local/make_all.nf": "2404f07898159790fbe27d33bbe05af7f6ab3ec59e3c8d80ddbc269bde1fdc78",
    "modules/local/multiqc.nf": "f181a556035ef24fd45561019f76ffd44cc2a0c7c4cfeaa55edb5898332c8e79",
    "modules/local/qc_gate.nf": "cc097f6608315812b96a32ef54d6d5b0fa74a2527f8cabb53c6981f2058069b3",
    "modules/local/validate_profile.nf": "d02c80b9d96077dda84425973b91d93ebe94b1de1c43a45517a402f6d7dad806",
    "subworkflows/local/annotate.nf": "ea66f39fa2bb18ebe8db9a1d5a2c76f77a703d935eef63fa6b7b67ca7c25d571",
    "bin/bam_counts.py": "1f7fea4e87764bb56d82f4936f82915d167e6fa2682679ae63b4b9e98c399a13",
    "bin/compare_profiles.py": "5157a53dd3c38301f0c47e2a7f9cdbc8799d98f8c52aee29527e6ecf5d68b712",
    "bin/dehost_split.sh": "65f019481de31a08db32f1f87f3ac3df5b532e7018c49b1ceaa169d9bb0b84ea",
    "bin/kraken2_bracken_to_profile.py": "039ee58a4f03af5473145ba4a41e91594c40e2d44f819ca6db7bf836c16a2138",
    "bin/metax_prepare.py": "2f8c91ecbea4a75f76b1cdd805c2f636dca0dc4adf309e097c236962f7195ff6",
    "bin/metax_to_profile.py": "92b5432a70e2c6aa90b7f653915a6bc47a77a4b2ff659edf92e0ac91ff84389b",
    "bin/profile_io.py": "1b0bc15eef42ad1ffaf7261afdb63d5c004d7057112ebd118add01ee295f829e",
    "bin/qc_gate.py": "ed660ce3fb4bd6a0930b1d9efbcde262fcadf4c7c34871dc4f4b69745c1e27ec",
    "bin/stub_to_profile.py": "58c4e6e5014d00f2e24c91afca45f6fc4e8f48e16de4faf7f01609f23803f98d",
    "bin/sum_counts.py": "a4dd651a36005187231e592aee7a5a2544155b24a98a6a109219f4c170354d00",
    "bin/validate_profile.py": "ab11490f897d69876a1a6b2aa4375894628ebb2371ebfaca8912c83cae8ecf8b",
    "assets/missing.txt": "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855",
    "assets/no_fastp.json": "9b6538038a97185201a9a43c2322c357b832adbf4d1a9aee38c4c7c200312934",
    "assets/no_panel.txt": "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855",
    "tests/insilico/toy/toy.config": "33b60d775014513771d25c1d1b2d6827b4da0dc73c3f9f3cd75b5fc81432377a",
    "tests/insilico/toy/target_panel.txt": "fc3272260bd8afadf76924e12ea4c70e50ecc74ee8d68e85d14187c95cb151c5",
    "tests/insilico/toy/taxonomy/nodes.dmp": "08ec7e1208f43e4269dcf94745ddc8b8862f3d943004d50f4883e29ebe540831",
    "tests/insilico/toy/taxonomy/names.dmp": "a7a3265f379e560e88fb0a56c83388dcb9a620dc0b7a8be909e271fcd0d60af8",
}
SAFE = re.compile(r"^[A-Za-z0-9_./:+-]+$")
MAX_FILE = 128 * 1024
MAX_BUNDLE = 2 * 1024 * 1024


def require(ok, message):
    if not ok:
        raise ValueError(message)


def safe_path(value):
    path = Path(value).absolute()
    require(not path.is_symlink(), f"Symlink root not allowed: {path}")
    path = path.resolve()
    # The runner uses argv, but the reviewed Nextflow shell tasks interpolate paths.
    require(SAFE.fullmatch(str(path)), f"Unsafe task-interpolated path: {path}")
    return path


def sha(path):
    h = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def canonical_sha(data):
    return hashlib.sha256(json.dumps(data, sort_keys=True, separators=(",", ":"),
                                    ensure_ascii=False).encode()).hexdigest()


def write_json(path, data):
    Path(path).write_text(json.dumps(data, ensure_ascii=False, indent=2) + "\n")


def read_json(path):
    return json.loads(Path(path).read_text())


def verify_hashes(mapping):
    for name, expected in mapping.items():
        path = Path(name)
        require(path.is_file(), f"Missing pinned file: {path}")
        require(sha(path) == expected, f"SHA256 mismatch: {path}; STOP, review changes")


def verify_sources(repo):
    verify_hashes({repo / name: digest for name, digest in SOURCE_PINS.items()})
    actual = {str(p.relative_to(repo)) for folder in ("modules", "subworkflows")
              for p in (repo / folder).rglob("*.nf")}
    actual |= {str(p.relative_to(repo)) for p in (repo / "bin").iterdir() if p.is_file()}
    expected = {n for n in SOURCE_PINS if n.startswith(("modules/", "subworkflows/", "bin/"))}
    require(actual == expected, f"Runtime file set differs: extra={actual-expected}, missing={expected-actual}")
    return dict(SOURCE_PINS)


def full_index(db):
    for name in ("metax_db.mxi", "metax_db.tsv.gz"):
        path = db / name
        require(path.is_file() and path.stat().st_size > 0, f"Missing/empty DB file: {path}")
    with gzip.open(db / "metax_db.tsv.gz", "rt") as handle:
        fractional = next((line.strip() for line in handle if line.startswith("#fractional\t")), None)
    require(fractional == "#fractional\tfalse", f"Full non-fractional index required: {db}")
    for name in ("nodes.dmp", "names.dmp", "merged.dmp", "delnodes.dmp"):
        require((db / "taxonomy" / name).is_file(), f"Missing taxonomy: {db}/{name}")


def load_manifests(args):
    manifests = {}
    for name, digest in MANIFEST_PINS.items():
        data = read_json(args.baseline / name)
        require(canonical_sha(data) == digest, f"Reviewed manifest mismatch: {name}")
        manifests[name] = data
    prep, build = manifests["prepare.json"], manifests["build_metax.json"]
    require(prep.get("prepare_completed") is True and
            prep.get("human_index_sequences_verified") is True, "Incomplete original preparation")
    require(prep["repo"] == str(args.repo) and prep["toy"] == str(args.toy) and
            prep["out"] == str(args.baseline),
            "Original server locations required; do not rewrite prepared manifests to relocate")
    require(build.get("index_mode") == "full" and build.get("fractional") is False,
            "Wrong MetaX index mode")
    verify_hashes(prep["input_sha256"])
    verify_hashes(prep["prepared_sha256"])
    verify_hashes(build["files_sha256"])
    verify_hashes({prep["metax"]: prep["metax_sha256"]})
    for role in ("full", "patho"):
        full_index(args.baseline / f"db_{role}")
    # Manifest hashes bind original reads, eight index pieces, taxonomy and resources.
    # Added runtime .metamcache files are legitimate and intentionally not rejected.
    index = safe_path(prep["hisat2_index"])
    parts = [Path(f"{index}.{i}.ht2") for i in range(1, 9)]
    require(all(str(p) in prep["input_sha256"] for p in parts), "Unbound human index")
    with open(args.baseline / "dbsheet_metax.csv", newline="") as handle:
        rows = list(csv.DictReader(handle))
    require(len(rows) == 2 and {r["role"] for r in rows} == {"full", "patho"}, "Wrong DB roles")
    for row in rows:
        require(row["annotator"] == "metax" and
                Path(row["path"]).resolve() == args.baseline / f'db_{row["role"]}',
                "Wrong annotator or DB location")
    for kind in ("pe", "se"):
        with open(args.toy / f"samplesheet_{kind}.csv", newline="") as handle:
            rows = list(csv.DictReader(handle))
        require(len(rows) == 1 and rows[0]["sample"] == f"SAMOVAR_{kind.upper()}",
                f"Unexpected sample: {kind}")
    return manifests


TOOL_SPECS = {
    "python": (["--version"], r"Python (\d+\.\d+\.\d+)"),
    "nextflow": (["-version"], r"version (\d+\.\d+\.\d+)"),
    "java": (["-version"], r'version "(\d+)'),
    "fastp": (["--version"], r"fastp (\d+\.\d+\.\d+)"),
    "hisat2": (["--version"], r"version (\d+\.\d+\.\d+)"),
    "samtools": (["--version"], r"samtools (\d+\.\d+(?:\.\d+)?)"),
    "multiqc": (["--version"], r"version (\d+\.\d+(?:\.\d+)?)"),
    "metax": (["--version"], r"^metax (\d+\.\d+\.\d+)\s*$"),
}


def version_value(name, text):
    match = re.search(TOOL_SPECS[name][1], text, flags=re.M)
    require(match, f"Cannot identify {name} version")
    return match.group(1)


def capture(command, cwd=None):
    done = subprocess.run([str(x) for x in command], cwd=cwd, stdout=subprocess.PIPE,
                          stderr=subprocess.STDOUT, timeout=180)
    text = done.stdout.decode("utf-8", errors="replace")
    require(done.returncode == 0, f"Version/Git command failed: {shlex.join(map(str,command))}\n{text[-4000:]}")
    return text


def verify_tools(args, prep):
    require(sys.version_info >= (3, 10), "Python >=3.10 required")
    tools = {}
    for name, (flags, _) in TOOL_SPECS.items():
        executable = prep["metax"] if name == "metax" else shutil.which(name)
        require(executable, f"Missing main-environment tool: {name}")
        text = capture([executable, *flags])
        (args.out / "logs" / f"version_{name}.log").write_text(text)
        observed = version_value(name, text)
        previous = version_value(name, (args.baseline / "logs" / f"version_{name}.log").read_text(errors="replace"))
        require(observed == previous, f"{name} version changed: {previous} -> {observed}")
        tools[name] = {"executable": str(executable), "version": observed}
    require(tools["metax"]["version"] == "2.0.3", "MetaX 2.0.3 required")
    require(tools["hisat2"]["version"] == "2.2.3", "HISAT2 2.2.3 required")
    require(tools["fastp"]["version"] == "1.3.7", "fastp 1.3.7 required")
    require(tools["python"]["version"] == ".".join(map(str, sys.version_info[:3])),
            "Invoked Python differs from main PATH Python")
    for name in ("python", "python3"):
        executable = shutil.which(name)
        require(executable and Path(executable).resolve() == Path(sys.executable).resolve(),
                f"{name} does not use the main Python; deactivate MetaX environment")
    for name in ("bash", "awk", "grep", "git"):
        require(shutil.which(name), f"Missing tool: {name}")
    return tools


def config_arg(args):
    # -C ignores ambient ~/.nextflow/config and launch-directory config.
    # The SAME reviewed project+execution configs are used, not new thresholds.
    return f"{args.repo / 'nextflow.config'},{args.baseline / 'execution.config'}"


def pipeline_command(args, prep, kind):
    require(kind in ("pe", "se"), "Unknown run kind")
    return ["nextflow", "-C", config_arg(args), "-log",
            str(args.out / "logs" / f"{kind}.nextflow.log"),
            "run", str(args.repo / "main.nf"), "-profile", "standard",
            "-ansi-log", "false", "-work-dir", str(args.out / f"work_{kind}"),
            "-with-trace", str(args.out / "logs" / f"{kind}.trace.tsv"),
            "--mode", "fastq", "--input", str(args.toy / f"samplesheet_{kind}.csv"),
            "--dbsheet", str(args.baseline / "dbsheet_metax.csv"),
            "--host_index", prep["hisat2_index"], "--target_panel", prep["target_panel"],
            "--bracken_read_len", "126" if kind == "pe" else "50",
            "--kraken_mem", "2 GB", "--outdir", str(args.out / f"results_{kind}")]


def metax_environment(args, prep):
    toolbin = args.out / "toolbin"
    toolbin.mkdir()
    (toolbin / "metax").symlink_to(prep["metax"])
    env = os.environ.copy()
    env["PATH"] = str(toolbin) + os.pathsep + env.get("PATH", "")
    return env


def check_outputs(args, kind):
    root = args.out / f"results_{kind}"
    sample = f"SAMOVAR_{kind.upper()}"
    profiles = list(root.rglob("*.profile.json"))
    require(len(profiles) == 4, f"{kind}: expected four profiles, got {len(profiles)}")
    actual = set()
    for path in profiles:
        data = read_json(path)
        require(data["annotator"] == "metax" and data["annotator_version"] == "2.0.3"
                and data["sample_id"] == sample, f"Wrong profile identity: {path}")
        actual.add((data["fraction"], data["role"]))
    require(actual == {(f,r) for f in ("host","dehost") for r in ("full","patho")},
            f"Incomplete profile grid: {kind}")
    trace = args.out / "logs" / f"{kind}.trace.tsv"
    rows = list(csv.DictReader(io.StringIO(trace.read_text()), delimiter="\t"))
    require(len(rows) == 17 and all(r["status"] == "COMPLETED" and r["exit"] == "0" for r in rows),
            f"{kind}: expected 17 completed, exit-zero tasks")
    qc = read_json(root / sample / "qc" / f"{sample}_qc_report.json")
    require(qc["sample_id"] == sample and qc["status"] in
            ("FAIL", "FALLBACK", "PASS_POSITIVE", "PASS_NEGATIVE"), "Invalid sample QC artifact")
    return {"workflow_execution": "completed", "tasks": len(rows),
            "qc_status": qc["status"], "qc_reasons": qc.get("fallback_reasons", []),
            "biological_accuracy": "NOT asserted; independent GT review required"}


def run_kind(args, prep, kind, env):
    # Stop after any execution/artifact error. A valid FALLBACK remains a valid
    # scientific result and does not pretend to be a shell execution failure.
    launch = args.out / f"launch_{kind}"
    launch.mkdir()
    command = pipeline_command(args, prep, kind)
    record = {"state": "running", "command": command, "cwd": str(launch)}
    write_json(args.out / f"run_{kind}.json", record)
    print("$ " + shlex.join(command), flush=True)
    start = time.monotonic()
    try:
        with open(args.out / "logs" / f"pipeline_{kind}.log", "w") as log:
            result = subprocess.run(command, cwd=launch, env=env, stdout=log, stderr=subprocess.STDOUT)
        record.update(exit_code=result.returncode, wall_seconds=time.monotonic()-start)
        require(result.returncode == 0, f"Nextflow {kind} failed, exit={result.returncode}; inspect logs")
        record.update(check_outputs(args, kind))
        record["state"] = "completed"
    except BaseException as exc:
        record.update(state="stopped", error=f"{type(exc).__name__}: {exc}")
        raise
    finally:
        write_json(args.out / f"run_{kind}.json", record)
    print(f'{kind.upper()}: workflow completed; sample QC={record["qc_status"]}; accuracy not asserted.', flush=True)
    return record


def bounded_text(path, tail=False):
    with open(path, "rb") as handle:
        size = path.stat().st_size
        if tail:
            handle.seek(max(0, size-MAX_FILE))
        content = handle.read(MAX_FILE)
    text = content.decode("utf-8", errors="replace")
    if tail:
        text = "\n".join(text.splitlines()[-120:])
    if size > MAX_FILE:
        text = f"[TRUNCATED {'TAIL' if tail else 'HEAD'}; original bytes={size}]\n" + text
    return text


def collect(args):
    for kind in ("pe", "se"):
        root = args.out / f"results_{kind}"
        summaries = []
        for path in sorted(root.glob("*/qc/fastp/*.fastp.json")):
            try:
                data = read_json(path)
                summaries.append({"file": str(path.relative_to(args.out)), "sha256": sha(path),
                                  "summary": data.get("summary"),
                                  "filtering_result": data.get("filtering_result")})
            except (OSError, ValueError) as exc:
                summaries.append({"file": str(path), "collection_error": str(exc)})
        if summaries:
            write_json(args.out / f"fastp_summary_{kind}.json", summaries)
    candidates = [(p, False) for p in sorted(args.out.glob("*.json"))]
    for name in ("prepare.json", "build_metax.json", "reference_manifest.json",
                 "execution.config", "dbsheet_metax.csv"):
        candidates.append((args.baseline / name, False))
    for kind in ("pe", "se"):
        for name in (f"{kind}.trace.tsv", f"{kind}.nextflow.log", f"pipeline_{kind}.log"):
            candidates.append((args.out / "logs" / name, name.endswith(".log")))
        root = args.out / f"results_{kind}"
        for pattern in ("*.counts.json", "*.profile.json", "*.profile.tsv", "*.metax.txt",
                        "*.compare.json", "*.compare.tsv", "*_qc_report.json", "*_qc_mqc.tsv",
                        "*.hisat2.summary.txt"):
            candidates.extend((p, False) for p in sorted(root.rglob(pattern)))
        # Only error/output tails, never .command.run (may contain environment).
        for pattern in (".command.err", ".command.out", ".exitcode"):
            candidates.extend((p, True) for p in sorted((args.out / f"work_{kind}").glob(f"*/*/{pattern}"))[:40])
    candidates.extend((p, True) for p in sorted((args.out / "logs").glob("version_*.log")))
    candidates.append((args.out / "logs" / "resolved.config", False))
    chunks = ["MetaX confirmation v1. Execution, sample QC and GT accuracy are distinct.\n"
              "Excludes reads, alignments, indices/cache, installation logs and environment dumps.\n"]
    receipt, seen, used = [], set(), 0
    for path, tail in candidates:
        if path in seen or not path.is_file() or path.is_symlink():
            continue
        seen.add(path)
        label = ("baseline/" + str(path.relative_to(args.baseline))) if path.is_relative_to(args.baseline) else str(path.relative_to(args.out))
        info = {"file": label, "bytes": path.stat().st_size, "sha256": sha(path)}
        content = bounded_text(path, tail)
        block = f"\n===== {label} =====\n{content.rstrip()}\n"
        size = len(block.encode())
        info["included"] = used + size <= MAX_BUNDLE
        if info["included"]:
            chunks.append(block)
            used += size
        receipt.append(info)
    write_json(args.out / "collection.json", {"files": receipt, "content_limit_bytes": MAX_BUNDLE})
    chunks.append("\n===== collection.json =====\n" +
                  json.dumps({"files": receipt}, ensure_ascii=False, indent=2) + "\n")
    destination = args.out / "review_minimal.txt"
    destination.write_text("".join(chunks))
    return destination


def execute(args, check_only=False):
    require(not args.out.exists(), f"Output exists; refusing overwrite: {args.out}")
    for old in (args.repo, args.toy, args.baseline):
        require(not old.is_relative_to(args.out), f"Output would contain existing inputs: {args.out}")
    require(not args.out.is_relative_to(args.toy) and not args.out.is_relative_to(args.baseline),
            "Output must not be inside old inputs/results")
    args.out.mkdir(parents=True, exist_ok=False)
    (args.out / "logs").mkdir()
    state = {"runner_version": VERSION, "state": "preflight", "check_only": check_only,
             "repo": str(args.repo), "out": str(args.out), "main_python": sys.executable,
             "runner_sha256": sha(__file__), "started_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
             "meaning": "Execution != QC PASS != GT accuracy != clinical diagnosis"}
    try:
        # Capture Git even when a subsequent source mismatch stops the run.
        # A dirty tree is recorded, not silently cleaned or auto-committed.
        try:
            state["git_commit"] = capture(["git", "-C", args.repo, "rev-parse", "HEAD"]).strip()
            state["git_status"] = capture(["git", "-C", args.repo, "status", "--short", "--branch"])
        except Exception as exc:
            state["git_error"] = f"{type(exc).__name__}: {exc}"
        state["source_sha256"] = verify_sources(args.repo)
        manifests = load_manifests(args)
        prep = manifests["prepare.json"]
        state["manifests_canonical_sha256"] = {k: canonical_sha(v) for k,v in manifests.items()}
        state["tools"] = verify_tools(args, prep)
        require("git_commit" in state and "git_status" in state,
                f"Cannot record Git provenance: {state.get('git_error')}")
        env = metax_environment(args, prep)
        resolved = capture(["nextflow", "-C", config_arg(args), "config", args.repo,
                            "-profile", "standard", "-flat"], cwd=args.out)
        (args.out / "logs" / "resolved.config").write_text(resolved)
        state["commands"] = {k: pipeline_command(args, prep, k) for k in ("pe", "se")}
        state["state"] = "preflight_passed"
        write_json(args.out / "provenance.json", state)
        if not check_only:
            state["runs"] = {}
            for kind in ("pe", "se"):
                state["runs"][kind] = run_kind(args, prep, kind, env)
            state["state"] = "execution_completed"
        return 0
    except (Exception, KeyboardInterrupt) as exc:
        state.update(state="stopped", error=f"{type(exc).__name__}: {exc}")
        print(f"STOP: {state['error']}", file=sys.stderr)
        return 2
    finally:
        write_json(args.out / "provenance.json", state)
        bundle = collect(args)
        print(f"Send only: {bundle}\nbytes={bundle.stat().st_size}; sha256={sha(bundle)}", flush=True)


def cli(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, default=Path.cwd())
    parser.add_argument("--toy", type=Path)
    parser.add_argument("--baseline", type=Path)
    parser.add_argument("--out", type=Path)
    parser.add_argument("--check-only", action="store_true",
                        help="Validate and resolve config only; consumes a NEW output directory")
    args = parser.parse_args(argv)
    try:
        args.repo = safe_path(args.repo)
        args.toy = safe_path(args.toy or args.repo / "insilico_runs/toy_995")
        args.baseline = safe_path(args.baseline or args.repo / "insilico_runs/metax_samovar_v1")
        args.out = safe_path(args.out or args.repo / "insilico_runs/metax_confirm_v2")
        return execute(args, args.check_only)
    except (OSError, ValueError) as exc:
        print(f"STOP before creating output: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(cli())
