#!/usr/bin/env python3
"""Pinned MetaX integration kit; standard library only, Linux, Python >=3.10."""
import argparse
import csv
import gzip
import hashlib
import io
import json
import os
from pathlib import Path
import re
import resource
import shlex
import shutil
import subprocess
import sys
import tarfile
import time

VERSION = "1"
THREADS = 4
SOURCE_SPECIES = {"9606": "9606", "2886930": "10847", "4932": "4932", "562": "562"}
SAFE = re.compile(r"^[A-Za-z0-9_./:+-]+$")


def require(condition, message):
    if not condition:
        raise ValueError(message)


def safe_path(path):
    path = Path(path).resolve()
    require(bool(SAFE.fullmatch(str(path))), f"Unsafe shell-interpolated path: {path}")
    return path


def file_ok(path):
    require(Path(path).is_file() and Path(path).stat().st_size > 0,
            f"Missing/empty file: {path}")
    return Path(path)


def sha(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def write_json(path, data):
    Path(path).write_text(json.dumps(data, ensure_ascii=False, indent=2) + "\n")


def read_json(path):
    return json.loads(Path(path).read_text())


def csv_write(path, fields, rows):
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)


def records(path):
    header, chunks = None, []
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                if header is not None:
                    require(bool(chunks), f"Empty sequence: {header}")
                    yield header, "".join(chunks)
                header, chunks = line[1:].strip(), []
                require(bool(header), f"Empty FASTA header: {path}")
            elif line.strip():
                require(header is not None, f"Sequence before header: {path}")
                sequence = line.strip()
                require(bool(re.fullmatch("[ACGTRYSWKMBDHVNacgtryswkmbdhvn]+", sequence)),
                        f"Non-IUPAC DNA: {path}")
                chunks.append(sequence)
    if header is not None:
        require(bool(chunks), f"Empty sequence: {header}")
        yield header, "".join(chunks)


def seq_map(path):
    result = {}
    for header, seq in records(path):
        name = header.split()[0]
        require(name not in result, f"Duplicate FASTA contig: {name}")
        result[name] = (len(seq), hashlib.sha256(seq.upper().encode()).hexdigest())
    require(bool(result), f"No FASTA records: {path}")
    return result


def taxonomy_check(tax):
    nodes, names = {}, set()
    for line in file_ok(tax / "nodes.dmp").read_text().splitlines():
        fields = [x.strip() for x in line.split("|")]
        require(len(fields) >= 3 and fields[0] not in nodes, "Invalid/duplicate taxon")
        nodes[fields[0]] = (fields[1], fields[2])
    for line in file_ok(tax / "names.dmp").read_text().splitlines():
        fields = [x.strip() for x in line.split("|")]
        if len(fields) >= 4 and fields[3] == "scientific name":
            names.add(fields[0])
    for taxid in nodes:
        require(taxid in names, f"No scientific name: {taxid}")
        current, seen = taxid, set()
        while current != "1":
            require(current not in seen and current in nodes, f"Broken lineage: {taxid}")
            seen.add(current)
            current = nodes[current][0]
    require(nodes.get("1") == ("1", "no rank"), "Invalid taxonomy root")
    for taxid in SOURCE_SPECIES.values():
        require(nodes.get(taxid, (None, None))[1] == "species", f"Not species: {taxid}")
    require(nodes.get("2886930") == ("10847", "no rank"), "Unexpected phiX lineage")
    return {"nodes": len(nodes), "species": sorted(set(SOURCE_SPECIES.values())),
            "phiX_source_to_species": {"2886930": "10847"}}


def reheader(refs, out):
    """Change headers only. Preserve sequence case, order and bases."""
    manifest = []
    with open(out / "reference_full.fasta", "w") as full, \
            open(out / "reference_patho.fasta", "w") as patho, \
            open(out / "reference_pilot.fasta", "w") as pilot:
        for source, path in refs.items():
            seen, n = set(), 0
            for n, (header, seq) in enumerate(records(path), 1):
                original_id = header.split()[0]
                require(original_id not in seen, f"Duplicate contig in {path}: {original_id}")
                seen.add(original_id)
                genome, contig = f"samovar_{source}", f"contig_{n:06d}"
                label = f"{genome}|{contig}|s:{SOURCE_SPECIES[source]}"
                text = f">{label}\n" + "\n".join(seq[i:i + 80] for i in range(0, len(seq), 80)) + "\n"
                full.write(text)
                if source != "9606":
                    patho.write(text)
                if source == "562":
                    pilot.write(text)
                manifest.append({"source_taxid": source, "species_taxid": SOURCE_SPECIES[source],
                                 "source_fasta": str(path), "source_header": header,
                                 "genome_id": genome, "contig_id": contig, "bases": len(seq),
                                 "sequence_sha256_case_sensitive": hashlib.sha256(seq.encode()).hexdigest()})
            require(n > 0, f"Empty reference: {path}")
    write_json(out / "reference_manifest.json", manifest)
    # Independent roundtrip check against source sequence hashes.
    for kind, expected in (("full", manifest),
                           ("patho", [m for m in manifest if m["source_taxid"] != "9606"]),
                           ("pilot", [m for m in manifest if m["source_taxid"] == "562"])):
        observed = list(records(out / f"reference_{kind}.fasta"))
        require(len(observed) == len(expected), f"Roundtrip record mismatch: {kind}")
        for (header, sequence), row in zip(observed, expected):
            require(header == f'{row["genome_id"]}|{row["contig_id"]}|s:{row["species_taxid"]}'
                    and len(sequence) == row["bases"]
                    and hashlib.sha256(sequence.encode()).hexdigest()
                    == row["sequence_sha256_case_sensitive"], f"Roundtrip mismatch: {kind}")
    return manifest


def full_index(db):
    file_ok(db / "metax_db.mxi")
    file_ok(db / "metax_db.tsv.gz")
    with gzip.open(db / "metax_db.tsv.gz", "rt") as f:
        value = next((line.strip() for line in f if line.startswith("#fractional\t")), None)
    require(value == "#fractional\tfalse", f"Full MetaX 2 index required: {db}, got {value}")


def measure(argv):
    """Separate interpreter: RUSAGE_CHILDREN refers to this one synchronous command."""
    metrics, command = Path(argv[0]), argv[1:]
    start = time.monotonic()
    result = subprocess.run(command)
    r = resource.getrusage(resource.RUSAGE_CHILDREN)
    write_json(metrics, {"command": command, "exit_code": result.returncode,
                         "wall_seconds": time.monotonic() - start,
                         "peak_rss_kib_linux": r.ru_maxrss,
                         "user_seconds": r.ru_utime, "system_seconds": r.ru_stime})
    return result.returncode


def run_logged(command, log, env=None, cwd=None, measured=False):
    command = [str(x) for x in command]
    print("$ " + shlex.join(command), flush=True)
    write_json(Path(str(log) + ".command.json"), {"command": command, "cwd": str(cwd or Path.cwd())})
    actual = ([sys.executable, str(Path(__file__).resolve()), "_measure",
               str(log) + ".resources.json"] + command) if measured else command
    with open(log, "w") as f:
        result = subprocess.run(actual, stdout=f, stderr=subprocess.STDOUT, env=env, cwd=cwd)
    print(f"exit={result.returncode}; log={log}", flush=True)
    require(result.returncode == 0, f"Command failed ({result.returncode}): {log}")


def check_tools():
    names = ("python", "nextflow", "java", "fastp", "hisat2", "hisat2-inspect",
             "samtools", "multiqc", "bash", "awk", "grep")
    found = {}
    for name in names:
        found[name] = shutil.which(name)
        require(found[name] is not None, f"Missing tool on main environment PATH: {name}")
    found["python_executable"] = sys.executable
    return found


def executable_metax(args):
    path = safe_path(args.metax or args.repo.parent / "metax_env_2_0_3/bin/metax")
    file_ok(path)
    text = subprocess.check_output([str(path), "--version"], text=True).strip()
    require(text == "metax 2.0.3", f"Unexpected MetaX: {text}")
    return path


def baseline_inputs(args):
    repo, toy, samovar = args.repo, args.toy, args.samovar
    require("version     = '1.10.0'" in file_ok(repo / "nextflow.config").read_text(),
            "Expected pipeline v1.10.0; inspect manifest version before using kit")
    module = file_ok(repo / "modules/local/annotate_metax.nf").read_text()
    require("metax=2.0.3" in module and "--host 9606" in module,
            "Expected pinned MetaX 2 adapter")
    for name in ("metax_prepare.py", "metax_to_profile.py"):
        file_ok(repo / "bin" / name)
    refs = {key: safe_path(toy / "refs/9606.fna" if key == "9606" else
                           samovar / f"data/test_genomes/meta/{key}.fna") for key in SOURCE_SPECIES}
    for path in refs.values():
        file_ok(path)
    generation, build = read_json(file_ok(toy / "generation.json")), read_json(file_ok(toy / "build.json"))
    require(generation.get("pairing_verified") is True, "Original generation did not verify pairing")
    reads = [toy / f"reads/{name}" for name in ("toy_PE_R1.fastq", "toy_PE_R2.fastq", "toy_SE50.fastq")]
    verified = {}
    for path in reads + [refs["9606"]]:
        path = safe_path(file_ok(path))
        key = str(path.relative_to(toy))
        expected = generation.get("files", {}).get(key)
        require(expected is not None, f"Missing original input hash: {key}")
        digest = sha(path)
        require(digest == expected, f"Modified Samovar input: {path}")
        verified[str(path)] = digest
    for kind, expected in (("pe", [reads[0], reads[1]]), ("se", [reads[2]])):
        sheet = file_ok(toy / f"samplesheet_{kind}.csv")
        with open(sheet, newline="") as f:
            reader = csv.DictReader(f)
            require(reader.fieldnames == ["sample", "fastq_1", "fastq_2", "single_end"],
                    f"Unexpected samplesheet columns: {sheet}")
            rows = list(reader)
        require(len(rows) == 1, f"Expected one original sample: {sheet}")
        row = rows[0]
        require(bool(re.fullmatch(r"[A-Za-z0-9_][A-Za-z0-9_.-]*", row["sample"])), "Unsafe sample name")
        require(row["single_end"].lower() == ("true" if kind == "se" else "false"),
                f"Wrong pairing mode: {sheet}")
        require(kind != "se" or row["fastq_2"] == "", "SE must have empty fastq_2")
        for column, want in zip(("fastq_1", "fastq_2"), expected):
            raw = Path(row[column])
            actual = safe_path(raw if raw.is_absolute() else sheet.parent / raw)
            require(actual == want.resolve(), f"Samplesheet does not reference expected toy read: {actual}")
    raw_index = Path(build["hisat2_index"])
    index = safe_path(raw_index if raw_index.is_absolute() else toy / raw_index)
    parts = [p for ext in ("ht2", "ht2l") for i in range(1, 9)
             if (p := Path(f"{index}.{i}.{ext}")).is_file()]
    require(len(parts) == 8 and len({p.suffix for p in parts}) == 1, "Expected one complete 8-part HISAT2 index")
    for path in parts:
        file_ok(path)
        hashes = build.get("files", {})
        expected = hashes.get(str(path))
        if expected is None:
            expected = hashes.get(str(path.relative_to(toy))) if path.is_relative_to(toy) else None
        require(expected is not None, f"No original HISAT2 hash: {path}")
        verified[str(path)] = sha(path)
        require(verified[str(path)] == expected, f"Modified HISAT2 part: {path}")
    tax = repo / "tests/insilico/toy/taxonomy"
    tax_info = taxonomy_check(tax)
    panel = repo / "tests/insilico/toy/target_panel.txt"
    require(set(file_ok(panel).read_text().split()) == {"10847", "4932", "562"}, "Unexpected target panel")
    config = file_ok(repo / "tests/insilico/toy/toy.config")
    return refs, index, verified, tax_info, config, panel


def configure(args, index, config, panel):
    out = args.out
    text = (
        f"includeConfig '{config}'\n"
        "process {\n"
        "  executor = 'local'\n"
        "  cpus = 4\n  memory = '8 GB'\n"
        "  resourceLimits = [cpus: 4, memory: 8.GB, time: 24.h]\n"
        "  withLabel: 'process_low' { cpus = 4; memory = '8 GB' }\n"
        "  withLabel: 'process_medium' { cpus = 4; memory = '8 GB' }\n"
        "  withLabel: 'process_kraken' { cpus = 4 }\n"
        "}\nexecutor.queueSize = 1\n"
        "conda.enabled = false\ndocker.enabled = false\nsingularity.enabled = false\n"
        "trace.overwrite = true\n"
    )
    (out / "execution.config").write_text(text)
    csv_write(out / "dbsheet_metax.csv", ["db_id", "annotator", "path", "role"],
              [{"db_id": f"samovar_metax_{role}", "annotator": "metax",
                "path": str(out / f"db_{role}"), "role": role} for role in ("full", "patho")])
    return {"hisat2_index": str(index), "target_panel": str(panel), "toy_config": str(config)}


def prepare(args):
    require(not args.out.exists(), f"Refusing to overwrite: {args.out}")
    require(args.out != args.toy and not args.out.is_relative_to(args.toy),
            "Output must not be inside original toy directory")
    tools = check_tools()
    metax = executable_metax(args)
    refs, index, verified, tax_info, config, panel = baseline_inputs(args)
    args.out.mkdir(parents=True)
    logs = args.out / "logs"
    logs.mkdir()
    # Reconstruct only to prove this SAME HISAT2 index represents the human FASTA.
    extracted = args.out / "human_from_hisat2.fasta"
    with open(extracted, "w") as f, open(logs / "hisat2_inspect.stderr.log", "w") as err:
        subprocess.run(["hisat2-inspect", str(index)], stdout=f, stderr=err, check=True)
    require(seq_map(extracted) == seq_map(refs["9606"]), "HISAT2 sequences differ from original human FASTA")
    extracted.unlink()
    contigs = reheader(refs, args.out)
    metadata = configure(args, index, config, panel)
    tracked = list(refs.values()) + [args.toy / "generation.json", args.toy / "build.json",
              args.toy / "samplesheet_pe.csv", args.toy / "samplesheet_se.csv", config, panel,
              args.repo / "tests/insilico/toy/taxonomy/nodes.dmp",
              args.repo / "tests/insilico/toy/taxonomy/names.dmp",
              args.repo / "nextflow.config", args.repo / "main.nf"]
    tracked += sorted((args.repo / "modules").rglob("*.nf"))
    tracked += sorted((args.repo / "subworkflows").rglob("*.nf"))
    tracked += sorted(p for p in (args.repo / "bin").iterdir() if p.is_file())
    fingerprints = {str(p): sha(p) for p in tracked}
    fingerprints.update(verified)
    metadata.update({"kit_version": VERSION, "repo": str(args.repo), "toy": str(args.toy),
                     "out": str(args.out), "metax": str(metax), "metax_sha256": sha(metax),
                     "tools": tools, "input_sha256": fingerprints, "taxonomy": tax_info,
                     "human_index_sequences_verified": True,
                     "microbial_source_identity": "Current Samovar FASTA captured, not a historical generation hash",
                     "reference_bases": sum(m["bases"] for m in contigs),
                     "prepared_sha256": {str(p): sha(p) for p in (
                         args.out / "execution.config", args.out / "dbsheet_metax.csv",
                         args.out / "reference_manifest.json",
                         *(args.out / f"reference_{k}.fasta" for k in ("full", "patho", "pilot")))}})
    for command, name in ((["python", "--version"], "python"), (["nextflow", "-version"], "nextflow"),
                          (["java", "-version"], "java"), ([metax, "--version"], "metax"),
                          (["fastp", "--version"], "fastp"), (["hisat2", "--version"], "hisat2"),
                          (["samtools", "--version"], "samtools"), (["multiqc", "--version"], "multiqc")):
        run_logged(command, logs / f"version_{name}.log")
    # Record the resolved, merged config without starting any pipeline task.
    run_logged(["nextflow", "-c", args.out / "execution.config", "config", args.repo,
                "-profile", "standard", "-flat"], logs / "nextflow_resolved.config", cwd=args.out)
    metadata["prepare_completed"] = True
    write_json(args.out / "prepare.json", metadata)
    print(f"PREPARED: {args.out}; {metadata['reference_bases']} reference bases")


def state(args):
    data = read_json(file_ok(args.out / "prepare.json"))
    require(data.get("prepare_completed") is True, "Preparation did not finish")
    require(data["repo"] == str(args.repo) and data["toy"] == str(args.toy) and data["out"] == str(args.out),
            "Different repo/toy/output from prepare.json")
    for path, digest in {**data["input_sha256"], **data["prepared_sha256"],
                         data["metax"]: data["metax_sha256"]}.items():
        require(sha(file_ok(path)) == digest, f"Changed since prepare: {path}")
    return data


def build_indexes(args):
    data = state(args)
    require(not (args.out / "build_metax.json").exists(), "Indexes already built; use run pe/se next")
    logs, metax = args.out / "logs", Path(data["metax"])
    # No auto-overwrite or implicit recovery after an interrupted build.
    for kind in ("pilot", "full", "patho"):
        require(not (args.out / f"db_{kind}").exists(), f"Existing index; do not overwrite: db_{kind}")
    def index_one(kind):
        db = args.out / f"db_{kind}"
        run_logged([metax, "index", args.out / f"reference_{kind}.fasta", "-o", db, "-t", THREADS],
                   logs / f"index_{kind}.log", measured=True)
        full_index(db)
        taxonomy = db / "taxonomy"
        taxonomy.mkdir()
        for name in ("nodes.dmp", "names.dmp"):
            shutil.copyfile(args.repo / "tests/insilico/toy/taxonomy" / name, taxonomy / name)
        for name in ("merged.dmp", "delnodes.dmp"):
            (taxonomy / name).write_text("")
        run_logged([metax, "info", "-d", db], logs / f"info_{kind}.log")
        return read_json(logs / f"index_{kind}.log.resources.json")
    pilot = index_one("pilot")
    manifest = read_json(args.out / "reference_manifest.json")
    pilot_bases = sum(m["bases"] for m in manifest if m["source_taxid"] == "562")
    scale = data["reference_bases"] / pilot_bases
    estimate = {"pilot": pilot, "pilot_bases": pilot_bases, "full_bases": data["reference_bases"],
                "scale": scale, "full_rss_projection_kib": pilot["peak_rss_kib_linux"] * scale * 1.5,
                "full_seconds_projection": pilot["wall_seconds"] * scale * 1.5,
                "note": "Screening projection with 1.5x margin, not a guarantee; fixed costs and mapping caches differ."}
    write_json(args.out / "pilot_estimate.json", estimate)
    require(estimate["full_rss_projection_kib"] < 5 * 1024**2,
            "Pilot suggests >5 GiB indexing RAM; stop and review before full build")
    results = {"pilot": pilot}
    for kind in ("full", "patho"):
        results[kind] = index_one(kind)
    hashes = {}
    for kind in ("full", "patho"):
        for p in sorted((args.out / f"db_{kind}").rglob("*")):
            if p.is_file():
                hashes[str(p)] = sha(p)
    write_json(args.out / "build_metax.json", {"resources": results, "files_sha256": hashes,
                                             "index_mode": "full", "fractional": False})
    print("BUILT: full and patho; original Kraken2 databases untouched")


def validate_db(args):
    build = read_json(file_ok(args.out / "build_metax.json"))
    for path, digest in build["files_sha256"].items():
        require(sha(Path(path)) == digest, f"Changed MetaX index/taxonomy: {path}")
    for kind in ("full", "patho"):
        full_index(args.out / f"db_{kind}")


def pipeline_command(args, data):
    kind = args.kind
    return ["nextflow", "-log", str(args.out / "logs" / f"{kind}.nextflow.log"),
            "run", str(args.repo / "main.nf"), "-profile", "standard",
            "-c", str(args.out / "execution.config"), "-ansi-log", "false",
            "-work-dir", str(args.out / f"work_{kind}"),
            "-with-trace", str(args.out / "logs" / f"{kind}.trace.tsv"),
            "--mode", "fastq", "--input", str(args.toy / f"samplesheet_{kind}.csv"),
            "--dbsheet", str(args.out / "dbsheet_metax.csv"),
            "--host_index", data["hisat2_index"], "--target_panel", data["target_panel"],
            "--bracken_read_len", "126" if kind == "pe" else "50",
            "--kraken_mem", "2 GB", "--outdir", str(args.out / f"results_{kind}"),
            *(["-resume"] if args.resume else [])]


def pipeline_run(args):
    data = state(args)
    check_tools()
    validate_db(args)
    kind, logs = args.kind, args.out / "logs"
    marker = args.out / f"run_{kind}.json"
    require(not marker.exists(), f"Already completed {kind}; preserve these results")
    exists = (args.out / f"work_{kind}").exists() or (args.out / f"results_{kind}").exists()
    require(not exists or args.resume, "Partial run exists: collect diagnostics, then explicit --resume")
    env = os.environ.copy()
    toolbin = args.out / "toolbin"
    toolbin.mkdir(exist_ok=True)
    link = toolbin / "metax"
    if link.is_symlink():
        require(link.resolve() == Path(data["metax"]), "Unexpected MetaX symlink")
    elif not link.exists():
        link.symlink_to(data["metax"])
    else:
        raise ValueError("Unexpected toolbin/metax file")
    env["PATH"] = str(toolbin) + os.pathsep + env["PATH"]
    # Do not activate MetaX's Python environment or replace main Python.
    command = pipeline_command(args, data)
    write_json(args.out / f"attempt_{kind}.json", {"state": "running", "command": command,
                                                  "main_python": sys.executable,
                                                  "PATH_metax": shutil.which("metax", path=env["PATH"])})
    launch = args.out / f"launch_{kind}"
    launch.mkdir(exist_ok=True)
    try:
        run_logged(command, logs / f"pipeline_{kind}.log", env=env, cwd=launch, measured=True)
        profiles = list((args.out / f"results_{kind}").rglob("*.profile.json"))
        require(len(profiles) == 4, f"Expected host/dehost x full/patho profiles; got {len(profiles)}")
        for path in profiles:
            require(read_json(path).get("annotator") == "metax", f"Wrong annotator: {path}")
    except Exception:
        write_json(args.out / f"attempt_{kind}.json", {"state": "failed", "command": command})
        raise
    write_json(marker, {"exit_code": 0, "command": command, "profile_json_count": len(profiles),
                        "meaning": "Workflow execution succeeded. QC and accuracy not asserted."})
    write_json(args.out / f"attempt_{kind}.json", {"state": "completed", "command": command})
    print(f"COMPLETED {kind}: execution only; inspect QC and profiles separately")


def collect(args):
    require(args.out.is_dir(), f"No output directory: {args.out}")
    selected = []
    for p in sorted(args.out.iterdir()):
        if p.is_file() and p.suffix in {".json", ".csv", ".config"}:
            selected.append(p)
    for folder in ("logs", "results_pe", "results_se"):
        if (args.out / folder).exists():
            for p in sorted((args.out / folder).rglob("*")):
                if p.is_file() and (
                        p.suffix in {".json", ".tsv", ".csv", ".txt", ".log", ".html", ".config"}
                        or p.name == "metax_environment.explicit.txt"):
                    selected.append(p)
    # Only small task command/error text; never reads, alignments, caches or indices.
    for kind in ("pe", "se"):
        work = args.out / f"work_{kind}"
        if work.exists():
            for pattern in (".command.sh", ".command.err", ".command.out", ".exitcode"):
                selected.extend(sorted(work.glob(f"*/*/{pattern}")))
    selected = list(dict.fromkeys(p for p in selected if p.is_file()))
    receipt = {"files": [{"path": str(p.relative_to(args.out)), "bytes": p.stat().st_size,
                           "sha256": sha(p)} for p in selected],
               "execution": {k: read_json(args.out / f"run_{k}.json") if
                             (args.out / f"run_{k}.json").exists() else "not completed"
                             for k in ("pe", "se")},
               "excluded": ["FASTQ", "FASTA", "BAM", "MetaX indices/cache", "full work directories"]}
    stamp = time.strftime("%Y%m%d_%H%M%S", time.gmtime()) + f"_{time.time_ns() % 10**9:09d}"
    archive = args.out / f"metax_samovar_results_{stamp}.tar.gz"
    with tarfile.open(archive, "x:gz") as tar:
        for p in selected:
            tar.add(p, arcname="metax_samovar_results/" + str(p.relative_to(args.out)), recursive=False)
        for name in ("install.log", "metax_environment.explicit.txt"):
            p = Path(__file__).parent / name
            if p.is_file():
                tar.add(p, arcname="metax_samovar_results/setup/" + name)
        content = (json.dumps(receipt, ensure_ascii=False, indent=2) + "\n").encode()
        info = tarfile.TarInfo("metax_samovar_results/collection.json")
        info.size = len(content)
        tar.addfile(info, io.BytesIO(content))
    require(archive.stat().st_size > 0, "Empty archive")
    print(f"SEND THIS FILE: {archive}\nbytes={archive.stat().st_size}; sha256={sha(archive)}")


def main():
    if len(sys.argv) > 1 and sys.argv[1] == "_measure":
        return measure(sys.argv[2:])
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, default=Path.cwd())
    parser.add_argument("--toy", type=Path)
    parser.add_argument("--samovar", type=Path)
    parser.add_argument("--out", type=Path)
    parser.add_argument("--metax", type=Path)
    sub = parser.add_subparsers(dest="action", required=True)
    sub.add_parser("prepare")
    sub.add_parser("index")
    run = sub.add_parser("run")
    run.add_argument("kind", choices=("pe", "se"))
    run.add_argument("--resume", action="store_true")
    sub.add_parser("collect")
    args = parser.parse_args()
    args.repo = safe_path(args.repo)
    args.toy = safe_path(args.toy or args.repo / "insilico_runs/toy_995")
    args.samovar = safe_path(args.samovar or args.repo.parent / "samovar")
    args.out = safe_path(args.out or args.repo / "insilico_runs/metax_samovar_v1")
    try:
        {"prepare": prepare, "index": build_indexes, "run": pipeline_run, "collect": collect}[args.action](args)
    except (OSError, ValueError, KeyError, subprocess.SubprocessError) as error:
        print(f"STOP: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
