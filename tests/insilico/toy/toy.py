#!/usr/bin/env python3
"""Real Samovar -> PE126/SE50 -> nipt_mngs_qc; no synthetic substitute."""
import argparse
from collections import Counter
import csv
import gzip
import hashlib
import importlib.metadata
import json
import math
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import sys
import time

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
SOURCE_TO_SPECIES = {"2886930": "10847", "4932": "4932", "562": "562"}
TAXIDS = {"9606", *SOURCE_TO_SPECIES}
PIN = "d6bb8246f9e31c7890f4393aa92bc7e0aa53d593"
MAX_RELATIVE_ERROR = .20 + 1e-12  # floating-point equality at the exact 20% boundary


def write_json(path, value):
    Path(path).write_text(json.dumps(value, ensure_ascii=False, indent=2) + "\n")


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(2**20), b""):
            h.update(block)
    return h.hexdigest()


def run(cmd, log, env=None, cwd=None):
    """Fail closed; completed files only, not append on object-backed mounts."""
    cmd = [str(x) for x in cmd]
    print("+", shlex.join(cmd), flush=True)
    start = time.monotonic()
    with open(log, "w") as out:
        result = subprocess.run(cmd, stdout=out, stderr=subprocess.STDOUT,
                                env=env, cwd=cwd)
    print(f"  exit={result.returncode}; {time.monotonic()-start:.1f}s; {log}",
          flush=True)
    if result.returncode:
        raise RuntimeError(f"Command failed ({result.returncode}); see {log}")


def text_capture(cmd, env=None):
    return subprocess.check_output([str(x) for x in cmd], text=True, env=env,
                                   stderr=subprocess.STDOUT).strip()


def fq(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as f:
        while True:
            h = f.readline().rstrip("\r\n")
            if not h:
                break
            s, plus, q = [f.readline().rstrip("\r\n") for _ in range(3)]
            if not h.startswith("@") or not plus.startswith("+") or len(s) != len(q):
                raise ValueError(f"Invalid FASTQ record in {path}: {h}")
            yield h, s, q


def read_id(header):
    return re.sub(r"/[12]$", "", header.lstrip("@").split()[0])


def taxid(header):
    # Microbial ISS IDs: taxid:562_<fragment>/1; host IDs: ...|taxid:9606/1.
    matches = re.findall(r"(?:^|[|@])taxid:(\d+)(?=[_/|\s]|$)", header)
    if len(matches) != 1 or matches[0] not in TAXIDS:
        raise ValueError(f"Unknown/ambiguous taxid in read ID: {header}")
    return matches[0]


def inspect_fastq(path, expected_length=None):
    counts, lengths, ids = Counter(), Counter(), set()
    id_hash = hashlib.sha256()
    for h, s, q in fq(path):
        ident = read_id(h)
        if ident in ids:
            raise ValueError(f"Duplicate ID in {path}: {ident}")
        ids.add(ident)
        counts[taxid(h)] += 1
        lengths[len(s)] += 1
        id_hash.update((ident + "\n").encode())
    if not ids or (expected_length and set(lengths) != {expected_length}):
        raise ValueError(f"Empty/wrong read lengths in {path}: {dict(lengths)}")
    return {"reads": sum(counts.values()), "taxids": dict(counts),
            "lengths": dict(lengths), "id_order_sha256": id_hash.hexdigest()}, ids


def check_pairs(r1, r2):
    import itertools
    n = 0
    for a, b in itertools.zip_longest(fq(r1), fq(r2)):
        if a is None or b is None or read_id(a[0]) != read_id(b[0]):
            raise ValueError(f"Mate mismatch at record {n+1}")
        if taxid(a[0]) != taxid(b[0]):
            raise ValueError("Mate taxid mismatch")
        n += 1
    return n


def fasta_records(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as f:
        name, chunks = None, []
        for line in f:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks).upper()
                name, chunks = line[1:].split()[0], []
            else:
                chunks.append(line.strip())
        if name is not None:
            yield name, "".join(chunks).upper()


def fasta_hashes(path):
    result = {}
    for name, seq in fasta_records(path):
        if name in result:
            raise ValueError(f"Duplicate FASTA contig {name}")
        result[name] = (len(seq), hashlib.sha256(seq.encode()).hexdigest())
    if not result:
        raise ValueError(f"Empty FASTA: {path}")
    return result


def write_csv(path, header, rows):
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(header)
        w.writerows(rows)


def parameters(args):
    model = Path(text_capture([args.samovar_python, "-c",
                               "import iss,pathlib; print(pathlib.Path(iss.__file__).parent/'profiles/HiSeq')"]))
    return {"schema": 1, "pairs": args.pairs, "seed": args.seed,
            "host_fraction": 0.995, "model": "hiseq", "generation_cores": 1,
            "host_fasta": str(args.host_fasta), "host_sha256": sha256(args.host_fasta),
            "microbes": {t: sha256(args.samovar_dir / f"data/test_genomes/meta/{t}.fna")
                        for t in sorted(SOURCE_TO_SPECIES)},
            "samovar_dir": str(args.samovar_dir),
            "samovar_commit": text_capture(["git", "-C", args.samovar_dir, "rev-parse", "HEAD"]),
            "samovar_dirty": bool(text_capture(["git", "-C", args.samovar_dir,
                                               "status", "--porcelain"])),
            "samovar_python": str(args.samovar_python),
            "iss_model_sha256": sha256(model),
            "taxonomy_sha256": {n: sha256(HERE/"taxonomy"/n) for n in ("nodes.dmp", "names.dmp")},
            "hisat2_index": str(args.hisat2_index) if args.hisat2_index else None,
            "source_to_species": SOURCE_TO_SPECIES}


def preflight(args):
    for value in (args.kraken_memory, args.task_memory):
        if not re.fullmatch(r"[1-9][0-9]*(?: (?:MB|GB|TB))", value):
            raise ValueError("Memory must be a positive integer with units, e.g. '2 GB'")
    for exe in ("seqkit", "fastp", "hisat2", "hisat2-build", "hisat2-inspect",
                "samtools", "kraken2", "kraken2-build", "bracken", "bracken-build",
                "nextflow", "multiqc"):
        if not shutil.which(exe):
            raise RuntimeError(f"Missing executable on pipeline PATH: {exe}")
    if args.pairs < 1 or args.threads < 1:
        raise ValueError("pairs and threads must be positive")
    for p in (args.host_fasta, args.samovar_python, args.samovar_dir/"bin/samovar",
              args.samovar_dir/"src/samovar/iss_config.py",
              args.samovar_dir/"workflow/iss_test/Snakefile"):
        if not p.is_file():
            raise FileNotFoundError(p)
    found = {p.stem for p in (args.samovar_dir/"data/test_genomes/meta").glob("*.fna")}
    if found != set(SOURCE_TO_SPECIES):
        raise ValueError(f"Expected exactly three toy genomes, found {sorted(found)}")
    pybin = args.samovar_python.parent
    if not all((pybin/e).is_file() for e in ("iss", "snakemake")):
        raise RuntimeError("ISS and Snakemake must be beside --samovar-python")
    version = text_capture([args.samovar_python, "-c",
                            "import sys; assert (3,10)<=sys.version_info[:2]<(3,14); "
                            "import samovar,numpy,pandas; print(sys.version)"])
    # Own run-local install config; never edit ../samovar/build or ~/.config.
    cfg = json.loads((args.samovar_dir/"build/config.json.example").read_text())
    cfg["root"] = str(args.samovar_dir)
    cfg["compilers"]["python"] = str(args.samovar_python)
    for name in ("python", "iss", "snakemake"):
        exe = args.samovar_python if name == "python" else pybin/name
        cfg["tools"][name] = ["", "bash", str(exe),
                              "reads_generator" if name == "iss" else "runtime"]
    cfg["genomes"]["processed"]["test"] = str(args.samovar_dir/"data/test_genomes")
    cfg_path = args.outdir/"samovar_config.json"
    write_json(cfg_path, cfg)
    env = os.environ.copy()
    env["SAMOVAR_CONFIG"] = str(cfg_path)
    env["PATH"] = str(pybin) + os.pathsep + env["PATH"]
    env["PYTHONPATH"] = str(args.samovar_dir/"src")
    args.samovar_env = env
    params = parameters(args)
    manifest = args.outdir/"manifest.json"
    if manifest.exists():
        if json.loads(manifest.read_text())["parameters"] != params:
            raise ValueError("Existing outdir has incompatible manifest; use a new --outdir")
    else:
        write_json(manifest, {"parameters": params, "python": version,
                             "created_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())})
    freeze = text_capture([args.samovar_python, "-c",
                           "import importlib.metadata as m; "
                           "print('\\n'.join(sorted(f'{d.metadata[\"Name\"]}=={d.version}' "
                           "for d in m.distributions())))"])
    (args.outdir/"samovar_packages.txt").write_text(freeze + "\n")
    versions = {}
    for name, flags in {"fastp": ["--version"], "hisat2": ["--version"],
                        "samtools": ["--version"], "kraken2": ["--version"],
                        "bracken": ["-v"], "multiqc": ["--version"],
                        "nextflow": ["-version"], "seqkit": ["version"]}.items():
        result = subprocess.run([name]+flags, text=True, errors="replace",
                                stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        versions[name] = {"executable": shutil.which(name), "version_output": result.stdout.strip(),
                          "exit_code": result.returncode}
    write_json(args.outdir/"tool_versions.json", versions)
    for action in ("generate", "prepare", "exec"):
        run([args.samovar_dir/"bin/samovar", action, "--help"],
            args.outdir/f"help_{action}.log", env)
    helptext = (args.outdir/"help_prepare.log").read_text()
    if "--startpoint" not in helptext or "--endpoint" not in helptext:
        raise RuntimeError("Samovar prepare lacks checkpoint slicing; use documented pinned commit")
    if shutil.disk_usage(args.outdir).free < 4*2**30:
        raise RuntimeError("Less than 4 GiB free disk; choose another working directory")
    print("Preflight OK", flush=True)


def generate(args):
    out = args.outdir
    marker = out/"generation.json"
    if marker.exists():
        for name, digest in json.loads(marker.read_text())["files"].items():
            if sha256(out/name) != digest:
                raise ValueError(f"Modified generated FASTQ: {name}")
        print("Generation verified; reusing existing reads", flush=True)
        return
    refs = out/"refs"
    refs.mkdir(exist_ok=True)
    host = refs/"9606.fna"
    # Copy plain sequence, retaining contig IDs and normalizing compression/wrapping.
    with open(host, "w") as f:
        for name, seq in fasta_records(args.host_fasta):
            f.write(f">{name}\n")
            for i in range(0, len(seq), 80):
                f.write(seq[i:i+80] + "\n")
    if fasta_hashes(host) != fasta_hashes(args.host_fasta):
        raise RuntimeError("Staged host FASTA differs")
    samovar_out = out/"samovar"
    cli = args.samovar_dir/"bin/samovar"
    run([cli, "generate", "--genome_dir", args.samovar_dir/"data/test_genomes/meta",
         "--host_genome", host, "--host_fraction", "0.995", "--n_samples", "1",
         "--total_reads", args.pairs, "--seed", args.seed, "--model", "hiseq",
         "--cores", "1", "--output_dir", samovar_out],
        out/"samovar_generate.log", args.samovar_env)
    run([cli, "prepare", "--startpoint", "setup_reads", "--endpoint", "setup_reads",
         "--cores", "1", "--output_dir", samovar_out],
        out/"samovar_prepare.log", args.samovar_env)
    run([cli, "exec", "--output_dir", samovar_out],
        out/"samovar_exec.log", args.samovar_env)
    raw1, raw2 = [samovar_out/f"initial/1_full_R{mate}.fastq" for mate in (1, 2)]
    raw, ids = inspect_fastq(raw1, 126)
    raw_other, other_ids = inspect_fastq(raw2, 126)
    if raw != raw_other or ids != other_ids or check_pairs(raw1, raw2) != args.pairs:
        raise ValueError("Raw mates/counts differ")
    expected_host = int(0.995*args.pairs)
    if raw["taxids"].get("9606", 0) != expected_host or raw["reads"] != args.pairs:
        raise ValueError(f"Samovar composition/count mismatch: {raw}")
    microbial = [raw["taxids"].get(t, 0) for t in SOURCE_TO_SPECIES]
    if max(microbial)-min(microbial) > 1:
        raise ValueError(f"Microbial inputs not equal: {microbial}")
    reads = out/"reads"
    reads.mkdir(exist_ok=True)
    pe1, pe2, se = [reads/n for n in ("toy_PE_R1.fastq", "toy_PE_R2.fastq",
                                     "toy_SE50.fastq")]
    for mate, dest in ((raw1, pe1), (raw2, pe2)):
        run(["seqkit", "shuffle", "--rand-seed", args.seed, "-o", dest, mate],
            out/f"shuffle_{dest.stem}.log")
    shuffled, shuffled_ids = inspect_fastq(pe1, 126)
    if ids != shuffled_ids or check_pairs(pe1, pe2) != args.pairs:
        raise ValueError("Shuffle lost/changed pairs")
    if raw["id_order_sha256"] == shuffled["id_order_sha256"]:
        raise ValueError("Shuffle did not change read order")
    with open(se, "w") as f:
        for h, s, q in fq(pe1):
            f.write(f"{h}\n{s[:50]}\n+\n{q[:50]}\n")
    se_stats, se_ids = inspect_fastq(se, 50)
    if se_ids != ids:
        raise ValueError("SE derivation changed IDs")
    write_csv(out/"samplesheet_pe.csv", ["sample", "fastq_1", "fastq_2", "single_end"],
              [["SAMOVAR_PE", pe1, pe2, "false"]])
    write_csv(out/"samplesheet_se.csv", ["sample", "fastq_1", "fastq_2", "single_end"],
              [["SAMOVAR_SE", se, "", "true"]])
    files = [pe1, pe2, se, raw1, raw2, host]
    write_json(marker, {"raw_R1": raw, "shuffled_R1": shuffled, "SE": se_stats,
                        "pairing_verified": True,
                        "files": {str(p.relative_to(out)): sha256(p) for p in files}})
    print(f"Generated: {raw['reads']} pairs; ground truth {raw['taxids']}", flush=True)


def build(args):
    out = args.outdir
    if not (out/"generation.json").exists():
        raise ValueError("Run generate first")
    marker = out/"build.json"
    if marker.exists():
        for name, digest in json.loads(marker.read_text())["files"].items():
            if sha256(Path(name)) != digest:
                raise ValueError(f"Modified index/database: {name}")
        print("Database/index hashes verified; reusing", flush=True)
        return
    host = out/"refs/9606.fna"
    idx = args.hisat2_index or out/"host_index/human"
    idx.parent.mkdir(exist_ok=True, parents=True)
    if not args.hisat2_index:
        run(["hisat2-build", "-p", args.threads, host, idx], out/"hisat2_build.log")
    extracted = out/"refs/index_extracted.fna"
    run(["hisat2-inspect", idx], extracted)
    if fasta_hashes(extracted) != fasta_hashes(host):
        raise ValueError("HISAT2 index does not match the simulated human FASTA")
    extracted.unlink()
    files = list(idx.parent.glob(idx.name+".*.ht2")) + list(idx.parent.glob(idx.name+".*.ht2l"))
    if len(files) != 8:
        raise ValueError(f"Expected exactly eight HISAT2 index files: {idx}")
    for role in ("full", "patho"):
        db = out/f"db_{role}"
        if db.exists() and any(db.iterdir()):
            raise ValueError(f"Partial database without build.json: {db}; "
                             "move both db_full/db_patho aside before retrying build")
        db.mkdir(exist_ok=True)
        tax = db/"taxonomy"
        tax.mkdir(exist_ok=True)
        for name in ("nodes.dmp", "names.dmp"):
            shutil.copyfile(HERE/"taxonomy"/name, tax/name)
        tids = sorted(SOURCE_TO_SPECIES) + (["9606"] if role == "full" else [])
        for tid in tids:
            original = host if tid == "9606" else args.samovar_dir/f"data/test_genomes/meta/{tid}.fna"
            lib = db/f"reference_{tid}.fa"
            with open(lib, "w") as f:
                for i, (name, seq) in enumerate(fasta_records(original), 1):
                    f.write(f">toy_{tid}_{i}|kraken:taxid|{tid}\n")
                    for n in range(0, len(seq), 80):
                        f.write(seq[n:n+80] + "\n")
            run(["kraken2-build", "--add-to-library", lib, "--db", db],
                out/f"kraken_add_{role}_{tid}.log")
        run(["kraken2-build", "--build", "--db", db, "--threads", args.threads],
            out/f"kraken_build_{role}.log")
        for length in (126, 50):
            run(["bracken-build", "-d", db, "-t", args.threads, "-k", "35", "-l", length],
                out/f"bracken_build_{role}_{length}.log")
            files.append(db/f"database{length}mers.kmer_distrib")
        files.extend(db/n for n in ("hash.k2d", "opts.k2d", "taxo.k2d"))
    write_csv(out/"dbsheet.csv", ["db_id", "annotator", "path", "role"],
              [[f"samovar_{r}", "kraken2_bracken", out/f"db_{r}", r]
               for r in ("full", "patho")])
    write_json(marker, {"hisat2_index": str(idx),
                        "files": {str(p): sha256(p) for p in files}})


def pipeline_fingerprint():
    paths = [REPO/"main.nf", REPO/"nextflow.config"]
    paths += list((REPO/"modules").rglob("*.nf")) + list((REPO/"subworkflows").rglob("*.nf"))
    paths += [p for p in (REPO/"bin").iterdir() if p.is_file()]
    paths += [HERE/"toy.py", HERE/"toy.config", HERE/"target_panel.txt"]
    return {str(p.relative_to(REPO)): sha256(p) for p in sorted(paths)}


def pipeline(args):
    out = args.outdir
    if not (out/"build.json").exists():
        raise ValueError("Run build first")
    idx = json.loads((out/"build.json").read_text())["hisat2_index"]
    # No previous green acceptance may survive an interrupted/failed rerun.
    write_json(out/"acceptance.json", {"pass": False, "state": "pipeline_running"})
    (out/"acceptance.tsv").unlink(missing_ok=True)
    resources = out/"execution_resources.config"
    resources.write_text(
        f"includeConfig '{HERE}/toy.config'\n"
        "process {\n"
        f"  withLabel: 'process_low' {{ cpus = {args.threads}; memory = '{args.task_memory}' }}\n"
        f"  withLabel: 'process_medium' {{ cpus = {args.threads}; memory = '{args.task_memory}' }}\n"
        f"  withLabel: 'process_kraken' {{ cpus = {args.threads} }}\n"
        "}\ntrace { overwrite = true }\n")
    for kind, length in (("pe", 126), ("se", 50)):
        done = out/f"run_{kind}.json"
        # Remove success sentinel BEFORE starting: interruption must not reuse old success.
        done.unlink(missing_ok=True)
        cmd = ["nextflow", "-log", out/f"{kind}.nextflow.log", "run", REPO/"main.nf",
               "-c", resources, "-ansi-log", "false", "-resume",
               "-work-dir", out/"work", "-with-trace", out/f"{kind}.trace.tsv",
               "--mode", "fastq", "--input", out/f"samplesheet_{kind}.csv",
               "--dbsheet", out/"dbsheet.csv", "--host_index", idx,
               "--target_panel", HERE/"target_panel.txt", "--bracken_read_len", length,
               "--kraken_mem", args.kraken_memory,
               "--outdir", out/f"results_{kind}"]
        run(cmd, out/f"pipeline_{kind}.log", cwd=out)
        write_json(done, {"exit_code": 0, "command": [str(x) for x in cmd],
                          "pipeline_files": pipeline_fingerprint(),
                          "execution_resources_sha256": sha256(resources),
                          "generation_sha256": sha256(out/"generation.json"),
                          "build_sha256": sha256(out/"build.json")})


def profile_metrics(tsv):
    with open(tsv) as f:
        taxa = list(csv.DictReader(f, delimiter="\t"))
    counts = {t["taxid"]: float(t["reads"]) for t in taxa}
    if len(counts) != len(taxa) or any(not math.isfinite(v) or v < 0 for v in counts.values()):
        raise ValueError(f"Duplicate taxid or invalid reads in {tsv}")
    expected = set(SOURCE_TO_SPECIES.values())
    missing = sorted(t for t in expected if counts.get(t, 0) < 10)
    extras = {t: v for t, v in counts.items() if t not in expected and v >= 10}
    total = sum(counts.get(t, 0) for t in expected)
    shares = {t: counts.get(t, 0)/total if total else 0 for t in sorted(expected)}
    errors = {t: abs(v-1/3)/(1/3) for t, v in shares.items()}
    return {"reads": counts, "missing": missing, "unexpected": extras,
            "shares": shares, "relative_errors": errors,
            "pass": not missing and not extras and all(v <= MAX_RELATIVE_ERROR for v in errors.values()),
            "contract_rows": taxa}


def bam_truth(bam, kind, scratch):
    """Foreground samtools; audit identities, not inference from abundance."""
    run(["samtools", "view", "-F", "0x900", bam], scratch)
    keys, counts = Counter(), Counter()
    with open(scratch) as f:
        for line in f:
            cols = line.split("\t")
            flag = int(cols[1])
            if kind == "pe":
                if not flag & 1 or bool(flag & 64) == bool(flag & 128):
                    raise ValueError("Invalid PE flags in validation BAM")
                mate = 1 if flag & 64 else 2
            else:
                if flag & 1:
                    raise ValueError("SE unexpectedly marked paired")
                mate = 0
            ident = read_id(cols[0])
            keys[(ident, mate)] += 1
            counts[taxid(cols[0])] += 1
    scratch.unlink()
    return keys, counts


def acceptance(args):
    out = args.outdir
    rows, details = [], {}

    def criterion(kind, name, ok, observed, expected):
        rows.append({"format": kind, "criterion": name, "pass": bool(ok),
                     "observed": observed, "expected": expected})

    try:
        gen = json.loads((out/"generation.json").read_text())
        raw = gen["raw_R1"]
        for rel, digest in gen["files"].items():
            if sha256(out/rel) != digest:
                raise ValueError(f"Changed generation artifact: {rel}")
        m = [raw["taxids"].get(t, 0) for t in SOURCE_TO_SPECIES]
        criterion("both", "raw_ground_truth",
                  raw["reads"] == args.pairs and raw["taxids"].get("9606") == int(.995*args.pairs)
                  and max(m)-min(m) <= 1,
                  raw["taxids"], "99.5% HOST; equal three microbes")
        for kind in ("pe", "se"):
            sample = "SAMOVAR_"+kind.upper()
            root = out/f"results_{kind}"
            sdir = root/sample
            state = json.loads((out/f"run_{kind}.json").read_text())
            fresh = (state["exit_code"] == 0 and state["pipeline_files"] == pipeline_fingerprint()
                     and state["execution_resources_sha256"] == sha256(out/"execution_resources.config")
                     and state["generation_sha256"] == sha256(out/"generation.json")
                     and state["build_sha256"] == sha256(out/"build.json"))
            criterion(kind, "pipeline_completed_current_inputs", fresh, fresh, True)
            required = [root/"multiqc/multiqc_report.html",
                        sdir/f"qc/{sample}_qc_report.json",
                        sdir/f"qc/fastp/{sample}.fastp.json",
                        sdir/f"qc/dehost/{sample}.hisat2.summary.txt"]
            for fraction in ("dehost", "host"):
                for role in ("full", "patho"):
                    required += [sdir/f"profiles/{sample}.{fraction}.{role}.profile.{ext}"
                                 for ext in ("tsv", "json")]
            for role in ("full", "patho"):
                required += [sdir/f"profiles/{sample}.{role}.dehost_vs_host.compare.{ext}"
                             for ext in ("tsv", "json")]
            required += [sdir/f"qc/counts/{sample}.{fr}.counts.json"
                         for fr in ("dehost", "host", "half")]
            missing = [str(p.relative_to(out)) for p in required
                       if not p.is_file() or p.stat().st_size == 0]
            criterion(kind, "required_reports", not missing, missing, "all present and nonempty")
            if missing:
                continue
            trimmed = sdir/"validation/trimmed"
            t1 = trimmed/"trim_R1.fastq.gz"
            info, t_ids = inspect_fastq(t1)
            raw_path = out/("reads/toy_PE_R1.fastq" if kind == "pe" else "reads/toy_SE50.fastq")
            original, orig_ids = inspect_fastq(raw_path, 126 if kind == "pe" else 50)
            criterion(kind, "trimmed_ids_subset", t_ids <= orig_ids,
                      len(t_ids - orig_ids), 0)
            trimmed_keys = Counter((i, 1 if kind == "pe" else 0) for i in t_ids)
            trim_microbes = sum(v for t, v in info["taxids"].items() if t != "9606")
            if kind == "pe":
                t2 = trimmed/"trim_R2.fastq.gz"
                check_pairs(t1, t2)
                trimmed_keys.update((i, 2) for i in t_ids)
                trim_microbes *= 2
            partitions, truth, nreads, combined = {}, {}, {}, Counter()
            for fr in ("host", "dehost", "half"):
                keys, counts = bam_truth(sdir/f"validation/partitions/{fr}.bam",
                                         kind, out/f"audit_{kind}_{fr}.sam")
                partitions[fr] = len(keys)
                truth[fr] = dict(counts)
                nreads[fr] = sum(counts.values())
                combined.update(keys)
                cj = json.loads((sdir/f"qc/counts/{sample}.{fr}.counts.json").read_text())
                criterion(kind, f"{fr}_counts_match_bam", cj["primary_reads"] == nreads[fr],
                          cj["primary_reads"], nreads[fr])
            exact = combined == trimmed_keys and all(v == 1 for v in combined.values())
            criterion(kind, "partition_ids_exact", exact,
                      {"partition_records": sum(nreads.values()), "trimmed_records": sum(trimmed_keys.values())},
                      "exact disjoint ID+mate partition")
            denom = sum(nreads.values())
            hostfrac = nreads["host"]/denom if denom else 0
            criterion(kind, "host_fraction", .993 <= hostfrac <= .997, hostfrac, "[0.993, 0.997]")
            if kind == "se":
                criterion(kind, "SE_half_empty", nreads["half"] == 0, nreads["half"], 0)
            dehost_microbes = sum(v for t, v in truth["dehost"].items() if t != "9606")
            criterion(kind, "trimmed_GT_microbial_reads", trim_microbes >= 1000, trim_microbes, ">=1000")
            criterion(kind, "DEHOST_GT_microbial_reads", dehost_microbes >= 1000, dehost_microbes, ">=1000")
            profiles = {}
            for role in ("full", "patho"):
                pr = profile_metrics(sdir/f"profiles/{sample}.dehost.{role}.profile.tsv")
                profiles[role] = pr
                criterion(kind, f"DEHOST_{role}_species", not pr["missing"] and not pr["unexpected"],
                          {"missing": pr["missing"], "unexpected": pr["unexpected"]}, "exact expected species at reads>=10")
                criterion(kind, f"DEHOST_{role}_composition",
                          all(v <= MAX_RELATIVE_ERROR for v in pr["relative_errors"].values()),
                          pr["shares"], "each 0.266667..0.400000")
            qc = json.loads((sdir/f"qc/{sample}_qc_report.json").read_text())
            criterion(kind, "QC_PASS", qc["status"] == "PASS_POSITIVE",
                      {"status": qc["status"], "reasons": qc.get("fallback_reasons")}, "PASS_POSITIVE")
            fastp = json.loads((sdir/f"qc/fastp/{sample}.fastp.json").read_text())
            criterion(kind, "fastp_record_units",
                      fastp["summary"]["after_filtering"]["total_reads"] == sum(nreads.values()),
                      fastp["summary"]["after_filtering"]["total_reads"], sum(nreads.values()))
            details[kind] = {"HOST_record_fraction": hostfrac, "counts": nreads,
                             "ground_truth": truth, "profiles": profiles, "QC": qc,
                             "fastp": fastp["summary"],
                             "trimmed_GT_microbes": trim_microbes,
                             "DEHOST_GT_microbes": dehost_microbes}
    except (OSError, ValueError, KeyError, RuntimeError) as e:
        criterion("both", "complete_audit", False, str(e), "all inputs available and valid")
    ok = bool(rows) and all(r["pass"] for r in rows) and set(details) == {"pe", "se"}
    write_json(out/"acceptance.json", {"pass": ok, "checks": rows, "details": details,
                                     "min_reads": 10, "relative_deviation": .20,
                                     "source_to_species": SOURCE_TO_SPECIES})
    write_csv(out/"acceptance.tsv", ["format", "criterion", "pass", "observed", "expected"],
              [[r[k] if not isinstance(r[k], (dict, list)) else json.dumps(r[k], ensure_ascii=False)
                for k in ("format", "criterion", "pass", "observed", "expected")] for r in rows])
    print("ACCEPTANCE", "PASS" if ok else "FAIL", flush=True)
    for r in rows:
        if not r["pass"]:
            print(f"  FAIL {r['format']}/{r['criterion']}: {r['observed']}", flush=True)
    return ok


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("stage", choices=["preflight", "generate", "build", "run", "check", "all"])
    p.add_argument("--samovar-dir", type=Path, default=REPO.parent/"samovar")
    p.add_argument("--samovar-python", type=Path,
                   help="Default: <samovar-dir>/.venv-nipt-toy/bin/python")
    p.add_argument("--host-fasta", type=Path)
    p.add_argument("--hisat2-index", type=Path)
    p.add_argument("--outdir", type=Path, required=True)
    p.add_argument("--pairs", type=int, default=600000)
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--threads", type=int, default=2)
    p.add_argument("--kraken-memory", default="2 GB",
                   help="Nextflow Kraken2 memory request; size to the actual database")
    p.add_argument("--task-memory", default="2 GB",
                   help="Nextflow fastp/HISAT2 and other non-Kraken task memory")
    args = p.parse_args()
    for name in ("samovar_dir", "samovar_python", "outdir", "hisat2_index", "host_fasta"):
        value = getattr(args, name)
        if value is not None:
            # Do not resolve Python symlink: its venv location is significant.
            setattr(args, name, Path(os.path.abspath(value)))
    args.host_fasta = args.host_fasta or args.samovar_dir/"data/test_genomes/host/9606.fna"
    args.samovar_python = args.samovar_python or args.samovar_dir/".venv-nipt-toy/bin/python"
    args.outdir.mkdir(parents=True, exist_ok=True)
    preflight(args)
    if args.stage in ("generate", "all"):
        generate(args)
    if args.stage in ("build", "all"):
        build(args)
    if args.stage in ("run", "all"):
        pipeline(args)
    if args.stage in ("check", "all"):
        if not acceptance(args):
            sys.exit(1)


if __name__ == "__main__":
    try:
        main()
    except (RuntimeError, ValueError, FileNotFoundError, subprocess.CalledProcessError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(1)
