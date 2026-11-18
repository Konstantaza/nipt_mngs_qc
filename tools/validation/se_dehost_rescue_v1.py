#!/usr/bin/env python3
"""One pinned toy-only SE dehosting experiment. No changes to pipeline/baseline.

Requires Python >=3.10, Linux, HISAT2 2.2.3 (including hisat2-inspect), samtools.
Run from the nipt_mngs_qc repository: python se_dehost_rescue_v1.py
The output directory must not already exist. Read summary.txt after any exit.
GT labels are used for validation/evaluation, never for alignment or removal.
"""
import argparse
from collections import Counter
import csv
import hashlib
import json
from pathlib import Path
import re
import resource
import shutil
import subprocess
import sys
import time

VERSION = "1.0.0"
EXPECTED = {"9606": 516, "562": 1000, "4932": 1000, "2886930": 1000}
TAXID = re.compile(r"(?:^|[|@])taxid:(\d+)(?=[_/|\s]|$)")
HUMAN_SHA = "3b00f25a9f6dd31c945a5a2891634ca84a0325b548f5a009143f397be0e3deac"
SE_SHA = "81cd49a32e2fad63497db5442a250acadd84ff68534d150ebe47bf84ae73d3dc"
INDEX_SHA = [
    "2ec75a1d7b70abf01ead6977abecb2028915a79310a050ea8f40fe2226928f26",
    "4e1644593c087af94379dc222971cadbfdefbe1c2bff1c3dc5122147a6944ae5",
    "52dd2eed2050e3f60993f32a6d5c6e405e85a5f300de87ddf58836b31fc50aca",
    "8f8ee6283d2595ed19f35f33bae8af2c00d677d4e1b6fb3b4f0a74bafeaf3842",
    "9780c4b048550c0304d7ce635632c19e0837e8e2845c85807b9ab55719edc392",
    "9d812e61110081da39b099a2387b43e2f3263cf431f425103df2b88f227dc422",
    "ca888f40c3caca805b37a5434c75de5550616e0795e7602fb91156f22dd90851",
    "7c9fa136d4413fa6173637e883b6998d32e1d675f88cddff9dcbcf331820f4b8",
]
RC = str.maketrans("ACGTRYSWKMBDHVN", "TGCAYRSWMKVHDBN")
PILOT_N = 256
THREADS = 2


def require(ok, message):
    if not ok:
        raise ValueError(message)


def sha(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1048576), b""):
            h.update(block)
    return h.hexdigest()


def write_json(path, obj):
    Path(path).write_text(json.dumps(obj, ensure_ascii=False, indent=2) + "\n")


def tool(name, override=None):
    found = shutil.which(override or name)
    require(found is not None, f"Executable not found: {override or name}")
    # Do not resolve symlinks: conda/dispatcher wrappers may use argv[0].
    return str(Path(found).absolute())


def source_taxid(name):
    hits = TAXID.findall(name)
    require(name.count("taxid:") == 1 and len(hits) == 1,
            f"Ambiguous/missing GT taxid in toy ID: {name}")
    require(hits[0] in EXPECTED, f"Unexpected source taxid: {name}")
    return hits[0]


def fasta_map(path):
    result, name, chunks = {}, None, []

    def store():
        if name is not None:
            seq = "".join(chunks).upper()
            require(seq and re.fullmatch("[ACGTRYSWKMBDHVN]+", seq),
                    f"Empty/non-IUPAC FASTA: {name}")
            require(name not in result, f"Duplicate FASTA contig: {name}")
            result[name] = [len(seq), hashlib.sha256(seq.encode()).hexdigest()]

    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                store()
                require(bool(line[1:].strip()), "Empty FASTA header")
                name, chunks = line[1:].split()[0], []
            elif line.strip():
                require(name is not None, "Sequence before FASTA header")
                chunks.append(line.strip())
        store()
    require(bool(result), "No FASTA sequences")
    return result


def fastq_records(path):
    with open(path) as f:
        while True:
            h = f.readline()
            if not h:
                break
            s, p, q = f.readline().rstrip("\n"), f.readline(), f.readline().rstrip("\n")
            require(h.startswith("@") and p.startswith("+") and len(s) == len(q) > 0,
                    f"Invalid/truncated FASTQ record in {path}")
            yield h[1:].split()[0], s, q


def parse_input_sam(text):
    records = {}
    for line in text.splitlines():
        c = line.split("\t")
        require(len(c) >= 11, "Malformed input SAM record")
        name, flag, seq, qual = c[0], int(c[1]), c[9], c[10]
        require(not flag & (0x100 | 0x800), f"Nonprimary input record: {name}")
        require(flag & 4 and not flag & (1 | 16 | 64 | 128),
                f"Expected unpaired, forward, unmapped SE input: {name}, flag={flag}")
        require(name not in records, f"Duplicate input QNAME: {name}")
        require(len(seq) == len(qual) == 50 and bool(re.fullmatch("[ACGTRYSWKMBDHVN]+", seq)),
                f"Expected 50nt DNA sequence and quality: {name}")
        require(all(33 <= ord(x) <= 126 for x in qual), f"Invalid quality: {name}")
        records[name] = {"seq": seq, "qual": qual, "taxid": source_taxid(name)}
    require(Counter(r["taxid"] for r in records.values()) == Counter(EXPECTED),
            "Unexpected toy GT balance; require 516 human and 1000 of each microbe")
    return records


def verify_original(records, path):
    matched, case_differences = set(), 0
    for name, seq, qual in fastq_records(path):
        if name not in records:
            continue
        require(name not in matched, f"Duplicate original QNAME: {name}")
        r = records[name]
        require(seq.upper() == r["seq"] and qual == r["qual"],
                f"BAM sequence/quality differs from original SE read: {name}")
        case_differences += seq != r["seq"]
        matched.add(name)
    require(matched == set(records), "Some DEHOST IDs missing from original FASTQ")
    return {"matched": len(matched), "raw_vs_bam_case_only_differences": case_differences,
            "note": "BAM encodes uppercase bases; exact qualities and case-insensitive bases checked."}


def write_fastq(path, records, names=None):
    with open(path, "w") as f:
        for name in records if names is None else names:
            r = records[name]
            f.write(f'@{name}\n{r["seq"]}\n+\n{r["qual"]}\n')


def input_diagnostics(records):
    rows = {}
    for tid in EXPECTED:
        rr = [r for r in records.values() if r["taxid"] == tid]
        qq = [ord(x) - 33 for r in rr for x in r["qual"]]
        rows[tid] = {
            "reads": len(rr), "length_counts": dict(Counter(len(r["seq"]) for r in rr)),
            "quality_min": min(qq), "quality_mean": sum(qq) / len(qq),
            "quality_max": max(qq), "Q30_base_fraction": sum(x >= 30 for x in qq) / len(qq),
            "N_bases": sum(r["seq"].count("N") for r in rr),
            "reads_with_N": sum("N" in r["seq"] for r in rr),
        }
    return rows


def preflight(args, out):
    hisat = tool("hisat2", args.hisat2)
    samtools = tool("samtools", args.samtools)
    inspect = tool("hisat2-inspect", args.inspect)
    version = subprocess.check_output([hisat, "--version"], text=True)
    require(re.search(r"\bversion 2\.2\.3\b", version) is not None,
            "This experiment requires HISAT2 2.2.3; no automatic installation")
    inspector_version = subprocess.check_output([inspect, "--version"], text=True)
    require(re.search(r"\bversion 2\.2\.3\b", inspector_version) is not None,
            "hisat2-inspect must also be 2.2.3")
    help_text = subprocess.check_output([hisat, "--help"], text=True)
    line = next((x for x in help_text.splitlines() if x.strip().startswith("--sensitive ")), "")
    require(line.split() == ["--sensitive", "--bowtie2-dp", "1", "-k", "30",
                             "--score-min", "L,0,-0.5"], "Unexpected sensitive preset expansion")
    toy = args.toy
    build = json.loads((toy / "build.json").read_text())
    index = Path(build["hisat2_index"])
    if not index.is_absolute():
        # Paths in build.json are interpreted relative to the repository cwd.
        index = args.repo / index
    original = toy / "reads/toy_SE50.fastq"
    human = toy / "refs/9606.fna"
    input_hashes = {}
    for path, expected in [(original, SE_SHA), (human, HUMAN_SHA)] + [
            (Path(f"{index}.{i}.ht2"), expected) for i, expected in enumerate(INDEX_SHA, 1)]:
        require(path.is_file(), f"Missing input: {path}")
        observed = sha(path)
        require(observed == expected, f"Pinned toy hash mismatch: {path}")
        input_hashes[str(path)] = observed
    reconstructed = out / "human_index_check.fasta"
    with open(reconstructed, "w") as f:
        subprocess.run([inspect, str(index)], stdout=f, check=True)
    require(fasta_map(reconstructed) == fasta_map(human), "Index/reference sequence mismatch")
    require(args.bam.is_file(), f"Missing BAM: {args.bam}")
    subprocess.run([samtools, "quickcheck", "-v", str(args.bam)], check=True)
    total = int(subprocess.check_output([samtools, "view", "-c", str(args.bam)], text=True))
    require(total == sum(EXPECTED.values()), f"Expected 3516 total BAM records, got {total}")
    text = subprocess.check_output([samtools, "view", str(args.bam)], text=True)
    records = parse_input_sam(text)
    original_check = verify_original(records, original)
    header = subprocess.check_output([samtools, "view", "-H", str(args.bam)], text=True)
    (out / "input_bam_header.sam").write_text(header)
    names = sorted(records)
    provenance = {
        "script_version": VERSION, "script_sha256": sha(Path(__file__)),
        "hisat2_version": version.strip(), "hisat2_inspect_version": inspector_version.strip(),
        "samtools_version": subprocess.check_output([samtools, "--version"], text=True, errors="replace"),
        "python_version": sys.version, "sensitive_expansion": line.strip(),
        "threads": THREADS, "bam": str(args.bam), "bam_sha256": sha(args.bam),
        "bam_header_sha256": hashlib.sha256(header.encode()).hexdigest(),
        "qname_set_sha256": hashlib.sha256(("\n".join(names) + "\n").encode()).hexdigest(),
        "record_content_sha256_sorted": hashlib.sha256(
            "".join(f'{n}\t{records[n]["seq"]}\t{records[n]["qual"]}\n' for n in names).encode()).hexdigest(),
        "pinned_input_sha256": input_hashes, "index_sequences_verified": True,
        "original_read_check": original_check,
        "GT_use": "Validation and post-alignment evaluation only; never a removal rule.",
    }
    write_json(out / "provenance.json", provenance)
    write_json(out / "input_diagnostics.json", input_diagnostics(records))
    write_fastq(out / "all_dehost.fastq", records)
    return hisat, index, records, provenance


def evaluate_sam(path, records):
    """Partition solely by the primary SAM mapped bit; join GT only afterwards."""
    decisions, secondary = {}, 0
    with open(path) as f:
        for line in f:
            if line.startswith("@"):
                continue
            c = line.rstrip("\n").split("\t")
            require(len(c) >= 11, "Malformed output SAM")
            name, flag = c[0], int(c[1])
            require(name in records, f"Unknown output ID: {name}")
            require(not flag & 1, f"Unexpected paired output: {name}")
            if flag & (0x100 | 0x800):
                secondary += 1
                continue
            require(name not in decisions, f"Duplicate primary output ID: {name}")
            seq, qual = c[9], c[10]
            if flag & 16:
                seq, qual = seq.translate(RC)[::-1], qual[::-1]
            require(seq == records[name]["seq"] and qual == records[name]["qual"],
                    f"Output sequence/quality not preserved: {name}")
            tags = {x.split(":", 2)[0]: x.split(":", 2)[2] for x in c[11:]}
            mapped = not bool(flag & 4)
            if mapped:
                require(c[2] != "*" and c[5] != "*" and "AS" in tags and "NM" in tags,
                        f"Incomplete mapped alignment: {name}")
            decisions[name] = {"mapped": mapped, "flag": flag, "rname": c[2], "pos": int(c[3]),
                               "mapq": int(c[4]), "cigar": c[5], "AS": tags.get("AS", ""),
                               "NM": tags.get("NM", "")}
    require(set(decisions) == set(records), "Missing primary output IDs")
    removed = {n for n, d in decisions.items() if d["mapped"]}
    retained = set(decisions) - removed
    require(not removed & retained and removed | retained == set(records), "Partition failure")
    counts = {}
    for tid in EXPECTED:
        all_names = {n for n in records if records[n]["taxid"] == tid}
        counts[tid] = {"input": len(all_names), "removed": len(all_names & removed),
                       "retained": len(all_names & retained)}
    human_left = counts["9606"]["retained"]
    total_left = len(retained)
    result = {
        "input_reads": len(records), "mapped_primary_removed": len(removed),
        "unmapped_primary_retained": total_left, "secondary_or_supplementary_records": secondary,
        "by_source_taxid": counts, "human_residual_reads": human_left,
        "human_fraction_residual": human_left / total_left if total_left else None,
        "all_ids_partition_once": True,
    }
    return result, decisions, retained


def measure_entry(argv):
    path, command = argv[0], argv[1:]
    start = time.monotonic()
    rc = subprocess.run(command).returncode
    r = resource.getrusage(resource.RUSAGE_CHILDREN)
    write_json(path, {"command": command, "exit_code": rc,
                     "wall_seconds": time.monotonic() - start,
                     "max_child_rss_kib_linux": r.ru_maxrss,
                     "rss_note": "Largest child-process high-water RSS, not summed concurrent RSS.",
                     "user_seconds": r.ru_utime, "system_seconds": r.ru_stime})
    return rc


def align(hisat, index, records, out, sensitive):
    out.mkdir()
    fastq = out / "input.fastq"
    write_fastq(fastq, records)
    sam = out / "aligned.sam"
    command = [hisat, "-x", str(index), "-U", str(fastq),
               "--no-spliced-alignment", "-p", str(THREADS),
               "--summary-file", str(out / "hisat2_summary.txt"), "-S", str(sam)]
    if sensitive:
        command.append("--sensitive")
    write_json(out / "command.json", command)
    with open(out / "stdout.log", "w") as stdout, open(out / "stderr.log", "w") as stderr:
        subprocess.run([sys.executable, str(Path(__file__).resolve()), "_measure",
                        str(out / "resources.json"), *command], stdout=stdout, stderr=stderr, check=True)
    metrics = json.loads((out / "resources.json").read_text())
    require(metrics["exit_code"] == 0, "Aligner failed")
    result, decisions, retained = evaluate_sam(sam, records)
    result["resources"] = metrics
    result["sam_bytes"] = sam.stat().st_size
    with open(out / "read_decisions.tsv", "w", newline="") as f:
        fields = ["qname", "source_taxid", "mapped", "flag", "rname", "pos", "mapq", "cigar", "AS", "NM"]
        writer = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for n, d in decisions.items():
            writer.writerow({"qname": n, "source_taxid": records[n]["taxid"], **d})
    write_fastq(out / "residual.fastq", records, [n for n in records if n in retained])
    write_json(out / "result.json", result)
    return result


def pilot_gate(pilot):
    scale = sum(EXPECTED.values()) / PILOT_N
    # Conservative screening, not a claim that max-child RSS equals aggregate RAM.
    projected_rss_gb = 4 * pilot["resources"]["max_child_rss_kib_linux"] / 1024**2
    projected_wall_s = pilot["resources"]["wall_seconds"] * scale
    projected_sam_mb = pilot["sam_bytes"] * scale / 1024**2
    require(projected_rss_gb <= 5, "Pilot memory projection >5GB; review resources before continuing")
    require(projected_wall_s <= 120, "Pilot projected runtime >120s; stop before full run")
    require(projected_sam_mb <= 100, "Pilot projected SAM >100MB; inspect before continuing")
    return {"rss_screen_4x_max_child_gb": projected_rss_gb,
            "linear_wall_projection_s": projected_wall_s, "linear_sam_projection_mb": projected_sam_mb}


def candidate_pass(result):
    frac = result["human_fraction_residual"]
    return (frac is not None and frac <= 0.05
            and all(result["by_source_taxid"][t]["retained"] == 1000
                    for t in EXPECTED if t != "9606"))


def summary_text(data):
    lines = [
        f"SE_DEHOST_RESCUE_V1 {VERSION}",
        f"status: {data['status']}",
        "experiment: toy-only, SE50, HISAT2 2.2.3, 2 threads",
        "baseline: 3516 reads; human=516; Ecoli=1000; yeast=1000; phiX=1000",
        f"baseline_GT_human_percent: {100 * 516 / 3516:.8f}",
        "rule: remove mapped primary SE reads; GT used for evaluation only",
        "criterion: residual GT human <=5% AND retain all 1000 reads of EACH microbe",
        "candidate: --sensitive = --bowtie2-dp 1 -k 30 --score-min L,0,-0.5",
    ]
    if "provenance" in data:
        p = data["provenance"]
        lines += [f"input_bam_sha256: {p['bam_sha256']}",
                  f"input_qname_set_sha256: {p['qname_set_sha256']}",
                  f"input_record_content_sha256_sorted: {p['record_content_sha256_sorted']}",
                  f"script_sha256: {p['script_sha256']}",
                  "pinned_inputs_and_index_sequences_verified: true",
                  f"original_read_validation: {json.dumps(p['original_read_check'])}"]
    for kind in ("control_pilot", "control_full", "candidate_pilot", "candidate_full"):
        if kind not in data:
            continue
        r = data[kind]
        frac = r["human_fraction_residual"]
        lines += ["", f"[{kind}]", f"input_reads: {r['input_reads']}",
                  f"mapped_removed: {r['mapped_primary_removed']}",
                  f"retained: {r['unmapped_primary_retained']}",
                  f"residual_GT_human_percent: {100 * frac:.8f}" if frac is not None
                  else "residual_GT_human_percent: undefined (empty residual)",
                  "source_taxid input removed retained"]
        lines += [f"{t} {c['input']} {c['removed']} {c['retained']}"
                  for t, c in r["by_source_taxid"].items()]
        lines += [f"wall_seconds: {r['resources']['wall_seconds']:.4f}",
                  f"max_child_rss_MiB: {r['resources']['max_child_rss_kib_linux']/1024:.3f}",
                  f"all_IDs_partition_once: {r['all_ids_partition_once']}"]
    if "candidate_full" in data:
        r = data["candidate_full"]
        ok = candidate_pass(r)
        lines += ["", f"predeclared_toy_criteria_met: {str(ok).lower()}",
                  f"human_removed_of_516: {r['by_source_taxid']['9606']['removed']}",
                  "recommendation: " + ("TOY_CANDIDATE_ONLY; no automatic pipeline change"
                                        if ok else "DO_NOT_DEPLOY; no further tuning in this experiment")]
    if data.get("error"):
        lines += ["", f"error: {data['error']}"]
    lines += ["", "Baseline files, QC thresholds and core pipeline were not changed.",
              "This is not validation on real samples, a full human genome or a production MetaX DB.",
              "A local run is not proof of current server BAM byte identity or independent replication.",
              "RAM metric is max child RSS, not summed concurrent RSS.", ""]
    return "\n".join(lines)


def cli():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, default=Path.cwd(), help="Repository root (default: cwd)")
    parser.add_argument("--toy", type=Path, help="Existing toy_995 directory")
    parser.add_argument("--bam", type=Path, help="Existing SE DEHOST BAM")
    parser.add_argument("--out", type=Path, help="NEW output directory; must not exist")
    parser.add_argument("--hisat2", help="HISAT2 executable (default: PATH)")
    parser.add_argument("--samtools", help="samtools executable (default: PATH)")
    parser.add_argument("--inspect", help="hisat2-inspect executable (default: PATH)")
    args = parser.parse_args()
    args.repo = args.repo.resolve()
    args.toy = (args.toy or args.repo / "insilico_runs/toy_995").resolve()
    args.bam = (args.bam or args.repo / (
        "insilico_runs/metax_samovar_v1/results_se/SAMOVAR_SE/validation/partitions/dehost.bam")).resolve()
    out = (args.out or args.repo / "insilico_runs/se_dehost_rescue_v1").resolve()
    # Refuse overwrite even after a previous failure. Use a new --out for an explicitly agreed retry.
    require(not out.exists(), f"Output exists; refusing overwrite: {out}")
    require(not args.bam.is_relative_to(out) and not args.toy.is_relative_to(out),
            "Output must not contain existing input locations")
    out.mkdir(parents=True, exist_ok=False)
    data = {"status": "STARTED"}

    def checkpoint():
        write_json(out / "experiment.json", data)
        (out / "summary.txt").write_text(summary_text(data))

    try:
        hisat, index, records, provenance = preflight(args, out)
        data["provenance"] = provenance
        checkpoint()
        pilot_records = dict(list(records.items())[:PILOT_N])
        data["control_pilot"] = align(hisat, index, pilot_records, out / "control_pilot", False)
        checkpoint()
        require(data["control_pilot"]["mapped_primary_removed"] == 0,
                "Control pilot maps original DEHOST reads: replay not reproduced; candidate NOT RUN")
        data["control_pilot_projection"] = pilot_gate(data["control_pilot"])
        data["control_full"] = align(hisat, index, records, out / "control_full", False)
        checkpoint()
        require(data["control_full"]["mapped_primary_removed"] == 0,
                "Full control maps original DEHOST reads: replay not reproduced; candidate NOT RUN")
        data["candidate_pilot"] = align(hisat, index, pilot_records, out / "candidate_pilot", True)
        checkpoint()
        data["candidate_pilot_projection"] = pilot_gate(data["candidate_pilot"])
        data["candidate_full"] = align(hisat, index, records, out / "candidate_full", True)
        data["status"] = "COMPLETED"
        data["predeclared_toy_criteria_met"] = candidate_pass(data["candidate_full"])
        checkpoint()
        print(summary_text(data))
        print(f"Send only: {out / 'summary.txt'}")
        return 0
    except Exception as exc:
        data["status"] = "STOPPED"
        data["error"] = f"{type(exc).__name__}: {exc}"
        checkpoint()
        print(summary_text(data), file=sys.stderr)
        print(f"Diagnostic summary: {out / 'summary.txt'}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "_measure":
        sys.exit(measure_entry(sys.argv[2:]))
    sys.exit(cli())
