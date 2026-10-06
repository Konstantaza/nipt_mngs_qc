#!/usr/bin/env python3
"""
make_toy_data.py — симулятор парных прочтений для smoke-test пайплайна.

Генерирует синтетические paired-end FASTQ (2x150) из заданных референсов
с заданными пропорциями. Фрагменты ~ N(170, 30) п.н. (подражание вкДНК),
ошибки секвенирования ~0.1%, качество Q36.

Использование:
  make_toy_data.py --ref human.fa:0.97 --ref myco.fa:0.02 ... \
      --n_pairs 100000 --out TOY_NORMAL --seed 42
"""

import argparse
import gzip
import random


def read_fasta(path):
    seqs = []
    name, buf = None, []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                if name:
                    seqs.append("".join(buf).upper())
                name, buf = line[1:], []
            else:
                buf.append(line)
    if name:
        seqs.append("".join(buf).upper())
    return seqs


RC = str.maketrans("ACGTN", "TGCAN")


def revcomp(s):
    return s.translate(RC)[::-1]


def simulate_pair(ref, rng, read_len, frag_mu=170, frag_sd=30, err=0.001):
    L = len(ref)
    flen = max(read_len + 1, int(rng.gauss(frag_mu, frag_sd)))
    flen = min(flen, L)
    start = rng.randrange(0, L - flen + 1)
    frag = ref[start:start + flen]
    r1 = list(frag[:read_len])
    r2 = list(revcomp(frag[-read_len:]))
    for read in (r1, r2):
        for i in range(len(read)):
            if rng.random() < err:
                read[i] = rng.choice("ACGT")
    return "".join(r1), "".join(r2)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ref", action="append", required=True,
                    help="fasta:пропорция (можно несколько)")
    ap.add_argument("--n_pairs", type=int, default=100000)
    ap.add_argument("--read_len", type=int, default=150)
    ap.add_argument("--out", required=True)
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()

    rng = random.Random(args.seed)
    refs, props = [], []
    for spec in args.ref:
        path, p = spec.rsplit(":", 1)
        refs.append(read_fasta(path))
        props.append(float(p))
    total = sum(props)
    props = [p / total for p in props]

    q = "I" * args.read_len  # Q36
    with gzip.open(f"{args.out}_R1.fastq.gz", "wt") as f1, \
         gzip.open(f"{args.out}_R2.fastq.gz", "wt") as f2:
        for i in range(args.n_pairs):
            k = rng.choices(range(len(refs)), weights=props)[0]
            ref = rng.choice(refs[k])
            r1, r2 = simulate_pair(ref, rng, args.read_len)
            f1.write(f"@{args.out}_{i}/1\n{r1}\n+\n{q}\n")
            f2.write(f"@{args.out}_{i}/2\n{r2}\n+\n{q}\n")
    print(f"{args.out}: {args.n_pairs} пар, пропорции {props}")


if __name__ == "__main__":
    main()
