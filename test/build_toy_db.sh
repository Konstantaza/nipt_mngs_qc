#!/usr/bin/env bash
# =============================================================================
# build_toy_db.sh — сборка toy-референсов, мини-БД Kraken2/Bracken и
# HISAT2-индекса «человека» для smoke-test пайплайна nipt_mngs_qc.
#
# Состав toy-мира:
#   human_5Mb.fa   — срез chr21 hg38 (5 Мб)        -> taxid 9606   (хост)
#   myco.fa        — Mycoplasma genitalium G37     -> taxid 2097   (бактерия)
#   b19.fa         — Human parvovirus B19          -> taxid 10798  (вирус)
#   ecuniculi.fa   — Encephalitozoon cuniculi chrI -> taxid 284813 (эук.-патоген)
#
# БД:
#   db_full  — все 4 организма (аналог standard Kraken2 DB)
#   db_patho — только E. cuniculi (аналог EukPathDB: эукариоты-патогены)
# =============================================================================
set -euo pipefail
export PATH="/workspace/.conda/envs/nipt_mngs/bin:$PATH"

WORK=/workspace/nipt_mngs_qc/test/toy
mkdir -p "$WORK/refs" "$WORK/db_full" "$WORK/db_patho" "$WORK/host_index"
cd "$WORK"

EFETCH="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"

# ---------- 1. Референсы ----------
if [ ! -s refs/myco.fa ]; then
  curl -s "${EFETCH}?db=nuccore&id=NC_000908.2&rettype=fasta&retmode=text" > refs/myco.fa
fi
if [ ! -s refs/b19.fa ]; then
  curl -s "${EFETCH}?db=nuccore&id=NC_000883.2&rettype=fasta&retmode=text" > refs/b19.fa
fi
if [ ! -s refs/ecuniculi.fa ]; then
  curl -s "${EFETCH}?db=nuccore&id=NC_003229.1&rettype=fasta&retmode=text" > refs/ecuniculi.fa
fi
if [ ! -s refs/human_5Mb.fa ]; then
  # срез chr21 (5 Мб) напрямую через NCBI efetch (seq_start/seq_stop)
  curl -s "${EFETCH}?db=nuccore&id=NC_000021.9&seq_start=10000001&seq_stop=15000000&rettype=fasta&retmode=text" \
    > refs/human_5Mb.fa
fi
ls -la refs/

# ---------- 2. Мини-таксономия (nodes.dmp / names.dmp) ----------
write_taxonomy () {
  local DB=$1
  mkdir -p "$DB/taxonomy"
  cat > "$DB/taxonomy/nodes.dmp" <<'EOF'
1	|	1	|	no rank	|
10239	|	1	|	superkingdom	|
10798	|	10239	|	species	|
131567	|	1	|	no rank	|
2	|	131567	|	superkingdom	|
2097	|	2	|	species	|
2759	|	131567	|	superkingdom	|
9606	|	2759	|	species	|
284813	|	2759	|	species	|
EOF
  cat > "$DB/taxonomy/names.dmp" <<'EOF'
1	|	root	|		|	scientific name	|
10239	|	Viruses	|		|	scientific name	|
10798	|	Human parvovirus B19	|		|	scientific name	|
131567	|	cellular organisms	|		|	scientific name	|
2	|	Bacteria	|		|	scientific name	|
2097	|	Mycoplasma genitalium	|		|	scientific name	|
2759	|	Eukaryota	|		|	scientific name	|
9606	|	Homo sapiens	|		|	scientific name	|
284813	|	Encephalitozoon cuniculi	|		|	scientific name	|
EOF
}

# ---------- 3. Библиотеки с kraken:taxid в заголовках ----------
make_lib () {  # $1=out fasta, $2=taxid, $3=in fasta
  python3 - "$1" "$2" "$3" <<'EOF'
import sys
out, taxid, inp = sys.argv[1], sys.argv[2], sys.argv[3]
with open(inp) as fh, open(out, "w") as oh:
    n = 0
    for line in fh:
        if line.startswith(">"):
            n += 1
            oh.write(f">seq{n}|kraken:taxid|{taxid}\n")
        else:
            oh.write(line)
EOF
}

if [ ! -s db_full/hash.k2d ]; then
  write_taxonomy db_full
  make_lib db_full/library_human.fa     9606   refs/human_5Mb.fa
  make_lib db_full/library_myco.fa      2097   refs/myco.fa
  make_lib db_full/library_b19.fa       10798  refs/b19.fa
  make_lib db_full/library_ecuniculi.fa 284813 refs/ecuniculi.fa
  for f in db_full/library_*.fa; do kraken2-build --add-to-library "$f" --db db_full; done
  kraken2-build --build --db db_full --threads 8
  bracken-build -d db_full -t 8 -k 35 -l 150
fi

if [ ! -s db_patho/hash.k2d ]; then
  write_taxonomy db_patho
  make_lib db_patho/library_ecuniculi.fa 284813 refs/ecuniculi.fa
  kraken2-build --add-to-library db_patho/library_ecuniculi.fa --db db_patho
  kraken2-build --build --db db_patho --threads 8
  bracken-build -d db_patho -t 8 -k 35 -l 150
fi

# ---------- 4. HISAT2-индекс «человека» ----------
if [ ! -s host_index/human.1.ht2 ]; then
  hisat2-build -p 8 refs/human_5Mb.fa host_index/human
fi

echo "=== Toy DBs и индекс готовы ==="
ls -la db_full db_patho host_index
