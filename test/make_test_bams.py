#!/usr/bin/env python3
"""
make_test_bams.py — генерация BAM-фикстур для критериев приёмки задачи 3.

TOY_SORTED (критерий 1): те же прочтения, что в TOY_NORMAL.{dehost,host}.bam,
но картированные hisat2 и отсортированные по координате — мейты перестают
быть соседними записями (модель координатно-сортированного BAM от NIPD).
Набор прочтений сохраняется, поэтому профиль TOY_SORTED обязан совпасть
с профилем TOY_NORMAL. Скрипт проверяет нетривиальность фикстуры: доля
соседних записей-мейтов должна быть мала (иначе тест ничего не проверяет).

TOY_SE (критерий 2): те же прочтения со сброшенными флагами пары
(0x1 PAIRED, 0x2 PROPER_PAIR, 0x8 MUNMAP, 0x40 READ1, 0x80 READ2) —
истинный single-end BAM.

Запуск (из корня репозитория, окружение с samtools и hisat2 на PATH):
  python3 test/make_test_bams.py
"""

import subprocess
import sys
import tempfile
from pathlib import Path

TOY = Path(__file__).resolve().parent / "toy"
BAMS = TOY / "bams"
REFS = TOY / "refs"
HOST_INDEX = TOY / "host_index" / "human"

# флаги, сбрасываемые при конверсии в single-end
CLEAR_FLAGS = 0x1 | 0x2 | 0x8 | 0x40 | 0x80


def run(cmd, **kw):
    print("+", " ".join(str(c) for c in cmd))
    res = subprocess.run([str(c) for c in cmd], **kw)
    if res.returncode != 0:
        sys.exit(f"команда завершилась с кодом {res.returncode}: {cmd}")


def extract_fastas_simple(bam, prefix, workdir):
    """BAM -> парные FASTA (collate|fasta, как в модуле BAM2FASTA)."""
    r1, r2 = workdir / f"{prefix}_1.fasta", workdir / f"{prefix}_2.fasta"
    cmd = (f"samtools collate -u -O {bam} | samtools fasta "
           f"-1 {r1} -2 {r2} -0 {workdir}/{prefix}_0.fasta "
           f"-s {workdir}/{prefix}_s.fasta -")
    run(["bash", "-c", cmd])
    return r1, r2


def align_and_sort(index_prefix, r1, r2, out_bam):
    """hisat2 (FASTA-вход) -> координатно-сортированный BAM.

    Набор прочтений сохраняется: hisat2 выдаёт по записи на прочтение
    (некартированные остаются в BAM как unmapped).
    """
    cmd = (f"hisat2 -x {index_prefix} -f -1 {r1} -2 {r2} -p 4 2>/dev/null | "
           f"samtools sort -O BAM -o {out_bam} -")
    run(["bash", "-c", cmd])


def strip_paired_flags(src_bam, dst_bam):
    """BAM -> single-end BAM: сброс флагов пары, RNEXT/PNEXT/TLEN в ноль.

    Имена делаются уникальными суффиксами /1,/2 (по исходным флагам READ1/
    READ2): в настоящем SE BAM каждое прочтение имеет своё имя, а дубли
    QNAME samtools fasta хэш-пейрингует обратно в «пары».
    """
    view = subprocess.Popen(["samtools", "view", "-h", str(src_bam)],
                            stdout=subprocess.PIPE, text=True)
    conv = subprocess.Popen(["samtools", "view", "-b", "-o", str(dst_bam), "-"],
                            stdin=subprocess.PIPE, text=True)
    for line in view.stdout:
        if line.startswith("@"):
            conv.stdin.write(line)
            continue
        f = line.rstrip("\n").split("\t")
        flag = int(f[1])
        if flag & 0x40:
            f[0] += "/1"
        elif flag & 0x80:
            f[0] += "/2"
        f[1] = str(flag & ~CLEAR_FLAGS)
        f[6], f[7], f[8] = "*", "0", "0"
        conv.stdin.write("\t".join(f) + "\n")
    conv.stdin.close()
    view.stdout.close()
    if conv.wait() != 0 or view.wait() != 0:
        sys.exit(f"strip_paired_flags failed for {src_bam}")


def fasta_pair_reads(bam, workdir, tag, use_collate):
    """Число прочтений в -1/-2 после samtools fasta (с collate или без).

    Нетривиальность TOY_SORTED: экстракция БЕЗ collate (старый код) обязана
    терять прочтения (samtools fasta пейрингует по хэшу имён с ограниченным
    кэшем; на координатно-сортированном входе часть пар превращается в
    «singletons» и не попадает в -1/-2), а С collate — восстанавливать все.
    """
    r1, r2 = workdir / f"{tag}_1.fa", workdir / f"{tag}_2.fa"
    se = workdir / f"{tag}_se.fa"
    if use_collate:
        cmd = (f"samtools collate -u -O {bam} | samtools fasta "
               f"-1 {r1} -2 {r2} -0 {se} -s {se}.s - 2>/dev/null; "
               f"cat {se}.s >> {se}")
    else:
        cmd = f"samtools fasta -1 {r1} -2 {r2} -0 {se} -s {se}.s {bam} 2>/dev/null"
    run(["bash", "-c", cmd])

    def n_fasta(p):
        with open(p) as fh:
            return sum(1 for ln in fh if ln.startswith(">"))
    return n_fasta(r1) + n_fasta(r2), n_fasta(se) + n_fasta(Path(str(se) + ".s"))


def paired_count(bam):
    res = subprocess.run(["samtools", "view", "-c", "-F", "0x900", "-f", "0x1",
                          str(bam)], capture_output=True, text=True)
    return int(res.stdout.strip())


def primary_count(bam):
    res = subprocess.run(["samtools", "view", "-c", "-F", "0x900", str(bam)],
                         capture_output=True, text=True)
    return int(res.stdout.strip())


def main():
    BAMS.mkdir(exist_ok=True)
    with tempfile.TemporaryDirectory() as td:
        work = Path(td)

        # комбинированный не-human индекс для dehost-фракции
        nonhuman = work / "nonhuman.fa"
        with open(nonhuman, "wb") as out:
            for ref in ("myco.fa", "b19.fa", "ecuniculi.fa"):
                out.write((REFS / ref).read_bytes())
        run(["hisat2-build", "-p", "4", nonhuman, work / "nonhuman_idx"])

        for frac, index in (("dehost", work / "nonhuman_idx"),
                            ("host", HOST_INDEX)):
            src = BAMS / f"TOY_NORMAL.{frac}.bam"
            r1, r2 = extract_fastas_simple(src, frac, work)
            align_and_sort(index, r1, r2, BAMS / f"TOY_SORTED.{frac}.bam")

        for frac in ("dehost", "host"):
            strip_paired_flags(BAMS / f"TOY_NORMAL.{frac}.bam",
                               BAMS / f"TOY_SE.{frac}.bam")

    # ---- валидация фикстур ----
    print("\n== валидация ==")
    with tempfile.TemporaryDirectory() as td2:
        work2 = Path(td2)
        for frac in ("dehost", "host"):
            bam = BAMS / f"TOY_SORTED.{frac}.bam"
            n_primary = primary_count(bam)
            n_new, _ = fasta_pair_reads(bam, work2, f"new_{frac}", True)
            n_old, _ = fasta_pair_reads(bam, work2, f"old_{frac}", False)
            print(f"TOY_SORTED.{frac}: primary={n_primary} "
                  f"extract(collate)={n_new} extract(без collate)={n_old} "
                  f"(потеря старого кода {100 * (n_new - n_old) / n_new:.1f}%)")
            if n_new != n_primary:
                sys.exit(f"ФИКСТУРА БИТАЯ: {bam} — collate|fasta вернул "
                         f"{n_new} вместо {n_primary} primary-прочтений")
            if n_old >= n_new:
                sys.exit(f"ФИКСТУРА ТРИВИАЛЬНА: {bam} — старый код ничего "
                         "не теряет, критерий 1 ничего не проверял бы")
        for frac in ("dehost", "host"):
            bam = BAMS / f"TOY_SE.{frac}.bam"
            if paired_count(bam) != 0:
                sys.exit(f"ФИКСТУРА БИТАЯ: {bam} — остались PAIRED-флаги")
            _, n_se = fasta_pair_reads(bam, work2, f"se_{frac}", True)
            print(f"TOY_SE.{frac}: primary={primary_count(bam)} "
                  f"paired=0 se_extracted={n_se}")
            if n_se != primary_count(bam):
                sys.exit(f"ФИКСТУРА БИТАЯ: {bam} — SE-экстракция вернула "
                         f"{n_se} вместо {primary_count(bam)}")
    print("фикстуры валидны")


if __name__ == "__main__":
    main()
