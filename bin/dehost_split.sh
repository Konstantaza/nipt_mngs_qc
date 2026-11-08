#!/usr/bin/env bash
# =============================================================================
# dehost_split.sh — разбиение primary-записей BAM на фракции host/dehost/half
#
# Вход: BAM с выравниваниями на геном человека (один проход выравнивателя,
# mode=fastq). Каждая PRIMARY-запись (-F 0x900: без secondary 0x100 и
# supplementary 0x800) попадает РОВНО в одну из трёх фракций — точная
# партиция входа (критерий приёмки 3: dehost + host + half == trimmed):
#
#   host   — оба мейта картированы    (PE: !unmap && !munmap; для PE-данных
#            эквивалентно samtools view -F 12 референсного Pipline1.sh)
#   dehost — оба мейта НЕ картированы (PE: unmap && munmap; эквивалент
#            samtools view -f 12)
#   half   — полукартированная пара: ровно один мейт картирован
#            (в hocort_split.nf такие пары молча терялись)
#
# SE-записи (без флага PAIRED) тоже покрыты партицией: mapped -> host,
# unmapped -> dehost (v1.10.0: mode=fastq поддерживает настоящий SE).
#
# Использование:
#   dehost_split.sh <in.bam> <out_prefix> [threads]
#   -> <out_prefix>.host.bam <out_prefix>.dehost.bam <out_prefix>.half.bam
# =============================================================================
set -euo pipefail

if [ "$#" -lt 2 ]; then
    echo "usage: dehost_split.sh <in.bam> <out_prefix> [threads]" >&2
    exit 2
fi

in_bam="$1"
out_prefix="$2"
threads="${3:-1}"

# host: оба мейта картированы (или SE mapped)
samtools view -@ "$threads" -b -F 0x900 \
    -e '(flag.paired && !flag.unmap && !flag.munmap) || (!flag.paired && !flag.unmap)' \
    -o "${out_prefix}.host.bam" "$in_bam"

# dehost: оба мейта некартированы (или SE unmapped)
samtools view -@ "$threads" -b -F 0x900 \
    -e '(flag.paired && flag.unmap && flag.munmap) || (!flag.paired && flag.unmap)' \
    -o "${out_prefix}.dehost.bam" "$in_bam"

# half: полукартированная пара — ровно один мейт картирован
samtools view -@ "$threads" -b -F 0x900 \
    -e 'flag.paired && ((flag.unmap && !flag.munmap) || (!flag.unmap && flag.munmap))' \
    -o "${out_prefix}.half.bam" "$in_bam"
