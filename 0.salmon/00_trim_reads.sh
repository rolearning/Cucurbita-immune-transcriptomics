#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob

project_root=${ANTI_PROJECT_ROOT:-"$(cd "$(dirname "$0")/.." && pwd)"}
cd "$project_root/0.salmon"

mkdir -p cleanFastq/Cm cleanFastq/Ls cleanFastq/Ta

for sp in Cm Ls Ta; do
    read1_files=(Fastq/"${sp}"_*_R1.fastq.gz)
    ((${#read1_files[@]} > 0)) || { echo "No R1 reads found for $sp in Fastq/" >&2; exit 1; }

    for r1 in "${read1_files[@]}"; do
        sample=$(basename "$r1" _R1.fastq.gz)
        r2="Fastq/${sample}_R2.fastq.gz"
        [ -f "$r2" ] || { echo "Missing paired read: $r2" >&2; exit 1; }

        echo "Trimming: $sample"
        trim_galore \
          -j 8 \
          -q 25 \
          --phred33 \
          --stringency 3 \
          --length 20 \
          --fastqc \
          --paired \
          "$r1" \
          "$r2" \
          --gzip \
          --output_dir "cleanFastq/${sp}"
    done
done
