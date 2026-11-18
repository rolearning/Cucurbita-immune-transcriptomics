#!/usr/bin/env bash
set -euo pipefail

project_root=${ANTI_PROJECT_ROOT:-"$(cd "$(dirname "$0")/.." && pwd)"}
cd "$project_root/0.salmon"

mkdir -p gentrome index_decoy

for sp in Cm Ls Ta; do
    mrna="../ref/${sp}.mrna"
    genome="../ref/${sp}.genome"
    [ -s "$mrna" ] || { echo "Missing transcript reference: $mrna" >&2; exit 1; }
    [ -s "$genome" ] || { echo "Missing genome reference: $genome" >&2; exit 1; }

    echo "Preparing $sp gentrome..."
    grep '^>' "$genome" | sed 's/^>//; s/ .*//' > "gentrome/${sp}.decoys.txt"
    [ -s "gentrome/${sp}.decoys.txt" ] || { echo "No decoy sequences found in $genome" >&2; exit 1; }
    cat "$mrna" "$genome" > "gentrome/${sp}.gentrome.fa"

    echo "Building $sp Salmon decoy-aware index..."
    salmon index \
      -t "gentrome/${sp}.gentrome.fa" \
      -d "gentrome/${sp}.decoys.txt" \
      -i "index_decoy/${sp}" \
      -p 16 \
      --keepDuplicates
done
