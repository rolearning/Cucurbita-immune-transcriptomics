#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob

project_root=${ANTI_PROJECT_ROOT:-"$(cd "$(dirname "$0")/.." && pwd)"}
cd "$project_root/0.salmon"

mkdir -p quant_decoy/Cm quant_decoy/Ls quant_decoy/Ta

# Validate annotations before starting potentially lengthy quantification.
for sp in Cm Ls Ta; do
    gff="../ref/${sp}.gff3"
    [ -s "$gff" ] || { echo "Missing annotation file: $gff" >&2; exit 1; }
done

# Quantify paired reads, retaining the existing completed-output skip rule.
for sp in Cm Ls Ta; do
    read1_files=(cleanFastq/"${sp}"/*_R1_val_1.fq.gz)
    ((${#read1_files[@]} > 0)) || { echo "No trimmed R1 reads found for $sp in cleanFastq/${sp}/" >&2; exit 1; }

    for r1 in "${read1_files[@]}"; do
        sample=$(basename "$r1" _R1_val_1.fq.gz)
        out="quant_decoy/${sp}/${sample}"
        if [ -s "$out/quant.sf" ]; then
            echo "Skip: ${sample}"
            continue
        fi

        r2="cleanFastq/${sp}/${sample}_R2_val_2.fq.gz"
        [ -f "$r2" ] || { echo "Missing paired read: $r2" >&2; exit 1; }
        index="index_decoy/${sp}"
        [ -d "$index" ] && [ -n "$(find "$index" -mindepth 1 -maxdepth 1 -print -quit)" ] || {
            echo "Missing or empty Salmon index: $index" >&2
            exit 1
        }

        echo "Running: ${sample}"
        salmon quant \
          -i "$index" \
          -l ISR \
          -1 "$r1" \
          -2 "$r2" \
          -p 16 \
          -o "$out"
    done
done

# Generate transcript-to-gene mappings from each species GFF3.
for sp in Cm Ls Ta; do
    awk -F'\t' '
    BEGIN{OFS="\t"; print "TXNAME","GENEID"}
    $3=="mRNA"{
        id=""; parent="";
        n=split($9,a,";");
        for(i=1;i<=n;i++){
            if(a[i]~/^ID=/){id=a[i]; sub(/^ID=/,"",id)}
            if(a[i]~/^Parent=/){parent=a[i]; sub(/^Parent=/,"",parent)}
        }
        if(id!="" && parent!="") print id,parent
    }' "../ref/${sp}.gff3" > "../ref/${sp}_tx2gene.tsv"
    [ "$(wc -l < "../ref/${sp}_tx2gene.tsv")" -gt 1 ] || {
        echo "No transcript-to-gene pairs found in ../ref/${sp}.gff3" >&2
        exit 1
    }
done

# Summarize mapping statistics from completed Salmon runs.
meta_files=(quant_decoy/*/*/aux_info/meta_info.json)
((${#meta_files[@]} > 0)) || { echo "No Salmon meta_info.json files found under quant_decoy/" >&2; exit 1; }
printf "Species\tSample\tProcessed\tMapped\tMapping_rate\n" > 02_salmon_mapping_rates.tsv

for f in "${meta_files[@]}"; do
    sample_dir=${f%/aux_info/meta_info.json}
    species_dir=${sample_dir%/*}
    species=${species_dir##*/}
    sample=${sample_dir##*/}

    processed=$(grep '"num_processed"' "$f" | awk -F': ' '{gsub(",","",$2); print $2}')
    mapped=$(grep '"num_mapped"' "$f" | awk -F': ' '{gsub(",","",$2); print $2}')
    rate=$(grep '"percent_mapped"' "$f" | awk -F': ' '{gsub(",","",$2); print $2}')
    [ -n "$processed" ] && [ -n "$mapped" ] && [ -n "$rate" ] || {
        echo "Missing mapping statistics in $f" >&2
        exit 1
    }

    printf "%s\t%s\t%s\t%s\t%s\n" "$species" "$sample" "$processed" "$mapped" "$rate" \
      >> 02_salmon_mapping_rates.tsv
done
