#!/bin/bash
# Merge the per-chunk PMBB Freeze4 PRS SNP extractions into one bed/bim/fam.
# Run after extract_pmbb_prs_chunks_loop.bsub (or the array version) has finished.
set -euo pipefail

module load plink/1.9-20210416

WDIR="/project/knathans_tecac/TGCT_PRS_jenny"
CHUNK_DIR="$WDIR/pmbb4/prs_chunk_bed"
OUT="$WDIR/pmbb4/pmbb4_prs_merged"

cd "$CHUNK_DIR"

# only chunks that actually produced output have a .bim here
prefixes=($(ls pmbb4_prs_*.bim 2>/dev/null | sed 's/\.bim$//'))

if [[ ${#prefixes[@]} -eq 0 ]]; then
    echo "No extracted chunk files found in $CHUNK_DIR"
    exit 1
fi

echo "Found ${#prefixes[@]} chunk files to merge"

if [[ ${#prefixes[@]} -eq 1 ]]; then
    plink --bfile "${prefixes[0]}" --make-bed --out "$OUT"
else
    first="${prefixes[0]}"
    merge_list="$CHUNK_DIR/merge_list.txt"
    > "$merge_list"
    for p in "${prefixes[@]:1}"; do
        echo "${p}.bed ${p}.bim ${p}.fam" >> "$merge_list"
    done

    plink --bfile "$first" \
          --merge-list "$merge_list" \
          --make-bed \
          --out "$OUT"
fi

echo "Merged bed/bim/fam written to ${OUT}.*"
wc -l "${OUT}.bim"
