#!/bin/bash
# Find manifest chunks that actually contain one of the target PRS SNP positions.
# Run once on the login node before submitting extract_pmbb_prs_chunks_loop.bsub

WDIR="/project/knathans_tecac/TGCT_PRS_jenny"
IMPUTED_DATA_DIR="/static/PMBB/PMBB-Release-2026-4.0/Imputed"
input_manifest="${IMPUTED_DATA_DIR}/metadata/imputed_variant_chunked_input_manifest.tsv"
snp_list="${WDIR}/snp-extract-list"

cd "$WDIR"

awk -F'_' '{print $1"_"$2}' "$snp_list" | sort -u > snp_chrpos.txt

> relevant_chunks.tsv
while IFS='_' read -r chr pos; do
  awk -v c="$chr" -v p="$pos" \
    'NR>1 && $2==c && p>=$4 && p<=$5 {print $0}' \
    "$input_manifest" >> relevant_chunks.tsv
done < snp_chrpos.txt

wc -l relevant_chunks.tsv
