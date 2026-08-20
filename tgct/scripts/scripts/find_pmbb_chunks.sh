
#!/bin/bash
# Find manifest chunks that actually contain one of the target PRS SNP positions.

WDIR="/home/jennyzli/TGCT_PRS_jenny2"
IMPUTED_DATA_DIR="/static/PMBB/PMBB-Release-2026-4.0/Imputed"
input_manifest="${IMPUTED_DATA_DIR}/metadata/imputed_variant_chunked_input_manifest.tsv"
snp_list="${WDIR}/gwas2021_prs_ids_pmbb.txt"

cd "$WDIR"

awk -F'_' '{print $1"_"$2}' "$snp_list" | sort -u > pmbb4/snp_chrpos.txt

> pmbb4/relevant_chunks.tsv
while IFS='_' read -r chr pos; do
  awk -v c="$chr" -v p="$pos" \
    'NR>1 && $2==c && p>=$4 && p<=$5 {print $0}' \
    "$input_manifest" >> pmbb4/relevant_chunks.tsv
done < pmbb4/snp_chrpos.txt

wc -l pmbb4/relevant_chunks.tsv
