#!/bin/bash
# Merge the 4 cohorts' common-SNP bfiles into one pooled bfile, then compute
# frequency once on the pooled sample (needed for makeScoreFile.R and for
# running --score a single time across everyone).
set -euo pipefail

module load plink/1.9-20210416

WDIR="/project/knathans_tecac/TGCT_PRS_jenny"
OUT="$WDIR/pooled/tecac_prs_pooled"
mkdir -p "$WDIR/pooled"

BFILES=(
  "$WDIR/discovery/penn_prs_common"
  "$WDIR/replication/rep_prs_common"
  "$WDIR/gencove/gencove_prs_common"
  "$WDIR/pmbb4/pmbb4_prs_merged"
)

first="${BFILES[0]}"
merge_list="$WDIR/pooled/merge_list.txt"
> "$merge_list"
for f in "${BFILES[@]:1}"; do
  echo "${f}.bed ${f}.bim ${f}.fam" >> "$merge_list"
done

plink --bfile "$first" --merge-list "$merge_list" --make-bed --allow-no-sex --out "$OUT"

# if plink reports conflicts, it writes ${OUT}-merge.missnp
if [[ -f "${OUT}-merge.missnp" ]]; then
    echo "Conflicts found, excluding and retrying"
    for f in "${BFILES[@]}"; do
        plink --bfile "$f" --exclude "${OUT}-merge.missnp" --make-bed --allow-no-sex --out "${f}_clean"
    done
    first="${BFILES[0]}_clean"
    > "$merge_list"
    for f in "${BFILES[@]:1}"; do
        echo "${f}_clean.bed ${f}_clean.bim ${f}_clean.fam" >> "$merge_list"
    done
    plink --bfile "$first" --merge-list "$merge_list" --make-bed --allow-no-sex --out "$OUT"
fi

plink --bfile "$OUT" --freq --out "$OUT"

score_file="$WDIR/gwas2021_rep_scores_b38.txt"
plink --bfile "$OUT" \
      --score "$score_file" 1 2 7 header sum \
      --out "${OUT}_score"

echo "Pooled bfile: ${OUT}.bed/.bim/.fam"
echo "Pooled freq: ${OUT}.frq"
echo "Pooled PRS: ${OUT}_score.profile"
wc -l "${OUT}.bim"
