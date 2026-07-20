#!/bin/bash
# Recode the merged PMBB PRS .bim SNP ID column from chr1_156233018_G_C
# to 1:156233018 (chr:pos, no alleles) to match the other 3 cohorts' convention.
set -euo pipefail

WDIR="/project/knathans_tecac/TGCT_PRS_jenny"
BFILE="$WDIR/pmbb4/pmbb4_prs_merged"

cp "${BFILE}.bim" "${BFILE}.bim.orig"

awk 'BEGIN{OFS="\t"} {
    n = split($2, a, "_")
    chr = a[1]
    sub(/^chr/, "", chr)
    $2 = chr":"a[2]
    print
}' "${BFILE}.bim.orig" > "${BFILE}.bim"

echo "Recoded $(wc -l < "${BFILE}.bim") variant IDs. Original saved as ${BFILE}.bim.orig"
head "${BFILE}.bim"
