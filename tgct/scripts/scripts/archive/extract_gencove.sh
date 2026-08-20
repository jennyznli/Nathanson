╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌
#!/bin/bash
set -euo pipefail

IMPUTED_DIR=/project/knathans_tecac/gencove/all/gencove_qc6.lifted.hg38.final
SNP_LIST=/project/knathans_tecac/TGCT_PRS_jenny/gwas2021_prs_ids.txt
OUT_DIR=/project/knathans_tecac/TGCT_PRS_jenny/gencove

plink --bfile "$IMPUTED_DIR" --extract "$SNP_LIST" --make-bed --out "$OUT_DIR/gencove_prs_merged"
