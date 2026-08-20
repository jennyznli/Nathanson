╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌
#!/bin/bash
set -euo pipefail

IMPUTED_DIR=/project/knathans_tecac/GWAS2024/pennwtcga/imputed
SNP_LIST=/project/knathans_tecac/TGCT_PRS_jenny/b38_snps
OUT_DIR=/project/knathans_tecac/TGCT_PRS_jenny/discovery

MERGE_LIST="$OUT_DIR/merge_list.txt"
> "$MERGE_LIST"

first=""
for bim in "$IMPUTED_DIR"/chr*-final-samp2.b38.bim; do
    prefix="${bim%.bim}"
    chr=$(basename "$prefix")
    out="$OUT_DIR/${chr}_extract"

    # some chromosomes have no matching SNPs -- plink exits nonzero in that
    # case, so don't let set -e kill the loop; just skip this chromosome
    plink --bfile "$prefix" --extract "$SNP_LIST" --make-bed --out "$out" || true

    [[ -s "${out}.bim" ]] || continue

    if [[ -z "$first" ]]; then
        first="$out"
    else
        echo "$out" >> "$MERGE_LIST"
    fi
done

plink --bfile "$first" --merge-list "$MERGE_LIST" --make-bed --out "$OUT_DIR/penn_prs_merged"
