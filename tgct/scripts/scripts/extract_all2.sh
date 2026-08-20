#!/bin/bash
set -euo pipefail

module load plink/1.9-20210416

WDIR="/home/jennyzli/TGCT_PRS_jenny2"
SNP_LIST="$WDIR/gwas2021_prs_ids.txt"
SNP_LIST23="$WDIR/gwas2021_prs_ids_23.txt"

mkdir -p "$WDIR/gencove" "$WDIR/replication/extract"

# extract gencove
GENCOVE=/project/knathans_tecac/gencove/all/gencove_qc6.lifted.hg38.final

plink --bfile "$GENCOVE" --extract "$SNP_LIST23" --make-bed --out "$WDIR/gencove/gencove_ext"

# extract replication
REP_DIR=/project/knathans_tecac/replication_vcfs
MERGE_LIST="$WDIR/replication/merge_list.txt"
> "$MERGE_LIST"

# each chromosome is its own subfolder containing a plink fileset, e.g.
# replication_vcfs/chr12/chr12_males.qc2.lifted.hg38.final.bim
shopt -s nullglob extglob
rep_dirs=("$REP_DIR"/chr@(+([0-9])|X|Y|M))
shopt -u nullglob extglob
[[ ${#rep_dirs[@]} -gt 0 ]] || { echo "ERROR: no directories matched $REP_DIR/chr<N>" >&2; exit 1; }

first=""
for dir in "${rep_dirs[@]}"; do
    chr=$(basename "$dir")
    out="$WDIR/replication/extract/${chr}_extract"
    prefix="$dir/${chr}_males.qc2.lifted.hg38.final"

    if [[ ! -f "${prefix}.bim" ]]; then
        echo "WARNING: missing ${prefix}.bim, skipping" >&2
        continue
    fi

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

[[ -n "$first" ]] || { echo "ERROR: no replication chromosome had matching SNPs" >&2; exit 1; }

if [[ -s "$MERGE_LIST" ]]; then
    plink --bfile "$first" --merge-list "$MERGE_LIST" --make-bed --out "$WDIR/replication/rep_ext"
else
    for ext in bed bim fam; do
        cp "${first}.${ext}" "$WDIR/replication/rep_ext.${ext}"
    done
fi


# discovery already done 

# pmbb4 in another script
