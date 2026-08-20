#!/bin/bash
# Usage: prs.sh <project_name> <score_file>
#
# score_file: tab-separated, header row required, columns:
#   1 MarkerName      chr:pos of SNP (must match .bim col2 after standardizing IDs)
#   2 EffectAllele    the allele the effect/weight applies to
#   3 AltAllele       the other allele -- used to validate that each of the 4
#                     genotype sources actually has the expected biallelic
#                     pair at that position before merging (catches things
#                     like MNP-split fragments or multiallelic sites sharing
#                     a position with the real SNP)
#   4 Effect          score value plink --score uses
#
# Each run creates/reuses $PRS_ROOT/<project_name>/ and does the full pipeline in there,
# so different score files (or the same file with a different SNP subset)
# can be run side by side as separate projects without clobbering each
# other, e.g.:
#   prs.sh gwas2021_z      gwas2021-hits-z.tsv
#   prs.sh gwas2021_effect gwas2021-hits-effect.tsv
#
# To submit as an LSF job (takes CLI args, so submit as `bsub ... command args`
# rather than `bsub < prs.sh` -- the latter can't pass positional args through):
#   mkdir -p /project/knathans_tecac/jenny/prs/logs   # one-time: LSF needs the
#                                                      # -o/-e directory to already
#                                                      # exist before job dispatch
#   bsub -J prs_gwas2021_z \
#        -o /project/knathans_tecac/jenny/prs/logs/gwas2021_z.%J.o \
#        -e /project/knathans_tecac/jenny/prs/logs/gwas2021_z.%J.e \
#        -n 1 -M 16000 -R "rusage[mem=16GB]" -W 8:00 \
#        bash prs.sh gwas2021_z gwas2021-hits-z.tsv
# LSF's .o/.e will duplicate what's already in the script's own timestamped
# log under $PRS_ROOT/logs/ (see LOGGING below) -- that's fine, the .o file
# also gets LSF's resource-usage summary appended, which the internal log doesn't have.
set -euo pipefail

module load plink/1.9-20210416
PLINK19="$(command -v plink)"
module load R/4.5

PRS_ROOT="/project/knathans_tecac/jenny/prs"

if [[ $# -ne 2 ]]; then
    echo "Usage: $0 <project_name> <score_file>" >&2
    exit 1
fi
PROJECT_NAME="$1"
SCORE_FILE="$2"

[[ -s "$SCORE_FILE" ]] || { echo "ERROR: score file not found or empty: $SCORE_FILE" >&2; exit 1; }
[[ "$(awk -F'\t' 'NR==1{print NF; exit}' "$SCORE_FILE")" -ge 4 ]] \
    || { echo "ERROR: $SCORE_FILE doesn't look tab-separated with >=4 columns (MarkerName, EffectAllele, AltAllele, Effect)" >&2; exit 1; }

# ------------------------------------------------------------
# LOGGING
# ------------------------------------------------------------
# one log file per run, named by project + timestamp so reruns of the same
# project don't clobber each other's logs; tee so it still prints live too
LOG_DIR="$PRS_ROOT/logs"
mkdir -p "$LOG_DIR"
LOG_FILE="$LOG_DIR/${PROJECT_NAME}_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "$LOG_FILE") 2>&1

log() {
    echo "[$(date +'%Y-%m-%d %H:%M:%S')] $*"
}
trap 'log "FAILED at line $LINENO (exit code $?)"' ERR

log "=== prs.sh start -- project=$PROJECT_NAME score_file=$SCORE_FILE ==="
log "Log file: $LOG_FILE"

# external source genotype locations
source "$PRS_ROOT/config/sources.sh"

WDIR="$PRS_ROOT/$PROJECT_NAME"
GENCOVE_WDIR="$WDIR/gencove"
REP_WDIR="$WDIR/replication"
REP_CHR_BED_DIR="$REP_WDIR/chr_bed"
PMBB_WDIR="$WDIR/pmbb"
PMBB_CHUNK_DIR="$PMBB_WDIR/chunk_bed"
DISCOVERY_WDIR="$WDIR/discovery"
POOLED_DIR="$WDIR/pooled"
OUT="$POOLED_DIR/$PROJECT_NAME"
# shared across all projects (unlike everything else above, which is per-project)
OUTPUT_DIR="$PRS_ROOT/output"

[[ -d "$WDIR" ]] && echo "NOTE: $WDIR already exists, outputs will be overwritten" >&2

mkdir -p "$GENCOVE_WDIR" \
    "$REP_WDIR" "$REP_CHR_BED_DIR" \
    "$PMBB_WDIR" "$PMBB_CHUNK_DIR" \
    "$DISCOVERY_WDIR" \
    "$POOLED_DIR" \
    "$OUTPUT_DIR"

cd "$WDIR"

# ------------------------------------------------------------
# SNP LISTS -- derived from the score file
# ------------------------------------------------------------
# naming: "." = ":"-separated, "_" = "_"-separated, chr first either way
SNP_LIST="$WDIR/snps_chr.pos.txt"          # chr:pos, e.g. 12:123456
SNP_LIST23="$WDIR/snps_chr.pos23.txt"      # chr:pos, chrX -> chr23 (gencove numbers X as 23)
SNP_LIST_CHRPOS="$WDIR/snps_chr_pos.txt"   # chr_pos, e.g. chr12_123456 (matches PMBB manifest format)
SNP_ALLELES="$WDIR/snps_alleles.txt"       # chr:pos<TAB>EffectAllele<TAB>AltAllele, one row per SNP

# cut first col of the score file
tail -n +2 "$SCORE_FILE" | cut -f1 > "$SNP_LIST"
# gencove's raw genotype files number the X chromosome as 23, not X
sed 's/^X:/23:/' "$SNP_LIST" > "$SNP_LIST23"
# PMBB manifest uses chr_pos (underscore, "chr" prefix) instead of chr:pos
awk -F':' '{print "chr"$1"_"$2}' "$SNP_LIST" | sort -u > "$SNP_LIST_CHRPOS"
# expected biallelic pair per position, used later to validate each source's alleles
tail -n +2 "$SCORE_FILE" | cut -f1,2,3 > "$SNP_ALLELES"

log "SNP lists built: $(wc -l < "$SNP_LIST") SNPs"

# ------------------------------------------------------------
# GENCOVE
# ------------------------------------------------------------
log "=== GENCOVE ==="
"$PLINK19" --bfile "$GENCOVE_BFILE" --extract "$SNP_LIST23" --make-bed --out "$GENCOVE_WDIR/gencove_ext"

# ------------------------------------------------------------
# REPLICATION
# ------------------------------------------------------------
log "=== REPLICATION ==="
# each replication chromosome is its own subfolder containing a plink fileset, e.g.
# replication_vcfs/chr12/chr12_males.qc2.lifted.hg38.final.bim
MERGE_LIST="$REP_WDIR/merge_list.txt"
> "$MERGE_LIST"

shopt -s nullglob extglob
rep_dirs=("$REP_DIR"/chr@(+([0-9])|X|Y|M))
shopt -u nullglob extglob
[[ ${#rep_dirs[@]} -gt 0 ]] || { echo "ERROR: no directories matched $REP_DIR/chr<N>" >&2; exit 1; }

first=""
for dir in "${rep_dirs[@]}"; do
    chr=$(basename "$dir")
    out="$REP_CHR_BED_DIR/${chr}_ext"
    prefix="$dir/${chr}_males.qc2.lifted.hg38.final"

    if [[ ! -f "${prefix}.bim" ]]; then
        echo "WARNING: missing ${prefix}.bim, skipping" >&2
        continue
    fi

    # some chromosomes have no matching SNPs -- plink exits nonzero in that
    # case, so don't let set -e kill the loop; just skip this chromosome
    "$PLINK19" --bfile "$prefix" --extract "$SNP_LIST" --make-bed --out "$out" || true

    [[ -s "${out}.bim" ]] || continue

    if [[ -z "$first" ]]; then
        first="$out"
    else
        echo "$out" >> "$MERGE_LIST"
    fi
done

[[ -n "$first" ]] || { echo "ERROR: no replication chromosome had matching SNPs" >&2; exit 1; }

if [[ -s "$MERGE_LIST" ]]; then
    "$PLINK19" --bfile "$first" --merge-list "$MERGE_LIST" --make-bed --out "$REP_WDIR/rep_ext"
else
    for ext in bed bim fam; do
        cp "${first}.${ext}" "$REP_WDIR/rep_ext.${ext}"
    done
fi
log "Replication done: $(wc -l < "$REP_WDIR/rep_ext.bim") SNPs"

# ------------------------------------------------------------
# DISCOVERY
# ------------------------------------------------------------
log "=== DISCOVERY ==="
# chr*-final-samp2.b38 plink filesets, one per chromosome, directly under
# DISCOVERY_DIR (flat files, not per-chromosome subfolders like REP_DIR).
# chrX doesn't follow that naming pattern at all -- it's just X.b38 -- so it
# never matches the glob below and is handled separately.
DISC_MERGE_LIST="$DISCOVERY_WDIR/merge_list.txt"
> "$DISC_MERGE_LIST"

first=""
disc_found=0
disc_total=0
disc_glob_hits=0
for bimfile in "$DISCOVERY_DIR"/chr*-final-samp2.b38.bim; do
    disc_glob_hits=$((disc_glob_hits + 1))
    [[ -f "$bimfile" ]] || continue  # nullglob is off, so a true no-match leaves the literal pattern here
    prefix="${bimfile%.bim}"
    log "  glob matched: $prefix"
    disc_total=$((disc_total + 1))

    chr=$(basename "$prefix" | sed -E 's/^chr([^-]+)-final-samp2\.b38$/\1/')
    out="$DISCOVERY_WDIR/chr${chr}_ext"

    # most chromosomes won't have any of our (small) SNP list on them -- that's
    # expected, not an error, so don't let set -e kill the loop when plink
    # exits nonzero for "no variants remaining after --extract"
    "$PLINK19" --bfile "$prefix" --extract "$SNP_LIST" --make-bed --out "$out" || true

    if [[ -s "${out}.bim" ]]; then
        disc_found=$((disc_found + 1))
        log "  chr${chr}: $(wc -l < "${out}.bim") matching SNP(s)"
        if [[ -z "$first" ]]; then
            first="$out"
        else
            echo "$out" >> "$DISC_MERGE_LIST"
        fi
    else
        log "  chr${chr}: no matching SNPs, skipped"
    fi
done
log "Glob \"$DISCOVERY_DIR/chr*-final-samp2.b38.bim\" matched $disc_total chromosome file(s)"

# --- chrX (stored as X.b38, unlike every other chromosome's naming) ---
# its .fam has FID=0 for every sample instead of FID==IID like the other
# chromosomes -- fix that with --update-ids before merging it with them
disc_chrx_src="$DISCOVERY_DIR/X.b38"
disc_chrx_idmap="$DISCOVERY_WDIR/chrX_ids.idmap"
disc_chrx_fixed="$DISCOVERY_WDIR/chrX_ids_fixed"

# idmap columns: old FID, old IID, new FID, new IID -- new FID := old IID
awk '{print $1, $2, $2, $2}' "${disc_chrx_src}.fam" > "$disc_chrx_idmap"
"$PLINK19" --bfile "$disc_chrx_src" --update-ids "$disc_chrx_idmap" --make-bed --out "$disc_chrx_fixed"

disc_chrx_out="$DISCOVERY_WDIR/chrX_ext"
"$PLINK19" --bfile "$disc_chrx_fixed" --extract "$SNP_LIST" --make-bed --out "$disc_chrx_out" || true
disc_total=$((disc_total + 1))

if [[ -s "${disc_chrx_out}.bim" ]]; then
    disc_found=$((disc_found + 1))
    log "  chrX: $(wc -l < "${disc_chrx_out}.bim") matching SNP(s)"
    if [[ -z "$first" ]]; then
        first="$disc_chrx_out"
    else
        echo "$disc_chrx_out" >> "$DISC_MERGE_LIST"
    fi
else
    log "  chrX: no matching SNPs, skipped"
fi

log "Discovery: $disc_found / $disc_total chromosome file(s) had at least one matching SNP"
[[ -n "$first" ]] || { echo "ERROR: no discovery chromosome had matching SNPs at all" >&2; exit 1; }

# --- subset every chromosome fileset to the common sample set before merging ---
# chrX doesn't have the same sample count as the autosomes, and plink's
# --merge-list assumes every fileset shares the same samples
disc_prefixes=("$first")
if [[ -s "$DISC_MERGE_LIST" ]]; then
    while read -r p; do disc_prefixes+=("$p"); done < "$DISC_MERGE_LIST"
fi

disc_common_samples="$DISCOVERY_WDIR/common_samples.txt"
for p in "${disc_prefixes[@]}"; do
    awk '{print $1, $2}' "${p}.fam"
done | sort | uniq -c | awk -v n="${#disc_prefixes[@]}" '$1==n{print $2, $3}' > "$disc_common_samples"

wc -l "$disc_common_samples"

for p in "${disc_prefixes[@]}"; do
    "$PLINK19" --bfile "$p" --keep "$disc_common_samples" --make-bed --out "${p}_common"
done

DISC_OUT="$DISCOVERY_WDIR/disc_ext"
first_common="${first}_common"
> "$DISC_MERGE_LIST"
for p in "${disc_prefixes[@]:1}"; do
    echo "${p}_common.bed ${p}_common.bim ${p}_common.fam" >> "$DISC_MERGE_LIST"
done

if [[ -s "$DISC_MERGE_LIST" ]]; then
    "$PLINK19" --bfile "$first_common" --merge-list "$DISC_MERGE_LIST" --make-bed --out "$DISC_OUT"
else
    for ext in bed bim fam; do
        cp "${first_common}.${ext}" "${DISC_OUT}.${ext}"
    done
fi

log "Discovery merged: $(wc -l < "${DISC_OUT}.bim") SNPs, written to ${DISC_OUT}.*"

# ------------------------------------------------------------
# PMBB
# ------------------------------------------------------------
log "=== PMBB ==="
relevant_chunks="${PMBB_WDIR}/relevant_chunks.tsv"

# make file with relevant chunks
> "$relevant_chunks"
while IFS='_' read -r chr pos; do
  awk -v c="$chr" -v p="$pos" \
    'NR>1 && $2==c && p>=$4 && p<=$5 {print $0}' \
    "$PMBB_MANIFEST" >> "$relevant_chunks"
done < "$SNP_LIST_CHRPOS"

wc -l "$relevant_chunks"

# extract chunks through loop with manifest
while read -r idx chrid chunk start_pos stop_pos idk; do
    input_prefix="${PMBB_DIR}/chunked_bed_files/PMBB-Release-2026-4.0_genetic_imputed.${idx}_${chrid}_chunk${chunk}_${start_pos}_${stop_pos}"
    output="${PMBB_CHUNK_DIR}/pmbb_prs_${chrid}_chunk${chunk}"
    temp_prs_snps="${PMBB_CHUNK_DIR}/prs_ids_${chrid}_chunk${chunk}.tmp"

    echo "Processing chunk ${chrid}_${chunk} (${start_pos}-${stop_pos})"

    if [[ ! -f "${input_prefix}.bim" ]]; then
        echo "WARNING: missing ${input_prefix}.bim, skipping chunk ${chrid}_${chunk}" >&2
        continue
    fi

    # grep exits 1 (not an error) when a chunk has none of our target SNPs --
    # that's expected and handled below, so don't let pipefail treat it as fatal
    grep -Ff "$SNP_LIST_CHRPOS" "${input_prefix}.bim" | awk '{print $2}' > "$temp_prs_snps" || true

    if [[ -s "$temp_prs_snps" ]]; then
        "$PLINK19" --bfile "$input_prefix" \
              --extract "$temp_prs_snps" \
              --make-bed \
              --out "$output"
    else
        echo "WARNING: no PRS SNPs found in chunk ${chrid}_${chunk} despite manifest overlap"
    fi
done < "$relevant_chunks"

# merge chunks
PMBB_OUT="$PMBB_WDIR/pmbb_ext"

cd "$PMBB_CHUNK_DIR"

# only chunks that actually produced output have a .bim here
prefixes=($(ls pmbb_prs_*.bim 2>/dev/null | sed 's/\.bim$//'))

if [[ ${#prefixes[@]} -eq 0 ]]; then
    echo "No extracted chunk files found in $PMBB_CHUNK_DIR"
    exit 1
fi

echo "Found ${#prefixes[@]} chunk files to merge"

if [[ ${#prefixes[@]} -eq 1 ]]; then
    "$PLINK19" --bfile "${prefixes[0]}" --make-bed --out "$PMBB_OUT"
else
    first="${prefixes[0]}"
    merge_list="$PMBB_CHUNK_DIR/merge_list.txt"
    > "$merge_list"
    for p in "${prefixes[@]:1}"; do
        echo "${p}.bed ${p}.bim ${p}.fam" >> "$merge_list"
    done

    "$PLINK19" --bfile "$first" \
          --merge-list "$merge_list" \
          --write-snplist \
          --make-bed \
          --out "$PMBB_OUT"
fi

log "PMBB merged: $(wc -l < "${PMBB_OUT}.bim") SNPs, written to ${PMBB_OUT}.*"

# ------------------------------------------------------------
# STANDARDIZE IDS + INTERSECT
# ------------------------------------------------------------
log "=== STANDARDIZE IDS + INTERSECT ==="
cd "$WDIR"

# standardize IDs to chr:pos format for all 4 datasets, validating alleles
# along the way, then intersect to keep only SNPs present in all 4
standardize_ids() {
    local prefix="$1"
    local out="$2"
    local tmp="${out}_uniqids"

    # bake alleles into the ID (chr:pos:A1:A2) so any position collision --
    # MNP-split fragments, true multiallelic sites, etc. -- becomes
    # distinguishable instead of silently sharing an ID (this is what broke
    # the pmbb merge previously: two rows both named "20:53580827"). done as
    # a plain text rewrite of the .bim's ID column, not a plink call -- .bed/
    # .fam don't reference the ID column at all (only row order/count
    # matters to them), and --update-name can't do this on its own since
    # it's a straight old-ID -> new-ID lookup and can't give two rows that
    # share an old ID two different new IDs
    awk '{chr = ($1 == 23) ? "X" : $1; $2 = chr":"$4":"$5":"$6; print}' "${prefix}.bim" > "${tmp}.bim"
    cp "${prefix}.bed" "${tmp}.bed"
    cp "${prefix}.fam" "${tmp}.fam"

    # keep only rows whose allele SET matches what the score file expects at
    # that position (order-independent -- A1/A2 order differences between
    # datasets are normal and plink's merge reconciles them automatically;
    # this is about the two alleles actually being the same two alleles)
    awk -v allele_file="$SNP_ALLELES" '
        BEGIN {
            while ((getline line < allele_file) > 0) {
                split(line, f, "\t")
                a = (f[2] < f[3]) ? f[2] f[3] : f[3] f[2]
                expected[f[1]] = a
            }
        }
        {
            chr = ($1 == 23) ? "X" : $1
            chrpos = chr":"$4
            obs = ($5 < $6) ? $5 $6 : $6 $5
            if ((chrpos in expected) && expected[chrpos] == obs) print $1":"$4":"$5":"$6
        }
    ' "${prefix}.bim" > "${out}.keep_ids"

    "$PLINK19" --bfile "$tmp" --extract "${out}.keep_ids" --make-bed --out "${out}_filtered" || true

    # collapse surviving unique IDs back down to plain chr:pos so this
    # dataset still lines up with the other 3 for the common-SNP intersection
    awk '{chr = ($1 == 23) ? "X" : $1; print $2, chr":"$4}' "${out}_filtered.bim" > "${out}.idmap"
    "$PLINK19" --bfile "${out}_filtered" --update-name "${out}.idmap" 2 1 --make-bed --out "$out"

    local n_before n_after
    n_before=$(wc -l < "${prefix}.bim")
    n_after=$(wc -l < "${out}.bim")
    log "  $(basename "$prefix"): $n_after / $n_before variant(s) kept after allele-match check"
}

standardize_ids "$PMBB_WDIR/pmbb_ext"        "$PMBB_WDIR/pmbb_ext_std"
standardize_ids "$REP_WDIR/rep_ext"           "$REP_WDIR/rep_ext_std"
standardize_ids "$GENCOVE_WDIR/gencove_ext"        "$GENCOVE_WDIR/gencove_ext_std"
standardize_ids "$DISCOVERY_WDIR/disc_ext"    "$DISCOVERY_WDIR/disc_ext_std"

awk '{print $2}' "$PMBB_WDIR/pmbb_ext_std.bim"      | sort -u > "$WDIR/pmbb_snps.sorted"
awk '{print $2}' "$REP_WDIR/rep_ext_std.bim"         | sort -u > "$WDIR/replication_snps.sorted"
awk '{print $2}' "$GENCOVE_WDIR/gencove_ext_std.bim"      | sort -u > "$WDIR/gencove_snps.sorted"
awk '{print $2}' "$DISCOVERY_WDIR/disc_ext_std.bim"  | sort -u > "$WDIR/discovery_snps.sorted"

# keep only SNPs present in all 4 (each list already deduped, so count==4 means "in every file")
cat "$WDIR/pmbb_snps.sorted" "$WDIR/replication_snps.sorted" "$WDIR/gencove_snps.sorted" "$WDIR/discovery_snps.sorted" \
  | sort | uniq -c | awk '$1==4{print $2}' > "$WDIR/common_snps.txt"

log "Common SNPs across all 4 sources: $(wc -l < "$WDIR/common_snps.txt")"

# subset each dataset to the common SNP list
"$PLINK19" --bfile "$PMBB_WDIR/pmbb_ext_std"       --extract "$WDIR/common_snps.txt" --make-bed --out "$PMBB_WDIR/pmbb_common"
"$PLINK19" --bfile "$REP_WDIR/rep_ext_std"          --extract "$WDIR/common_snps.txt" --make-bed --out "$REP_WDIR/rep_common"
"$PLINK19" --bfile "$GENCOVE_WDIR/gencove_ext_std"       --extract "$WDIR/common_snps.txt" --make-bed --out "$GENCOVE_WDIR/gencove_common"
"$PLINK19" --bfile "$DISCOVERY_WDIR/disc_ext_std"   --extract "$WDIR/common_snps.txt" --make-bed --out "$DISCOVERY_WDIR/disc_common"

# ------------------------------------------------------------
# MERGE + CALCULATE PRS
# ------------------------------------------------------------
log "=== MERGE + CALCULATE PRS ==="
BFILES=(
  "$REP_WDIR/rep_common"
  "$GENCOVE_WDIR/gencove_common"
  "$PMBB_WDIR/pmbb_common"
  "$DISCOVERY_WDIR/disc_common"
)

first="${BFILES[0]}"
merge_list="$POOLED_DIR/merge_list.txt"
> "$merge_list"
for f in "${BFILES[@]:1}"; do
  echo "${f}.bed ${f}.bim ${f}.fam" >> "$merge_list"
done

"$PLINK19" --bfile "$first" --merge-list "$merge_list" --make-bed --allow-no-sex --out "$OUT"

"$PLINK19" --bfile "$OUT" --freq --out "$OUT"

"$PLINK19" --bfile "$OUT" \
      --score "$SCORE_FILE" 1 2 4 header sum \
      --out "$OUT"

FINAL_OUTPUT="$OUTPUT_DIR/${PROJECT_NAME}.txt"
cp "${OUT}.profile" "$FINAL_OUTPUT"

log "=== Done. PRS output: ${OUT}.profile (copied to $FINAL_OUTPUT) ==="
