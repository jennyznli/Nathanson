#!/bin/bash
set -euo pipefail

# inputs
PREFIX="IDT-Target-HYB-PreCancer.hg38"
RAW_TARGETS="/home/jennyzli/jenny/gray/input/Targets-XGEN.FB037FC25D2F4BA7A4BBF676810B2ED6.bed"
RAW_PROBES="/home/jennyzli/jenny/gray/input/Probes-XGEN.FB037FC25D2F4BA7A4BBF676810B2ED6.bed"
REF_FASTA="/home/jennyzli/resources/Genomes/Human/hg38/fa/Homo_sapiens_assembly38.fasta"
REF_FAI="${REF_FASTA}.fai"
REF_DICT="${REF_FASTA%.fasta}.dict"
SLOP_BP=50
OUTDIR="/home/jennyzli/jenny/gray/bed"

# ---------------------

mkdir -p "$OUTDIR"

LOGBASE="$OUTDIR/${PREFIX}.process_targets_bed.$(date +%Y%m%d_%H%M%S)"
LOG_OUT="${LOGBASE}.o"
LOG_ERR="${LOGBASE}.e"
exec > >(tee -a "$LOG_OUT") 2> >(tee -a "$LOG_ERR" >&2)
echo "Logging stdout -> $LOG_OUT"
echo "Logging stderr -> $LOG_ERR" >&2

GENOME="$OUTDIR/hg38.genome"
awk '{print $1"\t"$2}' "$REF_FAI" > "$GENOME"

report () {
  # $1 = label, $2 = bed file
  local n_rows total_bp
  n_rows=$(grep -vc '^#' "$2" || true)
  total_bp=$(awk '!/^#/{sum+=$3-$2} END{print sum+0}' "$2")
  printf "%-28s rows=%-8s total_bp=%s\n" "$1" "$n_rows" "$total_bp"
}

check_contigs () {
  # $1 = bed file; warn about contigs not in the reference
  local missing
  missing=$(comm -23 <(cut -f1 "$1" | sort -u) <(cut -f1 "$REF_FAI" | sort -u))
  if [[ -n "$missing" ]]; then
    echo "WARNING: contigs in $(basename "$1") but not in reference:"
    echo "$missing"
  else
    echo "OK: all contigs in $(basename "$1") found in reference"
  fi
}

process_bed () {
  # $1 = label (targets/probes), $2 = raw BED, $3 = merge args
  # sets MERGED to the unpadded merged BED
  local label="$1" raw="$2" merge_args="$3"
  local chr="$OUTDIR/${PREFIX}.${label}.chr.bed"
  local sorted="$OUTDIR/${PREFIX}.${label}.chr.sorted.bed"
  MERGED="$OUTDIR/${PREFIX}.${label}.merged.bed"
  local ilist="$OUTDIR/${PREFIX}.${label}.merged.interval_list"

  echo
  echo "##### ${label} #####"
  echo "--- input contig naming ---"
  grep -v '^#\|^track\|^browser' "$raw" | cut -f1 | sort -u | head -5 || true
  report "${label} raw" "$raw"

  echo "--- 1. chr-prefix (drops #/track/browser header lines) ---"
  awk 'BEGIN{OFS="\t"} /^#|^track|^browser/{next} $1 !~ /^chr/{$1="chr"$1} {print}' "$raw" > "$chr"
  report "${label} chr-prefixed" "$chr"

  echo "--- 2. sort ---"
  bedtools sort -i "$chr" -g "$REF_FAI" > "$sorted"

  echo "--- 3. merge overlapping/adjacent intervals ---"
  # shellcheck disable=SC2086
  bedtools merge -i "$sorted" $merge_args > "$MERGED"
  report "${label} merged" "$MERGED"
  check_contigs "$MERGED"

  echo "--- 4. interval_list (BED + sequence dictionary header, for Picard) ---"
  gatk BedToIntervalList -I "$MERGED" -O "$ilist" -SD "$REF_DICT"
}

echo "Reference contig naming:"
grep -m3 '^>' "$REF_FASTA"

# Targets keep their names (col 4, comma-joined when merged); probes are merged to plain
# 3-col intervals since thousands of probe names per region aren't useful
process_bed targets "$RAW_TARGETS" "-c 4 -o distinct"
TARGETS_MERGED="$MERGED"
process_bed probes "$RAW_PROBES" ""

echo
echo "##### targets: slop +-${SLOP_BP}bp, then re-sort + re-merge (padding can create new overlaps) #####"
TARGETS_SLOP="$OUTDIR/${PREFIX}.targets.slop${SLOP_BP}.bed"
TARGETS_PADDED="$OUTDIR/${PREFIX}.targets.slop${SLOP_BP}.merged.bed"
bedtools slop -b "$SLOP_BP" -g "$GENOME" -i "$TARGETS_MERGED" > "$TARGETS_SLOP"
report "targets slop${SLOP_BP}" "$TARGETS_SLOP"
bedtools sort -i "$TARGETS_SLOP" -g "$REF_FAI" | bedtools merge -i - -c 4 -o distinct > "$TARGETS_PADDED"
report "targets slop${SLOP_BP} merged" "$TARGETS_PADDED"
check_contigs "$TARGETS_PADDED"