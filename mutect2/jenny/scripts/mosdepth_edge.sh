#!/bin/bash
# Usage: bash run_mosdepth_edge.sh <sample_id> <bam_or_cram>
set -euo pipefail

SAMPLE="$1"
BAM="$2"
OUTDIR="./mosdepth"   # all examination output lands here, relative to wherever this is run from

TARGETS_MERGED="/home/jennyzli/jenny/gray/bed/IDT-Target-HYB-PreCancer.hg38.targets.merged.bed"
GENOME="/home/jennyzli/jenny/gray/bed/hg38.genome"          # chrom<TAB>size, from process_targets_bed.sh step 5
WINDOW_BP=100                 # how far past each edge to look; widen if depth hasn't plateaued

mkdir -p "$OUTDIR"

# 1. Widen targets by WINDOW_BP just to define the extraction window (NOT the final padding).
bedtools slop -b "$WINDOW_BP" -g "$GENOME" -i "$TARGETS_MERGED" \
    | bedtools sort -i - -g "$GENOME" \
    | bedtools merge -i - > "$OUTDIR/windows.bed"

# 2. Per-base depth (mosdepth's default output includes {prefix}.per-base.bed.gz).
mosdepth --threads 4 "$OUTDIR/${SAMPLE}" "$BAM"

# 3. Keep only the bases actually near a target edge (keeps the file small for R).
zcat "$OUTDIR/${SAMPLE}.per-base.bed.gz" \
    | bedtools intersect -a - -b "$OUTDIR/windows.bed" \
    > "$OUTDIR/${SAMPLE}.edge_window.bed"

echo "Done: $OUTDIR/${SAMPLE}.edge_window.bed"
echo "Also copy/symlink $TARGETS_MERGED into $OUTDIR so the R script can find it."
