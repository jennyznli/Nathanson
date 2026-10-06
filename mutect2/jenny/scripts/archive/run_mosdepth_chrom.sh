#!/bin/bash
# Windowed depth across a single chromosome (fast -- mosdepth's own -c/--chrom
# restricts the scan, --by bins it, so no per-base expansion needed).
# Used by plot_chrom_coverage.R for the whole-chromosome overview.
#
# Usage: bash run_mosdepth_chrom.sh <sample_id> <bam_or_cram> <chrom> [window_bp]
set -euo pipefail

SAMPLE="$1"
BAM="$2"
CHROM="$3"
WINDOW_BP="${4:-5000}"
OUTDIR="./mosdepth"   # same dir run_mosdepth_edge.sh writes to

mkdir -p "$OUTDIR"

# Double-check -c/--chrom and --by against `mosdepth --help` on your install --
# not verified against your specific mosdepth version here.
mosdepth --threads 4 -c "$CHROM" --by "$WINDOW_BP" "$OUTDIR/${SAMPLE}.${CHROM}" "$BAM"

echo "Done: $OUTDIR/${SAMPLE}.${CHROM}.regions.bed.gz"
