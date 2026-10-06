#!/bin/bash

# Run mosdepth general coverage QC on all samples

set -euo pipefail

# -----------------------------
# Configuration
# -----------------------------

CASES="input/cases.tsv"
TARGETS="/home/jennyzli/jenny/gray/bed/IDT-Target-HYB-PreCancer.hg38.targets.merged.bed"
OUTDIR="output/mosdepth"
THREADS=4

mkdir -p "$OUTDIR"


# -----------------------------
# Run mosdepth
# -----------------------------

while IFS=$'\t' read -r SAMPLE BAM; do

    # Skip blank lines
    [[ -z "$SAMPLE" ]] && continue

    echo "========================================"
    echo "Sample: $SAMPLE"
    echo "BAM:    $BAM"
    echo "========================================"

    mosdepth \
        --threads "$THREADS" \
        --by "$TARGETS" \
        --thresholds 50,100,250,500,1000 \
        "$OUTDIR/${SAMPLE}" \
        "$BAM"

    echo "Finished: $SAMPLE"
    echo

done < "$CASES"

echo "========================================"
echo "All samples finished."
echo "Results: $OUTDIR"
echo "========================================"