#!/usr/bin/env bash
# Usage: bash collect_hs_metrics.sh <sample_id> <bam_or_cram>
set -euo pipefail

SAMPLE="$1"
BAM="$2"
OUTDIR="./hs_metrics"

REF="/home/jennyzli/resources/Genomes/Human/hg38/fa/Homo_sapiens_assembly38.fasta"
TARGET_INTERVALS="/home/jennyzli/jenny/gray/bed/IDT-Target-HYB-PreCancer.hg38.targets.merged.interval_list"
BAIT_INTERVALS="/home/jennyzli/jenny/gray/bed/IDT-Target-HYB-PreCancer.hg38.probes.merged.interval_list"

mkdir -p "$OUTDIR"

gatk CollectHsMetrics \
    -I "$BAM" \
    -O "$OUTDIR/${SAMPLE}.hs_metrics.txt" \
    -R "$REF" \
    --BAIT_INTERVALS "$BAIT_INTERVALS" \
    --TARGET_INTERVALS "$TARGET_INTERVALS" \
    --PER_TARGET_COVERAGE "$OUTDIR/${SAMPLE}.per_target_coverage.txt"

echo "Done: $OUTDIR/${SAMPLE}.hs_metrics.txt"
