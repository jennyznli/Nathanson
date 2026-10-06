#!/usr/bin/env bash
# Submits one lightweight bsub job per sample in input/cases.tsv.
# Run directly on the login node (this script itself is not a bsub job) --
# CollectHsMetrics is per-sample and cheap, unlike edge_effect_plot.R, so
# there's no combine step needed.
#
# Usage: bash submit_hs_metrics.sh
set -euo pipefail

mkdir -p logs hs_metrics

while IFS=$'\t' read -r sample bam; do
    bsub -J "hs_${sample}" \
         -o "logs/hs_metrics.${sample}.%J.out" \
         -e "logs/hs_metrics.${sample}.%J.err" \
         -n 1 -R "rusage[mem=4000]" -M 4000 -W 00:30 \
         bash scripts/collect_hs_metrics.sh "$sample" "$bam"
done < input/cases.tsv
