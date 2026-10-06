#!/bin/bash
# This script collates alignment metrics from FastQC and Picard using MultiQC.
# 
#BSUB -J multiqc_alignment_metrics
#BSUB -o logs/multiqc_alignment_metrics.%J.out
#BSUB -e logs/multiqc_alignment_metrics.%J.err
#BSUB -n 1
#BSUB -R "rusage[mem=1000]"
#BSUB -M 1000

set -euo pipefail

# Assign variables
FASTQC_DIRECTORY=/home/jennyzli/jenny/gray/fastq/fastqc
METRICS_DIRECTORY=/home/jennyzli/jenny/gray/metrics/GRCh38-broad
MARKDUP_DIRECTORY=/home/jennyzli/jenny/gray/bam_input/work
MOSDEPTH_DIRECTORY=/home/jennyzli/jenny/gray/mosdepth
OUTPUT_DIRECTORY=/home/jennyzli/jenny/gray/multiqc

# Create directory for output
mkdir -p "$OUTPUT_DIRECTORY"

multiqc \
  "$FASTQC_DIRECTORY" \
  "$METRICS_DIRECTORY" \
  "$MARKDUP_DIRECTORY" \
  "$MOSDEPTH_DIRECTORY" \
  --dirs \
  --outdir "$OUTPUT_DIRECTORY"