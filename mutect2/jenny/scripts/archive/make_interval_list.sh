#!/usr/bin/env bash
# One-time setup: convert the merged targets BED to a Picard .interval_list
# (needs a sequence dictionary header, unlike plain BED).
#
# Usage: bash make_interval_list.sh
set -euo pipefail

DICT="/home/jennyzli/resources/Genomes/Human/hg38/fa/Homo_sapiens_assembly38.dict"
TARGETS_BED="/home/jennyzli/jenny/gray/bed/IDT-Target-HYB-PreCancer.hg38.targets.merged.bed"
OUT="/home/jennyzli/jenny/gray/bed/IDT-Target-HYB-PreCancer.hg38.targets.merged.interval_list"

gatk BedToIntervalList \
    -I "$TARGETS_BED" \
    -O "$OUT" \
    -SD "$DICT"

echo "Wrote $OUT"
