#!/bin/bash
# External source genotype data -- fixed reference locations, don't change per project.
# Sourced by prs.sh from $PRS_ROOT/config/sources.sh. Not standalone: relies on
# prs.sh having already done `set -euo pipefail`.

DISCOVERY_DIR="/project/knathans_tecac/GWAS2024/pennwtcga/imputed"
GENCOVE_BFILE="/project/knathans_tecac/gencove/all/gencove_qc6.lifted.hg38.final"
REP_DIR="/project/knathans_tecac/replication_vcfs"
PMBB_DIR="/static/PMBB/PMBB-Release-2026-4.0/Imputed"
PMBB_MANIFEST="${PMBB_DIR}/metadata/imputed_variant_chunked_input_manifest.tsv"

