#!/usr/bin/env bash
# Usage: concat_vep_reports.sh <src_dir> <out_file> [glob_pattern]
set -euo pipefail

if [ $# -lt 2 ]; then
    echo "Usage: $0 <src_dir> <out_file> [glob_pattern]" >&2
    exit 1
fi

SRC_DIR="$1"
OUT_FILE="$2"
PATTERN="${3:-*.csv}"

shopt -s nullglob
files=("${SRC_DIR}"/${PATTERN})
shopt -u nullglob

if [ ${#files[@]} -eq 0 ]; then
    echo "No files matching ${PATTERN} found in ${SRC_DIR}" >&2
    exit 1
fi

ref_header=$(head -n1 "${files[0]}")
mismatched=()
for f in "${files[@]}"; do
    h=$(head -n1 "$f")
    if [ "$h" != "$ref_header" ]; then
        mismatched+=("$f")
    fi
done

if [ ${#mismatched[@]} -gt 0 ]; then
    echo "Header mismatch -- refusing to concatenate. Files with a different header than" >&2
    echo "${files[0]}:" >&2
    printf '  %s\n' "${mismatched[@]}" >&2
    exit 1
fi

mkdir -p "$(dirname "${OUT_FILE}")"
awk 'FNR==1 && NR!=1{next}{print}' "${files[@]}" > "${OUT_FILE}"

echo "Merged ${#files[@]} files -> ${OUT_FILE} ($(wc -l < "${OUT_FILE}") lines incl. header)"
