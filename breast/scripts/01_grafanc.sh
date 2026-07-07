#!/bin/bash

### BREAST CANCER ANCESTRY ANALYSIS ###

# Set up 
WDIR="/project/knathans_tecac/jenny/breast/analysis/grafanc"
INPUT_FILE="/project/knathans_tecac/jenny/grafanc/results/ancestry_final.txt"
ID_DIR="/project/knathans_tecac/jenny/breast/analysis/ids"
mkdir -p "$ID_DIR"
cd "$WDIR"

CASES_FILE="/project/knathans_tecac/jenny/phenotype/breast_cancer_filtered_patients_ids.txt"
CONTROLS_FILE="/project/knathans_tecac/jenny/phenotype/non_cancer_female_nocrep_ids.txt"

### INTERSECT GRAFANC, UNRELATED, CASE/CONTROLS ###
tail -n +2 "$INPUT_FILE" | cut -f1 | sort -u > "grafanc_ids.txt"
comm -12 "grafanc_ids.txt" $CASES_FILE | sort | comm -12 - <(sort PMBB-Release-2024-3.0_genetic_exome.3rd_degree_unrelated.txt) > case_ids.txt 
comm -12 "grafanc_ids.txt" $CONTROLS_FILE | sort | comm -12 - <(sort PMBB-Release-2024-3.0_genetic_exome.3rd_degree_unrelated.txt) > control_ids.txt

cat case_ids.txt control_ids.txt | sort -u > case_control_ids.txt
awk '{print "0\t" $1}' case_ids.txt > case_pids.txt
awk '{print "0\t" $1}' control_ids.txt > control_pids.txt
awk '{print "0\t" $1}' case_control_ids.txt > case_control_pids.txt

cp case_ids.txt "$ID_DIR/case_ids.txt"
cp control_ids.txt "$ID_DIR/control_ids.txt"
cp case_pids.txt "$ID_DIR/case_pids.txt"
cp control_pids.txt "$ID_DIR/control_pids.txt"
cp case_control_pids.txt "$ID_DIR/case_control_pids.txt"
cp case_control_ids.txt "$ID_DIR/case_control_ids.txt"

### ADD CASE/CONTROL TO ANCESTRY FILE ###
awk '{print $1 "\t" "Case"}' case_ids.txt > sample_status.tmp
awk '{print $1 "\t" "Control"}' control_ids.txt >> sample_status.tmp
awk -F'\t' 'BEGIN {OFS="\t"}
    NR==FNR {
        status[$1] = $2
        next
    }
    NR>FNR {
        if (FNR==1) {
            print $0, "Status"
        } else {
            if ($1 in status) {
                print $0, status[$1]
            }
        }
    }
' sample_status.tmp "$INPUT_FILE" > ancestry_breast.txt
rm sample_status.tmp

breast_samples=$(tail -n +2 ancestry_breast.txt | wc -l)
breast_cases=$(awk -F'\t' '$NF=="Case"' ancestry_breast.txt | wc -l)
breast_controls=$(awk -F'\t' '$NF=="Control"' ancestry_breast.txt | wc -l)

echo "Ancestry data with case/control status:"
echo "  Total samples: $breast_samples"
echo "  Cases: $breast_cases"
echo "  Controls: $breast_controls"
echo ""

### CREATE ANCESTRY IDS ### 
unique_ancestries=$(cut -f3 ancestry_breast.txt | tail -n +2 | sort -u)
echo "Found ancestries:"
echo "$unique_ancestries"
echo ""

for ancestry in $unique_ancestries; do
    ancestry_ids="$ID_DIR/${ancestry}_ids.txt"
    ancestry_pids="$ID_DIR/${ancestry}_pids.txt"
    ancestry_cases="$ID_DIR/${ancestry}_case_ids.txt"
    ancestry_controls="$ID_DIR/${ancestry}_control_ids.txt"
    ancestry_case_pids="$ID_DIR/${ancestry}_case_pids.txt"
    ancestry_control_pids="$ID_DIR/${ancestry}_control_pids.txt"
    
    # All samples
    awk -F'\t' -v anc="$ancestry" 'NR>1 && $3==anc {print $1}' ancestry_breast.txt > "$ancestry_ids"
    awk -F'\t' -v anc="$ancestry" 'NR>1 && $3==anc {print "0\t" $1}' ancestry_breast.txt > "$ancestry_pids"
    
    # Cases only
    awk -F'\t' -v anc="$ancestry" 'NR>1 && $3==anc && $NF=="Case" {print $1}' ancestry_breast.txt > "$ancestry_cases"
    awk -F'\t' -v anc="$ancestry" 'NR>1 && $3==anc && $NF=="Case" {print "0\t" $1}' ancestry_breast.txt > "$ancestry_case_pids"
    
    # Controls only
    awk -F'\t' -v anc="$ancestry" 'NR>1 && $3==anc && $NF=="Control" {print $1}' ancestry_breast.txt > "$ancestry_controls"
    awk -F'\t' -v anc="$ancestry" 'NR>1 && $3==anc && $NF=="Control" {print "0\t" $1}' ancestry_breast.txt > "$ancestry_control_pids"
done
echo ""

SUMMARY_FILE="breast_ancestry_summary.txt"

cat << EOF > "$SUMMARY_FILE"
=== BREAST CANCER ANCESTRY ANALYSIS SUMMARY ===
Total samples with ancestry and case/control data: $breast_samples
Cases: $breast_cases
Controls: $breast_controls

=== CASE-CONTROL RATIOS BY CONTINENTAL ANCESTRY ===
Ancestry          Cases  Controls  Total  Case%
EOF

for ancestry in $unique_ancestries; do
    cases_count=$(awk -F'\t' -v anc="$ancestry" '$3==anc && $NF=="Case"' ancestry_breast.txt | wc -l)
    controls_count=$(awk -F'\t' -v anc="$ancestry" '$3==anc && $NF=="Control"' ancestry_breast.txt | wc -l)
    total_count=$((cases_count + controls_count))
    
    if [[ $total_count -gt 0 ]]; then
        case_percent=$(awk "BEGIN {printf \"%.1f\", $cases_count * 100 / $total_count}")
        printf "%-16s  %5d  %8d  %5d  %5s%%\n" "$ancestry" "$cases_count" "$controls_count" "$total_count" "$case_percent" >> "$SUMMARY_FILE"
    fi
done

echo "=== ANALYSIS COMPLETE ==="
echo ""
echo "Output files created: $WDIR"
echo "  Main data:"
echo "    - ancestry_breast.txt (ancestry data with case/control status)"
echo "    - $SUMMARY_FILE (analysis summary)"
echo ""
echo "  ID files: $ID_DIR"
echo "    - case_ids.txt, control_ids.txt, case_control_ids.txt"
echo "    - case_pids.txt, control_pids.txt, case_control_pids.txt"
echo ""
