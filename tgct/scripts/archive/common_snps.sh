#!/bin/bash
set -euo pipefail

module load plink/1.9-20210416

# ------------------------------------------------------------
# standardize variant IDs to chr:pos before intersecting, since
# native IDs aren't consistent across the 4 datasets (rsID vs
# chr:pos vs whatever else) -- rebuilt from .bim col1 (chr) +
# col4 (bp position), with chr 23 remapped to X
# ------------------------------------------------------------
standardize_ids() {
    local prefix="$1"
    local out="$2"
    # idmap columns: 1=old ID (existing .bim col2), 2=new ID (chr:pos)
    awk '{chr = ($1 == 23) ? "X" : $1; print $2, chr":"$4}' "${prefix}.bim" > "${out}.idmap"
    # explicit column numbers (new=2, old=1) so this doesn't depend on
    # remembering plink's default --update-name column order
    plink --bfile "$prefix" --update-name "${out}.idmap" 2 1 --make-bed --out "$out"
}

standardize_ids ../pmbb4/pmbb4_ext         ../pmbb4/pmbb4_ext_std
standardize_ids ../replication/rep_ext     ../replication/rep_ext_std
standardize_ids ../discovery/discovery_prs ../discovery/discovery_prs_std
standardize_ids ../gencove/gencove_ext     ../gencove/gencove_ext_std

awk '{print $2}' ../pmbb4/pmbb4_ext_std.bim        | sort -u > pmbb4_snps.sorted
awk '{print $2}' ../replication/rep_ext_std.bim  | sort -u > replication_snps.sorted
# this one is the same
awk '{print $2}' ../discovery/discovery_prs_std.bim       | sort -u > discovery_snps.sorted
awk '{print $2}' ../gencove/gencove_ext_std.bim  | sort -u > gencove_snps.sorted

# 2. keep only SNPs present in all 4 (each list already deduped, so count==4 means "in every file")
cat pmbb4_snps.sorted replication_snps.sorted discovery_snps.sorted gencove_snps.sorted \
  | sort | uniq -c | awk '$1==4{print $2}' > common_snps.txt

wc -l common_snps.txt   # sanity check, should be ~60

# 3. subset each dataset to the common SNP list
plink --bfile ../pmbb4/pmbb4_ext_std       --extract common_snps.txt --make-bed --out ../pmbb4/pmbb4_common
plink --bfile ../replication/rep_ext_std --extract common_snps.txt --make-bed --out ../replication/rep_common
plink --bfile ../discovery/discovery_prs_std      --extract common_snps.txt --make-bed --out ../discovery/discovery_common
plink --bfile ../gencove/gencove_ext_std --extract common_snps.txt --make-bed --out ../gencove/gencove_common

