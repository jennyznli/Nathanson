
# cut -d'_' -f1,2 snp_list.txt > snp_list_chrpos.txt

module load R/4.5

Rscript makeScoreFile.R -a gwas2021_rep_scores_b38.txt -f pooled/tecac_prs_pooled.frq -o tgct_score_file1.txt
# need to make a case/control file or prep the file for comparequantiltomedia 
# can try for just replication 
awk 'NR==FNR {cases[$2]=1; next}
     {
       s=0
       if ($2 in cases) s=1
       print $1, $2, s
     }' tgct_prs_ids.txt pooled/tecac_prs_pooled.fam \
> tgct_case_control_combined.txt

compareQuantileToMedian.R 
