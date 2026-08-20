
# extract all relevant SNPs 
bsub -o extract_all.o -e extract_all.e bash extract_all2.sh

# extract all relevant SNPs from PMBB 
bash find_pmbb_chunks.sh
bash extract_pmbb_chunks.bsub
bash merge_pmbb_chunks.bsub

# standardize the IDS and  overlap these snps to have common list 
bash common_snps.sh

# calculate prs 
bash merge_pooled_prs.bsub

