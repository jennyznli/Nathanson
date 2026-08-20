#!/bin/bash

module load plink/1.9-20210416

# 4. recompute the PRS on each common-SNP subset
plink --bfile pmbb4/pmbb4_prs_common       --score tgct_score_file1.txt 1 2 7 header sum --out pmbb4/pmbb4_prs_adj
plink --bfile replication/rep_prs_common --score tgct_score_file1.txt 1 2 7 header sum --out replication/rep_prs_adj
plink --bfile discovery/penn_prs_common  --score tgct_score_file1.txt 1 2 7 header sum --out discovery/penn_prs_adj
plink --bfile gencove/gencove_prs_common --score tgct_score_file1.txt 1 2 7 header sum --out gencove/gencove_prs_adj
