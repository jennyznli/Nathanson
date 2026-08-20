#!/bin/bash

module load plink/1.9-20210416

plink --bfile pmbb4/pmbb4_prs_common	   --score gwas2021_rep_scores_b38.txt 1 2 7 header sum --out pmbb4/pmbb4_prs
plink --bfile replication/rep_prs_common --score gwas2021_rep_scores_b38.txt 1 2 7 header sum --out replication/rep_prs
plink --bfile discovery/penn_prs_common  --score gwas2021_rep_scores_b38.txt 1 2 7 header sum --out discovery/penn_prs
plink --bfile gencove/gencove_prs_common --score gwas2021_rep_scores_b38.txt 1 2 7 header sum --out gencove/gencove_prs
