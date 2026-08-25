# ========================
# RELAPSE
# ========================
library(here)
source(here("R", "config.R"))
setwd(here("tgct"))

# ========================
# LOAD DATA
# ========================
clin <- read_excel(here("tgct", "ss", "clinical data for cases for John.xls"), col_names = TRUE)
ss0 <- read_tsv(here("tgct", "ss", "pennsamples_masterslist_062226.txt"))
ss <- read_excel(here("tgct", "ss", "20260710_tgct_master.xlsx"))
ids <- read.table(here("tgct", "data", "20260710_clinical_ids.txt"))

cov <- fread(here("PMBB", "4.0", "PMBB-Release-2026-4.0_phenotype_covariates.txt"), header = TRUE)
flags <- fread(here("PMBB", "4.0", "rgcname_pmbbid_metadata_flags_freeze4.csv"), header = TRUE)

# ========================
# FIND SNPS in AR REGION
# ========================
ids <- ids %>% rename(Study = V1, Study_ID = V2)
# 9:124561497 and X:67140772

# extract patients that match clinical data
# extract SNP of interest


