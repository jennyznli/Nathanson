library(here)
setwd("/Users/jennyzli/Documents/Nathanson")
source(here("R", "config.R"))

# ========================
# LOAD DATA
# ========================
ss <- read_excel(file.path("ch", "ss", "20251110_aml_ngs.xlsx"))

up <- read.csv(here("simplexo", "data", "simplexo_up_map.csv"))
flags <- fread(here("PMBB", "3.0", "rgcname_pmbbid_metadata_flags.csv"))

up_pattern <- "^UPENN-PMBB_UP[0-9]+_UP[0-9]+$"
key <- flags %>% select(PMBB_ID, RGC_sample_name)
