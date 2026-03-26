# ========================
# Packages
# ========================
packages <- c("tidyr", "dplyr", "plotly", "readr", "readxl", "here",
              "stringr", "ggplot2",  "impute", "pals", "geneplotter"
)
purrr::walk(packages, ~ require(.x, character.only = TRUE))
here()

DATE <- format(Sys.Date(), "%Y%m%d")

# ============================================================
# PMBB CREP FLAGS
# ============================================================
# key <- read.csv(here("rgcname_pmbbid_metadata_flags.csv"))
# write_tsv(key, "crep_pmbb.tsv")
setwd("/Users/jennyzli/Documents/Nathanson/ch")

df <- read_excel(here("ch", "data", "brca_carriers_ch_freq_w_sampnum_20250716.xlsx"))
print(dim(df))
# 57170    13

pmbb_consent <- df %>% filter(PMBB_Consent == 1)
print(dim(pmbb_consent))
# 53062    13

## CREP only
crep <- df %>% filter(CREP == 1)
print(dim(crep))
# 2864   13

## CREP no PMBB
crep_no_pmbb <- crep %>% filter(PMBB_Consent == 0)
print(dim(crep_no_pmbb))
# 2193   13

## get PMBB and CREP consented for PMBB
crep_pmbb <- crep %>% filter(PMBB_Consent == 1)
print(dim(crep_pmbb))
# 671   13

## get no consent OR crep
crep_pmbb <- crep %>% filter(PMBB_Consent == 1)
print(dim(crep_pmbb))
# 671   13

# write IDs
write.table(pmbb_consent$PMBB_ID, "pmbb_consent_ids.txt", quote=FALSE, row.names = FALSE, col.names = FALSE)
write.table(crep$PMBB_ID, "crep_pmbb_ids.txt", quote=FALSE, row.names = FALSE, col.names = FALSE)
write.table(crep_no_pmbb$PMBB_ID, "crep_pmbb_noconsent_ids.txt", quote=FALSE, row.names = FALSE, col.names = FALSE)
write.table(crep_pmbb$PMBB_ID, "crep_pmbb_consent_ids.txt", quote=FALSE, row.names = FALSE, col.names = FALSE)

# ============================================================
# get matching PMBB ids?
# ============================================================
#


