# ========================
# Packages
# ========================
packages <- c("tidyr", "dplyr", "plotly", "readr", "readxl", "here",
              "stringr", "ggplot2",  "impute", "pals", "data.table", "lubridate"
)
purrr::walk(packages, ~ require(.x, character.only = TRUE))
here()


DATE <- format(Sys.Date(), "%Y%m%d")
TIMESTAMP <- format(Sys.time(), "%Y%m%d_%H%M%S")

setwd(here("tgct"))

# ========================
# LOAD DATA
# ========================
OCC <- here("PMBB", "3.0", "PMBB-Release-2024-3.0_phenotype_condition_occurrence.txt")
COV <- here("PMBB", "3.0", "PMBB-Release-2024-3.0_covariates.txt")
PER <- here("PMBB", "3.0", "PMBB-Release-2024-3.0_phenotype_person.txt")

cov <- fread(COV, header = TRUE)
person <- fread(PER, header = TRUE)
flags <- fread(here("PMBB", "3.0", "rgcname_pmbbid_metadata_flags.csv"), header = TRUE)
icds <- readLines((here("tgct", "data", "tgct_cancer_filtered_patients_ids.txt")))

# ========================
# MATCH PMBB TO KT IDS
# ========================
kt_pattern <- "^UPENN-PMBB_KT[0-9]+_KT[0-9]+$"

kt <- flags %>%
    mutate(match_col = apply(select(., RGC_sample_name, ID1, ID2), 1, function(row) {
        m <- row[grepl(kt_pattern, row)]
        if (length(m) > 0) m[1] else NA
    })) %>%
    filter(!is.na(match_col))
print(dim(kt))
# 900

kt <- select(kt, c("PMBB_ID",  "RGC_sample_name", "ID1", "ID2", "CREP", "match_col"))
kt$ID <- sub(".*(KT[0-9]{7}).*", "\\1", kt$match_col)

kt <- select(kt, c("PMBB_ID",  "match_col", "ID"))
kt$SampNum <- numbers <- as.numeric(gsub("KT", "", kt$ID))

write.csv(kt, here("PMBB", "3.0", "pmbb3_kt_key.csv"))

# ========================
# PROGENY - CASE PULL
# ========================
length(icds) # 779
# length(unique(icds))

intersect(icds, kt)

kt_not_icd <- setdiff(kt$PMBB_ID, icds)
length(kt_not_icd)
# 226 kt ids that are not detected by ICD, need to be checked

kt_icd <- intersect(kt$PMBB_ID, icds)
length(kt_icd)
# 674 kt ids that are identified by icd codes

icd_not_kt <- setdiff(icds, kt$PMBB_ID)
length(icd_not_kt)
# 105 icd codes not KT ides




