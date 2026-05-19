# ========================
# TGCT Cancer Selection
# ========================
library(here)
source(here("tgct", "scripts", "cancer_case_selection.R"))

tgct_results <- select_cancer_cases(
    cancer_name = "tgct",
    icd_codes = c("^C62", "^Z80.43", "^186", "^V10.47"),
    gender_filter = "Male",
    crep_filter = NULL,
    min_instances = 2,
    min_timespan = NULL,
    exclude = FALSE,
    age_filter = NULL,
    output_prefix = "tgct",
    pmbb_dir = here("PMBB"),
    data_dir = here("tgct", "data"),
    log_dir = here("tgct", "log")
)
