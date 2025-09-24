# ========================
# Breast Cancer Selection
# ========================
library(here)
source(here("simplexo", "scripts", "cancer_case_selection.R"))

breast_results <- select_cancer_cases(
    cancer_name = "breast_simplexo3_v2",
    icd_codes = c("^C50", "^Z85.4", "^174", "^V10.3", "^D05", "^Z86.000", "^233.0"),
    gender_filter = "Female",
    crep_filter = NULL,
    min_instances = 2,
    min_timespan = NULL,
    age_filter = NULL,
    exclude = FALSE,
    exclude_codes = c("^C(?!44)", "^Z85(?!\\.828)", "^(?!173)(?:1[4-9][0-9]|20[0-8]|209\\.[0-3])", "^V10(?!\\.83)"),
    min_exclude_instances = 0,
    min_exclude_timespan = NULL,
    output_prefix = NULL,
    pmbb_dir = here("PMBB"),
    data_dir = here("simplexo", "data"),
    log_dir = here("simplexo", "log")
)
