# ========================
# Non Cancer selection
# includes benign samples
# ========================
library(here)
source(here("simplexo", "scripts", "cancer_case_selection.R"))

non_cancer <- select_cancer_cases(
    cancer_name = "cancer_free",
    icd_codes = NULL,
    gender_filter = "Male",
    crep_filter = FALSE,
    min_instances = NULL,
    min_timespan = NULL,
    age_filter = NULL,
    exclude = TRUE,
    exclude_codes = c("^C(?!44)", "^D", "^Z85(?!\\.828)", "^Z86.0", "^(?!173)(?:1[4-9][0-9]|20[0-8]|209\\.[0-3])", "^V10(?!\\.83)"),
    min_exclude_instances = 0,
    min_exclude_timespan = NULL,
    output_prefix = NULL,
    pmbb_dir = here("PMBB"),
    data_dir = here("simplexo", "data"),
    log_dir = here("simplexo", "log")
)

non_cancer_no_benign <- select_cancer_cases(
    cancer_name = "cancer_free_benign",
    icd_codes = NULL,
    gender_filter = "Male",
    crep_filter = FALSE,
    min_instances = NULL,
    min_timespan = NULL,
    age_filter = NULL,
    exclude = TRUE,
    exclude_codes = c("^C(?!44)", "^D", "^Z85(?!\\.828)", "^Z86.0", "^(?!173)(?:1[4-9][0-9]|20[0-8]|209\\.[0-3])", "^V10(?!\\.83)"),
    min_exclude_instances = 0,
    min_exclude_timespan = NULL,
    output_prefix = NULL,
    pmbb_dir = here("PMBB"),
    data_dir = here("simplexo", "data"),
    log_dir = here("simplexo", "log")
)

