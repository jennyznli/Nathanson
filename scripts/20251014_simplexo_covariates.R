# ========================
# SIMPLEXO - AGE EXTRACTION & COVARIATE CREATION
# ========================

library(here)
setwd(here("simplexo"))

source(here("R", "load_packages.R"))
source(here("R", "sample_selection.R"))
source(here("R", "control_selection.R"))

# ========================
# FILE PATHS
# ========================

PMBB_DIR <- here("PMBB", "3.0")
OCC_FILE <- here(PMBB_DIR, "PMBB-Release-2024-3.0_phenotype_condition_occurrence.txt")
COV_FILE <- here(PMBB_DIR, "PMBB-Release-2024-3.0_covariates.txt")
PER_FILE <- here(PMBB_DIR, "PMBB-Release-2024-3.0_phenotype_person.txt")
EXOME_FILE <- here(PMBB_DIR, "PMBB-Release-2024-3.0_genetic_exome.norm.commonsnps.eigenvec")
IMP_FILE <- here(PMBB_DIR, "PMBB-Release-2024-3.0_genetic_imputed.eigenvec")

# ========================
# LOAD ALL DATA FILES
# ========================

cov <- fread(COV_FILE, header = TRUE)
person <- fread(PER_FILE, header = TRUE)
exome <- fread(EXOME_FILE, header = TRUE)[, 1:11]
imp <- fread(IMP_FILE, header = TRUE)[, 1:11]

all_ids <- unique(read.table(here("simplexo", "data", "simplexo_overall_v1_case_ids.txt"))$V1)
icd_cases <- fread(here("simplexo", "data", "simplexo_cancer_filtered_patients.txt"))
icd_case_ids <- readLines(here("simplexo", "data", "simplexo_cancer_filtered_patients_ids.txt"))

pmcr <- read.csv(here("simplexo", "ss", "pmbb_147_pmcrbreastage.csv"))
progeny_pmbb <- read.csv(here("simplexo", "data", "simplexo_progeny_pmbb.txt"))

# ========================
# EXTRACT DIAGNOSIS AGES FROM ICD CODES
# ========================
# Merge with date of birth
pmbb_dob <- person %>% select(person_id, birth_datetime)
dob_dx <- merge(icd_cases, pmbb_dob, by = "person_id", all.x = TRUE)

# Calculate diagnosis age (years)
dob_dx$CaDxAge_ICD <- as.numeric(difftime(
    dob_dx$first_date,
    dob_dx$birth_datetime,
    units = "days"
)) / 365.25

icd_ages <- dob_dx %>% select(person_id, CaDxAge_ICD)

cat("ICD ages extracted:", nrow(icd_ages), "records\n") # 3448

# ========================
# EXTRACT DIAGNOSIS AGES FROM PMCR (Cancer Registry)
# ========================

pmcr_case <- pmcr %>% filter(PmcrBreastFlag == 1)

pmcr_ages <- pmcr_case %>%
    select(PMBB_ID, PmcrBreastFirstDxDate_age) %>%
    rename(CaDxAge_PMCR = PmcrBreastFirstDxDate_age)

cat("PMCR cases:", nrow(pmcr_case), "records\n") # 2510

# ========================
# EXTRACT DIAGNOSIS AGES FROM PROGENY
# ========================
# Function to extract earliest age from semicolon-separated string
get_earliest_age <- function(age_str) {
    if (is.na(age_str) || is.null(age_str) || age_str == "" || age_str == "NA") {
        return(NA)
    }

    age_str <- as.character(age_str)
    age_parts <- strsplit(age_str, ";")[[1]]
    ages <- suppressWarnings(as.numeric(age_parts))
    valid_ages <- ages[!is.na(ages)]

    if (length(valid_ages) == 0) {
        warning(paste("No valid ages found in:", age_str))
        return(NA)
    }

    if (any(valid_ages < 0)) {
        warning(paste("Negative age found in:", age_str))
    }

    return(min(valid_ages))
}

progeny_pmbb$CaDxAge_Progeny <- sapply(progeny_pmbb$CaDxAge, get_earliest_age)
progeny_pmbb$CaDxAge_Progeny[is.infinite(progeny_pmbb$CaDxAge_Progeny)] <- NA

progeny_ages <- progeny_pmbb %>% select(PMBB_ID, CaDxAge_Progeny)

cat("Progeny records:", nrow(progeny_pmbb), "records\n") # 1520

# ========================
# MERGE ALL DIAGNOSIS AGES
# ========================

base_df <- data.frame(PMBB_ID = all_ids)

merged <- base_df %>%
    left_join(progeny_ages, by = "PMBB_ID") %>%
    left_join(pmcr_ages, by = "PMBB_ID") %>%
    left_join(icd_ages, by = c("PMBB_ID" = "person_id")) %>%
    mutate(
        # Priority order: Progeny > PMCR > ICD
        Age = coalesce(CaDxAge_Progeny, CaDxAge_PMCR, CaDxAge_ICD),
        # Find earliest diagnosis across all sources
        EarliestDiagnosis = pmin(CaDxAge_Progeny, CaDxAge_PMCR, CaDxAge_ICD, na.rm = TRUE),
        # Flag cases where selected age differs from earliest by >2 years
        AgeDifference = abs(Age - EarliestDiagnosis) > 2
    )
write.table(merged,
            here("simplexo", "data", "simplexo_overall_case_ages_merged.txt"),
            row.names = FALSE, quote = FALSE, sep = "\t")


# ========================
# SELECT FINAL CASE IDs WITH AGES
# ========================

merged_sel <- merged %>%
    select(PMBB_ID, Age) %>%
    filter(!is.na(Age)) %>%
    mutate(Age = round(Age, 2))

cat("Final cases with age:", nrow(merged_sel), "records\n") # 4059
cat("Cases without age:", nrow(base_df) - nrow(merged_sel), "records\n") # 4

write.table(sort(merged_sel$PMBB_ID),
            here("simplexo", "data", "simplexo_overall_case_ids.txt"),
            row.names = FALSE, col.names = FALSE, quote = FALSE)

# ========================
# CREATE FINAL COVARIATE FILES
# ========================

# Start with sample ages for all PMBB participants
pmbb_ages <- cov %>%
    select(person_id, Sample_age) %>%
    rename(PMBB_ID = person_id, Age = Sample_age)

# Prepare case ages
case_ages <- merged_sel %>%
    select(PMBB_ID, Age) %>%
    rename(Case_Age = Age)

# Prioritize case diagnosis age, fall back to sample age
all_ages <- pmbb_ages %>%
    left_join(case_ages, by = "PMBB_ID") %>%
    mutate(Age = coalesce(Case_Age, Age)) %>%
    select(PMBB_ID, Age)

# Add PCs (10)
# cov_imp_exome <- all_ages %>%
#     left_join(exome, by = c("PMBB_ID" = "person_id")) %>%
#     left_join(imp, by = c("PMBB_ID" = "person_id")) %>%
#     rename(IID = PMBB_ID)

cov_imp <- all_ages %>%
    inner_join(imp, by = c("PMBB_ID" = "person_id")) %>%
    rename(IID = PMBB_ID)

# cov_exome <- all_ages %>%
#     left_join(exome, by = c("PMBB_ID" = "person_id")) %>%
#     rename(IID = PMBB_ID)

# write.table(cov_imp_exome,
#             here("simplexo", "data", "simplexo_imp10_exome10_covariates.txt"),
#             row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
#
write.table(cov_imp,
            here("simplexo", "data", "simplexo_all_imp10_covariates.txt"),
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

# write.table(cov_exome,
#             here("simplexo", "data", "simplexo_all_exome10_covariates.txt"),
#             row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

# # Add PCs (5)
# cov_imp5_exome5 <- all_ages %>%
#     left_join(exome[, 1:6], by = c("PMBB_ID" = "person_id")) %>%
#     left_join(imp[, 1:6], by = c("PMBB_ID" = "person_id")) %>%
#     rename(IID = PMBB_ID)
#
# cov_imp5 <- all_ages %>%
#     left_join(imp[, 1:6], by = c("PMBB_ID" = "person_id")) %>%
#     rename(IID = PMBB_ID)
#
# cov_exome5 <- all_ages %>%
#     left_join(exome[, 1:6], by = c("PMBB_ID" = "person_id")) %>%
#     rename(IID = PMBB_ID)
#
# write.table(cov_imp5_exome5,
#             here("simplexo", "data", "simplexo_imp5_exome5_covariates.txt"),
#             row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
#
# write.table(cov_imp5,
#             here("simplexo", "data", "simplexo_imp5_covariates.txt"),
#             row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
#
# write.table(cov_exome5,
#             here("simplexo", "data", "simplexo_exome5_covariates.txt"),
#             row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")


# ========================
# CREATE STRATIFIED COVARIATE FILES
# ========================
imp <- imp %>% filter(!is.na(imputed_PC1))
dim(imp)

scase_50 <- readLines(here("simplexo", "data", "simplexo_50_case_ids.txt"))
case_malig <- readLines(here("simplexo", "data", "simplexo_malig_case_ids.txt"))
case_er_neg <- readLines(here("simplexo", "data", "simplexo_er_neg_case_ids.txt"))
case_er_pos <- readLines(here("simplexo", "data", "simplexo_er_pos_case_ids.txt"))
case_fhx <- readLines(here("simplexo", "data", "simplexo_fhx_case_ids.txt"))

cat("Loaded stratified case lists:\n")
cat("  Age ≤50:", length(case_50), "\n")
cat("  Invasive:", length(case_malig), "\n")
cat("  ER-:", length(case_er_neg), "\n")
cat("  ER+:", length(case_er_pos), "\n")
cat("  Family history:", length(case_fhx), "\n")

create_stratified_covariates <- function(case_ids, analysis_name) {
    # Filter case ages for this subset
    subset_case_ages <- merged_sel %>%
        filter(PMBB_ID %in% case_ids) %>%
        select(PMBB_ID, Age) %>%
        rename(Case_Age = Age)

    # Get all PMBB ages and prioritize case diagnosis age
    subset_all_ages <- pmbb_ages %>%
        left_join(subset_case_ages, by = "PMBB_ID") %>%
        mutate(Age = coalesce(Case_Age, Age)) %>%
        select(PMBB_ID, Age)

    # Add imputed PCs (10)
    subset_cov_imp <- subset_all_ages %>%
        inner_join(imp, by = c("PMBB_ID" = "person_id")) %>%
        rename(IID = PMBB_ID)

    write.table(subset_cov_imp,
                here("simplexo", "data", paste0("simplexo_", analysis_name, "_imp10_covariates.txt")),
                row.names = FALSE, col.names = TRUE, quote = FALSE, sep = " ")

    return(subset_cov_imp)
}

create_stratified_covariates(case_50, "50")
create_stratified_covariates(case_malig, "malig")
create_stratified_covariates(case_er_neg, "er_neg")
create_stratified_covariates(case_er_pos, "er_pos")
create_stratified_covariates(case_fhx, "fhx")

