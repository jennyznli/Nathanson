# ========================
# SIMPLEXO - CASE & CONTROL SELECTION
# ========================

library(here)
setwd(here("simplexo"))

source(here("R", "load_packages.R"))
source(here("R", "sample_selection.R"))
source(here("R", "control_selection.R"))


# ========================
# FUNCTIONS
# ========================
# USAGE: progeny_er_merged <- merge_duplicates(progeny_er, "SampNum")
merge_duplicates <- function(df, group_col) {
    df %>%
        group_by(across(all_of(group_col))) %>%
        summarise(
            across(-any_of(group_col), ~ paste(unique(as.character(.x)), collapse = ";")),
            .groups = "drop"
        )
}

# ========================
# BREAST CANCER ICD SELECTION
# ========================

breast_results <- select_samples(
    sample_name = "simplexo",
    icd_codes = c("^C50", "^Z85.3", "^D05", "^Z86.000",
                  "^174", "^V10.3", "^233.0"),
    gender_filter = "Female",
    crep_filter = NULL,
    min_instances = 3,
    min_timespan = NULL,
    age_filter = NULL,
    exclude = FALSE,
    pmbb_dir = here("PMBB"),
    data_dir = here("simplexo", "data"),
    log_dir = here("simplexo", "log")
)

# ========================
# MALIGNANT NEOPLASM CODES (FOR EXCLUSION)
# ========================

malig_neoplasms <- c(
    "^C(?!44)",           # Malignant neoplasms except skin cancer (C44)
    "^Z85(?!\\.828)",     # Personal history of malignant neoplasm except skin
    "^D05",               # Carcinoma in situ of breast
    "^Z86.000",           # Personal history of in-situ neoplasm of breast
    "^(?!173)(?:1[4-9][0-9]|20[0-8]|209\\.[0-3])",  # ICD-9 malignant neoplasms except 173
    "^233.0",             # ICD-9 carcinoma in situ of breast
    "^V10(?!\\.83)"       # ICD-9 personal history of malignant neoplasm except skin
)

# ========================
# CONTROL SELECTION
# ========================

controls <- select_controls(
    control_name = "simplexo",
    exclude_codes = malig_neoplasms,
    gender_filter = "Female",
    age_filter = NULL,
    crep_filter = FALSE,
    pmbb_dir = here("PMBB", "3.0"),
    data_dir = here("simplexo", "data"),
    log_dir = here("simplexo", "log")
)
breast_controls <- controls

# ========================
# LOAD DATA
# ========================

PMBB_DIR <- here("PMBB", "3.0")
OCC <- here(PMBB_DIR, "PMBB-Release-2024-3.0_phenotype_condition_occurrence.txt")
COV <- here(PMBB_DIR, "PMBB-Release-2024-3.0_covariates.txt")
PER <- here(PMBB_DIR, "PMBB-Release-2024-3.0_phenotype_person.txt")

cov <- fread(COV, header = TRUE)
person <- fread(PER, header = TRUE)
flags <- fread(here("PMBB", "3.0", "rgcname_pmbbid_metadata_flags.csv"))

progeny_case <- read_excel(here("simplexo", "ss", "br_pts_for_exwas_10022025.xlsx"))
pmcr <- read.csv(here("simplexo", "ss", "pmbb_147_pmcrbreastage.csv"))

# ========================
# MATCH PMBB TO UPXXXX IDS
# ========================

up_pattern <- "^UPENN-PMBB_UP[0-9]+_UP[0-9]+$"

# Find all cases that match pattern among ID columns
up <- flags %>%
    mutate(match_col = apply(
        select(., RGC_sample_name, ID1, ID2, ID3, ID4), 1,
        function(row) {
            m <- row[grepl(up_pattern, row)]
            if (length(m) > 0) m[1] else NA
        }
    )) %>%
    filter(!is.na(match_col))

print(paste("UP ID matches found:", nrow(up)))  # 2864

# get relevant columns
up <- select(up, PMBB_ID, RGC_sample_name, ID1, ID2, ID3, ID4, CREP, match_col)

# get VCFID
up$VCFID <- sub(".*(UP[0-9]{4}).*", "\\1", up$match_col)

# final mapping table
up <- select(up, PMBB_ID, match_col, VCFID)
up$SampNum <- as.numeric(gsub("[a-zA-Z]", "", up$VCFID))

print(paste("Total PMBB samples with UP ID:", nrow(up)))  # 2864

write.csv(up, here("simplexo", "data", "simplexo_up_map.csv"), row.names = FALSE)

# ========================
# PROGENY - CASE PULL
# ========================

progeny_case_merged <- merge_duplicates(progeny_case, "SampNum")

print(paste("Initial Progeny cases:", nrow(progeny_case)))  # 5046
print(paste("Duplicate entries:", sum(duplicated(progeny_case$SampNum))))
print(paste("Unique Progeny cases after merging:", nrow(progeny_case_merged)))  # 4273

# ========================
# FINAL CASES - MERGE ALL
# ========================
### left join Progeny with PMBB ###
progeny_pmbb <- merge(progeny_case_merged, up, by = "SampNum")
print(paste("Progeny samples with PMBB IDs:", nrow(progeny_pmbb)))  # 1561

# filter by females
progeny_pmbb <- progeny_pmbb %>% filter(Gender == "F")
print(paste("Female Progeny samples:", nrow(progeny_pmbb)))  # 1520

write.csv(progeny_pmbb, here("simplexo", "data", "simplexo_progeny_pmbb.txt"), row.names = FALSE)

progeny_pmbb_ids <- sort(progeny_pmbb$PMBB_ID)

### Merge with ICD cases ###
icd_cases <- fread(here("simplexo", "data", "simplexo_cancer_filtered_patients.txt"))
icd_case_ids <- readLines(here("simplexo", "data", "simplexo_cancer_filtered_patients_ids.txt"))
print(paste("ICD case IDs:", length(icd_case_ids)))  # 3448

all_ids <- sort(unique(c(progeny_pmbb_ids, icd_case_ids)))
print(paste("Total unique case IDs:", length(all_ids)))  # 4063

write.table(all_ids,
            here("simplexo", "data", "simplexo_overall_v1_case_ids.txt"),
            row.names = FALSE, col.names = FALSE, quote = FALSE)
# not finished since we need to ensure we have all age data, so run other covariates script

### RUN OTHER COVARIATE SCRIPT... ###

# ========================
# FINALIZE CONTROLS
# ========================

# Exclude cases from controls
controls <- breast_controls$final_controls %>%
    filter(!(person_id %in% all_ids)) %>%
    select(person_id, Sample_age, Sequenced_gender)

print(paste("Final control count:", nrow(controls)))  # 18278

write.table(controls$person_id,
            here("simplexo", "data", "simplexo_overall_control_ids.txt"),
            quote = FALSE, row.names = FALSE, col.names = FALSE)

# ========================================================================
# FAMILY HISTORY
# ========================================================================

all_ids <- readLines(here("simplexo", "data", "simplexo_overall_case_ids.txt"))
length(all_ids)
# 4059

###### PROGENY ######
progeny_fh <- read_excel(
    here("simplexo", "ss", "br_pts_for_exwas_10022025.xlsx")
)[, c("SampNum", "Gender", "NumFDR", "NumFDR_BC", "NumSDR", "NumSDR_BC")]

print(paste("Progeny family history entries:", nrow(progeny_fh)))  # 5046

progeny_fh <- merge_duplicates(progeny_fh, "SampNum")

print(paste("Unique Progeny family history cases:", nrow(progeny_fh)))  # 4273

# Convert 888 (unknown) to NA
progeny_fh$NumFDR <- as.numeric(progeny_fh$NumFDR)
progeny_fh$NumFDR_BC <- as.numeric(progeny_fh$NumFDR_BC)
progeny_fh$NumSDR <- as.numeric(progeny_fh$NumSDR)
progeny_fh$NumSDR_BC <- as.numeric(progeny_fh$NumSDR_BC)

progeny_fh <- progeny_fh %>%
    mutate(
        NumFDR = na_if(NumFDR, 888),
        NumFDR_BC = na_if(NumFDR_BC, 888),
        NumSDR = na_if(NumSDR, 888),
        NumSDR_BC = na_if(NumSDR_BC, 888)
    )

# Mark cases with family history
progeny_fh <- progeny_fh %>%
    mutate(
        Family_History = if_else(
            coalesce(NumFDR_BC, 0) > 0 | coalesce(NumSDR_BC, 0) > 0,
            1L, 0L
        )
    )

# merge with PMBB
progeny_pmbb_fh <- merge(progeny_fh, up, by = "SampNum") %>% filter(Gender == "F")
dim(progeny_pmbb_fh)

# 1520 total progeny in pmbb
progeny_pmbb_fh <- progeny_pmbb_fh %>% filter(Family_History == 1)
dim(progeny_pmbb_fh)
# 1100

progeny_pmbb_fh <- progeny_pmbb_fh %>% filter(PMBB_ID %in% all_ids)
dim(progeny_pmbb_fh)
# 1096
write.csv(progeny_pmbb_fh, here("simplexo", "data", "simplexo_progeny_family_history_df.csv"))

###### PMCR ######

hx <- read.csv(here("simplexo", "ss", "pmbb_1093_familyhx.csv"))

print(paste("Total family history records:", nrow(hx)))  # 41280
print(paste("Unique individuals with family history:", length(unique(hx$PMBB_ID))))  # 4285, original number

# Filter for breast cancer family history
breast_hx <- hx %>% filter(grepl("breast", MEDICAL_HX, ignore.case = TRUE))
breast_hx <- breast_hx %>% filter(!(RELATION %in% c("null", "Negative History")))
print(paste("Cases with breast cancer family history:", length(unique(breast_hx$PMBB_ID))))  # 1611

# Combine all family history IDs from Progeny and PMCR
pmcr_hx_ids <- unique(breast_hx$PMBB_ID)
hx_ids <- sort(unique(c(progeny_pmbb_fh$PMBB_ID, pmcr_hx_ids)))
print(paste("Total unique family history IDs:", length(hx_ids)))  # 2448

# Intersect with final case IDs (because initial list given to Colleen was w/ 3 instances)
final_hx_ids <- sort(intersect(all_ids, hx_ids))
print(paste("Final family history case IDs:", length(final_hx_ids)))  # 2382

write.table(final_hx_ids,
            here("simplexo", "data", "simplexo_fhx_case_ids.txt"),
            quote = FALSE, row.names = FALSE, col.names = FALSE)

###### ICD CODES ######
# family_history <- select_samples(
#     sample_name = "simplexo_fh",
#     icd_codes = c("^Z80.3", "^V16.3"),
#     gender_filter = "Female",
#     crep_filter = NULL,
#     min_instances = 1,
#     min_timespan = NULL,
#     age_filter = NULL,
#     exclude = FALSE,
#     pmbb_dir = here("PMBB"),
#     data_dir = here("simplexo", "data"),
#     log_dir = here("simplexo", "log")
# )
# result in 0 cases...

# ========================================================================
# ER STATUS
# ========================================================================
##### PROGENY #####
# take the unmerged version, filter out the Not Done, Not Reported, NA, Not Reported
progeny_pmbb_unmerged <- merge(progeny_case, up, by = "SampNum") %>% filter(Gender == "F")
print(paste("Female Progeny cases:", nrow(progeny_pmbb_unmerged)))  # 1841

# Filter to valid ER status
progeny_er <- progeny_pmbb_unmerged %>%
    filter(ERstatus %in% c("Positive", "Negative"))
cat("Progeny with ER status:", nrow(progeny_er), "\n") # 1023

# Merge duplicates
progeny_er_merged <- merge_duplicates(progeny_er, "SampNum")
cat("After merging duplicates:", nrow(progeny_er_merged), "\n") # 919

# Remove conflicting status (both positive and negative)
conflicting_progeny <- progeny_er_merged %>%
    filter(ERstatus %in% c("Positive;Negative", "Negative;Positive"))
cat("Progeny with conflicting ER status:", nrow(conflicting_progeny), "\n") # 21

progeny_er_clean <- progeny_er_merged %>%
    filter(!(ERstatus %in% c("Positive;Negative", "Negative;Positive")))

progeny_er_pos <- progeny_er_clean %>% filter(ERstatus == "Positive")
progeny_er_neg <- progeny_er_clean %>% filter(ERstatus == "Negative")
cat("Progeny ER+:", nrow(progeny_er_pos), "\n") # 593
cat("Progeny ER-:", nrow(progeny_er_neg), "\n") # 305

##### PMCR #####
er <- read.csv(here("simplexo", "ss", "pmbb_1093_brca_er_20251024.csv"))
print(paste("Total ER status records:", nrow(er)))  # 4563
print(paste("Unique individuals:", length(unique(er$PMBB_ID))))  # 4285

# Filter to valid ER status
er <- er %>% filter(EstrogenReceptorSummarry %in% c("Positive", "Negative"))
er_merged <- merge_duplicates(er, "PMBB_ID")
cat("After merging duplicates:", nrow(er_merged), "\n")
# 1849

# Remove conflicting status
conflicting_pmcr <- er_merged %>%
    filter(EstrogenReceptorSummarry %in% c("Positive;Negative", "Negative;Positive"))
cat("PMCR with conflicting ER status:", nrow(conflicting_pmcr), "\n") # 15

er_clean <- er_merged %>%
    filter(!(EstrogenReceptorSummarry %in% c("Positive;Negative", "Negative;Positive")))

er_pos <- er_clean %>% filter(EstrogenReceptorSummarry == "Positive")
er_neg <- er_clean %>% filter(EstrogenReceptorSummarry == "Negative")
cat("PMCR ER+:", nrow(er_pos), "\n") # 1468
cat("PMCR ER-:", nrow(er_neg), "\n") # 366

##### COMBINE #####
pos_ids_combined <- unique(c(er_pos$PMBB_ID, progeny_er_pos$PMBB_ID))
neg_ids_combined <- unique(c(er_neg$PMBB_ID, progeny_er_neg$PMBB_ID))

# Check for conflicts between sources
conflicts <- intersect(pos_ids_combined, neg_ids_combined)
cat("\nConflicts between sources:", length(conflicts), "\n") # 5

if (length(conflicts) > 0) {
    cat("WARNING: Some people have conflicting ER status between Progeny and PMCR\n")
    cat("Excluding these", length(conflicts), "individuals from both lists\n")
    pos_ids_combined <- setdiff(pos_ids_combined, conflicts)
    neg_ids_combined <- setdiff(neg_ids_combined, conflicts)
}

# Intersect with overall case list
pos_ids <- sort(intersect(pos_ids_combined, all_ids))
neg_ids <- sort(intersect(neg_ids_combined, all_ids))

cat("\nFinal counts:\n")
cat("ER+ cases:", length(pos_ids), "\n") #1871
cat("ER- cases:", length(neg_ids), "\n") #615
cat("Overlap (in both lists):", length(intersect(pos_ids, neg_ids)), "\n")

write.table(pos_ids,
            here("simplexo", "data", "simplexo_er_pos_case_ids.txt"),
            quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(neg_ids,
            here("simplexo", "data", "simplexo_er_neg_case_ids.txt"),
            quote = FALSE, row.names = FALSE, col.names = FALSE)

# ========================================================================
# INVASIVE ONLY (NO IN SITU)
# ========================================================================

##### ICD CODES #####
inv_only <- select_samples(
    sample_name = "simplexo_malig",
    icd_codes = c("^C50", "^Z85.4", "^174", "^V10.3"),
    gender_filter = "Female",
    crep_filter = NULL,
    min_instances = 3,
    min_timespan = NULL,
    age_filter = NULL,
    exclude = FALSE,
    pmbb_dir = here("PMBB"),
    data_dir = here("simplexo", "data"),
    log_dir = here("simplexo", "log")
)

##### PROGENY #####
progeny_pmbb_flag <- progeny_pmbb_unmerged %>%
    mutate(
        is_insitu_dx = str_detect(
            Diagnosis,
            regex("DCIS|LCIS|ductal carcinoma in situ|in situ|lobular carcinoma in situ",
                  ignore_case = TRUE)
        ),
        is_insitu_morph = str_detect(
            Morphology,
            regex("DCIS|LCIS|ductal carcinoma in situ|in situ|lobular carcinoma in situ",
                  ignore_case = TRUE)
        ),
        is_insitu = is_insitu_dx | is_insitu_morph,
        is_insitu = replace_na(is_insitu, FALSE)
    )

# Exclude in situ cases and benign...
progeny_invasive_records <- progeny_pmbb_flag %>%
    filter(!is_insitu) %>%
    filter(!(Diagnosis %in% c("Normal Benign Tissue", "99.Not Reported", "Not Reported"))) %>%
    filter(Gender == "F")

cat("Progeny invasive records (before merge):", nrow(progeny_invasive_records), "\n") # 1529
progeny_invasive <- merge_duplicates(progeny_invasive_records, "SampNum")
cat("Unique Progeny patients with invasive:", nrow(progeny_invasive), "\n") # 1320

##### COMBINE #####
icd_ids <- inv_only$filtered_patients$person_id
inv_only_ids <- unique(c(progeny_invasive$PMBB_ID, icd_ids))
final_inv_ids <- sort(intersect(inv_only_ids, all_ids))
cat("Final patients w/ invasive:", length(final_inv_ids), "\n") # 3791

write.table(final_inv_ids,
            here("simplexo", "data", "simplexo_malig_case_ids.txt"),
            quote = FALSE, row.names = FALSE, col.names = FALSE)

# ========================================================================
# UNDER 50
# ========================================================================
ages <- read.table(here("simplexo", "data", "simplexo_overall_case_ages_merged.txt"),
                   header = TRUE)

ages <- ages %>%
    select(PMBB_ID, Age) %>%
    filter(!is.na(Age)) %>%
    mutate(Age = round(Age, 2))

print(paste("Cases with age information:", nrow(ages)))

young <- ages %>% filter(Age <= 50) %>% filter(PMBB_ID %in% all_ids)

print(paste("Cases age <= 50:", length(unique(young$PMBB_ID))))  # 1900

write.table(sort(unique(young$PMBB_ID)),
            here("simplexo", "data", "simplexo_50_case_ids.txt"),
            quote = FALSE, row.names = FALSE, col.names = FALSE)

