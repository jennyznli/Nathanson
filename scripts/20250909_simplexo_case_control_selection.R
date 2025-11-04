# ========================
# SIMPLEXO - CASE & CONTROL SELECTION
# ========================

library(here)
setwd(here("simplexo"))

source(here("R", "load_packages.R"))
source(here("R", "sample_selection.R"))
source(here("R", "control_selection.R"))

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

print(paste("Initial Progeny cases:", nrow(progeny_case)))  # 5046

# Merge duplicate rows (separate entries for each breast)
print(paste("Duplicate entries:", sum(duplicated(progeny_case$SampNum))))

progeny_case_merged <- progeny_case %>%
    group_by(SampNum) %>%
    summarise(
        across(-any_of("SampNum"), ~ paste(unique(as.character(.x)), collapse = "; ")),
        .groups = "drop"
    )

print(paste("Unique Progeny cases after merging:", nrow(progeny_case_merged)))  # 4273

# ========================
# FINAL CASES - MERGE ALL
# ========================
### merge Progeny with PMBB ###
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

###### PROGENY ######
progeny_fh <- read_excel(
    here("simplexo", "ss", "br_pts_for_exwas_10022025.xlsx")
)[, c("SampNum", "Gender", "NumFDR", "NumFDR_BC", "NumSDR", "NumSDR_BC")]

print(paste("Progeny family history entries:", nrow(progeny_fh)))  # 5046

# Merge duplicate rows
progeny_fh <- progeny_fh %>%
    group_by(SampNum) %>%
    summarise(
        Gender = unique(Gender),
        NumFDR = unique(NumFDR),
        NumFDR_BC = unique(NumFDR_BC),
        NumSDR = unique(NumSDR),
        NumSDR_BC = unique(NumSDR_BC),
        .groups = "drop"
    )

print(paste("Unique Progeny family history cases:", nrow(progeny_fh)))  # 4273

# Convert 888 (unknown) to NA
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

sum(progeny_fh$Family_History == 1)
dim(progeny_fh)
# 2606/4273 in progeny w/ family history

# merge with PMBB
progeny_pmbb_fh <- merge(progeny_fh, up, by = "SampNum")
dim(progeny_pmbb)
# 1520 total progeny in pmbb
progeny_pmbb_fh <- progeny_pmbb_fh %>% filter(Family_History == 1)
dim(progeny_pmbb_fh)
# 1126
progeny_pmbb_fh <- progeny_pmbb_fh %>% filter(Gender == "F")
dim(progeny_pmbb_fh)
# 1100/1561

write.csv(progeny_pmbb_fh, here("simplexo", "data", "simplexo_progeny_family_history_df.csv"))

###### PMCR ######

hx <- read.csv(here("simplexo", "ss", "pmbb_1093_familyhx.csv"))

print(paste("Total family history records:", nrow(hx)))  # 41280
print(paste("Unique individuals with family history:", length(unique(hx$PMBB_ID))))  # 4285, original number

# Filter for breast cancer family history
breast_hx <- hx %>% filter(grepl("breast", MEDICAL_HX, ignore.case = TRUE))
pmcr_hx_ids <- unique(breast_hx$PMBB_ID)

print(paste("Cases with breast cancer family history:", length(pmcr_hx_ids)))  # 2228

# Combine all family history IDs from Progeny and PMCR
hx_ids <- sort(unique(c(progeny_pmbb_fh$PMBB_ID, pmcr_hx_ids)))
print(paste("Total unique family history IDs:", length(hx_ids)))  # 3063

# Intersect with final case IDs (because initial list given to Colleen was w/ 3 instances)
final_hx_ids <- sort(intersect(all_ids, hx_ids))
print(paste("Final family history case IDs:", length(final_hx_ids)))  # 2980

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
progeny_pmbb_unmerged <- merge(progeny_case, up, by = "SampNum")
progeny_pmbb_unmerged <- progeny_pmbb_unmerged %>% filter(Gender == "F")
print(paste("Female Progeny samples:", nrow(progeny_pmbb_unmerged)))  # 1841

progeny_er <- progeny_pmbb_unmerged %>% filter(ERstatus %in% c("Positive", "Negative"))
dim(progeny_er)
# 1023

progeny_er_merged <- progeny_er %>%
    group_by(SampNum) %>%
    summarise(
        across(-any_of("SampNum"), ~ paste(unique(as.character(.x)), collapse = ";")),
        .groups = "drop"
    )
dim(progeny_er_merged)
# 919

# filter out those w/ both positive and negative status
progeny_er_merged <- progeny_er_merged %>% filter(!(ERstatus %in% c("Positive;Negative","Negative;Positive")))
dim(progeny_er_merged)
# 898

progeny_er_neg <- progeny_er_merged %>% filter(ERstatus == "Negative")
progeny_er_pos <- progeny_er_merged %>% filter(ERstatus == "Positive")

print(paste("Progeny ER positive:", nrow(progeny_er_pos)))  # 593
print(paste("Progeny ER negative:", nrow(progeny_er_neg)))  # 305

##### PMCR #####

er <- read.csv(here("simplexo", "ss", "pmbb_1093_brca_er_20251024.csv"))

print(paste("Total ER status records:", nrow(er)))  # 4563
print(paste("Unique individuals:", length(unique(er$PMBB_ID))))  # 4285

# filter out indeterminates and nulls
er <- er %>% filter(EstrogenReceptorSummarry %in% c("Positive", "Negative"))

# Collapse data by PMBB_ID
er_merged <- er %>%
    group_by(PMBB_ID) %>%
    summarise(across(everything(), ~ {
        if (length(unique(.)) == 1) {
            unique(.)[1]
        } else {
            paste(unique(.), collapse = ";")
        }
    }))
dim(er_merged)
# 1849

# filter out those w/ both positive and negative status
er_merged <- er_merged %>% filter(!(EstrogenReceptorSummarry %in% c("Positive;Negative","Negative;Positive")))
dim(er_merged)
# 1834

er_pos <- er_merged %>% filter(EstrogenReceptorSummarry == "Positive")
er_neg <- er_merged %>% filter(EstrogenReceptorSummarry == "Negative")
length(er_pos$PMBB_ID) # 1468
length(er_neg$PMBB_ID) # 366

# Combine ER status from both sources but intersect w/ overall case list
pos_ids <- sort(intersect(unique(c(er_pos$PMBB_ID, progeny_er_pos$PMBB_ID)), all_ids))
neg_ids <- sort(intersect(unique(c(er_neg$PMBB_ID, progeny_er_neg$PMBB_ID)), all_ids))

print(paste("Final ER positive cases:", length(pos_ids))) # 1876
print(paste("Final ER negative cases:", length(neg_ids))) # 620

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
progeny_pmbb_insitu <- progeny_pmbb %>%
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

# Exclude in situ cases
progeny_no_insitu <- progeny_pmbb_insitu %>% filter(!is_insitu)

print(paste("Progeny invasive only:", nrow(progeny_no_insitu)))  # 1261

# Combine invasive cases
inv_only_ids <- unique(c(progeny_no_insitu$PMBB_ID, inv_only$filtered_patients$person_id))

print(paste("Total invasive only cases:", length(inv_only_ids)))  # 3780

write.table(sort(inv_only_ids),
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

young <- ages %>% filter(Age <= 50)

print(paste("Cases age <= 50:", nrow(young)))  # 1900

write.table(sort(young$PMBB_ID),
            here("simplexo", "data", "simplexo_50_case_ids.txt"),
            quote = FALSE, row.names = FALSE, col.names = FALSE)

