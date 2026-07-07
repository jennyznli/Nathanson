# ========================
# SIMPLEXO - F3 vs. F4
# ========================
library(here)
setwd(here("simplexo"))
source(here("R/config.R"))

pmbb4 <- here("PMBB", "4.0")

# covariates
cov    <- fread(file.path(pmbb4, "PMBB-Release-2026-4.0_phenotype_covariates.txt"),           header = TRUE)
person <- fread(file.path(pmbb4, "PMBB-Release-2026-4.0_phenotype_person.txt"),               header = TRUE)

# ICD
obs    <- fread(file.path(pmbb4, "PMBB-Release-2026-4.0_phenotype_observation.txt"),          header = TRUE)
cond   <- fread(file.path(pmbb4, "PMBB-Release-2026-4.0_phenotype_condition_occurrence.txt"), header = TRUE)

# pmcr
brca   <- fread(file.path(pmbb4, "PMBB-Release-2026-4.0_phenotype_cancer_brca.txt"),          header = TRUE)

# progeny
progeny_case <- read_excel(here("simplexo", "ss", "br_pts_for_exwas_10022025.xlsx"))

# mapping
flags <- read.csv(here("PMBB", "3.0", "rgcname_pmbbid_metadata_flags.csv"))
up <- read.csv(here("simplexo", "data", "simplexo_up_map.csv"))

# ========================
# PROGENY - PREPROCESS
# ========================
dim(progeny_case)
# 5046

progeny_case <- progeny_case %>% filter(Gender == "F") %>%
    merge_duplicates("SampNum") %>%
    inner_join(up, by = "SampNum")
dim(progeny_case)
# 1520

progeny_ids <- sort(progeny_pmbb$PMBB_ID)

# ============================================================
# PMCR - PREPROCESS
# ============================================================
pmcr4 <- fread(file.path(pmbb4, "PMBB-Release-2026-4.0_phenotype_cancer_pmcr.txt")) %>%
    filter(SeerSiteRecode2023ExpandedGroup == "Breast") %>%
    left_join(brca, by = "person_id") %>%
    left_join(cov,  by = "person_id") %>%
    filter(sequenced_gender == "Female")
dim(pmcr4)
# 3520

pmcr3 <- fread(here("PMBB", "3.0", "PMBB-Release-2024-3.1_phenotype_cancer_PMCR")) %>%
    filter(SeerSiteRecode2023ExpandedGroup == "Breast") %>%
    left_join(brca, by = "person_id") %>%
    left_join(cov,  by = "person_id") %>%
    filter(sequenced_gender == "Female")
dim(pmcr3)
# 2936

prop.table(table(pmcr3$ReportingHospital))
# CCH        HUP        LGH        PAH       PMPH       PPMC
# 0.03717599 0.76909959 0.03069577 0.13028649 0.01773533 0.01500682

prop.table(table(pmcr4$ReportingHospital))
# CCH        HUP        LGH        PAH       PMPH       PPMC
# 0.03779483 0.76499005 0.02728048 0.13384484 0.01648196 0.01960784

# merge duplicate tumors
pmcr3 <- pmcr3 %>% merge_duplicates("person_id")
pmcr4 <- pmcr4 %>% merge_duplicates("person_id")
dim(pmcr3)
# 2553
dim(pmcr4)
# 3090

# extract IDs
pmcr3_ids <- unique(as.character(pmcr3$person_id))
pmcr4_ids <- unique(as.character(pmcr4$person_id))

# ========================
# ICD SELECTION
# ========================
source(here("R", "sample_selection_f3.R"))
breast_results_f3 <- select_samples_f3(
    sample_name = "simplexo",
    icd_codes = c("^C50", "^Z85.3", # ICD10 - malignant
                  "^D05", "^Z86.000", # ICD10 - DCIS
                  "^174", "^V10.3", "^233.0"), #ICD9 - both
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
f3_df <- breast_results_f3$filtered_patients
f3_ids <- breast_results_f3$filtered_patients$person_id
dim(f3_df)
# 3448

source(here("R", "sample_selection_f4.R"))
breast_results_f4 <- select_samples_f4(
    sample_name = "simplexo",
    icd_codes = c("^C50", "^Z85.3", # ICD10 - malignant
                  "^D05", "^Z86.000", # ICD10 - DCIS
                  "^174", "^V10.3", "^233.0"), #ICD9 - both
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
f4_df <- breast_results_f4$filtered_patients
f4_ids <- breast_results_f4$filtered_patients$person_id
dim(f4_df)
# 45822

# ============================================================
# UPSET PLOT
# ============================================================
library(UpSetR)


id_sets <- list(
    PMCR3 = unique(as.character(pmcr3_ids)),
    PMCR4 = unique(as.character(pmcr4_ids)),
    ICD3    = unique(as.character(f3_ids)),
    ICD4    = unique(as.character(f4_ids))
)

png(file.path(here("simplexo", "figures", "upset_f3_f4.png")),
    width = 10, height = 6, units = "in", res = 300)

upset(
    fromList(id_sets),
    sets           = c("PMCR3", "PMCR4", "ICD3", "ICD4"),
    keep.order     = TRUE,
    order.by       = "freq",
    text.scale     = 1.4,
    mainbar.y.label = "Patients in intersection",
    sets.x.label    = "Patients per set"
)

dev.off()

id_sets <- list(
    PMCR3 = unique(as.character(pmcr3_ids)),
    PMCR4 = unique(as.character(pmcr4_ids)),
    ICD3    = unique(as.character(f3_ids)),
    ICD4    = unique(as.character(f4_ids)),
    Progeny = unique(as.character(progeny_ids))
)

png(file.path(here("simplexo", "figures", "upset_f3_f4_progeny.png")),
    width = 10, height = 6, units = "in", res = 300)

upset(
    fromList(id_sets),
    sets           = c("PMCR3", "PMCR4", "ICD3", "ICD4", "Progeny"),
    keep.order     = TRUE,
    order.by       = "freq",
    text.scale     = 1.4,
    mainbar.y.label = "Patients in intersection",
    sets.x.label    = "Patients per set"
)

dev.off()


id_sets <- list(
    PMCR4 = unique(as.character(pmcr4_ids)),
    ICD4    = unique(as.character(f4_ids)),
    Progeny = unique(as.character(progeny_ids))
)

png(file.path(here("simplexo", "figures", "upset_f4.png")),
    width = 7, height = 6, units = "in", res = 300)

upset(
    fromList(id_sets),
    sets           = c("PMCR4", "ICD4", "Progeny"),
    keep.order     = TRUE,
    order.by       = "freq",
    text.scale     = 1.4,
    mainbar.y.label = "Patients in intersection",
    sets.x.label    = "Patients per set"
)

dev.off()

# ============================================================
# EXAMINE - ICD MISSED BY PMCR
# ============================================================
# hypothesis - these will be personal history dominated?

icd_only <- membership %>% filter(ICD3 & ICD4 & !PMCR3 & !PMCR4) %>% pull(person_id)

f4_not_pmcr    <- setdiff(f4_ids, pmcr_ids)
f4_not_pmcr_df <- f4_df %>% filter(person_id %in% f4_not_pmcr)   # FIX: was f4_not_registry

# batch distribution of the extras vs the full f4 cohort
prop.table(table(f4_not_pmcr_df$batch))
prop.table(table(f4_df$batch))

# ============================================================
# EXAMINE - PMCR MISSED BY ICD
# ============================================================
icd_ids <- union(f3_ids, f4_ids)
pmcr4 <- pmcr4 %>% mutate(caught = person_id %in% icd_ids)

pmcr_missed <- setdiff(pmcr_ids, icd_ids)
pmcr_missed_df <- pmcr4 %>% filter(person_id %in% pmcr_missed)
dim(pmcr_missed_df)
#> 34

# ---- where are the missed cases reported? ----
table(pmcr_missed_df$ReportingHospital)
#>  CCH  HUP  LGH  PAH PMPH
#>    3   37   20    9    7

# ---- miss RATE per hospital (missed / total registry at that site) ----
# rate matters, not raw count: HUP dominates counts simply by volume.
pmcr4 %>%
    count(ReportingHospital, caught) %>%
    pivot_wider(names_from = caught, values_from = n, values_fill = 0) %>%
    rename(missed = `FALSE`, caught = `TRUE`) %>%
    mutate(miss_rate = missed / (missed + caught)) %>%
    arrange(desc(miss_rate))

# ---- how many qualifying ICD codes do the missed cases actually have? ----
# min_instances = 3, so a registry case is dropped for one of two reasons:
#   0 codes      -> genuinely ICD-invisible (outside care / path-only / sex filter)
#   1-2 codes    -> present but below the 3-instance threshold
all_patients <- breast_results_f4$all_patients
pmcr4_icd <- pmcr4 %>%
    left_join(all_patients, by = "person_id") %>%
    mutate(num_unique_codes = tidyr::replace_na(num_unique_codes, 0))

# whole registry cohort, by code bucket
pmcr4_icd %>%
    mutate(bucket = cut(num_unique_codes, c(-Inf, 0, 1, 2, Inf),
                        labels = c("0", "1", "2", "3+"))) %>%
    count(bucket)
#> 0     34
#> 1     37
#> 2    111
#> 3+  2943

# missed cases only, by code bucket
pmcr4_icd %>%
    filter(person_id %in% pmcr_missed) %>%
    mutate(bucket = cut(num_unique_codes, c(-Inf, 0, 1, 2, Inf),
                        labels = c("0", "1", "2", "3+"))) %>%
    count(bucket)
#> 0     34   <- no qualifying code at all
#> 1      2
#> 2      4
#> 3+    29   <- HAVE 3+ codes but still missed -> investigate (sex? code source?)


# ============================================================
# 7. WHAT THE f3 -> f4 METHODOLOGY CHANGE DID
# ============================================================
gained_f3_to_f4 <- setdiff(f4_ids, f3_ids)   #> gained 1138 (likely history codes)
lost_f3_to_f4   <- setdiff(f3_ids, f4_ids)   #> lost 4

gained_df <- f4_df %>%
    filter(person_id %in% gained_f3_to_f4) %>%
    filter(!(batch %in% c("3")))

# ---- restrict to freezes 1 & 2 to compare on a common population ----
# NOTE: f3_df uses `Batch` (capital), f4_df uses `batch` (lowercase).
f3_df_f12 <- f3_df %>% filter(Batch %in% c("1", "2"))
f4_df_f12 <- f4_df %>% filter(batch %in% c("1", "2"))
dim(f3_df_f12)   #> 3448
dim(f4_df_f12)   #> 3953  (still differs, likely the expanded ICD code set)

# who f4 picks up that f3 didn't, within freezes 1 & 2
diff_f4_vs_f3 <- anti_join(f4_df_f12, f3_df_f12, by = "person_id")

# ============================================================
# 8. ICD (f4) vs REGISTRY difference sets
# ============================================================
diff_f4_vs_pmcr <- anti_join(f4_df, pmcr4, by = "person_id")   # FIX: was pmcr_breast (undefined)




# ========================
# COMPARE ALL SOURCES
# ========================
# ICD selection
dim(breast_icd_df)
# 4582

# both have
sum(brca$person_id %in% breast_icd)
# 3481

# ones the ICD selection has that the PMBB sheet doesn't
icd_only <- setdiff(breast_icd, brca$person_id)
length(icd_only)
# 1527
icd_only_df <- breast_icd_df[!(breast_icd_df$person_id %in% brca$person_id), ]
table(icd_only_df$batch)
# 1   2   3
# 766 580 181

# ones the PMBB sheet has that the ICD selection doesn't
pmbb_only <- setdiff(brca$person_id, breast_icd)
length(pmbb_only)
# 69
pmbb_only_df <- brca[!(brca$person_id %in% breast_icd), ]
pmbb_only_df <- pmbb_only_df %>% left_join(cov, by = "person_id")
# 1  2  3
# 20 51  5
# only 5 from the most recent freeze?

# > table(pmbb_only_df$ReportingHospital)
#
# CCH  HUP  LGH  PAH PMPH
# 3   37   20    9    7

# > table(brca$ReportingHospital)
#
# CCH  HUP  LGH  PAH PMPH PPMC
# 135 2720   96  477   60   69
# so may be more non HUP?

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
source(here("R", "control_selection_f3.R"))
breast_controls_f3 <- select_controls_f3(
    control_name = "simplexo",
    exclude_codes = malig_neoplasms,
    gender_filter = "Female",
    age_filter = NULL,
    crep_filter = FALSE,
    pmbb_dir = here("PMBB"),
    data_dir = here("simplexo", "data"),
    log_dir = here("simplexo", "log")
)
f3c_df <- breast_controls_f3$excluded_icd_summary
f3c_ids <- breast_controls_f3$final_controls$person_id

source(here("R", "control_selection_f4.R"))
breast_controls_f4 <- select_controls_f4(
    control_name = "simplexo",
    exclude_codes = malig_neoplasms,
    gender_filter = "Female",
    age_filter = NULL,
    crep_filter = FALSE,
    pmbb_dir = here("PMBB", "4.0"),
    data_dir = here("simplexo", "data"),
    log_dir = here("simplexo", "log")
)
f4c_df <- breast_controls_f4$excluded_icd_summary
f4c_ids <- breast_controls_f4$final_controls$person_id

length(f3c_ids)
# 18278
length(f4c_ids)
# 18720

f4c_only <- f4c_ids[!(f4c_ids %in% f3c_ids)]
length(f4c_only)
# 4185

f4c_only_df <- f4c_df[f4c_only, ]

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

