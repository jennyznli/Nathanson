# ========================
# SIMPLEXO - F3 vs. F4
# ========================
library(here)
setwd(here("simplexo"))
source(here("R/config.R"))

pmbb4 <- here("PMBB", "4.0")

# covariates
cov <- fread(file.path(pmbb4, "PMBB-Release-2026-4.0_phenotype_covariates.txt"), header = TRUE)
person <- fread(file.path(pmbb4, "PMBB-Release-2026-4.0_phenotype_person.txt"), header = TRUE)

# ICD
obs <- fread(file.path(pmbb4, "PMBB-Release-2026-4.0_phenotype_observation.txt"), header = TRUE)
cond <- fread(file.path(pmbb4, "PMBB-Release-2026-4.0_phenotype_condition_occurrence.txt"), header = TRUE)

# breast cancer detailed info
brca <- fread(file.path(pmbb4, "PMBB-Release-2026-4.0_phenotype_cancer_brca.txt"), header = TRUE)
hx <- fread(file.path(pmbb4, "PMBB-Release-2026-4.0_phenotype_family_hx.txt"), header = TRUE)

# progeny
progeny_umerged <- read_excel(here("simplexo", "ss", "br_pts_for_exwas_10022025.xlsx"))
flags <- read.csv(here("PMBB", "3.0", "rgcname_pmbbid_metadata_flags.csv"))
up <- read.csv(here("simplexo", "data", "simplexo_up_map.csv"))

# id lists
exome_ids <- read.table(here("PMBB", "4.0", "PMBB-Release-2026-4.0_genetic_exome.sample_list.txt"), header = FALSE)$V1
imputed_ids <- read.csv(here("PMBB", "4.0", "PMBB-Release-2026-4.0_genetic_imputed.sample_list.txt"), header = FALSE)$V1
length(exome_ids)
# 70925
length(imputed_ids)
# 70493

seq_ids <- intersect(exome_ids, imputed_ids)
length(seq_ids)
# 70408

# ========================
# PROGENY - PREPROCESS
# ========================
dim(progeny_umerged)
# 5046

progeny_case <- progeny_umerged %>% filter(Gender == "F") %>%
    merge_duplicates("SampNum") %>%
    inner_join(up, by = "SampNum")
dim(progeny_case)
# 1520

progeny_ids <- sort(progeny_case$PMBB_ID)

# ============================================================
# PMCR - PREPROCESS
# ============================================================
pmcr3 <- fread(here("PMBB", "3.0", "PMBB-Release-2024-3.1_phenotype_cancer_PMCR")) %>%
    filter(SeerSiteRecode2023ExpandedGroup == "Breast") %>%
    left_join(brca, by = "person_id") %>%
    left_join(cov,  by = "person_id") %>%
    filter(sequenced_gender == "Female") %>%
    merge_duplicates("person_id")
dim(pmcr3)
# 2553

pmcr4 <- fread(file.path(pmbb4, "PMBB-Release-2026-4.0_phenotype_cancer_pmcr.txt")) %>%
    filter(SeerSiteRecode2023ExpandedGroup == "Breast") %>%
    left_join(brca, by = "person_id") %>%
    left_join(cov,  by = "person_id") %>%
    filter(sequenced_gender == "Female") %>%
    merge_duplicates("person_id")
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
# 4582

# ============================================================
# UPSET PLOT
# ============================================================
library(UpSetR)

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

# id_sets <- list(
#     PMCR3 = unique(as.character(pmcr3_ids)),
#     PMCR4 = unique(as.character(pmcr4_ids)),
#     ICD3    = unique(as.character(f3_ids)),
#     ICD4    = unique(as.character(f4_ids))
# )
#
# png(file.path(here("simplexo", "figures", "upset_f3_f4.png")),
#     width = 10, height = 6, units = "in", res = 300)
#
# upset(
#     fromList(id_sets),
#     sets           = c("PMCR3", "PMCR4", "ICD3", "ICD4"),
#     keep.order     = TRUE,
#     order.by       = "freq",
#     text.scale     = 1.4,
#     mainbar.y.label = "Patients in intersection",
#     sets.x.label    = "Patients per set"
# )
#
# dev.off()

# id_sets <- list(
#     PMCR4 = unique(as.character(pmcr4_ids)),
#     ICD4    = unique(as.character(f4_ids)),
#     Progeny = unique(as.character(progeny_ids))
# )
#
# png(file.path(here("simplexo", "figures", "upset_f4.png")),
#     width = 7, height = 6, units = "in", res = 300)
#
# upset(
#     fromList(id_sets),
#     sets           = c("PMCR4", "ICD4", "Progeny"),
#     keep.order     = TRUE,
#     order.by       = "freq",
#     text.scale     = 1.4,
#     mainbar.y.label = "Patients in intersection",
#     sets.x.label    = "Patients per set"
# )
#
# dev.off()

# ============================================================
# EXAMINE - PMCR MISSED BY ICD
# ============================================================
pmcr_missed <- setdiff(pmcr4_ids, f4_ids)
pmcr_missed_df <- pmcr4 %>% filter(person_id %in% pmcr_missed)
dim(pmcr_missed_df)
# 34

table(pmcr_missed_df$ReportingHospital)
# CCH  HUP  LGH  PAH PMPH
# 1    9   17    3    4
prop.table(table(pmcr_missed_df$ReportingHospital))
# CCH        HUP        LGH        PAH       PMPH
# 0.02941176 0.26470588 0.50000000 0.08823529 0.11764706

# compare proportions to normal distribution
pmcr42 <- fread(file.path(pmbb4, "PMBB-Release-2026-4.0_phenotype_cancer_pmcr.txt")) %>%
    filter(SeerSiteRecode2023ExpandedGroup == "Breast") %>%
    left_join(brca, by = "person_id") %>%
    left_join(cov,  by = "person_id") %>%
    filter(sequenced_gender == "Female")
table(pmcr42$ReportingHospital)
# CCH  HUP  LGH  PAH PMPH PPMC
# 133 2692   96  471   58   69
prop.table(table(pmcr42$ReportingHospital))
# CCH        HUP        LGH        PAH       PMPH       PPMC
# 0.03779483 0.76499005 0.02728048 0.13384484 0.01648196 0.01960784

# ============================================================
# EXAMINING F4 ONLY
# ============================================================
icd4_only_ids <- setdiff(f4_ids, union(pmcr4_ids, progeny_ids))
gained_f3_to_f4_ids <- setdiff(icd4_only_ids, f3_ids)
existing_by_f3_ids  <- f4_df %>% filter(batch %in% c("1", "2")) %>%
    pull(person_id) %>% unique()
newly_id_existing_ids <- intersect(gained_f3_to_f4_ids, existing_by_f3_ids)  # methodology effect
new_enrollee_ids      <- setdiff(gained_f3_to_f4_ids, existing_by_f3_ids)    # trivially new patients

cat("Gained f3 -> f4:", length(gained_f3_to_f4_ids), "\n") #1138
cat("  - already enrolled by freeze 3 (methodology effect):", length(newly_id_existing_ids), "\n") # 509
cat("  - new freeze-4 enrollees (batch 3, trivial):", length(new_enrollee_ids), "\n") #629

# f2-3
patient_code_type %>% filter(person_id %in% newly_id_existing_ids) %>%
    count(has_diagnosis, has_history) %>%
    mutate(pct = round(100 * n / sum(n), 1)) %>% print()
# has_diagnosis has_history     n   pct
# 1 FALSE         TRUE           65  17.5
# 2 TRUE          FALSE          21   5.6
# 3 TRUE          TRUE          286  76.9

# f4 only
patient_code_type %>% filter(person_id %in% new_enrollee_ids) %>%
    count(has_diagnosis, has_history) %>%
    mutate(pct = round(100 * n / sum(n), 1)) %>% print()
# has_diagnosis has_history     n   pct
# 1 FALSE         TRUE            2   1.1
# 2 TRUE          FALSE          24  13.3
# 3 TRUE          TRUE          155  85.6


# ========================
# CONTROL SELECTION
# ========================
# source(here("R", "control_selection_f3.R"))
# breast_controls_f3 <- select_controls_f3(
#     control_name = "simplexo",
#     exclude_codes = malig_neoplasms,
#     gender_filter = "Female",
#     age_filter = NULL,
#     crep_filter = FALSE,
#     pmbb_dir = here("PMBB"),
#     data_dir = here("simplexo", "data"),
#     log_dir = here("simplexo", "log")
# )
# f3c_df <- breast_controls_f3$excluded_icd_summary
# f3c_ids <- breast_controls_f3$final_controls$person_id

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
# f4c_df <- breast_controls_f4$excluded_icd_summary
# f4c_ids <- breast_controls_f4$final_controls$person_id

# length(f3c_ids)
# 18278
length(f4c_ids)
# 18720
