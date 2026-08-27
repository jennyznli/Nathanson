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
progeny_unmerged <- read_excel(here("simplexo", "ss", "br_pts_for_exwas_10022025.xlsx"))
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
dim(progeny_unmerged)
# 5046

progeny_case <- progeny_unmerged %>% filter(Gender == "F") %>%
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


# ========================================================================
# EARLIEST AGE BY CODE TYPE, within the merged ICD selection
# ========================================================================
# for everyone already in icd_ids (the combined-code selection), classify
# each of their individual matching occurrences as diagnosis-code vs
# history-code, then compute per-person earliest age separately by type:
#   - earliest age from a diagnosis code only
#   - earliest age from a personal-history code only
#   - earliest age overall (min across both -- same as CaDxAge_ICD)
# to see whether "age of dx" derived from a history code runs
# systematically younger/older than one derived from a true diagnosis code.
diagnosis_codes <- c("^C50", "^D05", "^174", "^233.0")  # ICD10 malignant/DCIS + ICD9 malignant/CIS
history_codes   <- c("^Z85.3", "^Z86.000", "^V10.3")    # ICD10 + ICD9 personal history of

occurrences <- breast_results_f4$matching_occurrences %>%
    filter(person_id %in% icd_ids) %>%
    mutate(
        is_diagnosis = grepl(paste(diagnosis_codes, collapse = "|"), source_value, perl = TRUE, ignore.case = TRUE),
        is_history   = grepl(paste(history_codes, collapse = "|"), source_value, perl = TRUE, ignore.case = TRUE)
    )

earliest_by_type <- occurrences %>%
    group_by(person_id) %>%
    summarise(
        earliest_diagnosis_date = if (any(is_diagnosis)) min(event_date[is_diagnosis], na.rm = TRUE) else as.Date(NA),
        earliest_history_date   = if (any(is_history))   min(event_date[is_history],   na.rm = TRUE) else as.Date(NA),
        earliest_overall_date   = min(event_date, na.rm = TRUE),
        .groups = "drop"
    ) %>%
    left_join(person %>% select(person_id, birth_datetime), by = "person_id") %>%
    mutate(
        Age_Diagnosis = as.numeric(difftime(earliest_diagnosis_date, birth_datetime, units = "days")) / 365.25,
        Age_History   = as.numeric(difftime(earliest_history_date,   birth_datetime, units = "days")) / 365.25,
        Age_Overall   = as.numeric(difftime(earliest_overall_date,   birth_datetime, units = "days")) / 365.25
    )

write.table(earliest_by_type, here("simplexo", "data", "simplexo4_earliest_age_by_codetype_df.txt"),
            row.names = FALSE, quote = FALSE, sep = "\t")

age_long <- earliest_by_type %>%
    select(person_id, Age_Diagnosis, Age_History, Age_Overall) %>%
    pivot_longer(cols = starts_with("Age_"), names_to = "Category", values_to = "Age") %>%
    filter(!is.na(Age)) %>%
    mutate(Category = recode(Category,
                             Age_Diagnosis = "Earliest Diagnosis Code",
                             Age_History   = "Earliest History Code",
                             Age_Overall   = "Earliest Overall"))

age_long %>%
    group_by(Category) %>%
    summarise(n = n(), mean_age = mean(Age), median_age = median(Age)) %>%
    print()

p_earliest_compare <- ggplot(age_long, aes(x = Age, fill = Category)) +
    geom_histogram(position = "identity", bins = 20, alpha = 0.4, color = "black") +
    scale_fill_manual(values = c(
        "Earliest Diagnosis Code" = "steelblue",
        "Earliest History Code"   = "pink",
        "Earliest Overall"        = "seagreen"
    )) +
    labs(title = "Earliest Age by ICD Code Type (within merged ICD selection)",
         x = "Age", y = "Count", fill = "Code Type") +
    theme_minimal()

ggsave(here("simplexo", "figures", "simplexo4_earliest_age_by_codetype_histogram.png"),
       plot = p_earliest_compare, width = 8, height = 4, dpi = 300)

# # ------------------------------------------------------------------------
# # diagnosis-only patients (never had a history code) vs. everyone with a
# # history code (overall) -- Age_Diagnosis vs Age_History distributions
# # ------------------------------------------------------------------------
# patient_code_flags <- occurrences %>%
#     group_by(person_id) %>%
#     summarise(
#         any_diagnosis = any(is_diagnosis),
#         any_history   = any(is_history),
#         .groups = "drop"
#     )
#
# diagnosis_only_ids <- patient_code_flags %>%
#     filter(any_diagnosis & !any_history) %>%
#     pull(person_id)
#
# cat("Patients with diagnosis codes only (never a history code):", length(diagnosis_only_ids), "\n")
# cat("As % of all ICD-identified patients:",
#     round(100 * length(diagnosis_only_ids) / length(icd_ids), 1), "%\n")
#
# history_pop_ids <- patient_code_flags %>%
#     filter(any_history) %>%
#     pull(person_id)
#
# cat("Patients with any history code:", length(history_pop_ids), "\n")
#
# diagnosis_only_vs_history_df <- bind_rows(
#     earliest_by_type %>%
#         filter(person_id %in% diagnosis_only_ids) %>%
#         transmute(person_id, Age = Age_Diagnosis, Category = "Diagnosis Only"),
#     earliest_by_type %>%
#         filter(person_id %in% history_pop_ids) %>%
#         transmute(person_id, Age = Age_History, Category = "History (Overall)")
# ) %>%
#     filter(!is.na(Age))
#
# diagnosis_only_vs_history_df %>%
#     group_by(Category) %>%
#     summarise(n = n(), mean_age = mean(Age), median_age = median(Age)) %>%
#     print()
#
# p_diag_only_vs_history <- ggplot(diagnosis_only_vs_history_df, aes(x = Age, fill = Category)) +
#     geom_histogram(position = "identity", bins = 20, alpha = 0.5, color = "black") +
#     scale_fill_manual(values = c(
#         "Diagnosis Only"     = "steelblue",
#         "History (Overall)"  = "pink"
#     )) +
#     labs(title = "Diagnosis-Only Patients vs. Overall History-Code Age",
#          x = "Age", y = "Count", fill = "Code Type") +
#     theme_minimal()
#
# ggsave(here("simplexo", "figures", "simplexo4_diagnosis_only_vs_history_histogram.png"),
#        plot = p_diag_only_vs_history, width = 8, height = 4, dpi = 300)

# ------------------------------------------------------------------------
# history-only patients (never had a diagnosis code) vs. everyone with a
# diagnosis code (overall) -- Age_History vs Age_Diagnosis distributions
# ------------------------------------------------------------------------
history_only_ids <- patient_code_flags %>%
    filter(any_history & !any_diagnosis) %>%
    pull(person_id)

cat("Patients with history codes only (never a diagnosis code):", length(history_only_ids), "\n")
cat("As % of all ICD-identified patients:",
    round(100 * length(history_only_ids) / length(icd_ids), 1), "%\n")
# 74, 1%

diagnosis_pop_ids <- patient_code_flags %>%
    filter(any_diagnosis) %>%
    pull(person_id)

cat("Patients with any diagnosis code:", length(diagnosis_pop_ids), "\n")
# 4508

history_only_vs_diagnosis_df <- bind_rows(
    earliest_by_type %>%
        filter(person_id %in% history_only_ids) %>%
        transmute(person_id, Age = Age_History, Category = "History Only"),
    earliest_by_type %>%
        filter(person_id %in% diagnosis_pop_ids) %>%
        transmute(person_id, Age = Age_Diagnosis, Category = "Diagnosis (Overall)")
) %>%
    filter(!is.na(Age))

history_only_vs_diagnosis_df %>%
    group_by(Category) %>%
    summarise(n = n(), mean_age = mean(Age), median_age = median(Age)) %>%
    print()
# Category                n mean_age median_age
# <chr>               <int>    <dbl>      <dbl>
# 1 Diagnosis (Overall)  4508     55.4       55.1
# 2 History Only           74     54.5       51.2

p_hist_only_vs_diagnosis <- ggplot(history_only_vs_diagnosis_df, aes(x = Age, fill = Category)) +
    geom_histogram(position = "identity", bins = 30, alpha = 0.5, color = "black") +
    scale_fill_manual(values = c(
        "History Only"          = "pink",
        "Diagnosis (Overall)"   = "steelblue"
    )) +
    labs(title = "History-Only Patients vs. Overall Diagnosis-Code Age",
         x = "Age", y = "Count", fill = "Code Type") +
    theme_minimal()

ggsave(here("simplexo", "figures", "simplexo4_history_only_vs_diagnosis_histogram.png"),
       plot = p_hist_only_vs_diagnosis, width = 8, height = 5, dpi = 300)

# ------------------------------------------------------------------------
# patients with BOTH a diagnosis code and a history code -- within this
# subgroup, compare their own Age_Diagnosis vs. Age_History (paired by
# person) to see whether the two code types disagree on age even when
# both are present for the same patient
# ------------------------------------------------------------------------
both_ids <- patient_code_flags %>%
    filter(any_diagnosis & any_history) %>%
    pull(person_id)

cat("Patients with BOTH a diagnosis code and a history code:", length(both_ids), "\n")
cat("As % of all ICD-identified patients:",
    round(100 * length(both_ids) / length(icd_ids), 1), "%\n")

both_ages_df <- earliest_by_type %>%
    filter(person_id %in% both_ids) %>%
    select(person_id, Age_Diagnosis, Age_History, Age_Overall) %>%
    mutate(
        Age_Diff = Age_History - Age_Diagnosis,
        # which code type actually supplied the final/overall (earliest) age
        Final_Age_Source = case_when(
            Age_Diagnosis < Age_History ~ "Diagnosis",
            Age_Diagnosis > Age_History ~ "History",
            TRUE ~ "Tie"
        )
    )

cat("Median Age_Diagnosis (both group):", median(both_ages_df$Age_Diagnosis, na.rm = TRUE), "\n")
# 54.9
cat("Median Age_History (both group):", median(both_ages_df$Age_History, na.rm = TRUE), "\n")
# 56.9
cat("Median Age_Overall (both group):", median(both_ages_df$Age_Overall, na.rm = TRUE), "\n")
# 54.7
cat("Median difference (History - Diagnosis):", median(both_ages_df$Age_Diff, na.rm = TRUE), "\n")
# 0.61

# how m
# how many/what proportion of the "both" group has their final (earliest
# overall) age coming from the history code vs. the diagnosis code
both_ages_df %>%
    count(Final_Age_Source) %>%
    mutate(pct = round(100 * n / sum(n), 1)) %>%
    print()
# Final_Age_Source     n   pct
# <chr>            <int> <dbl>
# 1 Diagnosis         3600  88.2
# 2 History            327   8
# 3 Tie                153   3.8

both_long <- both_ages_df %>%
    select(person_id, Age_Diagnosis, Age_History, Age_Overall) %>%
    pivot_longer(cols = c(Age_Diagnosis, Age_History, Age_Overall), names_to = "Category", values_to = "Age") %>%
    mutate(Category = recode(Category,
                             Age_Diagnosis = "Diagnosis Code",
                             Age_History   = "History Code",
                             Age_Overall   = "Final (Earliest of Both)"))

both_long %>%
    group_by(Category) %>%
    summarise(n = n(), mean_age = mean(Age), median_age = median(Age)) %>%
    print()

# Category                     n mean_age median_age
# <chr>                    <int>    <dbl>      <dbl>
# 1 Diagnosis Code            4080     55.2       55.0
# 2 Final (Earliest of Both)  4080     54.9       54.7
# 3 History Code              4080     56.9       56.9

# pairwise plots instead of one 3-way overlay (triple overlay was hard
# to distinguish) -- Diagnosis vs History, Diagnosis vs Final, History vs Final
both_colors <- c(
    "Diagnosis Code"            = "steelblue",
    "History Code"              = "pink",
    "Final (Earliest of Both)"  = "seagreen"
)

pairwise_hist <- function(df, cat1, cat2) {
    ggplot(df %>% filter(Category %in% c(cat1, cat2)), aes(x = Age, fill = Category)) +
        geom_histogram(position = "identity", bins = 20, alpha = 0.5, color = "black") +
        scale_fill_manual(values = both_colors) +
        labs(title = paste0(cat1, " vs. ", cat2),
             x = "Age", y = "Count", fill = "Code Type") +
        theme_minimal()
}

p_both_diag_vs_hist <- pairwise_hist(both_long, "Diagnosis Code", "History Code")
p_both_diag_vs_final <- pairwise_hist(both_long, "Diagnosis Code", "Final (Earliest of Both)")
p_both_hist_vs_final <- pairwise_hist(both_long, "History Code", "Final (Earliest of Both)")

ggsave(here("simplexo", "figures", "simplexo4_both_diag_vs_history_histogram.png"),
       plot = p_both_diag_vs_hist, width = 8, height = 5, dpi = 300)
ggsave(here("simplexo", "figures", "simplexo4_both_diag_vs_final_histogram.png"),
       plot = p_both_diag_vs_final, width = 8, height = 5, dpi = 300)
ggsave(here("simplexo", "figures", "simplexo4_both_history_vs_final_histogram.png"),
       plot = p_both_hist_vs_final, width = 8, height = 5, dpi = 300)

