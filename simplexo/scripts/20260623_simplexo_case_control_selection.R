# ========================
# SIMPLEXO F4
# CASE CONTROL SELECTION
# ========================
library(here)
setwd(here())
source(here("R/config.R"))

pmbb4 <- here("PMBB", "4.0")

# covariates
cov <- fread(file.path(pmbb4, "PMBB-Release-2026-4.0_phenotype_covariates.txt"), header = TRUE)
person <- fread(file.path(pmbb4, "PMBB-Release-2026-4.0_phenotype_person.txt"), header = TRUE)

# ICD
obs <- fread(file.path(pmbb4, "PMBB-Release-2026-4.0_phenotype_observation.txt"), header = TRUE)
cond <- fread(file.path(pmbb4, "PMBB-Release-2026-4.0_phenotype_condition_occurrence.txt"), header = TRUE)

# PMCR
pmcr_all <- fread(file.path(pmbb4, "PMBB-Release-2026-4.0_phenotype_cancer_pmcr.txt"), header = TRUE)
brca <- fread(file.path(pmbb4, "PMBB-Release-2026-4.0_phenotype_cancer_brca.txt"), header = TRUE)
hx <- fread(file.path(pmbb4, "PMBB-Release-2026-4.0_phenotype_family_hx.txt"), header = TRUE)

# progeny
progeny_unmerged <- read_excel(here("simplexo", "ss", "br_pts_for_exwas_10022025.xlsx"))

flags <- read.csv(here("PMBB", "3.0", "rgcname_pmbbid_metadata_flags.csv"))
up <- read.csv(here("simplexo", "data", "simplexo_up_map.csv"))

# id lists
exome_ids <- read.table(here("PMBB", "4.0", "PMBB-Release-2026-4.0_genetic_exome.sample_list.txt"), header = FALSE)$V1
imputed_ids <- read.csv(here("PMBB", "4.0", "PMBB-Release-2026-4.0_genetic_imputed.sample_list.txt"), header = FALSE)$V1

exome_no_pcs <- read.table(here("PMBB", "4.0", "PMBB-Release-2026-4.0_genetic_exome.no1KG.eigenvectors.tsv"), header = TRUE)
exome_with_pcs <- read.table(here("PMBB", "4.0", "PMBB-Release-2026-4.0_genetic_exome.with1KG.eigenvectors.tsv"), header = TRUE)
imp_no_pcs <- read.table(here("PMBB", "4.0", "PMBB-Release-2026-4.0_genetic_imputed.with1KG.eigenvectors.tsv"), header = TRUE)
imp_with_pcs <- read.table(here("PMBB", "4.0", "PMBB-Release-2026-4.0_genetic_imputed.no1KG.eigenvectors.tsv"), header = TRUE)

# ========================
# AVAILABLE SEQUENCING
# ========================
length(exome_ids)
# 70925
length(imputed_ids)
# 70493

seq_ids <- intersect(exome_ids, imputed_ids)
length(seq_ids)
# 70408

# ========================
# PROGENY
# ========================
# filter out benign tissues, restrict to females
progeny_clean <- progeny_unmerged %>%
    inner_join(up, by = "SampNum") %>%
    filter(Diagnosis != "Normal Benign Tissue" | is.na(Diagnosis)) %>%
    filter(Gender == "F") %>%
    rename(person_id = PMBB_ID)

# merge to individual level
progeny <- progeny_clean %>%
    merge_duplicates("SampNum")

# function to extract earliest age from semicolon-separated string
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

progeny$CaDxAge_Progeny <- sapply(progeny$CaDxAge, get_earliest_age)
progeny$CaDxAge_Progeny[is.infinite(progeny$CaDxAge_Progeny)] <- NA
progeny_ages <- progeny %>% select(person_id, CaDxAge_Progeny)

progeny_ids <- sort(unique(progeny$person_id))
length(progeny_ids)
# 1517

# ============================================================
# PMCR
# ============================================================
# restrict to breast and females
pmcr_unmerged <- pmcr_all %>%
    filter(SeerSiteRecode2023ExpandedGroup == "Breast") %>%
    left_join(brca, by = "person_id") %>%
    left_join(cov,  by = "person_id") %>%
    filter(sequenced_gender == "Female")

# calc age of dx
pmcr <- pmcr_unmerged %>%
    merge_duplicates("person_id") %>%
    left_join(person %>% select("person_id", "birth_datetime"), by = "person_id")

pmcr$CaDxAge_PMCR <- as.numeric(difftime(
    pmcr$FirstDxDate,
    pmcr$birth_datetime,
    units = "days"
)) / 365.25

pmcr_ages <- pmcr %>% select(person_id, CaDxAge_PMCR)
pmcr_ids <- sort(unique(pmcr$person_id))

length(pmcr_ids)
# 3090

# ========================
# ICD SELECTION
# ========================
source(here("R", "sample_selection_f4.R"))
breast_results_f4 <- select_samples_f4(
    sample_name = "simplexo4",
    icd_codes = c("^C50", "^Z85.3",            # ICD10 - malignant
                  "^D05", "^Z86.000",          # ICD10 - DCIS
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

# add birth date for age of dx calc
icd <- breast_results_f4$filtered_patients %>%
    left_join(person %>% select(person_id, birth_datetime), by = "person_id")
icd$CaDxAge_ICD <- as.numeric(difftime(
    icd$first_date,
    icd$birth_datetime,
    units = "days"
)) / 365.25
icd_ages <- icd %>% select(person_id, CaDxAge_ICD)

icd_ids <- sort(unique(breast_results_f4$filtered_patients$person_id))

length(icd_ids)
# 4582


# ============================================================
# INITIAL CASE LIST
# ============================================================
case_ids <- sort(unique(c(icd_ids, progeny_ids, pmcr_ids)))
length(case_ids)
# 5149

write.csv(case_ids, file.path(here("simplexo", "data", "simplexo4_case_ids_temp.csv")))

# ============================================================
# FIND AGE of DX
# ============================================================
case_seq_ids <- intersect(case_ids, seq_ids)
length(case_seq_ids)
# 5120

base_df <- data.frame(person_id = case_seq_ids)
age_df <- base_df %>%
    left_join(progeny_ages, by = "person_id") %>%
    left_join(pmcr_ages, by = "person_id") %>%
    left_join(icd_ages, by = "person_id") %>%
    mutate(
        # priority order: progeny, pmcr, ICD
        CaDxAge = coalesce(CaDxAge_Progeny, CaDxAge_PMCR, CaDxAge_ICD),
        # earliest diagnosis
        EarliestDiagnosis = pmin(CaDxAge_Progeny, CaDxAge_PMCR, CaDxAge_ICD, na.rm = TRUE),
        EarliestDiagnosis = if_else(is.infinite(EarliestDiagnosis), NA_real_, EarliestDiagnosis),
        AgeDxDiscrepancyFlag = abs(CaDxAge - EarliestDiagnosis) > 2
    )

age_df_sel <- age_df %>%
    select(person_id, CaDxAge) %>%
    filter(!is.na(CaDxAge)) %>%
    mutate(Age = round(CaDxAge, 2))

dim(age_df_sel)
# 5117

final_ids <- sort(unique(age_df_sel$person_id))
final_df <- cov %>% filter(person_id %in% final_ids) %>%
    left_join(age_df_sel, by = "person_id")

write.table(age_df_sel, here("simplexo", "data", "simplexo4_overall_case_age_df.txt"),
            row.names = FALSE, quote = FALSE, sep = "\t")

write.table(final_df, here("simplexo", "data", "simplexo4_overall_cov_df.txt"),
            row.names = FALSE, quote = FALSE, sep = "\t")

write.table(final_ids, here("simplexo", "data", "simplexo4_overall_case_ids.txt"),
            row.names = FALSE, col.names = FALSE, quote = FALSE)

# ========================
# QUICK CHECK
# ========================
table(final_df$batch)
# 1    2    3
# 2095 2393  629

table(final_df$crep_highrisk_flag)
# 0    1
# 3549 1568

summary(final_df$Age)
sd(final_df$Age, na.rm = TRUE)

table(final_df$AgeDxDiscrepancyFlag, useNA = "ifany")

# ========================
# MALIGNANT NEOPLASM CODES
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
source(here("R", "control_selection_f4.R"))
dim(cov)
# 70925

breast_controls_f4 <- select_controls_f4(
    control_name = "simplexo4",
    exclude_codes = malig_neoplasms,
    gender_filter = "Female",
    age_filter = NULL,
    crep_filter = FALSE,
    pmbb_dir = here("PMBB", "4.0"),
    data_dir = here("simplexo", "data"),
    log_dir = here("simplexo", "log")
)

controls <- breast_controls_f4$final_controls %>%
    filter(!(person_id %in% case_ids)) %>%
    filter(!(person_id %in% pmcr_all$person_id)) %>%  # should not be in any PMCR since they are all cancer
    filter(!(person_id %in% progeny$person_id)) %>%   # should not be in progeny but may be redundnt
    filter(person_id %in% seq_ids) %>%
    select(person_id, sample_age, sequenced_gender)
dim(controls)
# 18209

control_ids <- sort(unique(controls$person_id))

write.table(controls, here("simplexo", "data", "simplexo4_overall_control_df.txt"),
            quote = FALSE, row.names = FALSE, col.names = FALSE)

# ========================
# FINAL COVARIATEs
# ========================
all_ids <- unique(c(control_ids, final_ids))

pmbb_ages <- cov %>%
    select(person_id, sample_age)

case_ages <- age_df_sel %>%
    select(person_id, CaDxAge)

# merge
all_ages <- pmbb_ages %>%
    filter(person_id %in% all_ids) %>%
    left_join(case_ages, by = "person_id") %>%
    # select dx age, then sample age
    mutate(Age = coalesce(CaDxAge, sample_age)) %>%
    select(person_id, Age)
dim(all_ages)
# 23326

all_cov <- cov %>% filter(person_id %in% all_ids) %>%
    left_join(all_ages, by = "person_id") %>%
    mutate(Group = if_else(person_id %in% final_ids, "1", "0")) %>%
    select(person_id, batch, Group, Age) %>%
    filter(!is.na(Age))
case_control_ids <- sort(all_cov$person_id)
length(case_control_ids)
# 23323

table(all_cov$Group)
final_case_control_ids <- unique(all_cov$person_id)

## FINAL IDS
final_case_ids <- sort(intersect(final_ids, case_control_ids))
final_control_ids <- sort(intersect(control_ids, case_control_ids))
length(final_case_ids)
length(final_control_ids)
# 5117
# 18206

## PIDs
case_control_pids <- as.data.frame(case_control_ids)
case_control_pids$x <- case_control_pids$case_control_ids

final_case_pids <- as.data.frame(final_case_ids)
final_case_pids$x <- final_case_pids$case_ids

final_control_pids <- as.data.frame(final_control_ids)
final_control_pids$x <- final_control_pids$control_ids

write.table(all_cov, here("simplexo", "cov", "simplexo4_overall_case_control_cov_df.txt"),
            row.names = FALSE, quote = FALSE, sep = "\t")

write.table(case_control_ids, here("simplexo", "cov", "simplexo4_overall_case_control_ids.txt"),
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")

write.table(case_control_pids, here("simplexo", "cov", "simplexo4_overall_case_control_pids.txt"),
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")

write.table(final_case_ids, here("simplexo", "cov", "simplexo4_overall_case_ids.txt"),
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")

write.table(final_case_pids, here("simplexo", "cov", "simplexo4_overall_case_pids.txt"),
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")

write.table(final_control_ids, here("simplexo", "cov", "simplexo4_overall_control_ids.txt"),
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")

write.table(final_control_pids, here("simplexo", "cov", "simplexo4_overall_control_pids.txt"),
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")

write.table(all_cov, here("simplexo", "cov", "simplexo4_overall_case_control_age_df.txt"),
            row.names = FALSE, quote = FALSE, sep = "\t")

### sex files
sex_file <- cov %>%
    mutate(
        FID = person_id,
        IID = person_id,
        SEX = case_when(
            sequenced_gender == "Male"   ~ 1,
            sequenced_gender == "Female" ~ 2,
            TRUE ~ NA_real_
        )
    ) %>%
    select(FID, IID, SEX)

write.table(sex_file, file.path(here("simplexo", "cov", "sex.txt")),
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

# ========================
# COVARIATE FILES
# ========================
### exome no 1kg
exome_no_cov <- all_cov %>%
    mutate(
        FID = person_id,
        IID = person_id
    ) %>%
    select(-Group) %>%
    left_join(exome_no_pcs, by = c("FID" = "person_id")) %>%
    select(FID, IID, Age, batch, num_range("exome_PC", 1:6))

write.table(exome_no_cov, file.path(here("simplexo", "cov", "simplexo4_all_no1kg_exome6_covariates.txt")),
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

### update file
update <- sex_file %>% mutate(
    a = FID,
    b = FID
) %>% select(FID, IID, a, b)
write.table(update, file.path(here("simplexo", "cov", "array_update_ids.txt")),
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)


### exome 1kg
exome_with_cov <- all_cov %>%
    mutate(
        FID = person_id,
        IID = person_id
    ) %>%
    select(-Group) %>%
    left_join(exome_with_pcs, by = c("FID" = "person_id")) %>%
    select(FID, IID, Age, batch, num_range("exome_PC", 1:6))

write.table(exome_with_cov, file.path(here("simplexo", "cov", "simplexo4_all_with1kg_exome6_covariates.txt")),
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

### imputed no 1kg
imp_no_cov <- all_cov %>%
    mutate(
        FID = person_id,
        IID = person_id
    ) %>%
    select(-Group) %>%
    left_join(imp_no_pcs, by = c("FID" = "person_id")) %>%
    select(FID, IID, Age, batch, num_range("imputed_PC", 1:6))

write.table(exome_no_cov, file.path(here("simplexo", "cov", "simplexo4_all_no1kg_imp6_covariates.txt")),
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

### imputed no 1kg
imp_no_cov <- all_cov %>%
    mutate(
        FID = person_id,
        IID = person_id
    ) %>%
    select(-Group) %>%
    left_join(imp_with_pcs, by = c("FID" = "person_id")) %>%
    select(FID, IID, Age, batch, num_range("imputed_PC", 1:6))

write.table(exome_no_cov, file.path(here("simplexo", "cov", "simplexo4_all_no1kg_imp6_covariates.txt")),
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)


# ========================
# DIAGNOSTIC PLOTS
# ========================
all_cov$batch <- as.factor(all_cov$batch)
all_cov$Group <- as.factor(all_cov$Group)

### stacked bar plots
group_totals <- all_cov %>% count(Group)
group_totals

p_stack <- ggplot(all_cov, aes(x = Group, fill = batch)) +
    geom_bar(position = "stack", color = "darkgray") +
    geom_text(stat = "count", aes(label = after_stat(count)), position = position_stack(vjust = 0.5), size = 3) +
    geom_text(data = group_totals, aes(x = Group, y = n, label = n), inherit.aes = FALSE, vjust = -0.5) +
    scale_fill_brewer(type = "qual", palette = "Set3") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    labs(title = paste("Case vs Control - Batch"), x = "Group", y = "Count", fill = all_cov$batch) +
    theme_minimal()
ggsave(here("simplexo", "figures", paste0("simplexo4_stacked_bar_batch.png")),
       plot = p_stack, width = 4, height = 4, dpi = 300)

### histogram by age
p_age_overlay <- ggplot(all_cov, aes(x = Age, fill = Group)) +
    geom_histogram(bins = 30, alpha = 0.6, position = "identity", color = "darkgray") +
    scale_fill_manual(values = c("1" = "purple", "0" = "darkgray")) +
    labs(title = "Age Distribution - Case vs Control", x = "Age", y = "Count") +
    theme_minimal()
ggsave(here("simplexo", "figures", "simplexo4_age_histogram_overlay.png"),
       plot = p_age_overlay, width = 6, height = 4, dpi = 300)

# ========================================================================
# FAMILY HISTORY
# ========================================================================
### PROGENY
progeny_fh <- progeny[, c("person_id", "SampNum", "Gender", "NumFDR", "NumFDR_BC", "NumSDR", "NumSDR_BC")]

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

progeny_pmbb_fh <- progeny_fh %>%
    filter(Family_History == 1) %>%
    filter(person_id %in% case_seq_ids)

progeny_fh_ids <- unique(progeny_pmbb_fh$person_id)
length(progeny_fh_ids)
# 1096

### PMCR
breast_hx <- hx %>%
    filter(grepl("breast", MEDICAL_HX, ignore.case = TRUE)) %>%
    filter(!(RELATION %in% c("null", "", "Negative History"))) %>%
    left_join(cov %>% select(person_id, sequenced_gender), by = "person_id") %>%
    filter(sequenced_gender == "Female")

pmcr_fh_ids <- unique(breast_hx$person_id)
length(pmcr_fh_ids)
# 11010

### ICD CODES
family_history <- select_samples_f4(
    sample_name = "simplexo4_fh",
    icd_codes = c("^Z80.3", "^V16.3"),
    gender_filter = "Female",
    crep_filter = NULL,
    min_instances = 1,
    min_timespan = NULL,
    age_filter = NULL,
    exclude = FALSE,
    pmbb_dir = here("PMBB"),
    data_dir = here("simplexo", "data"),
    log_dir = here("simplexo", "log")
)

icd_fh_ids <- unique(family_history$filtered_patients$person_id)
length(icd_fh_ids)
# 4357

### COMBINE
all_fh_ids <- unique(c(progeny_fh_ids, pmcr_fh_ids, icd_fh_ids))
all_fh_ids <- sort(intersect(final_ids, all_fh_ids))
length(all_fh_ids)
# 3223

write.table(all_fh_ids, here("simplexo", "data", "simplexo4_fhx_case_ids.txt"),
            quote = FALSE, row.names = FALSE, col.names = FALSE)

# ========================================================================
# ER STATUS
# ========================================================================
### PROGENY
progeny_er <- progeny_clean %>%
    filter(ERstatus %in% c("Positive", "Negative")) %>%
    merge_duplicates("person_id")
dim(progeny_er)
# 918

progeny_erp <- progeny_er %>% filter(ERstatus == "Positive")
progeny_ern <- progeny_er %>% filter(ERstatus == "Negative")
cat("Progeny ER+:", nrow(progeny_erp), "\n")
# 593
cat("Progeny ER-:", nrow(progeny_ern), "\n")
# 304

### PMCR
pmcr_er <- pmcr_unmerged %>%
    filter(EstrogenReceptorSummary %in% c("Positive", "Negative")) %>%
    merge_duplicates("person_id")
cat("After merging duplicates:", nrow(pmcr_er), "\n")
# 2727

pmcr_erp <- pmcr_er %>% filter(EstrogenReceptorSummary == "Positive")
pmcr_ern <- pmcr_er %>% filter(EstrogenReceptorSummary == "Negative")
cat("PMCR ER+:", nrow(pmcr_erp), "\n")
# 2076
cat("PMCR ER-:", nrow(pmcr_ern), "\n")
# 620

###### ICD CODES ######
# nothing...

###  COMBINE
all_erp_ids <- unique(c(pmcr_erp$person_id, progeny_erp$person_id))
all_ern_ids <- unique(c(pmcr_ern$person_id, progeny_ern$person_id))

# conflicts between sources
er_conflicts <- intersect(all_erp_ids, all_ern_ids)
length(er_conflicts)
# 16

# remove conflicts and make sure they have sequencing
all_erp_ids <- sort(intersect(setdiff(all_erp_ids, er_conflicts), final_ids))
all_ern_ids <- sort(intersect(setdiff(all_ern_ids, er_conflicts), final_ids))
cat("ER+ cases:", length(all_erp_ids), "\n") #2290
cat("ER- cases:", length(all_ern_ids), "\n") #789

write.table(all_erp_ids, here("simplexo", "data", "simplexo4_er_pos_case_ids.txt"),
            quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(all_ern_ids, here("simplexo", "data", "simplexo4_er_neg_case_ids.txt"),
            quote = FALSE, row.names = FALSE, col.names = FALSE)

# ========================================================================
# INVASIVE ONLY (NO IN SITU)
# ========================================================================
### ICD CODES
inv_only <- select_samples_f4(
    sample_name = "simplexo4_malig",
    icd_codes = c("^C50", "^Z85.4", # ICD10
                  "^174", "^V10.3"), # ICD9
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

icd_malig_ids <- unique(inv_only$filtered_patients$person_id)
length(icd_malig_ids)
# 4598

###  PROGENY
progeny_clean <- progeny_clean %>%
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

progeny_malig <- progeny_clean %>%
    filter(!is_insitu) %>%
    merge_duplicates("SampNum")

progeny_malig_ids <- unique(progeny_malig$person_id)
length(progeny_malig_ids)
# 1340

### PMCR
pmcr_malig <- pmcr_unmerged %>%
    filter(!(grepl("situ", HistologicType, ignore.case = TRUE))) %>%
    merge_duplicates("person_id")

pmcr_malig_ids <- unique(pmcr_malig$person_id)
length(pmcr_malig_ids)
# 2968

### COMBINE
all_malig_ids <- sort(intersect(unique(c(icd_malig_ids, progeny_malig_ids, pmcr_malig_ids)), final_ids))
length(all_malig_ids)
# 4783

write.table(all_malig_ids, here("simplexo", "data", "simplexo4_malig_case_ids.txt"),
            quote = FALSE, row.names = FALSE, col.names = FALSE)

# ========================================================================
# UNDER 50
# ========================================================================
young <- age_df_sel %>% filter(Age <= 50)
young_ids <- unique(young$person_id)
length(young_ids)
# 2200

write.table(young_ids, here("simplexo", "data", "simplexo4_50_case_ids.txt"),
            quote = FALSE, row.names = FALSE, col.names = FALSE)


# ========================================================================
# MULTIPLE PHENO FILE
# ========================================================================
make_binary_pheno <- function(person_id, case_ids, control_ids) {
    case_when(
        person_id %in% case_ids    ~ 1,
        person_id %in% control_ids ~ 0,
        TRUE                       ~ NA_real_
    )
}

pheno <- all_cov %>%
    mutate(
        FID            = person_id,
        IID            = person_id,
        Overall        = make_binary_pheno(person_id, final_ids,     control_ids),
        Under50        = make_binary_pheno(person_id, young_ids,     control_ids),
        FamilyHistory  = make_binary_pheno(person_id, all_fh_ids,    control_ids),
        Malignant      = make_binary_pheno(person_id, all_malig_ids, control_ids),
        ERPos          = make_binary_pheno(person_id, all_erp_ids,   control_ids),
        ERNeg          = make_binary_pheno(person_id, all_ern_ids,   control_ids)
    ) %>%
    # FID/IID + phenotype columns only -- preprocess_regenie.py --pheno-file treats every
    # non-FID/IID column as a phenotype to test, so batch/Age/Group/person_id must not ride along
    select(FID, IID, Overall, Under50, FamilyHistory, Malignant, ERPos, ERNeg)

(colSums(is.na(pheno)))
# FID           IID       Overall       Under50 FamilyHistory     Malignant         ERPos         ERNeg
# 0             0             0          2917          1894           334          2827          4328


write.table(pheno, here("simplexo", "data", "simplexo4_pheno.txt"),
            quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")

write_xlsx(pheno, here("simplexo", "data", "simplexo4_pheno.xlsx"))

