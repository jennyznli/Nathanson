# ========================
# IDENTIFYING CANCER & AGE OF DX
# ========================
library(here)
setwd("/Users/jennyzli/Documents/Nathanson")
source(here("R", "config.R"))

# ========================
# LOAD DATA
# ========================
### BRCA1/2 carrier pull ###
# master
progeny <- read_excel(here("ch", "ss", "brca_carriers_ch_freq_w_seen_in_crep_20251020.xlsx"), sheet = "Data_from_master_table")
progeny <- progeny %>% filter(DNA == "D")
progeny$SampNum <- as.numeric(progeny$SampNum)

# tumor
progeny_tumor <- read_excel(here("ch", "ss", "brca_carriers_ch_freq_w_seen_in_crep_20251020.xlsx"), sheet = "Tumor")

# ICD codes
icd_hem <- read_excel(here("ch", "data", "icd_hematologic.xlsx"))

# key
up <- read.csv(here("simplexo", "data", "simplexo_up_map.csv"))

### breast cancer pull ###
progeny_breast <- read_excel(here("simplexo", "ss", "br_pts_for_exwas_10022025.xlsx")) %>% filter(DNA == "D")

### PMBB ###
PMBB_DIR = here("PMBB", "3.0")
COV <- here(PMBB_DIR, "PMBB-Release-2024-3.0_covariates.txt")
PER <- here(PMBB_DIR, "PMBB-Release-2024-3.0_phenotype_person.txt")
cov <- fread(COV, header = TRUE)
person <- fread(PER, header = TRUE)

# ======================
# SELECT CANCER
# ======================
malig_insitu <- c(
    "^C(?!44)",           # Malignant neoplasms except skin cancer (C44)
    "^Z85(?!\\.828)",     # Personal history of malignant neoplasm except skin
    "^D0(?!4)[0-9]",      # Benign/in situ neoplasms except skin cancer (DO4)
    "^Z86.00",            # Personal history of in-situ neoplasm
    "^(?!173)(?:1[4-9][0-9]|2[0-3][0-9])",  # ICD-9 malignant/insitu neoplasms except 173
    "^V10(?!\\.83)"       # ICD-9 personal history of malignant neoplasm except skin
)

hem <- paste0("^", icd_hem$code)
cancer_codes <- c(malig_insitu, hem)
# 69

cancers <- select_samples(
    sample_name = "hem_solid",
    icd_codes = cancer_codes,
    min_instances = 1,
    exclude = FALSE,
    data_dir = here("ch", "data"),
    pmbb_dir = here("PMBB"),
    log_dir = here("ch", "log")
)

cancer_ids <- sort(unique(cancers$filtered_patients$person_id))
length(cancer_ids)
# 29051

write.csv(cancers$filtered_patients, here("ch", "data", "pmbb_hem_solid_cancer_df.csv"))
write.csv(cancer_ids, here("ch", "data", "pmbb_hem_solid_cancer_ids.csv"))

# ========================
# PROGENY - BRCA12
# ========================
### TUMOR ###
# collapse to individual
progeny_tumor_merged <- merge_duplicates(progeny_tumor, "globalid")

# get earliest dx age
progeny_tumor_merged$CaDxAge_Progeny <- sapply(progeny_tumor_merged$CaDxAge, get_earliest_age)

# get earliest dx date
progeny_tumor_merged$CaDxDate_Progeny <- sapply(
    progeny_tumor_merged$CaDxDate,
    get_earliest_date,
    USE.NAMES = FALSE
)
progeny_tumor_merged$CaDxDate_Progeny <- as.Date(
    progeny_tumor_merged$CaDxDate_Progeny,
    origin = "1970-01-01"
)

### MASTER ###
progeny_merged <- merge_duplicates(progeny, "SampNum")

# filter down to PMBB only
progeny_pmbb <- progeny_merged %>%
    left_join(up, by = "SampNum") %>%
    filter(!is.na(PMBB_ID))
length(unique(progeny_pmbb$PMBB_ID))
# 1278

# merge tumor and master
progeny_pmbb$globalid <- as.numeric(progeny_pmbb$globalid)
progeny_pmbb_tumor <- inner_join(progeny_pmbb, progeny_tumor_merged, by = "globalid")

progeny_b12_ids <- sort(unique(progeny_pmbb_tumor$PMBB_ID))
length(progeny_b12_ids)
# 724

# ========================
# PROGENY - BREAST CANCER
# ========================
progeny_breast_merged <- merge_duplicates(progeny_breast, "SampNum")

progeny_breast_merged$CaDxAge_Progeny <- sapply(progeny_breast_merged$CaDxAge, get_earliest_age)

progeny_breast_merged$CaDxDate_Progeny <- sapply(
    progeny_breast_merged$CaDxDate,
    get_earliest_date,
    USE.NAMES = FALSE
)
progeny_breast_merged$CaDxDate_Progeny <- as.Date(
    progeny_breast_merged$CaDxDate_Progeny,
    origin = "1970-01-01"
)

progeny_breast_pmbb <- progeny_breast_merged %>%
    left_join(up, by = "SampNum") %>%
    filter(!is.na(PMBB_ID))

progeny_breast_pmbb_ids <- sort(progeny_breast_pmbb$PMBB_ID)
length(progeny_breast_pmbb_ids)
# 1561

# ========================
# MERGE ALL SOURCES
# ========================
final_cancer_ids <- sort(unique(c(cancer_ids, progeny_breast_pmbb_ids, progeny_b12_ids)))
length(final_cancer_ids)
# 29606

write.csv(final_cancer_ids, here("ch", "data", "pmbb_ch_cancer_ids.csv"), row.names = FALSE)

# ========================
# COMPILE DX AGES
# ========================
###  ICD CODES ###
dx_date <- cancers$all_patients %>% dplyr::select("person_id", "first_date")

# infer age of cancer diagnosis with first date of diagnosis and birth date
dx_age <- dx_date %>%
    left_join(person[, c("birth_datetime", "person_id")], by = "person_id") %>%
    mutate(Dx_Age = as.numeric(difftime(first_date, birth_datetime, units = "days")) / 365.25)
colnames(dx_age) <- c("PMBB_ID", "CaDxDate_ICD", "Birth_Date", "CaDxAge_ICD")

### PROGENY - BREAST ###
simplexo_age <- read.table(here("simplexo", "data", "simplexo_overall_case_ages_merged.txt"), header = TRUE)
dim(simplexo_age)
# 4063
colnames(simplexo_age)

### PROGENY - BRCA12 ###
tumor_age <- progeny_pmbb_tumor %>% dplyr::select(PMBB_ID, CaDxDate_Progeny, CaDxAge_Progeny)
dim(tumor_age)
# 724

# create merged dataframe
df <- data.frame(PMBB_ID = final_cancer_ids)

# diagnosis age based on priority: Progeny breast > Progeny BRCA1/2 > ICD ages
df$CaDxAge <- coalesce(
    simplexo_age$Age[match(df$PMBB_ID, simplexo_age$PMBB_ID)],
    tumor_age$CaDxAge_Progeny[match(df$PMBB_ID, tumor_age$PMBB_ID)],
    dx_age$CaDxAge_ICD[match(df$PMBB_ID, dx_age$PMBB_ID)]
)

# ========================
# ADD SAMPLE AGES
# ========================
df$Sample_Age <- cov$Sample_age[match(df$PMBB_ID, cov$person_id)]

df <- df %>% filter(!is.na(Sample_Age), !is.na(CaDxAge))
dim(df)
# 29590

write.csv(df, here("ch", "data", "pmbb_cancer_age.csv"))

