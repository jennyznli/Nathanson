# ========================
# SIMPLEXO F4
# CASE CONTROL SELECTION
# ========================
library(here)
setwd(here("tgct"))
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
hx <- fread(file.path(pmbb4, "PMBB-Release-2026-4.0_phenotype_family_hx.txt"), header = TRUE)

flags <- read.csv(here("PMBB", "3.0", "rgcname_pmbbid_metadata_flags.csv"))
up <- read.csv(here("simplexo", "data", "simplexo_up_map.csv"))

# id lists
exome_ids <- read.table(here("PMBB", "4.0", "PMBB-Release-2026-4.0_genetic_exome.sample_list.txt"), header = FALSE)$V1
imputed_ids <- read.csv(here("PMBB", "4.0", "PMBB-Release-2026-4.0_genetic_imputed.sample_list.txt"), header = FALSE)$V1
genotype_ids <- read.csv(here("PMBB", "4.0", "PMBB-Release-2026-4.0_genetic_genotype.sample_list.txt"), header = FALSE)$V1

# others
brad <- fread(here("tgct", "data", "geno_map_20250828.csv"), header = TRUE)
mrn <- read_excel(here("tgct", "data", "KT_IDS_MRNs_UIDs.xlsx"))

x <- read_tsv(here("tgct", "data", "tgct_final_controls.txt"))
# DiscoveryID     GencoveID       PMBBID2 ReplicationID
# 65            81           921            66

# mixed models


# ========================
# AVAILABLE SEQUENCING
# ========================
length(exome_ids)
# 70925
length(imputed_ids)
# 70493
length(genotype_ids)
# 27535 in most recent...

seq_ids <- intersect(exome_ids, imputed_ids)
length(seq_ids)
# 70408

# ========================
# TGCT Cancer Selection
# ========================
library(here)
source(here("R", "sample_selection_f4.R"))

tgct_results <- select_samples_f4(
    sample_name = "tgct4",
    icd_codes = c("^C62", "^Z85.47",  # ICD10
                  "^186", "^V10.47"), # ICD9
    gender_filter = "Male",
    crep_filter = NULL,
    min_instances = 2,
    min_timespan = NULL,
    exclude = FALSE,
    age_filter = NULL,
    pmbb_dir = here("PMBB"),
    data_dir = here("tgct", "data"),
    log_dir = here("tgct", "log")
)

tgct_res <- tgct_results$filtered_patients
dim(tgct_res)
# 1   2   3
# 170 646  66


### should chekc not in PMCR

# ========================
# ...
# ========================
# 3 KT ids in the f4


# ========================
# Non Cancer selection
# ========================
source(here("R", "control_selection_f4.R"))

malig_neoplasms <- c(
    "^C(?!44)",           # Malignant neoplasms except skin cancer (C44)
    "^Z85(?!\\.828)",     # Personal history of malignant neoplasm except skin
    "^D05",               # Carcinoma in situ of breast
    "^Z86.000",           # Personal history of in-situ neoplasm of breast
    "^(?!173)(?:1[4-9][0-9]|20[0-8]|209\\.[0-3])",  # ICD-9 malignant neoplasms except 173
    "^233.0",             # ICD-9 carcinoma in situ of breast
    "^V10(?!\\.83)"       # ICD-9 personal history of malignant neoplasm except skin
)

malig_neoplasms <- c("^C(?!44)", "^Z85(?!\\.828)", "^(?!173)(?:1[4-9][0-9]|20[0-8]|209\\.[0-3])", "^V10(?!\\.83)")
malign_benign_neoplasms <- c("^C(?!44)", "^D", "^(?!173)(?:1[4-9][0-9]|2[0-3][0-9])", "^Z85", "^V10(?!\\.83)")

tgct_controls <- select_controls(
    control_name = "tgct4",
    exclude_codes = malig_neoplasms,
    gender_filter = "Male",
    age_filter = NULL,
    crep_filter = FALSE,
    pmbb_dir = here("PMBB", "3.0"),
    data_dir = here("tgct", "data"),
    log_dir = here("tgct", "log")
)

# ========================
# MATCH PMBB TO KT IDS
# ========================
kt_pattern <- "^UPENN-PMBB_KT[0-9]+_KT[0-9]+$"

kt <- flags %>%
    mutate(match_col = apply(select(., RGC_sample_name, ID1, ID2, ID3, ID4), 1, function(row) {
        m <- row[grepl(kt_pattern, row)]
        if (length(m) > 0) m[1] else NA
    })) %>%
    filter(!is.na(match_col))
print(dim(kt))
# 904

kt <- select(kt, c("PMBB_ID",  "RGC_sample_name", "ID1", "ID2", "ID3", "ID4", "CREP", "match_col"))
kt$ID <- sub(".*(KT[0-9]{7}).*", "\\1", kt$match_col)

kt <- select(kt, c("PMBB_ID",  "match_col", "ID"))
# kt$SampNum <- numbers <- as.numeric(gsub("KT", "", kt$ID))
kt$SampNum <- as.numeric(gsub("[a-zA-Z]", "", kt$ID))

## add in MRN
dim(mrn) #1165    4

# clean out 0 from MRN col
mrn$`2_hup_mrn`[mrn$`2_hup_mrn` == 0] <- NA
mrn <- mrn %>%
    filter(!is.na(`2_ktid`), !is.na(`2_hup_mrn`))
dim(mrn) #739   3
length(unique(mrn$`2_ktid`))

# take out KT and take out numeric
mrn$SampNum <- as.numeric(gsub("[a-zA-Z]", "", mrn$`2_ktid`))

# KTR IDs
# which(is.na(mrn$SampNum)) #249 352 442
# mrn[249,]
# mrn[352,]

#  1165    3
# colnames(mrn) <- c("2_ktid", "MRN")
kt_mrn <- left_join(kt, mrn, by = "SampNum")

write.csv(kt_mrn, here("PMBB", "3.0", "pmbb3_kt_key.csv"))

# ========================
# COMPARE W/ ICD CODES
# ========================
icd_df <- fread(here("tgct", "data", "tgct4_filtered_patients.txt"), header = TRUE)
dim(icd_df)

icds <- fread(here("tgct", "data", "tgct4_filtered_patients_ids.txt"), header = FALSE)$V1
# 882

kt_not_icd <- setdiff(kt_mrn$PMBB_ID, icds)
length(kt_not_icd)
# 207 kt ids that are not detected by ICD, need to be checked

kt_icd <- intersect(kt_mrn$PMBB_ID, icds)
length(kt_icd)
# 697 kt ids that are identified by icd codes

icd_not_kt <- setdiff(icds, kt_mrn$PMBB_ID)
length(icd_not_kt)
# 185 icd codes not KT ides

# left join the kt not icd codes w/ their kt ids
kt_not_icd_df <- kt_mrn %>% filter(PMBB_ID %in% kt_not_icd)
dim(kt_not_icd_df)
# 207

# join all 779 samples
# or maybe merge the kt_mrn with the icd_df
x <- left_join(icd_df[,1], kt_mrn, by = c("person_id" = "PMBB_ID"))
dim(x)
write.csv(x, here("tgct", "data", "tgct_icd_mrn_df.csv"), row.names = FALSE, quote = FALSE)

x <- kt_mrn %>% filter(PMBB_ID %in% icd_df$person_id)
dim(icd_df)

# sum(is.na(kt_not_icd_df$`2_hup_mrn`))
kt_icd_df <- kt_mrn %>% filter(PMBB_ID %in% kt_icd)
kt_icd_df <- left_join(
    kt_icd_df,
    icd_df,
    by = c("PMBB_ID" = "person_id")
)

# save results!
write.table(kt_not_icd, here("tgct", "data", "tgct_kt_not_icd.txt"),
            col.names = FALSE, row.names = FALSE, quote = FALSE)
write.csv(kt_not_icd_df, here("tgct", "data", "tgct_kt_not_icd_df.csv"), row.names = FALSE, quote = FALSE)

write.table(kt_icd, here("tgct", "data", "tgct_kt_icd.txt"),
            col.names = FALSE, row.names = FALSE, quote = FALSE)
write.csv(kt_icd_df, here("tgct", "data", "tgct_kt_icd_df.csv"), row.names = FALSE, quote = FALSE)




