# ========================
# TGCT CLINICAL DATA
# ========================

library(here)
source(here("R", "load_packages.R"))
setwd(here("tgct"))

# ========================
# LOAD DATA
# ========================
OCC <- here("PMBB", "3.0", "PMBB-Release-2024-3.0_phenotype_condition_occurrence.txt")
COV <- here("PMBB", "3.0", "PMBB-Release-2024-3.0_covariates.txt")
PER <- here("PMBB", "3.0", "PMBB-Release-2024-3.0_phenotype_person.txt")

cov <- fread(COV, header = TRUE)
person <- fread(PER, header = TRUE)
flags <- fread(here("PMBB", "3.0", "rgcname_pmbbid_metadata_flags.csv"), header = TRUE)
icds <- readLines((here("tgct", "data", "tgct_cancer_filtered_patients_ids.txt")))
brad <- fread(here("tgct", "data", "geno_map_20250828.csv"), header = TRUE)
mrn <- read_excel(here("tgct", "data", "KT_IDS_MRNs_UIDs.xlsx"))

# ========================
# TGCT Cancer Selection
# ========================
library(here)
source(here("R", "sample_selection.R"))

tgct_results <- select_samples(
    sample_name = "tgct",
    icd_codes = c("^C62", "^Z80.43", "^186", "^V10.47"),
    gender_filter = "Male",
    crep_filter = NULL,
    min_instances = 2,
    min_timespan = NULL,
    exclude = FALSE,
    age_filter = NULL,
    output_prefix = "tgct",
    pmbb_dir = here("PMBB"),
    data_dir = here("tgct", "data"),
    log_dir = here("tgct", "log")
)

# ========================
# Non Cancer selection
# ========================

library(here)
source(here("R", "control_selection.R"))

malig_neoplasms <- c("^C(?!44)", "^Z85(?!\\.828)", "^(?!173)(?:1[4-9][0-9]|20[0-8]|209\\.[0-3])", "^V10(?!\\.83)")
malign_benign_neoplasms <- c("^C(?!44)", "^D", "^(?!173)(?:1[4-9][0-9]|2[0-3][0-9])", "^Z85", "^V10(?!\\.83)")

tgct_controls <- select_controls(
    control_name = "tgct",
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
#
icd_df <- fread(here("tgct", "data", "tgct_cancer_filtered_patients.txt"), header = TRUE)

length(icds) # 779
# length(unique(icds))

kt_not_icd <- setdiff(kt_mrn$PMBB_ID, icds)
length(kt_not_icd)
# 226 kt ids that are not detected by ICD, need to be checked

kt_icd <- intersect(kt_mrn$PMBB_ID, icds)
length(kt_icd)
# 678 kt ids that are identified by icd codes

icd_not_kt <- setdiff(icds, kt_mrn$PMBB_ID)
length(icd_not_kt)
# 101 icd codes not KT ides

# left join the kt not icd codes w/ their kt ids
kt_not_icd_df <- kt_mrn %>% filter(PMBB_ID %in% kt_not_icd)
dim(kt_not_icd_df)

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



