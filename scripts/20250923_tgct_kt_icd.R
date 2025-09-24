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
print(dim(kt2))

kt <- select(kt, c("PMBB_ID",  "RGC_sample_name", "ID1", "ID2", "ID3", "ID4", "CREP", "match_col"))
kt$ID <- sub(".*(KT[0-9]{7}).*", "\\1", kt$match_col)

kt <- select(kt, c("PMBB_ID",  "match_col", "ID"))
kt$SampNum <- numbers <- as.numeric(gsub("KT", "", kt$ID))

write.csv(kt, here("PMBB", "3.0", "pmbb3_kt_key.csv"))

# ========================
# COMPARE W/ ICD CODES
# ========================
length(icds) # 779
# length(unique(icds))

kt_not_icd <- setdiff(kt$PMBB_ID, icds)
length(kt_not_icd)
# 226 kt ids that are not detected by ICD, need to be checked

kt_icd <- intersect(kt$PMBB_ID, icds)
length(kt_icd)
# 678 kt ids that are identified by icd codes

icd_not_kt <- setdiff(icds, kt$PMBB_ID)
length(icd_not_kt)
# 101 icd codes not KT ides


# left join the kt not icd codes w/ their kt ids
kt_not_icd_df <- kt %>% filter(PMBB_ID %in% kt_not_icd)
colnames(kt_not_icd_df) <- c("PMBB_ID", "KT_ID", "ID", "SampNum")
write.table(kt_not_icd, here("tgct", "data", "tgct_kt_not_icd.txt"),
            col.names = FALSE, row.names = FALSE, quote = FALSE)
write.csv(kt_not_icd_df, here("tgct", "data", "tgct_kt_not_icd_df.csv"), row.names = FALSE, quote = FALSE)

dim(kt_not_icd_df)
length(kt_not_icd)

# ========================
# BRAD'S FILE
# ========================
dim(brad)
# 59624    10 -- all of PMBB??
table(brad$seq_site)
# CAG   NCI  PMBB REGEN  TCGA UPENN
# 61  1896   295 53223  3961   150    38

# kt_pattern <- "^UPENN-PMBB_KT[0-9]+_KT[0-9]+$"
brad_kt <- brad[apply(brad, 1, function(row) any(grepl("KT", row, fixed = TRUE))), ]
