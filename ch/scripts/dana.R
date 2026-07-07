# ========================
# CHECKING DANA's BRCA1/2 carrier list
# ========================
library(here)
setwd("/Users/jennyzli/Documents/Nathanson")
source(here("R", "config.R"))

# ========================
# READ DATA
# ========================
# case lists
b1_case <- read.csv(file.path("ch", "data", "pmbb_brca1_case_ids.csv"))$x
b2_case <- read.csv(file.path("ch", "data", "pmbb_brca2_case_ids.csv"))$x
b12_case <- read.csv(file.path("ch", "data", "pmbb_brca12_case_ids.csv"))$x
b12_both_case <- read.csv(file.path("ch", "data", "pmbb_brca12_both_case_ids.csv"))$x
b12_df <- read.csv(file.path("ch", "data", "pmbb_brca12_cov.csv"))
b12_df2 <- read.csv(file.path("ch", "data", "pmbb_brca12_cov_df.csv"))

up <- read.csv(here("simplexo", "data", "simplexo_up_map.csv"))

# dana's list
dana <- read.csv(here("brca", "pmbb_1099_breastcancercases_wildtypesamples_pmbb_20260416_pmbbid.csv"))
# 703

# PMBB
cov <- fread(here("PMBB", "4.0", "PMBB-Release-2026-4.0_phenotype_covariates.txt"), header = TRUE)
person <- fread(here("PMBB", "4.0", "PMBB-Release-2026-4.0_phenotype_person.txt"), header = TRUE)
brca <- fread(here("PMBB", "4.0", "PMBB-Release-2026-4.0_phenotype_cancer_brca.txt"), header = TRUE)

obs <- fread(here("PMBB", "4.0", "PMBB-Release-2026-4.0_phenotype_observation.txt"), header = TRUE)
cond <- fread(here("PMBB", "4.0", "PMBB-Release-2026-4.0_phenotype_condition_occurrence.txt"), header = TRUE)

# ========================
# CHECK
# ========================
dana <- dana %>% mutate(
    Dana_BRCA12 = ifelse(grepl("BRCA", PositiveBiomarkerNames), 1, 0)
)

dana <- dana %>% left_join(b12_df %>% select(person_id, BRCA12_Case, BRCA1_Case, BRCA2_Case, BRCA1_BRCA2_Case),
                   by = c("PMBB_ID" = "person_id"))

sum(dana$Dana_BRCA12, na.rm = TRUE)
# 69
sum(dana$BRCA12_Case, na.rm = TRUE)
# 72

dana <- dana %>% mutate(
    mismatch = (BRCA12_Case != Dana_BRCA12)
)
mismatch <- dana %>% filter(mismatch == TRUE)
jenny_detect <- mismatch %>% filter(BRCA12_Case == 1)
dana_detect <- mismatch %>% filter(Dana_BRCA12 == 1)
missing <- dana %>%
    filter(is.na(mismatch)) %>%
    left_join(cov %>% select(person_id, batch),
                              by = c("PMBB_ID" = "person_id"))
dim(missing)
# 51 in freeze 4.0

excel_sheets <- list(
    "Complete"               = dana,
    "Jenny_Found_Dana_None"  = jenny_detect,
    "Dana_Found_Jenny_None"  = dana_detect,
    "Freeze_4.0"             = missing
)

write_xlsx(
    x = excel_sheets,
    path = here("ch", "data", "pmbb_1099_breastcancercases_wildtypesamples_pmbb_20260623_pmbbid_jenny.xlsx")
)

# ========================
# COMPILE
# ========================
# > length(unique(cond$person_id))
# [1] 68634
# > length(unique(obs$person_id))
# [1] 67828

# PMBB2802414598396 was found in dana's but not mine
# same with PMBB7467803196897

f23_b2_ids <- read.csv(file.path("ch", "data", "pmbb_brca2_prev.csv"))
"PMBB2802414598396" %in% f23_b2_ids
"PMBB7467803196897" %in% f23_b2_ids

b1_mine <- read.csv(file.path("ch", "data", "brca1_all_filtered.csv"))
b2_mine <- read.csv(file.path("ch", "data", "brca2_all_filtered.csv"))

"PMBB2802414598396" %in% b1_mine$Sample.ID
"PMBB2802414598396" %in% b2_mine$Sample.ID
# these both had low LoF

"PMBB7467803196897" %in% b1_mine$Sample.ID
"PMBB7467803196897" %in% b2_mine$Sample.ID







