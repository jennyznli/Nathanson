# ========================
# SIMPLEXO CLINICAL DATA
# ========================

library(here)
source(here("R", "load_packages.R"))
setwd(here("simplexo"))

# ========================
# LOAD DATA
# ========================
PMBB_DIR = here("PMBB", "3.0")
OCC <- here(PMBB_DIR, "PMBB-Release-2024-3.0_phenotype_condition_occurrence.txt")
COV <- here(PMBB_DIR, "PMBB-Release-2024-3.0_covariates.txt")
PER <- here(PMBB_DIR, "PMBB-Release-2024-3.0_phenotype_person.txt")
cov <- fread(COV, header = TRUE)
person <- fread(PER, header = TRUE)

flags <- fread(here("ss", "rgcname_pmbbid_metadata_flags.csv"))
progeny_case <- read_excel(here("ss", "br_pts_for_exwas_08292025.xlsx"))
pmcr <- read.csv(here("ss", "pmbb_147_pmcrbreastage.csv"))

# ========================
# MATCH PMBB TO UPXXXX IDS
# ========================
up_pattern <- "^UPENN-PMBB_UP[0-9]+_UP[0-9]+$"

# find all cases that match pattern among three columns, filter for those w/ this ID
up <- flags %>%
    mutate(match_col = apply(select(., RGC_sample_name, ID1, ID2), 1, function(row) {
        m <- row[grepl(up_pattern, row)]
        if (length(m) > 0) m[1] else NA
    })) %>%
    filter(!is.na(match_col))
print(dim(up))
# 2864   14

up <- select(up, c("PMBB_ID",  "RGC_sample_name", "ID1", "ID2", "CREP", "match_col"))
up$VCFID <- sub(".*(UP[0-9]{4}).*", "\\1", up$match_col)

up <- select(up, c("PMBB_ID",  "match_col", "VCFID"))
up$SampNum <- numbers <- as.numeric(gsub("UP", "", up$VCFID))
# 2864 samples in PMBB that have UP ID

# ========================
# PROGENY - CASE PULL
# ========================
dim(progeny_case)
# 5046   25

# there are duplicates in progeny bc of separate entries for each breast
sum(duplicated(progeny_case$SampNum))
dups <- progeny_case[progeny_case$SampNum %in% progeny_case$SampNum[duplicated(progeny_case$SampNum)], ]
dim(dups)
# around 773 dups

### MERGE DUPLICATE ROWS
progeny_case <- progeny_case %>%
    group_by(SampNum) %>%
    summarise(
        across(-any_of("SampNum"), ~ paste(unique(as.character(.x)), collapse = "; ")),
        .groups = "drop"
    )
dim(progeny_case) #4273

# ========================
# FINAL CASES - MERGE ALL
# ========================
# from icd code selection
icds <- unique(read.table(here("data", "breast_simplexo3_v2_cancer_filtered_patients_ids.txt"))$V1)
print(length(unique(icds))) # 3729
icds_df <- as.data.frame(icds) # for manual search

# PROGENY IN THE PMBB
progeny_pmbb <- merge(progeny_case, up, by = "SampNum")
dim(progeny_pmbb) # 1561 samples

# merge w/ covariates to check the gender/sex
check <- merge(progeny_pmbb, cov, by.x = "PMBB_ID", by.y = "person_id")
dim(check) # 1561   40

# check that genders match
# check$Gender2 <- ifelse(check$Gender == "F", "Female", "Male")
# (check$Gender2 == check$Sequenced_gender) # they all match!

progeny_pmbb <- progeny_pmbb %>% filter(Gender == "F")
progeny_pmbb_ids <- sort(progeny_pmbb$PMBB_ID)
write.table(progeny_pmbb_ids, here("data", "progeny_ids.txt"), row.names = FALSE, col.names = FALSE, quote = FALSE)
length(progeny_pmbb_ids)
# 1520

# ========================
# ICD AGE INFERENCE
# ========================
# read summary results from previous ICD analysis for cases w/ first and last occurrence of diagnosis
icd_cases <- fread(here("breast_simplexo3_v2_cancer_filtered_patients.txt"))
icd_case_ids <- readLines(here("breast_simplexo3_v2_cancer_filtered_patients_ids.txt"))

all_ids <- sort(unique(c(progeny_pmbb_ids, icd_case_ids)))
write.table(all_ids, here("data", "simplexo_pmbb_ids.txt"), row.names = FALSE, col.names = FALSE, quote = FALSE)

# get DOB from person file
pmbb_dob <- person %>% select(person_id, birth_datetime)
dob_dx <- merge(icd_cases, pmbb_dob, by = "person_id")

# calculate diagnosis age: first date of cancer occurrence - DOB
dob_dx$CaDxAge_ICD <- as.numeric(difftime(dob_dx$first_date, dob_dx$birth_datetime, units = "days")) / 365
# manual check - seems right
# select(dob_dx, c("first_date", "birth_datetime", "CaDxAge"))

# ========================
# AGE MERGING
# ========================
### GET AGES FROM COLLEEN'S PULL
dim(pmcr)
pmcr_case <- pmcr %>% filter((pmcr$PmcrBreastFlag) == 1) #2510    5
dim(pmcr_case)
pmcr_ages <- pmcr_case %>% select(c("PMBB_ID", "PmcrBreastFirstDxDate_age")) # 2510    5
colnames(pmcr_ages) <- c("PMBB_ID", "CaDxAge_Colleen")

### GET AGES FROM PROGENY
# Function to get the earliest date from a comma-separated string
get_earliest_age <- function(age_str) {
    if (!is.na(age_str) && age_str != "") {
        # Split by comma, convert to numeric, and return the minimum age
        ages <- as.numeric(strsplit(age_str, "; ")[[1]])
        return(min(ages, na.rm = TRUE))
    } else {
        return(NA)  # Return NA if the string is empty or NA
    }
}

# apply the function to extract the earliest age from CaDxAge in progeny_pmbb
dim(progeny_pmbb) #1520   29
progeny_pmbb$CaDxAge_Progeny <- sapply(progeny_pmbb$CaDxAge, get_earliest_age)
progeny_pmbb$CaDxAge_Progeny[is.infinite(progeny_pmbb$CaDxAge_Progeny)] <- NA
progeny_ages <- progeny_pmbb %>% select(PMBB_ID, CaDxAge_Progeny)

### GET AGES FROM ICD CODES
icd_ages <- dob_dx %>% select(person_id, CaDxAge_ICD) #3729    2
dim(icd_ages)

### MERGE
base_df <- data.frame(PMBB_ID = all_ids)
merged <- base_df %>%
    left_join(progeny_ages, by = "PMBB_ID") %>%  # first try to merge progeny ages
    left_join(pmcr_ages, by = "PMBB_ID") %>%   # then merge colleen's ages
    left_join(icd_ages, by = c("PMBB_ID" = "person_id")) %>%  # then merge icd ages
    mutate(
        # Calculate final Age using coalesce
        Age = coalesce(CaDxAge_Progeny, CaDxAge_Colleen, CaDxAge_ICD),

        # Find the earliest diagnosis age
        EarliestDiagnosis = pmin(CaDxAge_Progeny, CaDxAge_Colleen, CaDxAge_ICD, na.rm = TRUE),

        # Check if the final Age is different from the earliest diagnosis
        AgeDifference = abs(Age - EarliestDiagnosis) > 2
    )
diff <- merged %>% filter(AgeDifference == TRUE)

write.table(merged_data, here("data", "simplexo3_case_ages_merged.txt"))

merged_sel <- merged %>%
    select(PMBB_ID, Age) %>%
    filter(!is.na(Age))
dim(merged_sel)
write.table(merged_sel, here("data", "simplexo3_case_ages.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

# some with empty ages
# 1766 2982 3055
merged_sel[1766,]
merged_sel[2982,]
merged_sel[3055,]

# > merged_sel[1766,]
# PMBB_ID Age
# 1766 PMBB4770145685120  NA
# > merged_sel[2982,]
# PMBB_ID Age
# 2982 PMBB7373064082967  NA
# > merged_sel[3055,]
# PMBB_ID Age
# 3055 PMBB7499670247051  NA



