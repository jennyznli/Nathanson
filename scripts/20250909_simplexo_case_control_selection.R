# ========================
# SIMPLEXO - CASE & CONTROL SELECTION
# ========================
library(here)
setwd(here("simplexo"))

source(here("R", "load_packages.R"))
source(here("R", "sample_selection.R"))
source(here("R", "control_selection.R"))

# ========================
# Breast Cancer ICD Selection
# ========================
# no filtering, all
breast_results <- select_samples(
    sample_name = "simplexo",
    icd_codes = c("^C50", "^Z85.3", "^D05","^Z86.000",
                  "^174", "^V10.3", "^233.0"),
    gender_filter = "Female",
    crep_filter = NULL,
    min_instances = 2,
    min_timespan = NULL,
    age_filter = NULL,
    exclude = FALSE,
    pmbb_dir = here("PMBB"),
    data_dir = here("simplexo", "data"),
    log_dir = here("simplexo", "log")
)

# ========================
# Non Cancer ICD selection
# ========================
# malig_neoplasms <- c("^C(?!44)", "^Z85(?!\\.828)", "^D05", "^Z86.000",
#                      "^(?!173)(?:1[4-9][0-9]|20[0-8]|209\\.[0-3])","^233.0",
#                      "^V10(?!\\.83)")
#
# controls <- select_controls(
#     control_name = "simplexo",
#     exclude_codes = malig_neoplasms,
#     gender_filter = "Female",
#     age_filter = NULL,
#     crep_filter = FALSE,
#     pmbb_dir = here("PMBB", "3.0"),
#     data_dir = here("simplexo", "data"),
#     log_dir = here("simplexo", "log")
# )

# ========================
# LOAD DATA
# ========================
PMBB_DIR = here("PMBB", "3.0")
OCC <- here(PMBB_DIR, "PMBB-Release-2024-3.0_phenotype_condition_occurrence.txt")
COV <- here(PMBB_DIR, "PMBB-Release-2024-3.0_covariates.txt")
PER <- here(PMBB_DIR, "PMBB-Release-2024-3.0_phenotype_person.txt")
cov <- fread(COV, header = TRUE)
person <- fread(PER, header = TRUE)

flags <- fread(here("PMBB", "3.0", "rgcname_pmbbid_metadata_flags.csv"))
progeny_case <- read_excel(here("simplexo", "ss", "br_pts_for_exwas_08292025.xlsx"))
pmcr <- read.csv(here("simplexo", "ss", "pmbb_147_pmcrbreastage.csv"))
icds <- unique(read.table(here("simplexo", "data", "simplexo_cancer_filtered_patients_ids.txt"))$V1)

# ========================
# MATCH PMBB TO UPXXXX IDS
# ========================
up_pattern <- "^UPENN-PMBB_UP[0-9]+_UP[0-9]+$"

# find all cases that match pattern among three columns, filter for those w/ this ID
up <- flags %>%
    mutate(match_col = apply(select(., RGC_sample_name, ID1, ID2, ID3, ID4), 1, function(row) {
        m <- row[grepl(up_pattern, row)]
        if (length(m) > 0) m[1] else NA
    })) %>%
    filter(!is.na(match_col))
print(dim(up))
# 2864   14

up <- select(up, c("PMBB_ID",  "RGC_sample_name", "ID1", "ID2", "ID3", "ID4", "CREP", "match_col"))
up$VCFID <- sub(".*(UP[0-9]{4}).*", "\\1", up$match_col)

up <- select(up, c("PMBB_ID",  "match_col", "VCFID"))
up$SampNum <- as.numeric(gsub("[a-zA-Z]", "", up$VCFID))
dim(up)
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
print(length(unique(icds))) # 3729
icds_df <- as.data.frame(icds) # for manual search

# PROGENY IN THE PMBB
progeny_pmbb <- merge(progeny_case, up, by = "SampNum") #inner join
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
# MERGE CASE LISTS
# ========================
# read summary results from previous ICD analysis for cases w/ first and last occurrence of diagnosis
icd_cases <- fread(here("simplexo", "data", "simplexo_cancer_filtered_patients.txt"))
icd_case_ids <- readLines(here("simplexo", "data", "simplexo_cancer_filtered_patients_ids.txt"))

all_ids <- sort(unique(c(progeny_pmbb_ids, icd_case_ids)))
length(all_ids)
# write.table(all_ids, here("simplexo", "data", "simplexo_overall_case_ids.txt"), row.names = FALSE, col.names = FALSE, quote = FALSE)
# check if we have age info for all... so don't write yet

# ========================
# INFER AGES FROM COVARIATE
# ========================
# get DOB from person file
pmbb_dob <- person %>% select(person_id, birth_datetime)
dob_dx <- merge(icd_cases, pmbb_dob, by = "person_id", all.x=TRUE)

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
    # Handle NA, NULL, or empty values
    if (is.na(age_str) || is.null(age_str) || age_str == "" || age_str == "NA") {
        return(NA)
    }
    age_str <- as.character(age_str)
    age_parts <- strsplit(age_str, "; ")[[1]]
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
    left_join(progeny_ages, by = "PMBB_ID") %>%  # merge progeny ages
    left_join(pmcr_ages, by = "PMBB_ID") %>%   # merge colleen's ages (cancer registry)
    left_join(icd_ages, by = c("PMBB_ID" = "person_id")) %>%  # merge icd ages
    mutate(
        # Calculate final Age in order from progeny, colleen, icd
        Age = coalesce(CaDxAge_Progeny, CaDxAge_Colleen, CaDxAge_ICD),
        # Find the earliest diagnosis age
        EarliestDiagnosis = pmin(CaDxAge_Progeny, CaDxAge_Colleen, CaDxAge_ICD, na.rm = TRUE),
        # Check if the final Age is different from the earliest diagnosis
        AgeDifference = abs(Age - EarliestDiagnosis) > 2
    )
# diff <- merged %>% filter(AgeDifference == TRUE)

write.table(merged, here("simplexo", "data", "simplexo_overall_case_ages_merged.txt"))

# x <- merged %>% filter(is.na(Age))
merged_sel <- merged %>%
    select(PMBB_ID, Age) %>%
    filter(!is.na(Age)) %>%
    mutate(Age = round(Age, 2))

dim(merged_sel)
# 3 didn't have any ages
length(sort(merged_sel$PMBB_ID))
write.table(sort(merged_sel$PMBB_ID), here("simplexo", "data", "simplexo_overall_case_ids.txt"), row.names = FALSE, col.names = FALSE, quote = FALSE)

# ========================
# GET 3 DIFFERENT SAMPLE SETS
# ========================
# merged_cov <- left_join(merged_sel, cov, by = c("PMBB_ID" = "person_id")) #4285   14
# merged_crep <- merged_cov %>% filter(CREP_HighRisk_Flag == 1) %>% select("PMBB_ID", "Age", "Sequenced_gender")
# merged_no_crep <- merged_cov %>% filter(CREP_HighRisk_Flag == 0) %>% select("PMBB_ID", "Age", "Sequenced_gender")
# merged_all <- merged_cov %>% select("PMBB_ID", "Age", "Sequenced_gender")
# dim(merged_crep)
# dim(merged_no_crep)
# dim(merged_cov)
# # unique(merged_cov$Sequenced_gender)
#
# write.table(merged_crep, here("simplexo", "data", "simplexo_case_crep_ages.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
# write.table(merged_no_crep, here("simplexo", "data", "simplexo_case_nocrep_ages.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
# write.table(merged_all, here("simplexo", "data", "simplexo_case_all_ages.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

# ========================
# Non Cancer selection
# ========================
# malig_neoplasms <- c("^C(?!44)", "^Z85(?!\\.828)", "^(?!173)(?:1[4-9][0-9]|20[0-8]|209\\.[0-3])", "^V10(?!\\.83)")
# malign_benign_neoplasms <- c("^C(?!44)", "^D", "^(?!173)(?:1[4-9][0-9]|2[0-3][0-9])", "^Z85", "^V10(?!\\.83)")

malig_neoplasms <- c("^C(?!44)", "^Z85(?!\\.828)", "^D05", "^Z86.000",
                     "^(?!173)(?:1[4-9][0-9]|20[0-8]|209\\.[0-3])","^233.0",
                     "^V10(?!\\.83)")

controls <- select_controls(
    control_name = "simplexo",
    exclude_codes = malig_neoplasms,
    gender_filter = "Female",
    age_filter = NULL,
    crep_filter = FALSE,
    pmbb_dir = here("PMBB", "3.0"),
    data_dir = here("simplexo", "data"),
    log_dir = here("simplexo", "log")
)

# breast_controls <- controls
# ========================
# CONTROLs - EXCLUDE CASES FROM EHR NOT CAUGHT IN ICD
# ========================
controls <- breast_controls$final_controls %>% filter(!(person_id %in% all_ids)) %>% select("person_id", "Sample_age", "Sequenced_gender")
dim(controls)
# 18278
write.table(controls, here("simplexo", "data", "simplexo_overall_control.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

# ========================
# FAMILY HISTORY
# ========================
progeny_fh <- read_excel(here("simplexo", "ss", "br_pts_for_exwas_10022025.xlsx"))
progeny_fh <- progeny_fh[,c("SampNum", "Gender", "NumFDR", "NumFDR_BC", "NumSDR", "NumSDR_BC")]
dim(progeny_fh)
# 5046    4

# NumFDR - Number of first degree relatives
# NumFDR_BC - Number of first degree relatives with breast cancer
# NumSDR - Number of second degree relatives
# NumSDR_BC - Number of second degree relatives with breast cancer
# 888 - unknown
sum(duplicated(progeny_fh$SampNum))
dups <- progeny_fh[progeny_fh$SampNum %in% progeny_fh$SampNum[duplicated(progeny_fh$SampNum)], ]
dim(dups)
# 1484    5

### MERGE DUPLICATE ROWS
progeny_fh <- progeny_fh %>%
    group_by(SampNum) %>%
    summarise(
        Gender = first(unique(Gender)),
        NumFDR = first(unique(NumFDR)),
        NumFDR_BC = first(unique(NumFDR_BC)),
        NumSDR = first(unique(NumSDR)),
        NumSDR_BC = first(unique(NumSDR_BC)),
        .groups = "drop"
    )
dim(progeny_fh)
# 4273    5

# convert 888 -> NA in the family history columns
progeny_fh <- progeny_fh %>%
    mutate(
        NumFDR    = na_if(NumFDR, 888),
        NumFDR_BC = na_if(NumFDR_BC, 888),
        NumSDR    = na_if(NumSDR, 888),
        NumSDR_BC = na_if(NumSDR_BC, 888)
    )

progeny_fh <- progeny_fh %>%
    mutate(
        Family_History = if_else(
            coalesce(NumFDR_BC, 0) > 0 | coalesce(NumSDR_BC, 0) > 0,
            1L, 0L
        )
    )

dim(progeny_fh)
# 4273    6
sum(progeny_fh$Family_History == 1)
# 2606 w/ family history total

# GET ONLY THOSE W/ PMBB IDS
progeny_pmbb_fh <- merge(progeny_fh, up, by = "SampNum")
dim(progeny_pmbb_fh) # 1561 samples
progeny_pmbb_fh <- progeny_pmbb_fh %>% filter(Family_History == 1)  %>% filter(Gender == "F")
dim(progeny_pmbb_fh)
# 1100 left w/ pmbb id, female, and family history

write.csv(progeny_pmbb_fh, here("simplexo", "data", "simplexo_progeny_family_history.csv"))
progeny_pmbb_fh_ids <- sort(progeny_pmbb_fh$PMBB_ID)
length(progeny_pmbb_fh_ids)
write.table(progeny_pmbb_fh_ids, here("simplexo", "data", "simplexo_progeny_fh_ids.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

# STILL NEED COLLEEN'S INFORMATION ??

# ========================
# ER STATUS
# ========================
progeny_case <- read_excel(here("simplexo", "ss", "br_pts_for_exwas_08292025.xlsx"))
progeny_pmbb <- merge(progeny_case, up, by = "SampNum")
dim(progeny_pmbb)
unique(factor(progeny_pmbb$ERstatus))
# Levels: Conflicting Indeterminate Negative Not Applicable Not Done Not Reported Positive
write.csv(table(progeny_pmbb$ERstatus), here("simplexo", "data", "er_status.csv"))
er_pos <- progeny_pmbb %>% filter(ERstatus == "Positive") #dim(er_pos)
er_neg <- progeny_pmbb %>% filter(ERstatus == "Negative") #dim(er_neg)
dim(er_pos)
dim(er_neg)

# 695
# 349

# ========================
# DCIS + MALIGNANT
# ========================
# invasive + dcis
dcis <- select_samples(
    sample_name = "simplexo_malig_dcis",
    icd_codes = c("^C50", "^Z85.4", "^174", "^V10.3", "^D05", "^Z86.000", "^233.0"),
    gender_filter = "Female",
    crep_filter = NULL,
    min_instances = 2,
    min_timespan = NULL,
    age_filter = NULL,
    exclude = FALSE,
    pmbb_dir = here("PMBB"),
    data_dir = here("simplexo", "data"),
    log_dir = here("simplexo", "log")
)

# invasive only
inv_only <- select_samples(
    sample_name = "simplexo_malig",
    icd_codes = c("^C50", "^Z85.4", "^174", "^V10.3"),
    gender_filter = "Female",
    crep_filter = NULL,
    min_instances = 2,
    min_timespan = NULL,
    age_filter = NULL,
    exclude = FALSE,
    pmbb_dir = here("PMBB"),
    data_dir = here("simplexo", "data"),
    log_dir = here("simplexo", "log")
)

# should be same as before?
malig_dcis_neoplasms <- c("^C(?!44)", "^Z85(?!\\.828)", "^(?!173)(?:1[4-9][0-9]|20[0-8]|209\\.[0-3])", "^V10(?!\\.83)",  "^D05", "^233.0")
breast_50_controls <- select_controls(
    control_name = "simplexo",
    exclude_codes = malig_dcis_neoplasms,
    gender_filter = "Female",
    age_filter = NULL,
    crep_filter = FALSE,
    pmbb_dir = here("PMBB", "3.0"),
    data_dir = here("simplexo", "data"),
    log_dir = here("simplexo", "log")
)

## get it from progeny...
progeny_pmbb
# "DCIS, NOS (C50._)"
# "Ductal carcinama in situ, NOS (C50._)"
# # look for DCIS, in situ in the Diagnosis column, then

dcis_rows <- progeny_pmbb %>%
    filter(str_detect(Diagnosis, regex("DCIS|ductal carcinoma in situ|in situ", ignore_case = TRUE)))

dcis_rows %>% select(Diagnosis) %>% distinct()
# 1                                           DCIS, NOS (C50._)
# 2                           Lobular carcinoma in situ (C50._)
# 3                       Ductal carcinoma in situ, NOS (C50._)
# 4                                      Carcinoma in situ, NOS
# 5               Ductal carcinoma in situ, comedo type (C50._)
# 6 Intraductal carcinoma and lobular carcinoma in situ (C50._)
# 7                 Ductal carcinoma in situ, papillary (C50._)
# 8                                   DCIS, comedo type (C50._)
# 9                                     DCIS, papillary (C50._)

# exclude any patients with only these diagnoses...
# merge w/ other file..

# ========================
# Non Cancer - OVER 50
# ========================
malig_dcis_neoplasms <- c("^C(?!44)", "^Z85(?!\\.828)", "^(?!173)(?:1[4-9][0-9]|20[0-8]|209\\.[0-3])", "^V10(?!\\.83)",  "^D05", "^233.0")
breast_50_controls <- select_controls(
    control_name = "simplexo_50",
    exclude_codes = malig_dcis_neoplasms,
    gender_filter = "Female",
    age_filter = 50,
    crep_filter = FALSE,
    pmbb_dir = here("PMBB", "3.0"),
    data_dir = here("simplexo", "data"),
    log_dir = here("simplexo", "log")
)

# ========================
# FAMILY HISTORY
# ========================
family_history <- select_samples(
    sample_name = "simplexo_fh",
    icd_codes = c("^Z80.3", "^V16.3"),
    gender_filter = "Female",
    crep_filter = NULL,
    min_instances = 2,
    min_timespan = NULL,
    age_filter = NULL,
    exclude = FALSE,
    pmbb_dir = here("PMBB"),
    data_dir = here("simplexo", "data"),
    log_dir = here("simplexo", "log")
)

# "^Z80.3", "^V16.3"
malig_dcis_fh_neoplasms <- c("^C(?!44)", "^Z85(?!\\.828)", "^(?!173)(?:1[4-9][0-9]|20[0-8]|209\\.[0-3])",
                          "^V10(?!\\.83)",  "^D05", "^233.0", "^Z80.3", "^V16.3")
breast_fh_controls <- select_controls(
    control_name = "simplexo_fh",
    exclude_codes = malig_dcis_fh_neoplasms,
    gender_filter = "Female",
    age_filter = NULL,
    crep_filter = FALSE,
    pmbb_dir = here("PMBB", "3.0"),
    data_dir = here("simplexo", "data"),
    log_dir = here("simplexo", "log")
)
