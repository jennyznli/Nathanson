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
breast_results <- select_samples(
    sample_name = "simplexo",
    icd_codes = c("^C50", "^Z85.3", "^D05","^Z86.000",
                  "^174", "^V10.3", "^233.0"),
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

# ========================
# Non Cancer ICD selection
# ========================
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
progeny_case <- read_excel(here("simplexo", "ss", "br_pts_for_exwas_10022025.xlsx"))
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
write.csv(up, here("simplexo", "data", "simplexo_up_map.csv"))

# ========================
# PROGENY - CASE PULL
# ========================
dim(progeny_case)
# 5046   25

# there are duplicates in progeny bc of separate entries for each breast
sum(duplicated(progeny_case$SampNum))
dups <- progeny_case[progeny_case$SampNum %in% progeny_case$SampNum[duplicated(progeny_case$SampNum)], ]
dim(dups)
length(unique(dups$SampNum))
# around 711 ppl w more than one case of breast cancer

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
### PROGENY CASES ###
progeny_pmbb <- merge(progeny_case, up, by = "SampNum")
dim(progeny_pmbb) # 1561 samples

# filter by females
progeny_pmbb <- progeny_pmbb %>% filter(Gender == "F")
write.csv(progeny_pmbb, here("simplexo", "data", "simplexo_progeny_pmbb.txt"))

progeny_pmbb_ids <- sort(progeny_pmbb$PMBB_ID)
length(progeny_pmbb_ids)
# 1520
write.table(progeny_pmbb_ids, here("data", "progeny_breast_case_ids.txt"), row.names = FALSE, col.names = FALSE, quote = FALSE)

### MERGE ICD CASES ###
icd_cases <- fread(here("simplexo", "data", "simplexo_cancer_filtered_patients.txt"))
icd_case_ids <- readLines(here("simplexo", "data", "simplexo_cancer_filtered_patients_ids.txt"))
length(icd_case_ids)
# 3448

all_ids <- sort(unique(c(progeny_pmbb_ids, icd_case_ids)))
length(all_ids)
# 4063

write.table(all_ids, here("simplexo", "data", "simplexo_overall_v1_case_ids.txt"), row.names = FALSE, col.names = FALSE, quote = FALSE)
# check if we have age info for all... so don't write yet

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

breast_controls <- controls
# ========================
# CONTROLs - EXCLUDE CASES FROM EHR NOT CAUGHT IN ICD
# ========================
controls <- breast_controls$final_controls %>% filter(!(person_id %in% all_ids)) %>% select("person_id", "Sample_age", "Sequenced_gender")
dim(controls)
# 18278
write.table(controls, here("simplexo", "data", "simplexo_overall_control.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
write.table(controls$person_id, here("simplexo", "data", "simplexo_overall_control_ids.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

# ========================================================================
# FAMILY HISTORY
# ========================================================================

###### PROGENY ######
progeny_fh <- read_excel(here("simplexo", "ss", "br_pts_for_exwas_10022025.xlsx"))[,c("SampNum", "Gender", "NumFDR", "NumFDR_BC", "NumSDR", "NumSDR_BC")]
dim(progeny_fh)
# 5046 samples

# NumFDR - Number of first degree relatives
# NumFDR_BC - Number of first degree relatives with breast cancer
# NumSDR - Number of second degree relatives
# NumSDR_BC - Number of second degree relatives with breast cancer
# 888 - unknown

# again duplicate samples from each breast
sum(duplicated(progeny_fh$SampNum))
dups <- progeny_fh[progeny_fh$SampNum %in% progeny_fh$SampNum[duplicated(progeny_fh$SampNum)], ]
dim(dups)
# 1484

### MERGE DUPLICATE ROWS
progeny_fh <- progeny_fh %>%
    group_by(SampNum) %>%
    summarise(
        Gender = (unique(Gender)),
        NumFDR = (unique(NumFDR)),
        NumFDR_BC = (unique(NumFDR_BC)),
        NumSDR = (unique(NumSDR)),
        NumSDR_BC = (unique(NumSDR_BC)),
        .groups = "drop"
    )
sum(duplicated(progeny_fh$SampNum))
dim(progeny_fh)
# 4273 total breast cancer cases

# convert 888 -> NA in the family history columns
progeny_fh <- progeny_fh %>%
    mutate(
        NumFDR    = na_if(NumFDR, 888),
        NumFDR_BC = na_if(NumFDR_BC, 888),
        NumSDR    = na_if(NumSDR, 888),
        NumSDR_BC = na_if(NumSDR_BC, 888)
    )

# mark with 1 if they have any breast cancer in first or second degree relatives
progeny_fh <- progeny_fh %>%
    mutate(
        Family_History = if_else(
            coalesce(NumFDR_BC, 0) > 0 | coalesce(NumSDR_BC, 0) > 0,
            1L, 0L
        )
    )
sum(progeny_fh$Family_History == 1)
dim(progeny_fh)
# 2606/4273 in progeny w/ family history

# GET ONLY THOSE W/ PMBB IDS
progeny_pmbb_fh <- merge(progeny_fh, up, by = "SampNum")
dim(progeny_pmbb)
# 1561 total progeny in pmbb
progeny_pmbb_fh <- progeny_pmbb_fh %>% filter(Family_History == 1)
dim(progeny_pmbb_fh)
# 1126
progeny_pmbb_fh <- progeny_pmbb_fh %>% filter(Gender == "F")
dim(progeny_pmbb_fh)
# 1100/1561 in PMBB and progeny who are female and have family history

write.csv(progeny_pmbb_fh, here("simplexo", "data", "simplexo_progeny_family_history_df.csv"))

###### CANCER REGISTRY ######
hx <- read.csv(file.path("pmbb_1093_familyhx.csv"))
all_ids <- readLines(here("simplexo", "data", "simplexo_overall_case_ids.txt"))

dim(hx)
# 41280
length(unique(hx$PMBB_ID)) #4285
breast_hx <- hx %>% filter(grepl("breast", MEDICAL_HX, ignore.case = TRUE))
length(unique(breast_hx$PMBB_ID))
# 2228 have family history / breast history
# have to check the negative history condition !!

pmcr_hx_ids <- unique(breast_hx$PMBB_ID)
length(pmcr_hx_ids)
# 2228

hx_ids <- sort(unique(c(progeny_pmbb_fh$PMBB_ID, pmcr_hx_ids)))
length(hx_ids)
# 3445

# have to intersect w/ the newly updated final_ids because the 4285 we gave colleen the first time that
# she pulled fhx was 2 instance of breast cancer, but now we're doing 3 instances
final_hx_ids <- sort(intersect(sort(all_ids), sort(hx_ids)))
length(final_hx_ids)
# 3317
write.table(final_hx_ids, here("simplexo", "data", "simplexo_fhx_case_ids.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

###### ICD CODES ###### aactually result in no cases...
# family_history <- select_samples(
#     sample_name = "simplexo_fh",
#     icd_codes = c("^Z80.3", "^V16.3"),
#     gender_filter = "Female",
#     crep_filter = NULL,
#     min_instances = 1,
#     min_timespan = NULL,
#     age_filter = NULL,
#     exclude = FALSE,
#     pmbb_dir = here("PMBB"),
#     data_dir = here("simplexo", "data"),
#     log_dir = here("simplexo", "log")
# )

# ========================================================================
# ER STATUS
# ========================================================================
#### PROGENY
progeny_er_pos <- progeny_pmbb %>% filter(ERstatus == "Positive") #dim(er_pos)
progeny_er_neg <- progeny_pmbb %>% filter(ERstatus == "Negative") #dim(er_neg)
dim(progeny_er_pos)
# 529
dim(progeny_er_neg)
# 268

### PMCR
er <- read.csv(file.path("data", "pmbb_1093_brca_er_20251024.csv"))
dim(er)
# 4563 ids
length(unique(er$PMBB_ID)) #4285 unique individuals which is what i gave her
# Collapse the data by PMBB_ID
er <- er %>%
    group_by(PMBB_ID) %>%
    summarise(across(everything(),
                     ~ ifelse(length(unique(.)) == 1, unique(.)[1], paste(unique(.), collapse = ";"))))
unique(er$EstrogenReceptorSummarry)
er <- er %>%
    mutate(FinalER = case_when(
        EstrogenReceptorSummarry %in% c("Positive", "Positive;Unknown", "Unknown;Positive") ~ "Positive",
        EstrogenReceptorSummarry %in% c("Negative", "Negative;Unknown", "Unknown;Negative") ~ "Negative",
        EstrogenReceptorSummarry %in% c("Unknown", "null", "Positive;Negative", "Negative;Positive") ~ "Unkown",
        TRUE ~ EstrogenReceptorSummarry
    ))
unique(er$FinalER)

er_pos <- er %>% filter(FinalER == "Positive")
dim(er_pos)
# 1468
er_neg <- er %>% filter(FinalER == "Negative")
dim(er_neg)
# 366

pos <- sort(intersect(unique(c(er_pos$PMBB_ID, progeny_er_pos$PMBB_ID)), all_ids))
length(pos)
neg <- sort(intersect(unique(c(er_neg$PMBB_ID, progeny_er_neg$PMBB_ID)), all_ids))
length(neg)

write.table(neg, here("simplexo", "data", "simplexo_er_neg_case_ids.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
write.table(pos, here("simplexo", "data", "simplexo_er_pos_case_ids.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

# ========================================================================
# INVASIVE ONLY
# ========================================================================
inv_only <- select_samples(
    sample_name = "simplexo_malig",
    icd_codes = c("^C50", "^Z85.4", "^174", "^V10.3"),
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

progeny_pmbb1 <- progeny_pmbb %>%
    mutate(
        is_insitu_dx = str_detect(Diagnosis, regex("DCIS|LCIS|ductal carcinoma in situ|in situ|lobular carcinoma in situ", ignore_case = TRUE)),
        is_insitu_morph = str_detect(Morphology, regex("DCIS|LCIS|ductal carcinoma in situ|in situ|lobular carcinoma in situ", ignore_case = TRUE)),
        is_insitu = is_insitu_dx | is_insitu_morph,
        is_insitu = replace_na(is_insitu, FALSE)
    )

progeny_no_insitu <- progeny_pmbb1 %>%
    filter(!is_insitu)
dim(progeny_no_insitu)
# 1261

inv_only_ids <- unique(c(progeny_no_insitu$PMBB_ID, inv_only$filtered_patients$person_id))
length(inv_only_ids)
# 3780 cases
write.table(sort(inv_only_ids), here("simplexo", "data", "simplexo_malig_case_ids.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

# ========================
# CASES - UNDER 50
# ========================
ages <- read.table(here("simplexo", "data", "simplexo_overall_case_ages_merged.txt"), header = TRUE)

ages <- ages %>%
    select(PMBB_ID, Age) %>%
    filter(!is.na(Age)) %>%
    mutate(Age = round(Age, 2))

dim(ages)
young <- ages %>% filter(Age <= 50)
dim(young)
# 1897
write.table(sort(ages$PMBB_ID), here("simplexo", "data", "simplexo_50_case_ids.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")


