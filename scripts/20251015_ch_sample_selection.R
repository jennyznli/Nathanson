# ========================
# Packages
# ========================

library(here)
setwd("/Users/jennyzli/Documents/Nathanson")
setwd( "/Users/jennyzli/Documents/Nathanson")
source(here("R", "load_packages.R"))
source(here("R", "sample_selection.R"))
source(here("R", "control_selection.R"))
library(VennDiagram)

# ========================
# LOAD DATA
# ========================
var <- read.csv(here("ch", "data", "all_genes_filtered_variants.csv"))
x <- read.csv(here("ch", "data", "all_genes_with_filter_info.csv"))
progeny <- read_excel(here("ch", "ss", "brca_carriers_ch_freq_w_seen_in_crep_20251020.xlsx"), sheet = "Data_from_master_table")
progeny <- progeny %>% filter(DNA == "D")
progeny_tumor <- read_excel(here("ch", "ss", "brca_carriers_ch_freq_w_seen_in_crep_20251020.xlsx"), sheet = "Tumor")
progeny$SampNum <- as.numeric(progeny$SampNum)

up <- read.csv(here("simplexo", "data", "simplexo_up_map.csv"))
brad <- read_excel(here("ch", "ss", "BW-Regeneron_Workbook_20250627.xlsx"), sheet = "CREP_workbook")
old_brca1 <- read_excel(here("ch", "ss", "BRCA1.BRCA2 P.LP_04.11.21.xlsx"), sheet = "BRCA1")
old_brca2 <- read_excel(here("ch", "ss", "BRCA1.BRCA2 P.LP_04.11.21.xlsx"), sheet = "BRCA2")

PMBB_DIR = here("PMBB", "3.0")
OCC <- here(PMBB_DIR, "PMBB-Release-2024-3.0_phenotype_condition_occurrence.txt")
COV <- here(PMBB_DIR, "PMBB-Release-2024-3.0_covariates.txt")
PER <- here(PMBB_DIR, "PMBB-Release-2024-3.0_phenotype_person.txt")
cov <- fread(COV, header = TRUE)
person <- fread(PER, header = TRUE)

# ========================
# BRCA GERMLINE MUTATION CARRIERS
# ========================
### PROGENY ###
dim(progeny) #4336
length(unique(progeny$globalid))
# 2816

table(progeny$Seen_In_CREP)
# get rid of the no, unkown, NA
# [1] "Yes"                       "No"
# [3] NA                          "Unknown"
# [5] "HUP CREP"                  "Pennsylvania Hospital"
# [7] "Chester County Hospital"   "Penn Medicine Cherry Hill"
# [9] "Penn Medicine Princeton"

progeny_merged <- merge_duplicates(progeny, "SampNum")

progeny_pmbb <- progeny_merged %>%
    left_join(up, by = "SampNum") %>%
    filter(!is.na(PMBB_ID))
dim(progeny_pmbb)
# 1278 are in the PMBB

brca1_progeny_pmbb <- progeny_pmbb %>% filter(BRCA1_Presence_of_Mutation %in% c("Yes", "Obligate Carrier", "Yes-unconfirmed"))
dim(brca1_progeny_pmbb) #681
brca2_progeny_pmbb <- progeny_pmbb %>% filter(BRCA2_Presence_of_Mutation %in% c("Yes", "Obligate Carrier", "Yes-unconfirmed"))
dim(brca2_progeny_pmbb) #606

length(intersect(brca1_progeny_pmbb$globalid, brca2_progeny_pmbb$globalid))
# 9 with both and in PMBB

##### SEEN IN CREP #####
progeny_crep <- progeny %>% filter(!(Seen_In_CREP %in% c("No", "Unknown", NA)))
dim(progeny_crep)
# 1648
progeny_crep_merged <- merge_duplicates(progeny_crep, "SampNum")

progeny_crep_pmbb <- progeny_crep_merged %>%
    left_join(up, by = "SampNum") %>%
    filter(!is.na(PMBB_ID))
dim(progeny_crep_merged)
# 783 are in the PMBB

brca1_progeny_crep_pmbb <- progeny_crep_pmbb %>% filter(BRCA1_Presence_of_Mutation %in% c("Yes", "Obligate Carrier", "Yes-unconfirmed"))
dim(brca1_progeny_crep_pmbb) #409
brca2_progeny_crep_pmbb <- progeny_crep_pmbb %>% filter(BRCA2_Presence_of_Mutation %in% c("Yes", "Obligate Carrier", "Yes-unconfirmed"))
dim(brca2_progeny_crep_pmbb) #380

### BRAD'S CALLING ###
# merge w/ map to PMBB ID...
length(unique(brad$Participant)) #2869

brad_pmbb <- left_join(brad, up, by = "VCFID")
dim(brad_pmbb) # 2869

brca1_brad <- brad_pmbb %>% filter(Mutation_Gene1 %in% c("BRCA1")) %>% filter(!(is.na(Participant)))
brca2_brad <- brad_pmbb %>% filter(Mutation_Gene1 %in% c("BRCA2")) %>% filter(!(is.na(Participant)))
length(unique(brca1_brad$Participant)) #678
length(unique(brca2_brad$Participant)) #592

### OLD SHEET ###
dim(old_brca1) # 153
dim(old_brca2) # 263
length(unique(old_brca1$SampleID)) # 152
length(unique(old_brca2$SampleID)) # 263
old_brca1 <- merge_duplicates(old_brca1, "SampleID")
dim(old_brca1) # 152

##### PLOT INTERSECTIONS #####
brca1_sources <- list(
    Progeny = unique(brca1_progeny_crep_pmbb$PMBB_ID),
    Brad_Workbook = unique(brca1_brad$PMBB_ID),
    P.LP = unique(old_brca1$SampleID)
)

brca2_sources <- list(
    Progeny = unique(brca2_progeny_crep_pmbb$PMBB_ID),
    Brad_Workbook = unique(brca2_brad$PMBB_ID),
    P.LP = unique(old_brca2$SampleID)
)

##### BRCA1 VENN #####
brca1_venn <- ggVennDiagram(brca1_sources,
                            category.names = c("Progeny", "Brad's Sheet", "P.LP")) +
    scale_fill_gradient(low = "white", high = "darkgray") +
    labs(title = "BRCA1 Carrier Sources Overlap") +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
          legend.position = "none")

ggsave(here("ch", "figures", "brca1_sources_venn.png"),
       plot = brca1_venn,
       width = 9, height = 6,
       dpi = 300, bg = "white")

##### BRCA2 VENN #####
brca2_venn <- ggVennDiagram(brca2_sources,
                            category.names = c("Progeny", "Brad's Sheet", "P.LP")) +
    scale_fill_gradient(low = "white", high = "darkgray") +
    labs(title = "BRCA2 Carrier Sources Overlap") +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
          legend.position = "none")

ggsave(here("ch", "figures", "brca2_sources_venn.png"),
       plot = brca2_venn,
       width = 9, height = 6,
       dpi = 300, bg = "white")

##### INTERSECTION OF F3 CREP + F2 OLD IDS #####
length(unique(c(brca2_progeny_pmbb$PMBB_ID, brca2_brad$PMBB_ID))) # 608
length(brca2_brad$PMBB_ID) # 592 so mostly the same, so just take the intersection of brad's and the progeny crep sheets

brca2_progeny_brad <- intersect(brca2_progeny_crep_pmbb$PMBB_ID, brca2_brad$PMBB_ID)
brca1_progeny_brad <- intersect(brca1_progeny_crep_pmbb$PMBB_ID, brca1_brad$PMBB_ID)

all_brca2_crep_ids <- unique(c(brca2_progeny_brad, old_brca2$SampleID))
length(all_brca2_crep_ids)
# 592

all_brca1_crep_ids <- unique(c(brca1_progeny_brad, old_brca1$SampleID))
length(all_brca1_crep_ids)
# 555

all_brca12_crep_ids <- sort(unique(c(all_brca1_crep_ids, all_brca2_crep_ids)))
length(all_brca12_crep_ids)
# 1412

##### INTERSECTION INCLUDING NON CREP - FOR EXCLUSION/CONTROLS #####

all_brca2_ids <- unique(c(
    brca2_progeny_pmbb$PMBB_ID,
    brca2_brad$PMBB_ID,
    old_brca2$SampleID
))
length(all_brca2_ids)
# 870

all_brca1_ids <- unique(c(
    brca1_progeny_pmbb$PMBB_ID,
    brca1_brad$PMBB_ID,
    old_brca1$SampleID
))
length(all_brca1_ids)
# 834

all_brca12_ids <- sort(unique(c(all_brca1_ids, all_brca2_ids)))
length(all_brca12_ids) # 1692

# use this to find non brca1/2 carriers as controls...
all_ids <- sort(unique(cov$person_id))
non_brca12 <- setdiff(all_ids, unique(c(all_brca1_ids, all_brca2_ids)))
length(non_brca12)
# 55481 total nonBRCA1/2
# need to separate this into no cancer/pretreatemtn and cancer

# ========================================================================
# STRATA 1: BRCA CARRIERS W/O ANY TYPE OF CANCER...
# ========================================================================
malig_insitu_neoplasms <- c(
    "^C(?!44)",           # Malignant neoplasms except skin cancer (C44)
    "^Z85(?!\\.828)",     # Personal history of malignant neoplasm except skin
    "^D0(?!4)[0-9]",  # Benign/in situ neoplasms except skin cancer (DO4)
    "^Z86.00",             # Personal history of in-situ neoplasm
    "^(?!173)(?:1[4-9][0-9]|2[0-3][0-9])",  # ICD-9 malignant/insitu neoplasms except 173
    "^V10(?!\\.83)"       # ICD-9 personal history of malignant neoplasm except skin
)

myeloid_neoplasms <- c(
    # Myeloproliferative/myelodysplastic disorders:
    "^D47.3",               # ICD-10: Essential thrombocythemia
    "^D46",                  # ICD-10: Myelodysplastic syndromes
    "^D47.4",                # ICD-10: Osteomyelofibrosis
    "^D75.81",               # ICD-10: Myelofibrosis
    "^D45",                  # ICD-10: Polycythemia vera
    "^238.71",  # Essential thrombocythemia
    "^238.72",  # Myelodysplastic syndrome, low grade
    "^238.73",  # Myelodysplastic syndrome, high grade
    "^238.74",  # Myelodysplastic syndrome with 5q deletion
    "^238.75",  # Myelodysplastic syndrome, unspecified
    "^289.83",  # Myelofibrosis
    "^238.4"    # Polycythemia vera
)
cancers <- select_samples(
      sample_name = "CH_cancer",              # Name of the type as string
      icd_codes = c(myeloid_neoplasms, malig_insitu_neoplasms),     # ICD9/10 code patterns (regex-compatible), or NULL if taking all ICD codes
      gender_filter = NULL,                # "Male", "Female", or NULL
      crep_filter = NULL,                  # TRUE - include only CREP samples, FALSE - exclude CREP samples, NULL - no filter for CREP
      age_filter = NULL,                   # int - minimum age, vector - range of edges, or NULL for no filter
      min_instances = 2,                   # Minimum instances of ICD codes
      min_timespan = NULL,                 # Minimum timespan between first and last diagnosis occurrence (days)
      exclude = FALSE,                     # Whether to exclude patients with certain ICD codes (e. g. those w/ history of other cancer)
      data_dir = here("ch", "data"),                     # Directory for outputs
      pmbb_dir = here("PMBB"),                     # Directory for PMBB input w/ covariates, condition_occurrences, demographic tables
      log_dir = here("ch", "log"),                      # Directory for log files
  )
dim(cancers$filtered_patients)
# 24756

##### ALL CANCER IDS #####
progeny_tumor_merged <- merge_duplicates(progeny_tumor, "globalid")
# 2376 ppl in tumor
# i should get the youngest age?and the earliest date?
# Function to extract earliest age from semicolon-separated string

progeny_tumor_merged$CaDxAge_Progeny <- sapply(progeny_tumor_merged$CaDxAge, get_earliest_age)
progeny_tumor_merged$CaDxDate_Progeny <- sapply(progeny_tumor_merged$CaDxDate, get_earliest_date)
progeny_tumor_merged$CaDxDate_Progeny <- as.Date(progeny_tumor_merged$CaDxDate_Progeny, origin = "1970-01-01")
progeny_tumor_merged <- progeny_tumor_merged %>%
    filter(!is.na(CaDxAge_Progeny))
dim(progeny_tumor_merged)
# 2306   20

dim(progeny_pmbb)
# 1278
progeny_pmbb$globalid <- as.numeric(progeny_pmbb$globalid)
progeny_pmbb_tumor <- inner_join(progeny_pmbb, progeny_tumor_merged, by = "globalid")
dim(progeny_pmbb_tumor)
#  714   54

cancer_ids <- sort(unique(c(progeny_pmbb_tumor$PMBB_ID, cancers$all_patients$person_id)))
# cancer_ids <- sort(unique(c(progeny_pmbb$PMBB_ID, cancers$filtered_patients$person_id)))
# not sure abt this...

length(cancer_ids)
# 29143 - not filtered

##### LOAD IN BREAST CANCER #####
progeny_breast <- read_excel(here("simplexo", "ss", "br_pts_for_exwas_10022025.xlsx"))
progeny_breast_merged <- merge_duplicates(progeny_breast, "SampNum")
progeny_breast_pmbb <- merge(progeny_breast_merged, up, by = "SampNum")

progeny_breast_pmbb_ids <- sort(progeny_breast_pmbb$PMBB_ID)
length(progeny_breast_pmbb_ids)
# 1561 (non gender filter)

length(unique(c(cancer_ids, progeny_breast_pmbb_ids)))
# 29432 - with just one instance, breast progeny adds some...

# ========================
# MAKE OVERALL DF FOR BRCA1/2 CARRIERS ONLY...
# ========================
brca12_df <- data.frame(
    PMBB_ID = all_brca12_crep_ids,
    Mutation = NA,
    Cancer = 0,
    Dx_Age = NA,
    Sample_Age = NA
)

### MUTATION INFO ###
brca12_df$Mutation <- ifelse(brca12_df$PMBB_ID %in% all_brca1_crep_ids, "BRCA1",
                              ifelse(brca12_df$PMBB_ID %in% all_brca2_crep_ids, "BRCA2", NA))

### CANCER STATUS ###
brca12_df$Cancer <- ifelse(brca12_df$PMBB_ID %in% unique(c(cancer_ids, progeny_breast_pmbb_ids)), 1, 0)

### DIAGNOSIS AGE AND DATE ###
# ICD CODES
dx_date <- cancers$all_patients %>%
    select("person_id", "first_date")
# 28883

dx_age <- dx_date %>%
    left_join(person[, c("birth_datetime", "person_id")], by = "person_id") %>%
    mutate(Dx_Age = as.numeric(difftime(first_date, birth_datetime, units = "days")) / 365.25)

head(dx_age)
colnames(dx_age) <- c("PMBB_ID", "Dx_Date", "Birth_Date", "Dx_Age")

# SIMPLEXO
simplexo_age <- read.table(file.path("simplexo", "data", "simplexo_overall_case_ages_merged.txt"), header = TRUE)

# PROGENY
tumor_age <- progeny_pmbb_tumor %>%
    select(PMBB_ID,
           CaDxDate = CaDxDate_Progeny,
           CaDxAge = CaDxAge_Progeny)

# Add Dx_Age based on priority: simplexo_age > tumor_age > dx_age
brca12_df$Dx_Age <- ifelse(
    brca12_df$Cancer == 1,
    coalesce(
        simplexo_age$Age[match(brca12_df$PMBB_ID, simplexo_age$PMBB_ID)],
        tumor_age$CaDxAge[match(brca12_df$PMBB_ID, tumor_age$PMBB_ID)],
        dx_age$Dx_Age[match(brca12_df$PMBB_ID, dx_age$PMBB_ID)]
    ),
    NA
)

### SAMPLE AGE ###
brca12_df$Sample_Age <- cov$Sample_age[match(brca12_df$PMBB_ID, cov$person_id)]

head(brca12_df)

# See which source provides ages
cat("Simplexo ages:", sum(!is.na(simplexo_age$Age[match(brca12_df$PMBB_ID, simplexo_age$PMBB_ID)])), "\n")
cat("Tumor ages:", sum(!is.na(tumor_age$CaDxAge[match(brca12_df$PMBB_ID, tumor_age$PMBB_ID)])), "\n")
cat("ICD ages:", sum(!is.na(dx_age$Dx_Age[match(brca12_df$PMBB_ID, dx_age$PMBB_ID)])), "\n")

# ========================
# MAKE OVERALL DF FOR NON-BRCA1/2 CARRIERS ONLY...
# ========================
# this is from cancer
non_brca12_df <- data.frame(
    PMBB_ID = non_brca12,
    Mutation = NA,
    Cancer = 0,
    Dx_Age = NA,
    Sample_Age = NA
)

### CANCER STATUS ###
non_brca12_df$Cancer <- ifelse(non_brca12_df$PMBB_ID %in% unique(c(cancer_ids, progeny_breast_pmbb_ids)), 1, 0)

# Add Dx_Age based on priority: simplexo_age > dx_age
non_brca12_df$Dx_Age <- ifelse(
    non_brca12_df$Cancer == 1,
    coalesce(
        simplexo_age$Age[match(non_brca12_df$PMBB_ID, simplexo_age$PMBB_ID)],
        dx_age$Dx_Age[match(non_brca12_df$PMBB_ID, dx_age$PMBB_ID)]
    ),
    NA
)

### SAMPLE AGE ###
non_brca12_df$Sample_Age <- cov$Sample_age[match(non_brca12_df$PMBB_ID, cov$person_id)]

head(non_brca12_df)

cat("Simplexo ages:", sum(!is.na(simplexo_age$Age[match(non_brca12_df$PMBB_ID, simplexo_age$PMBB_ID)])), "\n")
cat("ICD ages:", sum(!is.na(dx_age$Dx_Age[match(non_brca12_df$PMBB_ID, dx_age$PMBB_ID)])), "\n")

# ========================================================================
# STRATA 1: PRE DIAGNOSIS & NON CANCER
# ========================================================================
# only keep if sample age is BEFORE the cancer diagnosis age
case_strata1 <- brca12_df %>%
    filter(Sample_Age < Dx_Age | Cancer == 0)
dim(case_strata1)
# 570

control_strata1 <- non_brca12_df %>%
    filter(Sample_Age < Dx_Age | Cancer == 0)
dim(control_strata1)
# 32808

# ========================================================================
# STRATA 2: AFTER CANCER
# ========================================================================
case_strata2 <- brca12_df %>%
    filter(Sample_Age >= Dx_Age & Cancer == 1)
dim(case_strata2)
# 615

control_strata2 <- non_brca12_df %>%
    filter(Sample_Age >= Dx_Age & Cancer == 1)
dim(control_strata2)
# 22664

# ========================================================================
# NOW SEE HOW MANY BRCA CARRIERS HAVE A CH MUTATION?
# ========================================================================
### these calls are only CREP though...

length(unique(var$Sample.ID))
# 2335
length(unique(up$PMBB_ID))
# 2864

ch_ids <- as.data.frame(sort(unique(var$Sample.ID)))
colnames(ch_ids) <- "Sample.ID"
# map over
ch_id_pmbb <- ch_ids %>%
    left_join(up, by = c("Sample.ID" = "match_col")) %>%
    drop_na()
dim(ch_id_pmbb)
# 2332

case_strata1 <- case_strata1 %>%
    mutate(CH_Variant = PMBB_ID %in% ch_id_pmbb$PMBB_ID)
sum(case_strata1$CH_Variant)
dim(case_strata1)
# 315 / 570

case_strata2 <- case_strata2 %>%
    mutate(CH_Variant = PMBB_ID %in% ch_id_pmbb$PMBB_ID)
sum(case_strata2$CH_Variant)
dim(case_strata2)
# 311 / 615

# expect nothing lol
control_strata1 <- control_strata1 %>%
    mutate(CH_Variant = PMBB_ID %in% ch_id_pmbb$PMBB_ID)
sum(control_strata1$CH_Variant)
dim(control_strata1)
# 368 / 32808

control_strata2 <- control_strata2 %>%
    mutate(CH_Variant = PMBB_ID %in% ch_id_pmbb$PMBB_ID)
sum(control_strata2$CH_Variant)
dim(control_strata2)
# 908 / 22664








