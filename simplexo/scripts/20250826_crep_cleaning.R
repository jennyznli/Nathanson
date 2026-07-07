# ========================
# Packages
# ========================

library(here)
setwd(here("simplexo"))

source(here("R", "load_packages.R"))
source(here("R", "sample_selection.R"))
source(here("R", "control_selection.R"))

# ========================
# Analysis
# ========================
PMBB_DIR = here("PMBB", "3.0")
OCC <- file.path(PMBB_DIR, "PMBB-Release-2024-3.0_phenotype_condition_occurrence.txt")
COV <- file.path(PMBB_DIR, "PMBB-Release-2024-3.0_covariates.txt")
PER <- file.path(PMBB_DIR, "PMBB-Release-2024-3.0_phenotype_person.txt")

cov <- fread(COV, header = TRUE)
flags <- fread(file.path("ss", "rgcname_pmbbid_metadata_flags.csv"))
key <- read_excel(file.path("ss", "BW-Regeneron_Workbook_20250627.xlsx"), sheet = "CREP_workbook")
ss <- read_excel(file.path("ss", "br_pts_for_exwas_08292025.xlsx"))

# key$UP_ID <- sprintf("UP%04d", as.integer(key$Participant))
# sum(key$UP_ID == key$VCFID)
# it's all the same

### ANALYSIS W/ FLAGS SS ONLY ###
up_pattern <- "^UPENN-PMBB_UP[0-9]+_UP[0-9]+$"
kt_pattern <- "^UPENN-PMBB_KT[0-9]+_KT[0-9]+$"

# find all cases that match pattern among three columns, filter for those w/ this ID
up <- flags %>%
    mutate(match_col = apply(select(., RGC_sample_name, ID1, ID2), 1, function(row) {
        m <- row[grepl(up_pattern, row)]
        if (length(m) > 0) m[1] else NA
    })) %>%
    filter(!is.na(match_col))
print(dim(up))
# 2864   14

kt <- flags %>%
    mutate(match_col = apply(select(., RGC_sample_name, ID1, ID2), 1, function(row) {
        m <- row[grepl(kt_pattern, row)]
        if (length(m) > 0) m[1] else NA
    })) %>%
    filter(!is.na(match_col))
print(dim(kt))
# 900

up2 <- select(up, c("PMBB_ID",  "RGC_sample_name", "ID1", "ID2", "CREP", "match_col"))
up2$VCFID <- sub(".*(UP[0-9]{4}).*", "\\1", up2$match_col)
# flags2$VCFID <- sub(".*(KT[0-9]{4}).*", "\\1", flags2$match_col)

up3 <- select(up2, c("PMBB_ID",  "match_col", "VCFID"))
up3$SampNum <- numbers <- as.numeric(gsub("UP", "", up3$VCFID))
# 2864 samples in PMBB that have UP ID

# ========================
# OVERALL CREP ANALYSIS - can skip this...
# ========================
### PMBB COVARIATES FILE ###
crep <- cov %>% filter(CREP_HighRisk_Flag == 1)
print(dim(crep))

## are crep samples from progeny and pmbb the same?
crep_pmbb <- crep$person_id
crep_progeny <- flags3$PMBB_ID
length(crep_pmbb) #2833
length(crep_progeny) #2864

sum(crep_pmbb %in% crep_progeny) #2833, so all of the pmbb are in progeny
sum(!(crep_progeny %in% crep_pmbb)) #31 in progeny but not in pmbb
sum(!(crep_pmbb %in% crep_progeny)) #0 in pmbb but not in progeny

# venn.plot <- draw.pairwise.venn(
#     area1 = length(crep_pmbb),
#     area2 = length(crep_progeny),
#     cross.area = length(intersect(crep_pmbb, crep_progeny)),
#     category = c("PMBB CREP", "Flags CREP"),
#     fill = c("skyblue", "pink"),
#     lty = "blank",
#     cex = 1.5,
#     cat.cex = 1.5,
#     cat.pos = c(-20, 20)
# )
# ggsave(file.path("plots", "20250829_venn_crep_pmbb_progeny.png"), venn.plot)

progeny_not_pmbb <- setdiff(crep_progeny, crep_pmbb)
print(progeny_not_pmbb)

# [1] "PMBB1128910863353" "PMBB1203302743184" "PMBB1230383977275" "PMBB1422030088258" "PMBB1534550830439"
# [6] "PMBB1970453501303" "PMBB2025997391657" "PMBB2027093284064" "PMBB2388005759286" "PMBB2786896176432"
# [11] "PMBB3662861030833" "PMBB3704586042027" "PMBB3917572865526" "PMBB4195811438392" "PMBB4244298090288"
# [16] "PMBB4425729784182" "PMBB4811155336791" "PMBB4880769461608" "PMBB4984126299852" "PMBB5162255907663"
# [21] "PMBB6394380060914" "PMBB6447824256436" "PMBB7130575788188" "PMBB7152905566278" "PMBB7458355573554"
# [26] "PMBB7990632775668" "PMBB8516410259866" "PMBB9098182828989" "PMBB9210220031010" "PMBB9520609830734"
# [31] "PMBB9731498790294"

all_pmbb <- cov$person_id
sum(progeny_not_pmbb %in% all_pmbb)
# all of these are def in the PMBB

# get a list of those i'm pretty sure are in progeny? 2833
crep_sure <- intersect(crep_pmbb, crep_progeny)
length(crep_sure)

# ========================
# PROGENY ANALYSIS
# ========================
# No - not in EHR
# Yes - HUP, PennMed, EHR - ICD codes
# but NOT in ICD codes...

### ICD CODE ANALYSIS ###
icds <- unique(read.table(file.path("data", "breast_simplexo3_v2_cancer_filtered_patients_ids.txt"))$V1)
print(length(unique(icds))) # 3729
icds_df <- as.data.frame(icds) # for manual search

### EXAMINE DUPLICATES IN PROGENY
sum(duplicated(ss$SampNum))
dups <- ss[ss$SampNum %in% ss$SampNum[duplicated(ss$SampNum)], ]
# has entry for each breast so that's why

### MERGE DUPLICATE ROWS
ss2 <- ss %>%
    group_by(SampNum) %>%
    summarise(
        across(-any_of("SampNum"), ~ paste(unique(as.character(.x)), collapse = "; ")),
        .groups = "drop"
    )
dim(ss2) #4273 PATIENTS

# ========================
# FINAL CASE LIST
# ========================

# PROGENY IN THE PMBB
progeny_pmbb <- merge(ss2, up3, by = "SampNum")
# merge w/ covariates to check the gender/sex
check <- merge(progeny_pmbb, cov, by.x = "PMBB_ID", by.y = "person_id")
dim(check)
# 1561   40

check$Gender2 <- ifelse(check$Gender == "F", "Female", "Male")
(check$Gender2 == check$Sequenced_gender) # they all match!

progeny_pmbb1 <- check %>% filter(Gender == "F")
progeny_pmbb2 <- check %>% filter(Sequenced_gender == "Female")
sum(progeny_pmbb$PMBB_ID %in% progeny_pmbb2$PMBB_ID)

final_ids <- sort(progeny_pmbb2$PMBB_ID)
ages <- progeny_pmbb2 %>% select(PMBB_ID, CaDxAge)
write.table(final_ids, file.path("data", "progeny_ids.txt"), row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(ages, file.path("data", "pmbb_progeny_ages.txt"), row.names = FALSE, col.names = FALSE, quote = FALSE)

# get the case







# GET WENTING'S NOT DETECTED BY ICD
progeny_not_icd <- progeny_pmbb %>% filter(!(PMBB_ID %in% icds))
dim(progeny_not_icd) #600
table(progeny_not_icd$Seen_In_CREP)
# HUP CREP       NA       No  Unknown      Yes
# 2       32      390       21      155

write_xlsx(progeny_not_icd, file.path("ss", "20250901_progeny_not_icd.xlsx"))

# GET WENTING'S DETECTED BY PROGENY
progeny_icd <- progeny_pmbb %>% filter(PMBB_ID %in% icds)
dim(progeny_icd) #961
table(progeny_icd$Seen_In_CREP)
# Chester County Hospital                   HUP CREP
# 2                         35
# NA                         No
# 14                        216
# Penn Medicine Valley Forge      Pennsylvania Hospital
# 1                          4
# Unknown                        Yes
# 1                        688


# just merge my ICD visits onto hers
# f <- fread(file.path("data", "breast_simplexo3_v2_cancer_all_patients.txt"))
# f2 <- f %>% select(c("person_id", "num_visits", "CREP_HighRisk_Flag", "Batch", "all_matching_codes", "Class", "Sample_age", "Sample_date"))
# colnames(f2) <- c("PMBB_ID", "ICD_Num_Instances", "PMBB_CREP_HighRisk_Flag", "Freeze", "ICD_codes", "Genetic_Ancestry", "Sample_Age", "Sample_Date")
# final <- merge(progeny_pmbb, f2, by = "PMBB_ID")
# dim(final)
# final2 <- final %>% filter(ICD_Num_Instances >= 2)
# dim(final2)

progeny_pmbb_ehr <- progeny_pmbb %>% filter(!is.na(Seen_In_CREP)) %>%
    filter(Seen_In_CREP != "No") %>%
    filter(Seen_In_CREP != "NA") %>%
    filter(Seen_In_CREP != "Unknown")
print(dim(progeny_pmbb_ehr)) #877 in the EHR (CREP, Penn Med, etc.)
unique(progeny_pmbb_ehr$Seen_In_CREP)

progeny_pmbb_ehr_not_icd <- progeny_pmbb_ehr %>% filter(!(PMBB_ID %in% icds))
unique(progeny_pmbb_ehr_not_icd$Seen_In_CREP)
print(dim(progeny_pmbb_ehr_not_icd))
write_xlsx(progeny_pmbb_ehr_not_icd, file.path("ss", "20250901_progeny_ehr_not_icd.xlsx"))

icds <- unique(read.table(file.path("data", "breast_simplexo3_v2_cancer_filtered_patients_ids.txt"))$V1)



# GET ICD W/ RGC ID BUT NOT IN PROGENY
icd_only <- flags3 %>% filter(PMBB_ID %in% icd_not_progeny)
# 961 icds in Progeny_PMBB

# GET ICD
progeny_pmbb_not_ehr <- progeny_pmbb %>%
    filter(is.na(Seen_In_CREP) | Seen_In_CREP %in% c("No", "NA", "Unknown"))
dim(progeny_pmbb_not_ehr) #674  28

progeny_pmbb_not_ehr_not_icd <- progeny_pmbb_not_ehr %>% filter(!(PMBB_ID %in% icds))
dim(progeny_pmbb_not_ehr_not_icd) #443  28

progeny_pmbb_not_ehr_icd <- progeny_pmbb_not_ehr %>% filter((PMBB_ID %in% icds))
dim(progeny_pmbb_not_ehr_icd) #231  28
#

# sum(icd_not_progeny %in% crep_pmbb) #47
# sum(icd_not_progeny %in% crep_progeny) #47
# sum(icd_not_progeny %in% f2_ids) #811

# f2 <- cov %>% filter(Batch == 2)
# f2_ids <- f2$person_id

# filter for those in freeze 2.0?
## get anybody who can't be detected by ICD code in freeze 2.0? or are in Seen_In_CREP

dim(df)
progeny_case <- unique(df$PMBB_ID)
length(progeny_case)

# filter for those only with Yes in CREP, if it's NA or NO, get rid of them
# df3 <- df2 %>% filter(Seen_In_CREP != "Yes") %>% filter(Seen_In_CREP != "HUP CREP")

# how many at only penn medicine sites aren't covered by the ICD codes?
ss_pm <- ss2 %>% filter(!is.na(Seen_In_CREP)) %>% filter(Seen_In_CREP != "No") %>% filter(Seen_In_CREP != "Unknown")
# %>% filter(Seen_In_CREP != "Yes") %>% filter(Seen_In_CREP != "HUP CREP")
df <- merge(ss_pm, flags3, by = "SampNum")
unique(df$Seen_In_CREP)
sum(df$PMBB_ID %in% icds)
icd_not_progeny <- setdiff(df$PMBB_ID, icds)
length(icd_not_progeny)
# 189 detected by icd codes but not in progeny


unique(ss$Seen_In_CREP)
progeny_case_crep <- unique(df2$PMBB_ID)
length(progeny_case_crep)


length(icds) #4112
# overlap icd with crep
pmbb_crep_icds <- unique(intersect(icds, crep_pmbb))
length(pmbb_crep_icds)
# 1103 crep cases w/ breast cancer identified with icd codes

length(progeny_case) # 1561
sum(progeny_case %in% pmbb_crep_icds) #1016 are found in progeny and are also found in ICD codes, specifically CREP
sum(!(progeny_case %in% pmbb_crep_icds)) #545 are found in progeny but not in ICD codes, specifically CREP

sum(progeny_case %in% icds) #1042 are found in progeny and are also found in ICD codes, not just CREP flag so maybe older?
sum(!(progeny_case %in% icds)) #519 are found in progeny and not in ICD codes, not just CREP flag so maybe older?
# examine these differences??

progeny_not_icd <- setdiff(progeny_case, pmbb_crep_icds)
icd_not_progeny <- setdiff(pmbb_crep_icds, progeny_case)
length(progeny_not_icd) #545
length(icd_not_progeny) #87

venn.plot <- draw.pairwise.venn(
    area1 = length(progeny_case),
    area2 = length(pmbb_crep_icds),
    cross.area = length(intersect(progeny_case, pmbb_crep_icds)),
    category = c("Progeny", "ICD"),
    fill = c("skyblue", "pink"),
    lty = "blank",
    cex = 1.5,
    cat.cex = 1.5,
    cat.pos = c(-20, 20)
)

### EXAMINE THOSE FOUND IN ICD BUT NOT PROGENY
icd_codes <- c("^C50", "^Z85.4", "^174", "^V10.3")
pattern <- paste(icd_codes, collapse = "|")
occ2 <- occ %>% filter(person_id %in% icd_not_progeny)
occ3 <- occ2[
    grepl(pattern, condition_source_value, perl = TRUE, ignore.case = TRUE)
]

# are these even in the progeny db at all?
head(icd_not_progeny)


