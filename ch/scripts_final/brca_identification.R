# ========================
# IDENTIFYING GBRCA12 CARRIERS IN PMBB
# ========================
library(here)
setwd("/Users/jennyzli/Documents/Nathanson")
source(here("R", "config.R"))
source(here("ch", "scripts", "brca_qc.R"))
library(gridExtra)
library(grid)

# ========================
# READ DATA
# ========================
### newest BRCA1/2 annotations
all_b1 <- read.csv(file.path("ch", "data", "brca_vep3.BRCA1.vep.report.csv"))
all_b2 <- read.csv(file.path("ch", "data", "brca_vep3.BRCA2.vep.report.csv"))

COV <- here("PMBB", "3.0", "PMBB-Release-2024-3.0_covariates.txt")
PER <- here("PMBB", "3.0", "PMBB-Release-2024-3.0_phenotype_person.txt")
cov <- fread(COV, header = TRUE)
person <- fread(PER, header = TRUE)

flags <- read_excel(file.path("PMBB", "pmbb_flag_tables.xlsx"))
consent <- flags %>% filter(PMBB_Consent == 1)
consent_ids <- consent$PMBB_ID
length(consent_ids)
# 53062

up <- read.csv(here("simplexo", "data", "simplexo_up_map.csv"))
up$SampNum <- as.numeric(up$SampNum)

rgc <- read.csv(here("PMBB", "rgcname_pmbbid_metadata_flags.csv"))
crep <- rgc %>% filter(CREP == 1)
crep_ids <- crep$PMBB_ID

# ========================
# FREEZE 2
# ========================
f2_brca1 <- read_excel(here("ch", "ss", "BRCA1.BRCA2 P.LP_04.11.21.xlsx"), sheet = "BRCA1")
f2_brca2 <- read_excel(here("ch", "ss", "BRCA1.BRCA2 P.LP_04.11.21.xlsx"), sheet = "BRCA2")

### INITIAL STATS
length(unique(f2_brca1$SampleID))
length(unique(f2_brca2$SampleID))
# 152
# 263

range(f2_brca1$ALT_AlleleFrac)
range(f2_brca2$ALT_AlleleFrac)
# > range(f2_brca1$ALT_AlleleFrac)
# [1] 0.259 0.684
# > range(f2_brca2$ALT_AlleleFrac)
# [1] 0.185 0.633

range(f2_brca1$ALT_AlleleDepth)
range(f2_brca2$ALT_AlleleDepth)
# > range(f2_brca1$ALT_AlleleDepth)
# [1]  6 75
# > range(f2_brca2$ALT_AlleleDepth)
# [1]  4 60

range(f2_brca1$Total_Depth)
range(f2_brca2$Total_Depth)
# > range(f2_brca1$Total_Depth)
# [1]  17 146
# > range(f2_brca2$Total_Depth)
# [1]  19 111

### THOSE WITH BOTH
length(intersect(f2_brca1$SampleID, f2_brca2$SampleID))
# 1 with both

### TOTAL SAMPLES
length(union(f2_brca1$SampleID, f2_brca2$SampleID))
# 414 total samples

# ========================
# FREEZE 3.0
# ========================
### FOUND ###
fnd_vars <- read_excel(here("ch", "ss", "BW-Regeneron_Workbook_20250627.xlsx"), sheet = "CREP.found_variants") %>%
    filter(Mutation_Gene1 %in% c("BRCA1", "BRCA2")) %>%
    filter(!is.na(Sample.ID)) %>%
    left_join(up, by = c("Sample.ID" = "VCFID")) %>%
    mutate(variant_key = paste0("('", PMBB_ID, "', '", HGVSc, "')"))

fnd_lines <- read_excel(here("ch", "ss", "BW-Regeneron_Workbook_20250627.xlsx"), sheet = "CREP.found_lines") %>%
    filter(Gene %in% c("BRCA1", "BRCA2")) %>%
    filter(!is.na(Sample.ID)) %>%
    left_join(up, by = c("Sample.ID" = "match_col")) %>%
    mutate(variant_key = paste0("('", PMBB_ID, "', '", HGVSc, "')"))
f3_b12_fnd <- fnd_lines %>% filter(variant_key %in% fnd_vars$variant_key) %>%
    filter(Variant.LoF_level != 4)
dim(f3_b12_fnd)
# 1122

table(f3_b12_fnd$Gene)
# BRCA1 BRCA2
# 578   544

### LIKELY ###
lik_vars <- read_excel(here("ch", "ss", "BW-Regeneron_Workbook_20250627.xlsx"), sheet = "CREP.likely_variants") %>%
    filter(Mutation_Gene1 %in% c("BRCA1", "BRCA2")) %>%
    filter(!is.na(Sample.ID)) %>%
    left_join(up, by = c("Sample.ID" = "VCFID")) %>%
    mutate(variant_key = paste0("('", PMBB_ID, "', '", HGVSc, "')"))
lik_lines <- read_excel(here("ch", "ss", "BW-Regeneron_Workbook_20250627.xlsx"), sheet = "CREP.likely_lines") %>%
    filter(Gene %in% c("BRCA1", "BRCA2")) %>%
    filter(!is.na(Sample.ID)) %>%
    left_join(up, by = c("Sample.ID" = "match_col")) %>%
    mutate(variant_key = paste0("('", PMBB_ID, "', '", HGVSc, "')"))
f3_b12_lik <- lik_lines %>% filter(variant_key %in% lik_vars$variant_key) %>%
    filter(Variant.LoF_level != 4)
dim(f3_b12_lik)
# 32

table(f3_b12_lik$Gene)
# BRCA1 BRCA2
# 16 17

### UNDETECTABLE ###
f3_b12_und <- read_excel(here("ch", "ss", "BW-Regeneron_Workbook_20250627.xlsx"), sheet = "CREP.undetectable_variants") %>%
    filter(Mutation_Gene1 %in% c("BRCA1", "BRCA2")) %>%
    filter(!is.na(Sample.ID)) %>%
    left_join(up, c("Sample.ID" =  "VCFID")) %>%
    filter(!is.na(PMBB_ID))
dim(f3_b12_und)
# 71

table(f3_b12_und$Mutation_Gene1)
# BRCA1 BRCA2
# 65     6

### POSSIBLE / PUTATIVE -> RESOLVED  ###
f3_b12_pos_mm <- read_excel(file.path("ch", "ss", "brca12_palb2_progenymatch_070125_w_family_source.xlsx"), sheet = "Mutation Match") %>%
    select("Sample.ID", "WSL ID", "Mutation_GeneX", "HGVScX", "HGVSpX") %>%
    filter(Mutation_GeneX %in% c("BRCA1", "BRCA2"))
dim(f3_b12_pos_mm)
# 31

f3_b12_pos_ml <- read_excel(file.path("ch", "ss", "brca12_palb2_progenymatch_070125_w_family_source.xlsx"), sheet = "Most Likely a Match") %>%
    select("Sample.ID", "WSL ID", "Mutation_GeneX", "HGVScX", "HGVSpX") %>%
    filter(Mutation_GeneX %in% c("BRCA1", "BRCA2"))
dim(f3_b12_pos_ml)
# 4

f3_b12_pos_gm <- read_excel(file.path("ch", "ss", "brca12_palb2_progenymatch_070125_w_family_source.xlsx"), sheet = "Gene Match -Yes-Unconfirmed") %>%
    select("Sample.ID", "WSL ID", "Mutation_GeneX", "HGVScX", "HGVSpX") %>%
    filter(Mutation_GeneX %in% c("BRCA1", "BRCA2"))
dim(f3_b12_pos_gm)
# 104

f3_b12_pos <- rbind(f3_b12_pos_mm, f3_b12_pos_ml, f3_b12_pos_gm) %>%
    mutate(variant_key = paste0("('", Sample.ID, "', '", HGVScX, "')"))
dim(f3_b12_pos)
# 139

# attach variant annotation
f3_b12_pos_anno <- read_excel(here("ch", "ss", "BW-Regeneron_Workbook_20250627.xlsx"), sheet = "CREP.possible_lines") %>%
    filter(Gene %in% c("BRCA1", "BRCA2")) %>%
    filter(!is.na(Sample.ID)) %>%
    left_join(up, c("Sample.ID" =  "match_col")) %>%
    filter(!is.na(PMBB_ID)) %>%
    filter(Variant.LoF_level != 4)
dim(f3_b12_pos_anno) # 196
colnames(f3_b12_pos_anno)
f3_b12_res <- f3_b12_pos_anno %>% filter(variant_key %in% f3_b12_pos$variant_key)
dim(f3_b12_res)
# 136

table(f3_b12_res$Gene)
# BRCA1 BRCA2
# 61    75

### COMBINE ALL FREEZE 3
# first merge the ones with annotations
f3_b12_lik_fnd_pos <- rbind(f3_b12_lik, f3_b12_fnd, f3_b12_res) %>%
    filter(Variant.LoF_level != 4)
dim(f3_b12_lik_fnd_pos) # 1290
table(f3_b12_lik_fnd_pos$Gene)
# BRCA1 BRCA2
# 655   635

### QC STATS
range(f3_b12_lik_fnd_pos$Sample.Depth) # 14 117
range(f3_b12_lik_fnd_pos$Sample.AltDepth) # 4 59
range(f3_b12_lik_fnd_pos$Sample.AltFrac) # 0.179 0.792
length(unique(f3_b12_lik_fnd_pos$Sample.ID)) # 1283
(unique(f3_b12_lik_fnd_pos$Variant.Consequence))
table(f3_b12_lik_fnd_pos$Variant.LoF_level)
# 1    2
# 1274   16

write.csv(f3_b12_lik_fnd_pos, file.path("ch", "data", "crep_brca12_lik_fnd_pos_df.csv"))

# combine annotations with undetected
f3_b1_lik_pos_fnd <- f3_b12_lik_fnd_pos %>% filter(Gene == "BRCA1")
f3_b2_lik_pos_fnd <- f3_b12_lik_fnd_pos %>% filter(Gene == "BRCA2")

f3_b1_und <- f3_b12_und %>% filter(Mutation_Gene1 == "BRCA1")
f3_b2_und <- f3_b12_und %>% filter(Mutation_Gene1 == "BRCA2")

f3_b1_ids <- sort(unique(c(f3_b1_lik_pos_fnd$PMBB_ID, f3_b1_und$PMBB_ID)))
f3_b2_ids <- sort(unique(c(f3_b2_lik_pos_fnd$PMBB_ID, f3_b2_und$PMBB_ID)))
length(f3_b1_ids) # 719
length(f3_b2_ids) # 641

length(intersect(f3_b1_ids, f3_b2_ids))
# 8 with both

length(unique(union(f3_b1_ids, f3_b2_ids)))
# 1352 total

# ========================
# F2 + F3
# ========================
f23_b1_ids <- sort(unique(c(f3_b1_ids, f2_brca1$SampleID)))
f23_b2_ids <- sort(unique(c(f3_b2_ids, f2_brca2$SampleID)))

write.csv(f23_b1_ids, file.path("ch", "data", "pmbb_brca1_prev.csv"), row.names = FALSE)
write.csv(f23_b2_ids, file.path("ch", "data", "pmbb_brca2_prev.csv"), row.names = FALSE)

length(f23_b1_ids)
length(f23_b2_ids)
# > length(f23_b1_ids)
# [1] 871
# > length(f23_b2_ids)
# [1] 904

length(union(f23_b2_ids, f23_b1_ids))
# 1766
length(intersect(f23_b2_ids, f23_b1_ids))
# 9

# ========================
# MY ANNOTATIONS
# ========================
all_b1o <- process_brca_data(all_b1, "BRCA1", "all", 0.179, 0.8, 14, 4, 10)
all_b2o <- process_brca_data(all_b2, "BRCA2", "all", 0.179, 0.8, 14, 4, 10)

### DIAGNOSTIC PLOTS
create_diagnostic_pdf(
    all_b1o$step1, all_b1o$step2, all_b1o$step3,
    "BRCA1 (All PMBB)",
    file.path("ch", "figures", "brca1_all_qc_histograms.pdf")
)

create_diagnostic_pdf(
    all_b2o$step1, all_b2o$step2, all_b2o$step3,
    "BRCA2 (All PMBB)",
    file.path("ch", "figures", "brca2_all_qc_histograms.pdf")
)

write.csv(all_b1o$step3, file.path("ch", "data", "brca1_all_filtered.csv"), row.names = FALSE)
write.csv(all_b2o$step3, file.path("ch", "data", "brca2_all_filtered.csv"), row.names = FALSE)

all_b1o_lof12 <- all_b1o$step3 %>% filter(Variant.LoF_level %in% c(1, 2), AutoGVP %in% c("Pathogenic", "Likely_pathogenic"))
all_b2o_lof12 <- all_b2o$step3 %>% filter(Variant.LoF_level %in% c(1, 2), AutoGVP %in% c("Pathogenic", "Likely_pathogenic"))

all_b1_lof12_ids <- unique(all_b1o_lof12$Sample.ID)
all_b2_lof12_ids <- unique(all_b2o_lof12$Sample.ID)

length(all_b1_lof12_ids)
# 890
length(all_b2_lof12_ids)
# 1020

length(union(all_b1_lof12_ids, all_b2_lof12_ids))
# 1899

length(intersect(all_b1_lof12_ids, all_b2_lof12_ids))
# 11

# ========================
# COMBINE WITH EVERYTHING
# ========================
final_b1_ids <- sort(union(f23_b1_ids, all_b1_lof12_ids))
final_b2_ids <- sort(union(f23_b2_ids, all_b2_lof12_ids))
final_b12_ids <- sort(union(final_b1_ids, final_b2_ids))
final_b12_both_ids <- sort(intersect(final_b1_ids, final_b2_ids))

length(final_b1_ids)
# 965
length(final_b2_ids)
# 1039
length(final_b12_ids)
# 1992
length(final_b12_both_ids)
# 12

write.csv(final_b1_ids, file.path("ch", "data", "pmbb_brca1_case_ids.csv"), row.names = FALSE)
write.csv(final_b2_ids, file.path("ch", "data", "pmbb_brca2_case_ids.csv"), row.names = FALSE)
write.csv(final_b12_ids, file.path("ch", "data", "pmbb_brca12_case_ids.csv"), row.names = FALSE)

# ========================
# COMPARISON WITH OTHER CALLING
# ========================
### INTERSECTIONS
length(f23_b1_ids) #871
length(intersect(f23_b1_ids, all_b1_lof12_ids)) # 796

length(f23_b2_ids) # 904
length(intersect(f23_b2_ids, all_b2_lof12_ids)) # 885

### DIFFERENCES (w/ LOF 1 or 2)
undet_b1_ids <- setdiff(f23_b1_ids, all_b1_lof12_ids)
undet_b2_ids <- setdiff(f23_b2_ids, all_b2_lof12_ids)
undet_b1 <- all_b1o$step3 %>%
    filter(Sample.ID %in% undet_b1_ids, !(Sample.ID %in% na.omit(f3_b12_und$PMBB_ID)))
undet_b2 <- all_b2o$step3 %>%
    filter(Sample.ID %in% undet_b2_ids, !(Sample.ID %in% f3_b12_und$PMBB_ID))
dim(undet_b1)
# 26
dim(undet_b2)
# 26??

new_b1_ids <- setdiff(all_b1_lof12_ids, f23_b1_ids)
new_b2_ids <- setdiff(all_b2_lof12_ids, f23_b2_ids)
length(new_b1_ids) # 94
length(new_b2_ids) # 135

new_b1 <- all_b1o_lof12 %>%
    filter(Sample.ID %in% new_b1_ids) %>%
    left_join(cov %>% select(person_id, Batch), by = c("Sample.ID" = "person_id")) %>%
    left_join(up, by = c("Sample.ID" = "PMBB_ID")) %>%
    mutate(is_consented = Sample.ID %in% consent_ids) %>%
    mutate(crep = Sample.ID %in% crep_ids)
new_b2 <- all_b2o_lof12 %>%
    filter(Sample.ID %in% new_b2_ids) %>%
    left_join(cov %>% select(person_id, Batch), by = c("Sample.ID" = "person_id")) %>%
    left_join(up, by = c("Sample.ID" = "PMBB_ID")) %>%
    mutate(is_consented = Sample.ID %in% consent_ids) %>%
    mutate(crep = Sample.ID %in% crep_ids)
new <- rbind(new_b2, new_b1)

dim(new_b1)
# 94
dim(new_b2)
# 135

length(intersect(f23_b1_ids, all_b1_lof12_ids)) # 796

write.csv(new_b1, file.path("ch", "data", "brca1_undiscovered_df.csv"), row.names = FALSE)
write.csv(new_b2, file.path("ch", "data", "brca2_undiscovered_df.csv"), row.names = FALSE)
write.csv(new, file.path("ch", "data", "brca12_undiscovered_df.csv"), row.names = FALSE)

# Quick summary
table(new$Batch)
table(new$is_consented)
table(new$crep)

# 1   2
# 27 202
# > table(new$is_consented)
#
# FALSE  TRUE
# 98   131
# > table(new$crep)
#
# FALSE  TRUE
# 154    75



