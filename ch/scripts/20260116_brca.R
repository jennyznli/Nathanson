library(here)
setwd("/Users/jennyzli/Documents/Nathanson")
source(here("R", "config.R"))
source(here("ch", "scripts", "brca_qc2.R"))

library(ggplot2)
library(gridExtra)
library(grid)

# ========================
# READ DATA
# ========================
all_b1 <- read.csv(file.path("ch", "data", "brca_vep3.BRCA1.vep.report.csv"))
all_b2 <- read.csv(file.path("ch", "data", "brca_vep3.BRCA2.vep.report.csv"))

# ========================
# INITIAL EDA
# ========================
## BRCA1
dim(all_b1)
# [1] 7267528      62

table(all_b1$MANE_SELECT)
# NM_000988.5 NM_001330230.2 NM_001363254.2    NM_001661.4    NM_005440.5    NM_005899.5    NM_006373.4
# 103412            156          43915           1400         513801         388024         937441
# NM_007294.4    NM_145041.4
# 2413614        2865765

x <- all_b1 %>% filter(Gene == "BRCA1")
table(x$MANE_SELECT)
# NM_007294.4
# 2413614

###  BRCA2
dim(all_b2)
# [1] 11498648       62

table(all_b2$MANE_SELECT)
# .    NM_000059.4 NM_001136571.2    NM_014887.3    NM_023037.3    NM_052818.3
# 431255         814542          36101        2692967         221151        7302632

y <- all_b2 %>% filter(Gene == "BRCA2")
table(y$MANE_SELECT)
# NM_000059.4

# ========================
# PROCESS DATA
# ========================
all_b1o <- process_brca_data(all_b1, "BRCA1", "all3", 0.3, 0.7, 20, 3, 20)
all_b2o <- process_brca_data(all_b2, "BRCA2", "all3", 0.3, 0.7, 20, 3, 20)

dim(all_b1o$step3)
# 119718     62
dim(all_b2o$step3)
# 78870    62
# ========================
# GENERATE DIAGNOSTIC PLOTS
# ========================
create_diagnostic_pdf(
    all_b1o$step1, all_b1o$step2, all_b1o$step3,
    "BRCA1 (All PMBB)",
    file.path("ch", "figures", "brca1_all3_qc_histograms.pdf")
)

create_diagnostic_pdf(
    all_b2o$step1, all_b2o$step2, all_b2o$step3,
    "BRCA2 (All PMBB)",
    file.path("ch", "figures", "brca2_all3_qc_histograms.pdf")
)

# ========================
# SAVE FILTERED DATA
# ========================
write.csv(all_b1o$step3, file.path("ch", "data", "brca1_all3_filtered.csv"), row.names = FALSE)
write.csv(all_b2o$step3, file.path("ch", "data", "brca2_all3_filtered.csv"), row.names = FALSE)

# all_b1o_path <- all_b1o$step3 %>% filter(!grepl("benign", ClinVar.SIG, ignore.case = TRUE), AM.class != "likely_benign") %>% select(Sample.ID)
# all_b2o_path <- all_b2o$step3 %>% filter(!grepl("benign", ClinVar.SIG, ignore.case = TRUE), AM.class != "likely_benign") %>% select(Sample.ID)

all_b1o_path <- all_b1o$step3 %>% filter(Variant.LoF_level == 1) %>% select(Sample.ID)
all_b2o_path <- all_b2o$step3 %>% filter(Variant.LoF_level == 1) %>% select(Sample.ID)

dim(all_b1o_path) #836  58
dim(all_b2o_path) #983 58
all_b1o_ids <- all_b1o_path$Sample.ID
all_b2o_ids <- all_b2o_path$Sample.ID
# #
# write.csv(all_b1o$step3, file.path("ch", "data", "brca1_all3_filtered.csv"), row.names = FALSE)
# write.csv(all_b2o$step3, file.path("ch", "data", "brca2_all3_filtered.csv"), row.names = FALSE)

# ========================
# COMPARE W PREVIOUS FREEZE 2
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

table(f2_brca1$ExonicFunc.refGene)
# .  frameshift_deletion frameshift_insertion    nonsynonymous_SNV             stopgain
# 9                   77                   34                   22                   11

table(f2_brca1$GenomicRegion.refGene)
# exonic splicing
# 144        9

table(f2_brca1$CLNSIG)
# likely_patho, patho, or both

### FILTERING
f2_brca1_fil <- f2_brca1 %>% filter(Total_Depth > 20, ALT_AlleleDepth > 5, ALT_AlleleFrac > 0.3, ALT_AlleleFrac < 0.7)
f2_brca2_fil <- f2_brca2 %>% filter(Total_Depth > 20, ALT_AlleleDepth > 5, ALT_AlleleFrac > 0.3, ALT_AlleleFrac < 0.7)

f2_b1_ids <- unique(f2_brca1_fil$SampleID)
f2_b2_ids <- unique(f2_brca2_fil$SampleID)
length(f2_b1_ids)
length(f2_b2_ids)
# 148
# 256

# intersection w/o any filtering
length(intersect(all_b1$Sample.ID, f2_b1_ids))
length(intersect(all_b2$Sample.ID, f2_b2_ids))
# 147 ...
# 255 ...

# intersection @ step 1
length(intersect(all_b1o$step1$Sample.ID, f2_b1_ids))
length(intersect(all_b1o$step1$Sample.ID, f2_b2_ids))
# 147
# 221
# so did they not get rid of all the same ones? wihtout GQ...

# intersection @ step 2
length(intersect(all_b1o$step2$Sample.ID, f2_b1_ids))
length(intersect(all_b1o$step2$Sample.ID, f2_b2_ids))
# 146
# 187

# intersection @ step 2
length(intersect(all_b1o$step3$Sample.ID, f2_b1_ids))
length(intersect(all_b1o$step3$Sample.ID, f2_b2_ids))
# 146
# 163

# ========================
# COMPARE WITH BRAD'S f3
# ========================
f3_df <- read_excel(here("ch", "ss", "BW-Regeneron_Workbook_20250627.xlsx"), sheet = "CREP_workbook") %>% filter(Mutation_Gene1 %in% c("BRCA1", "BRCA2"))
up <- read.csv(here("simplexo", "data", "simplexo_up_map.csv"))
f3_df1 <- f3_df %>% left_join(up, by = "VCFID") %>% filter(!is.na(PMBB_ID))
f3_df <- f3_df1
f3_ids <- f3_df1$PMBB_ID
length(f3_ids)
# 1265

table(f3_df1$Mutation_Gene1)
# BRCA1 BRCA2
# 674   591

f3_fnd <- read_excel(here("ch", "ss", "BW-Regeneron_Workbook_20250627.xlsx"), sheet = "CREP.found_lines") %>% filter(Gene %in% c("BRCA1", "BRCA2")) %>% filter(!is.na(Sample.ID))
f3_possible <- read_excel(here("ch", "ss", "BW-Regeneron_Workbook_20250627.xlsx"), sheet = "CREP.possible_lines") %>% filter(Gene %in% c("BRCA1", "BRCA2")) %>% filter(!is.na(Sample.ID))

f3_fnd <- f3_fnd %>% filter(Sample.ID %in% f3_df$match_col)
f3_possible <- f3_possible %>% filter(Sample.ID %in% f3_df$match_col) %>% select(-variant_key)
f3 <- rbind(f3_fnd, f3_possible)
f3 <- f3 %>% left_join(up, by = c("Sample.ID" = "match_col"))
# "stop_gained"                                        "frameshift_variant&splice_region_variant"
# [3] "splice_donor_variant"                               "frameshift_variant"
# [5] "missense_variant"                                   "splice_polypyrimidine_tract_variant&intron_variant"
# [7] "missense_variant&splice_region_variant"             "inframe_deletion"
# [9] "splice_donor_region_variant&intron_variant"         "splice_acceptor_variant"
# [11] "stop_gained&splice_region_variant"                  "splice_donor_variant&coding_sequence_variant"
# [13] "start_lost"                                         "splice_donor_5th_base_variant&intron_variant"
# [15] "synonymous_variant"                                 "stop_gained&inframe_insertion"
# [17] "intron_variant"                                     "splice_region_variant&synonymous_variant"
# [19] "stop_gained&frameshift_variant"
# #
### INITIAL STATS
range(f3$Sample.Depth) # 14 117
range(f3$Sample.AltDepth) # 4 55
range(f3$Sample.AltFrac) # 0.179 0.792
length(unique(f3$Sample.ID)) # 602
length(unique(f3 %>% filter(Gene == "BRCA1"))$Sample.ID) # 604
length(unique(f3 %>% filter(Gene == "BRCA2"))$Sample.ID) # 586

### AFTER FILTERING
f3_fil <- f3 %>% filter(Sample.Depth >= 15, Sample.AltDepth >= 4, Sample.AltFrac >= 0.2, Sample.AltFrac <= 0.8)
f3_fil <- f3_fil %>% left_join(up, by = c("Sample.ID" = "match_col"))

range(f3_fil$Sample.Depth) # 15 117
range(f3_fil$Sample.AltDepth) # 4 55
range(f3_fil$Sample.AltFrac) # 0.207 0.792

f3_b1_fil <- f3_fil %>% filter(Gene == "BRCA1")
f3_b2_fil <- f3_fil %>% filter(Gene == "BRCA2")
f3_b1_ids <- unique(f3_b1_fil$PMBB_ID.x)
f3_b2_ids <- unique(f3_b2_fil$PMBB_ID.x)

# INTERSECTING...
# with master spreadsheet
length(intersect(unique(f3_df %>% filter(Mutation_Gene1 == "BRCA1"))$PMBB_ID, all_b1o$step1$Sample.ID))
length(intersect(unique(f3_df %>% filter(Mutation_Gene1 == "BRCA2"))$PMBB_ID, all_b2o$step1$Sample.ID))
# 661
# 591

# with original
length(intersect(unique(f3 %>% filter(Gene == "BRCA1"))$PMBB_ID, all_b1$Sample.ID))
length(intersect(unique(f3 %>% filter(Gene == "BRCA2"))$PMBB_ID, all_b2$Sample.ID))
# 603
# 581

# with step 1
length(intersect(unique(f3 %>% filter(Gene == "BRCA1"))$PMBB_ID, all_b1o$step1$Sample.ID))
length(intersect(unique(f3 %>% filter(Gene == "BRCA2"))$PMBB_ID, all_b2o$step1$Sample.ID))
# 603
# 581

# with step 2
length(intersect(unique(f3 %>% filter(Gene == "BRCA1"))$PMBB_ID, all_b1o$step2$Sample.ID))
length(intersect(unique(f3 %>% filter(Gene == "BRCA2"))$PMBB_ID, all_b2o$step2$Sample.ID))
# 579
# 579

# with step 3
length(intersect(unique(f3 %>% filter(Gene == "BRCA1"))$PMBB_ID, all_b1o$step3$Sample.ID))
length(intersect(unique(f3 %>% filter(Gene == "BRCA2"))$PMBB_ID, all_b2o$step3$Sample.ID))
# 591
# 578

# length(setdiff(unique(f3 %>% filter(Gene == "BRCA1"))$PMBB_ID, all_b1o$step3$Sample.ID))
# length(setdiff(unique(f3 %>% filter(Gene == "BRCA2"))$PMBB_ID, all_b2o$step3$Sample.ID))
# # 98
# # 38
# # i guess we can see if these are not in the consented?

# ========================
# PREVIOUSLY UNIDENTIFIED SAMPLES
# ========================
COV <- here("PMBB", "3.0", "PMBB-Release-2024-3.0_covariates.txt")
PER <- here("PMBB", "3.0", "PMBB-Release-2024-3.0_phenotype_person.txt")
cov <- fread(COV, header = TRUE)
person <- fread(PER, header = TRUE)
flags <- read_excel(file.path("PMBB", "pmbb_flag_tables.xlsx"))
consent <- flags %>% filter(PMBB_Consent == 1)
consent_ids <- consent$PMBB_ID
length(consent_ids)
# 53062

prev_b1_ids <- union(f2_b1_ids, f3_b1_ids)
prev_b2_ids <- union(f2_b2_ids, f3_b2_ids)

newid_b1 <- setdiff(all_b1o_ids, prev_b1_ids)
newid_b2 <- setdiff(all_b2o_ids, prev_b2_ids)

all_b1o_path <- all_b1o$step3 %>% filter(Sample.ID %in% newid_b1, Variant.LoF_level == 1)
all_b2o_path <- all_b2o$step3 %>% filter(Sample.ID %in% newid_b2, Variant.LoF_level == 1)

write.csv(all_b1o_path, file.path("ch", "data", "brca1_all3_undiscovered_df.csv"), row.names = FALSE)
write.csv(all_b2o_path, file.path("ch", "data", "brca2_all3_undiscovered_df.csv"), row.names = FALSE)

# those that may be previously unidentified?
sum(!(newid_b1 %in% consent_ids)) # 96
sum(!(newid_b2 %in% consent_ids)) # 86
sum((newid_b1 %in% consent_ids)) # 38
sum((newid_b2 %in% consent_ids)) # 93

new_b1_df <- cov %>% filter(person_id %in% newid_b1)
new_b2_df <- cov %>% filter(person_id %in% newid_b2)

write.csv(newid_b1, file.path("ch", "data", "brca1_all3_undiscovered.csv"), row.names = FALSE)
write.csv(newid_b2, file.path("ch", "data", "brca2_all3_undiscovered.csv"), row.names = FALSE)

table(new_b1_df$Batch)
table(new_b2_df$Batch)
# 1  2
# 88 94
# > table(new_b2_df$Batch)
#
# 1   2
# 105 161

# hmmm okay..
