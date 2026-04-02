# ========================
# DETERMINING GERMLINE BRCA12 CARRIERS IN PMBB
# ========================

library(here)
setwd("/Users/jennyzli/Documents/Nathanson")
source(here("R", "config.R"))
source(here("ch", "scripts", "brca_qc2.R"))

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

table(f2_brca1$ExonicFunc.refGene)
# .  frameshift_deletion frameshift_insertion    nonsynonymous_SNV             stopgain
# 9                   77                   34                   22                   11

table(f2_brca1$GenomicRegion.refGene)
# exonic splicing
# 144        9

table(f2_brca1$CLNSIG)
# likely_patho, patho, or both

# ========================
# FREEZE 3.0
# ========================
f3_b12 <- read_excel(here("ch", "ss", "BW-Regeneron_Workbook_20250627.xlsx"), sheet = "CREP_workbook") %>% filter(Mutation_Gene1 %in% c("BRCA1", "BRCA2"))  %>%
    left_join(up, "VCFID") %>%
    filter(!is.na(PMBB_ID))
f3_b12_fnd <- read_excel(here("ch", "ss", "BW-Regeneron_Workbook_20250627.xlsx"), sheet = "CREP.found_lines") %>% filter(Gene %in% c("BRCA1", "BRCA2")) %>% filter(!is.na(Sample.ID)) %>%
    left_join(up, c("Sample.ID" =  "match_col")) %>%
    filter(!is.na(PMBB_ID))
f3_b12_pos <- read_excel(here("ch", "ss", "BW-Regeneron_Workbook_20250627.xlsx"), sheet = "CREP.possible_lines") %>% filter(Gene %in% c("BRCA1", "BRCA2")) %>% filter(!is.na(Sample.ID)) %>%
    left_join(up, c("Sample.ID" =  "match_col")) %>%
    filter(!is.na(PMBB_ID))
f3_b12_out <- read_excel(here("ch", "ss", "BW-Regeneron_Workbook_20250627.xlsx"), sheet = "CREP.filtered_out_variant_lines") %>% filter(Gene %in% c("BRCA1", "BRCA2")) %>% filter(!is.na(Sample.ID)) %>%
    left_join(up, c("Sample.ID" =  "match_col")) %>%
    filter(!is.na(PMBB_ID))
f3_b12_und <- read_excel(here("ch", "ss", "BW-Regeneron_Workbook_20250627.xlsx"), sheet = "CREP.undetectable_variants") %>% filter(Mutation_Gene1 %in% c("BRCA1", "BRCA2"))%>%
    filter(!is.na(Sample.ID)) %>%
    left_join(up, c("Sample.ID" =  "VCFID")) %>%
    filter(!is.na(PMBB_ID))
f3_b12_mis <- read_excel(here("ch", "ss", "BW-Regeneron_Workbook_20250627.xlsx"), sheet = "CREP.missing_variants") %>% filter(Mutation_Gene1 %in% c("BRCA1", "BRCA2")) %>%
    filter(!is.na(Sample.ID)) %>%
    left_join(up, c("Sample.ID" =  "VCFID")) %>%
    filter(!is.na(PMBB_ID))

dim(f3_b12) # 1278
f3_b12_df <- rbind(f3_b12_pos %>% select(-variant_key), f3_b12_fnd)

# those with variant annotations
f3_b12_anno <- f3_b12_df %>% filter(VCFID %in% f3_b12$VCFID) # 1190
# those without variant annotations
f3_b12_noanno <- f3_b12 %>% filter(!(VCFID %in% f3_b12_df$VCFID)) # 89

set_f3_b12_fnd <- unique(f3_b12_fnd$PMBB_ID)
set_f3_b12_pos <- unique(f3_b12_pos$PMBB_ID)

f3_b12_anno <- f3_b12_anno %>%
    mutate(
        detection_status = case_when(
            PMBB_ID %in% set_f3_b12_fnd ~ "found",
            PMBB_ID %in% set_f3_b12_pos ~ "possible",
            TRUE ~ "filtered_out"  # These are all from f3_b12_out
        )
    )
table(f3_b12_anno$detection_status)
# found possible
# 1172       18

f3_b12_noanno <- f3_b12_noanno %>%
    mutate(
        detection_status = case_when(
            VCFID %in% f3_b12_und$Sample.ID ~ "undetected",  # Fixed column name
            VCFID %in% f3_b12_mis$Sample.ID ~ "missing",     # Fixed column name
            TRUE ~ "other"
        )
    )
table(f3_b12_noanno$detection_status)
# missing      other undetected
# 19          1         69

write.csv(f3_b12_anno, file.path("ch", "data", "crep_workbook_brca12_annotated.csv"))

### INITIAL STATS
range(f3_b12_anno$Sample.Depth) # 14 117
range(f3_b12_anno$Sample.AltDepth) # 4 55
range(f3_b12_anno$Sample.AltFrac) # 0.179 0.792
length(unique(f3_b12_anno$Sample.ID)) # 1176
length(unique(f3_b12_anno %>% filter(Gene == "BRCA1"))$Sample.ID) # 604
length(unique(f3_b12_anno %>% filter(Gene == "BRCA2"))$Sample.ID) # 586

unique(f3_b12_anno$Bio.type)
# protein_coding

unique(f3_b12_anno$Variant.Class)
# "SNV"          "insertion"    "deletion"     "substitution"

unique(f3_b12_anno$Variant.Consequence)
# [1] "stop_gained"
# [2] "frameshift_variant"
# [3] "stop_gained&frameshift_variant"
# [4] "splice_acceptor_variant"
# [5] "missense_variant"
# [6] "frameshift_variant&splice_region_variant"
# [7] "splice_donor_variant"
# [8] "splice_polypyrimidine_tract_variant&intron_variant"
# [9] "missense_variant&splice_region_variant"
# [10] "inframe_deletion"
# [11] "splice_donor_region_variant&intron_variant"
# [12] "stop_gained&splice_region_variant"
# [13] "splice_donor_variant&coding_sequence_variant"
# [14] "start_lost"
# [15] "splice_donor_5th_base_variant&intron_variant"
# [16] "synonymous_variant" -- should be filtered out?
# [17] "stop_gained&inframe_insertion"
# [18] "intron_variant"
# [19] "splice_region_variant&synonymous_variant"

unique(f3_b12_anno$Variant.LoF_level)
# 1 2 3 4

unique(f3_b12_anno$Variant.LoF_level)
# 1 2 3 4

unique(f3_b12_anno$AutoGVP)
# [1] "Pathogenic"             "."                      "Likely_pathogenic"
# [4] "Uncertain_significance" "Likely_benign"          "Benign"

### SAVE PATHOGENICS
f3_b12_lof1 <- f3_b12_anno %>% filter(f3_b12_anno$Variant.LoF_level == 1)
f3_b12_lof2 <- f3_b12_anno %>% filter(f3_b12_anno$Variant.LoF_level == 2)
dim(f3_b12_lof1)
# 1166
dim(f3_b12_lof2)
# 16

write.csv(f3_b12_lof1, file.path("ch", "data", "crep_workbook_brca12_anno_lof1.csv"))
write.csv(f3_b12_lof2, file.path("ch", "data", "crep_workbook_brca12_anno_lof2.csv"))
# these VUS are all found in CREP, so I'll keep them ...

# may need to add additional QC
f3_b12_lof12 <- f3_b12_anno %>% filter(f3_b12_anno$Variant.LoF_level %in% c("1", "2"))
#  1182   62
range(f3_b12_lof12$Sample.Depth) # 14 117
range(f3_b12_lof12$Sample.AltDepth) # 5 55
range(f3_b12_lof12$Sample.AltFrac) # 0.179 0.792
length(unique(f3_b12_lof12$Sample.ID)) # 1176

length(unique(f3_b12_lof12 %>% filter(Gene == "BRCA1"))$Sample.ID) # 604
length(unique(f3_b12_lof12 %>% filter(Gene == "BRCA2"))$Sample.ID) # 586

write.csv(f3_b12_lof12, file.path("ch", "data", "f3_b12_lof12.csv"))

f3_b1_lof12o <- process_brca_data(f3_b12_lof12, "BRCA1", "crep", 0.18, 0.7, 17, 4, 0)
f3_b2_lof12o <- process_brca_data(f3_b12_lof12, "BRCA2", "crep", 0.18, 0.7, 17, 4, 0)

write.csv(f3_b1_lof12o$step3, file.path("ch", "data", "f3_brca1_path_step3.csv"))
write.csv(f3_b2_lof12o$step3, file.path("ch", "data", "f3_brca2_path_step3.csv"))

f3_b1_lof12_qc <- f3_b1_lof12o$step3
f3_b2_lof12_qc <- f3_b2_lof12o$step3
dim(f3_b1_lof12_qc) # 596
dim(f3_b2_lof12_qc) # 576

create_diagnostic_pdf(
    f3_b1_lof12o$step1, f3_b1_lof12o$step2, f3_b1_lof12o$step3,
    "BRCA1 (CREP)",
    file.path("ch", "figures", "brca1_crep_qc_histograms.pdf")
)
create_diagnostic_pdf(
    f3_b2_lof12o$step1, f3_b2_lof12o$step2, f3_b2_lof12o$step3,
    "BRCA2 (CREP)",
    file.path("ch", "figures", "brca2_crep_qc_histograms.pdf")
)

### ADD THOSE MISSING FROM PROGENY
f3_b1_noanno <- f3_b12_noanno %>% filter(Mutation_Gene1 == "BRCA1")
f3_b2_noanno <- f3_b12_noanno %>% filter(Mutation_Gene1 == "BRCA2")

f3_b1_ids <- sort(unique(c(f3_b1_lof12o$step3$PMBB_ID, f3_b1_noanno$PMBB_ID)))
f3_b2_ids <- sort(unique(c(f3_b2_lof12o$step3$PMBB_ID, f3_b2_noanno$PMBB_ID)))
length(f3_b1_ids)
length(f3_b2_ids)
# > length(f3_b1_ids)
# [1] 667
# > length(f3_b2_ids)
# [1] 594

# ========================
# TRIPLE CHECK WITH PROGENY
# ========================
length(intersect(f3_b12_noanno$PMBB_ID, progeny_pmbb$PMBB_ID))
# all of them are in progeny so that's good

progeny <- read_excel(here("ch", "ss", "brca_carriers_ch_freq_w_seen_in_crep_20251020.xlsx"), sheet = "Data_from_master_table") %>% filter(DNA == "D")
dim(progeny)
# 2816

progeny_merged <- merge_duplicates(progeny, "SampNum")
progeny_merged$SampNum <- as.numeric(progeny_merged$SampNum)

progeny_pmbb <- progeny_merged %>%
    left_join(up, by = "SampNum") %>%
    filter(!is.na(PMBB_ID))
dim(progeny_pmbb)
# 1278

progeny_pmbb_cov <- progeny_pmbb %>% left_join(cov, by = c("PMBB_ID" = "person_id"))
table(progeny_pmbb_cov$Batch)
# 1    2
# 1 1277
# one is not in CREP in freeze 2

b1_progeny_pmbb <- progeny_pmbb %>% filter(BRCA1_Presence_of_Mutation %in% c("Yes", "Obligate Carrier"))
length(unique(b1_progeny_pmbb$PMBB_ID)) #681
b2_progeny_pmbb <- progeny_pmbb %>% filter(BRCA2_Presence_of_Mutation %in% c("Yes", "Obligate Carrier", "VUS Positive"))
length(unique(b2_progeny_pmbb$PMBB_ID)) #607
# VUS is only one PMBB6960038827474, which I manually checked has BRCA2

### EXAMINING THE UNKNOWN
# b1_progeny_pmbb_unk <- progeny_pmbb %>% filter(BRCA1_Presence_of_Mutation %in% c("Unknown", "Not Tested"))
# dim(b1_progeny_pmbb_unk) #203
# b2_progeny_pmbb_unk <- progeny_pmbb %>% filter(BRCA2_Presence_of_Mutation %in% c("Unknown", "Not Tested"))
# dim(b2_progeny_pmbb_unk) #209
#
# sum(b1_progeny_pmbb_unk$PMBB_ID %in% f3_b12_lof12$PMBB_ID) # 193
# sum(b2_progeny_pmbb_unk$PMBB_ID %in% f3_b12_lof12$PMBB_ID) # 178
#
# sum(b1_progeny_pmbb_unk$PMBB_ID %in% f3_b1_lof12_qc$PMBB_ID) # 0
# sum(b2_progeny_pmbb_unk$PMBB_ID %in% f3_b2_lof12_qc$PMBB_ID) # 0
#
# # SO THESE ARE MOSTLY JUST NOT CARRIERS of the other BRCA gene, EXCEPT FOR LIKE A FEW
# examine_b1 <- b1_progeny_pmbb_unk %>% filter(!(PMBB_ID %in% f3_b2_lof12_qc$PMBB_ID))
# examine_b2 <- b2_progeny_pmbb_unk %>% filter(!(PMBB_ID %in% f3_b1_lof12_qc$PMBB_ID))

# so i think these can be excluded as well...
#
# write.csv(examine_b1, file.path("ch", "data", "examine_b1.csv"))
# write.csv(examine_b2, file.path("ch", "data", "examine_b2.csv"))

# ========================
# TAKE THE UNION
# ========================
f23_b1_ids <- sort(unique(c(b1_progeny_pmbb$PMBB_ID, f2_brca1$SampleID, f3_b1_ids)))
f23_b2_ids <- sort(unique(c(b2_progeny_pmbb$PMBB_ID, f2_brca2$SampleID, f3_b2_ids)))

length(f23_b1_progeny_ids)
length(f23_b2_progeny_ids)
# > length(f23_b1_ids)
# [1] 834
# > length(f23_b2_ids)
# [1] 870

# ========================
# MY ANNOTATIONS
# ========================
all_b1o <- process_brca_data(all_b1, "BRCA1", "all3", 0.18, 0.7, 17, 4, 20)
all_b2o <- process_brca_data(all_b2, "BRCA2", "all3", 0.18, 0.7, 17, 4, 20)

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

all_b1o_lof1 <- all_b1o$step3 %>% filter(Variant.LoF_level == 1)
all_b2o_lof1 <- all_b2o$step3 %>% filter(Variant.LoF_level == 1)
length(unique(all_b1o_lof1$Sample.ID)) # 859
length(unique(all_b2o_lof1$Sample.ID)) # 1007

table(all_b1o_lof2$AutoGVP)
# .                 Benign          Likely_benign
# 38                 117119                   1553
# Likely_pathogenic             Pathogenic Uncertain_significance
# 11                    870

all_b1o_lof2 <- all_b1o$step3 %>% filter(Variant.LoF_level == 2)
all_b2o_lof2 <- all_b2o$step3 %>% filter(Variant.LoF_level == 2)
length(unique(all_b1o_lof2$Sample.ID)) # 58
length(unique(all_b2o_lof2$Sample.ID)) # 92

all_b1o_path <- all_b1o$step3 %>% filter(Variant.LoF_level %in% c(1, 2))
all_b2o_path <- all_b2o$step3 %>% filter(Variant.LoF_level %in% c(1, 2))

# ========================
# DOES IT HAVE OLD SAMPLES
# ========================
length(f23_b1_ids) #834
sum(all_b1o_lof1$Sample.ID %in% f23_b1_ids)
# 737

length(f23_b2_ids)
sum(f23_b2_ids %in% all_b2o_ids)

# ========================
# PREVIOUSLY UNIDENTIFIED SAMPLES
# ========================
# like for sure??
new_b1_ids <- setdiff(all_b1o_lof1$Sample.ID, f23_b1_ids)
new_b2_ids <- setdiff(all_b2o_lof1$Sample.ID, f23_b2_ids)
length(new_b1_ids) # 129
length(new_b2_ids) # 183

new_b1o_path <- all_b1o_path %>% filter(Sample.ID %in% new_b1_ids)
new_b2o_path <- all_b2o_path %>% filter(Sample.ID %in% new_b2_ids)

write.csv(all_b1o_path, file.path("ch", "data", "brca1_all3_undiscovered_df.csv"), row.names = FALSE)
write.csv(all_b2o_path, file.path("ch", "data", "brca2_all3_undiscovered_df.csv"), row.names = FALSE)

# check which freeze...
new_b1_df <- cov %>% filter(person_id %in% new_b1_ids)
new_b2_df <- cov %>% filter(person_id %in% new_b2_ids)

table(new_b1_df$Batch)
table(new_b2_df$Batch)
# 1  2
# 30 142
# > table(new_b2_df$Batch)
#
# 1   2
# 91 182

# so mostly in freeze 3 but some in freeze 2...

### IN CONSENTED??
sum(new_b1_ids %in% consent$PMBB_ID)
sum(new_b2_ids %in% consent$PMBB_ID)

