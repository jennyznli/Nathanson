library(here)
setwd("/Users/jennyzli/Documents/Nathanson")
source(here("R", "config.R"))
source(here("ch", "scripts", "brca_qc.R"))

library(ggplot2)
library(gridExtra)
library(grid)

# ========================
# READ DATA
# ========================
all_brca1 <- read.csv(file.path("ch", "data", "brca_vep.BRCA1.vep.report.csv"))
all_brca2 <- read.csv(file.path("ch", "data", "brca_vep.BRCA2.vep.report.csv"))
cons_brca1 <- read.csv(here("ch", "data", "PMBB-Release-2024-3.0_genetic_exome_BRCA1_NF.consented_pmbb-only.norm.vep.report.csv"))
cons_brca2 <- read.csv(here("ch", "data", "PMBB-Release-2024-3.0_genetic_exome_BRCA2_NF.consented_pmbb-only.norm.vep.report.csv"))
#
# all2_brca1 <- read.csv(file.path("ch", "data", "brca_vep2.BRCA1.vep.report.csv"))
# all2_brca2 <- read.csv(file.path("ch", "data", "brca_vep2.BRCA2.vep.report.csv"))

# ========================
# PROCESS DATA
# ========================
all_b1 <- process_brca_data(all_brca1, "BRCA1", "all")
all_b2 <- process_brca_data(all_brca2, "BRCA2", "all")
cons_b1 <- process_brca_data(cons_brca1, "BRCA1", "consented")
cons_b2 <- process_brca_data(cons_brca2, "BRCA2", "consented")
#
# all2_b1 <- process_brca_data(all2_brca1, "BRCA1", "all2")
# all2_b2 <- process_brca_data(all2_brca2, "BRCA2", "all2")

# ========================
# GENERATE DIAGNOSTIC PLOTS
# ========================
create_diagnostic_pdf(
    all_b1$step1, all_b1$step2, all_b1$step3,
    "BRCA1 (All PMBB)",
    file.path("ch", "figures", "brca1_all_qc_histograms.pdf")
)

create_diagnostic_pdf(
    all_b2$step1, all_b2$step2, all_b2$step3,
    "BRCA2 (All PMBB)",
    file.path("ch", "figures", "brca2_all_qc_histograms.pdf")
)

create_diagnostic_pdf(
    cons_b1$step1, cons_b1$step2, cons_b1$step3,
    "BRCA1 (Consented)",
    file.path("ch", "figures", "brca1_cons_qc_histograms.pdf")
)

create_diagnostic_pdf(
    cons_b2$step1, cons_b2$step2, cons_b2$step3,
    "BRCA2 (Consented)",
    file.path("ch", "figures", "brca2_cons_qc_histograms.pdf")
)


# create_diagnostic_pdf(
#     all2_b1$step1, all2_b1$step2, all2_b1$step3,
#     "BRCA2 (Consented)",
#     file.path("ch", "figures", "all2_brca2_cons_qc_histograms.pdf")
# )
#
#
# create_diagnostic_pdf(
#     all2_b2$step1, all2_b2$step2, all2_b2$step3,
#     "BRCA2 (Consented)",
#     file.path("ch", "figures", "all2_brca2_cons_qc_histograms.pdf")
# )

# ========================
# SAVE FILTERED DATA
# ========================
write.csv(all_b1$step3, file.path("ch", "data", "brca1_all_filtered.csv"), row.names = FALSE)
write.csv(all_b2$step3, file.path("ch", "data", "brca2_all_filtered.csv"), row.names = FALSE)
write.csv(cons_b1$step3, file.path("ch", "data", "brca1_cons_filtered.csv"), row.names = FALSE)
write.csv(cons_b2$step3, file.path("ch", "data", "brca2_cons_filtered.csv"), row.names = FALSE)

# save individuals with pathogenic variants and see how it overlaps with others
all_b1_path <- all_b1$step3 %>% filter(Variant.LoF_level %in% c(1, 2)) %>% select(Sample.ID)
all_b2_path <- all_b2$step3 %>% filter(Variant.LoF_level %in% c(1, 2)) %>% select(Sample.ID)
dim(all_b1_path) #399  58
dim(all_b2_path) #137  58

# all_b1_ids <- all_b1$step3 %>% filter(AutoGVP %in% c("Pathogenic", "Likely_pathogenic")) %>% select(Sample.ID)
# all_b2_ids <- all_b2$step3 %>% filter(AutoGVP %in% c("Pathogenic", "Likely_pathogenic")) %>% select(Sample.ID)
# dim(all_b1_ids) #383  58
# dim(all_b2_ids) #122  58

all_b1_ids <- all_b1_path$Sample.ID
all_b2_ids <- all_b2_path$Sample.ID

write.csv(all_b1$step3, file.path("ch", "data", "brca1_all_filtered.csv"), row.names = FALSE)
write.csv(all_b2$step3, file.path("ch", "data", "brca2_all_filtered.csv"), row.names = FALSE)
write.csv(cons_b1$step3, file.path("ch", "data", "brca1_cons_filtered.csv"), row.names = FALSE)
write.csv(cons_b2$step3, file.path("ch", "data", "brca2_cons_filtered.csv"), row.names = FALSE)

# ========================
# COMPARE W PREVIOUS FREEZE 2
# ========================
f2_brca1 <- read_excel(here("ch", "ss", "BRCA1.BRCA2 P.LP_04.11.21.xlsx"), sheet = "BRCA1")
f2_brca2 <- read_excel(here("ch", "ss", "BRCA1.BRCA2 P.LP_04.11.21.xlsx"), sheet = "BRCA2")
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

f2_brca1_fil <- f2_brca1 %>% filter(Total_Depth > 20, ALT_AlleleDepth > 5, ALT_AlleleFrac > 0.3, ALT_AlleleFrac < 0.7)
f2_brca2_fil <- f2_brca2 %>% filter(Total_Depth > 20, ALT_AlleleDepth > 5, ALT_AlleleFrac > 0.3, ALT_AlleleFrac < 0.7)

# ORIGINAL INDIVIDUALS
length(unique(f2_brca1_fil$SampleID))
length(unique(f2_brca2_fil$SampleID))
# > length(unique(f2_brca1_fil$SampleID))
# [1] 148
# > length(unique(f2_brca2_fil$SampleID))
# [1] 256

f2_b1_ids <- unique(f2_brca1_fil$SampleID)
f2_b2_ids <- unique(f2_brca2_fil$SampleID)
length(f2_b1_ids)
length(f2_b2_ids)
# 148
# 256

length(intersect(all_b1_ids, f2_b1_ids))
length(intersect(all_b2_ids, f2_b2_ids))
# 64 ...
# 32...

# ========================
# COMPARE WITH BRAD'S f3
# ========================
f3_df <- read_excel(here("ch", "ss", "BW-Regeneron_Workbook_20250627.xlsx"), sheet = "CREP_workbook") %>% filter(Mutation_Gene1 %in% c("BRCA1", "BRCA2"))
up <- read.csv(here("simplexo", "data", "simplexo_up_map.csv"))
f3_df <- f3_df %>% left_join(up, by = "VCFID") %>% filter(!is.na(PMBB_ID))
length(f3_ids)
# 1265



f3_fnd <- read_excel(here("ch", "ss", "BW-Regeneron_Workbook_20250627.xlsx"), sheet = "CREP.found_lines") %>% filter(Gene %in% c("BRCA1", "BRCA2")) %>% filter(!is.na(Sample.ID))
f3_possible <- read_excel(here("ch", "ss", "BW-Regeneron_Workbook_20250627.xlsx"), sheet = "CREP.possible_lines") %>% filter(Gene %in% c("BRCA1", "BRCA2")) %>% filter(!is.na(Sample.ID))

f3_fnd_fil <- f3_fnd %>% filter(Sample.ID %in% f3_df$match_col)
f3_pos_fil <- f3_possible %>% filter(Sample.ID %in% f3_df$match_col) %>% select(-variant_key)

f3_fil <- rbind(f3_fnd_fil, f3_pos_fil)
range(f3_fil$Sample.Depth) # 21 117
range(f3_fil$Sample.AltDepth) # 8
range(f3_fil$Sample.AltFrac) # 0.302 0.692

f3_fil <- f3_fil %>% filter(Sample.Depth > 20, Sample.AltDepth > 5, Sample.AltFrac > 0.3, Sample.AltFrac < 0.7)
f3_fil <- f3_fil %>% left_join(up, by = c("Sample.ID" = "match_col"))

# 1139   58
f3_b1_fil <- f3_fil %>% filter(Gene == "BRCA1")
f3_b2_fil <- f3_fil %>% filter(Gene == "BRCA2")

# ORIGINAL NUMBER
length(unique(f3_b1_fil$Sample.ID))
length(unique(f3_b2_fil$Sample.ID))
# 674
# 591

f3_b1_ids <- unique(f3_b1_fil$PMBB_ID)
f3_b2_ids <- unique(f3_b2_fil$PMBB_ID)

# INTERSECTING...
length(intersect(all_b1_ids, f3_b1_ids))
length(intersect(all_b2_ids, f3_b2_ids))
# 259
# 65

# ========================
# PREVIOUSLY UNIDENTIFIED SAMPLES
# ========================
prev_b1_ids <- union(f2_b1_ids, f3_b1_ids)
prev_b2_ids <- union(f2_b2_ids, f3_b2_ids)

newid_b1 <- setdiff(all_b1_ids, prev_b1_ids)
newid_b2 <- setdiff(all_b2_ids, prev_b2_ids)
length(newid_b1)
length(newid_b2)
# 75
# 38 ...


