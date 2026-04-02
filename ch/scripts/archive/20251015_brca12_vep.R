# ========================
# Packages
# ========================

library(here)
setwd("/Users/jennyzli/Documents/Nathanson")
source(here("R", "load_packages.R"))
source(here("R", "sample_selection.R"))
source(here("R", "control_selection.R"))
library(VennDiagram)

# ========================
# LOAD DATA
# ========================
brca1 <- read.csv(here("ch", "data", "PMBB-Release-2024-3.0_genetic_exome_BRCA1_NF.consented_pmbb-only.norm.vep.report.csv"))
brca2 <- read.csv(here("ch", "data", "PMBB-Release-2024-3.0_genetic_exome_BRCA2_NF.consented_pmbb-only.norm.vep.report.csv"))
brca1 <- brca1 %>% filter(Gene == "BRCA1")
brca2 <- brca2 %>% filter(Gene == "BRCA2")

found <- read_excel(here("ch", "ss", "BW-Regeneron_Workbook_20250627.xlsx"), sheet = "CREP.found_lines") %>% filter(Gene == "BRCA1")
likely <- read_excel(here("ch", "ss", "BW-Regeneron_Workbook_20250627.xlsx"), sheet = "CREP.likely_lines") %>% filter(Gene == "BRCA1")
possible <- read_excel(here("ch", "ss", "BW-Regeneron_Workbook_20250627.xlsx"), sheet = "CREP.possible_lines") %>% filter(Gene == "BRCA1")
out <- read_excel(here("ch", "ss", "BW-Regeneron_Workbook_20250627.xlsx"), sheet = "CREP.filtered_out_variant_lines") %>% filter(Gene == "BRCA1")

unique(found$ClinVar.SIG)
unique(found$AutoGVP)

unique(likely$ClinVar.SIG)
unique(likely$AutoGVP)

unique(possible$ClinVar.SIG)
unique(possible$AutoGVP)

unique(out$ClinVar.SIG)
unique(out$AutoGVP)

dim(brca1)
# 496899     58
dim(brca2)
# 538838     58

length(unique(brca1$Sample.ID))
# 45459
length(unique(brca2$Sample.ID))
# 51553

### FILTER FOR PATHOGENICITY ####
table(brca1$Variant.LoF_level)
# 1      2      3      4
# 305    209   8765 574025

brca1_path <- brca1 %>% filter(Variant.LoF_level == '1')
brca1_vus <- brca1 %>% filter(Variant.LoF_level == '2')

table(brca2$Variant.LoF_level)
# 1      2      3      4
# 541    199   3836 575799

brca2_path <- brca2 %>% filter(Variant.LoF_level == '1')
brca2_vus <- brca2 %>% filter(Variant.LoF_level == '2')







