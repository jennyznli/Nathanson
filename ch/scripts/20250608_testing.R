# ========================
# PACKAGES
# ========================
packages <- c("tidyr", "dplyr", "plotly", "readr", "readxl", "here",
              "stringr", "ggplot2",  "impute", "pals", "geneplotter"
)
purrr::walk(packages, ~ require(.x, character.only = TRUE))
here()

DATE <- format(Sys.Date(), "%Y%m%d")

# ========================
# GET CH GENES
# ========================
genes <- as.data.frame(read_excel(here("ch", "genes.xlsx"), col_names = FALSE))
colnames(genes) <- "ch_genes"
write.csv(genes, here("ch", "ch_genes.csv"), row.names = FALSE)

# ============================================================
# READ IN SS
# ============================================================
workbook <- read_excel(here("ch", "BW-Regeneron_Workbook_20250602.xlsx"), sheet = "CREP_workbook")
length(unique(workbook$Participant))
# 2869 unique participants from CREP
unique(workbook$Mutation_Gene1)
# [1] NA       "BRCA1"  "BRCA2"  "MSH2"   "TP53"   "PMS2"   "DSP"    "CHEK2"  "RAD51D" "CDKN2A" "APC"    "ATM"    "MLH1"   "PALB2"
# [15] "STK11"  "MSH6"   "NBN"    "MUTYH"  "SDHD"   "EPHB4"  "CDH1"   "BARD1"  "PTEN"   "BRIP1"  "NTHL1"  "RAD51C" "TSC1"   "EPCAM"
all_participants <- (unique(workbook$Participant))

up <- read_excel(here("ch", "BW-Regeneron_Workbook_20250602.xlsx"), sheet = "UP")
# 2869 in PMBB
pp <- read_excel(here("ch", "BW-Regeneron_Workbook_20250602.xlsx"), sheet = "PP")
# 305?
tb <- read_excel(here("ch", "BW-Regeneron_Workbook_20250602.xlsx"), sheet = "TB")
# 1000
kt <- read_excel(here("ch", "BW-Regeneron_Workbook_20250602.xlsx"), sheet = "KT")
# kate 831

# ============================================================
# READ IN CSV
# ============================================================
ch_genes <- read.csv(here("ch", "ch_genes.csv"), header = TRUE)
genes <- as.factor(ch_genes$ch_genes)

found_vars <- read.csv(here("ch", "csv", "found_variants.csv"))
likely_vars <- read.csv(here("ch", "csv", "likely_variants.csv"))

# possible_lines <- read.csv(here("ch", "csv", "possible_variant_lines.csv"))
# dim(possible_lines) # 403 ppl
#
# possible_vars <- read.csv(here("ch", "csv", "possible_variants.csv"))
# dim(possible_vars) # 407 ppl

found_vars <- read.csv(here("ch", "csv", "found_variants.csv"))
found_lines <- read.csv(here("ch", "csv", "found_variant_lines.csv"))

# 1290 so all variants in the found lines
sum(found_vars$VCFID %in% found_lines$Sample.ID)

# lines excluded? 53
excluded <- found_lines %>% filter(!(Sample.ID %in% found_vars$VCFID))
x <- !(found_lines$Sample.ID %in% found_vars$VCFID)

# ============================================================
# FIND BRCA1/2 CARRIERS
# ============================================================
found_vars <- read.csv(here("ch", "csv", "found_variants.csv"))
found_lines <- read.csv(here("ch", "csv", "found_variant_lines.csv"))
found_lines$Sample_ID <- sub(".*(UP\\d+).*", "\\1", found_lines$Sample.ID)

# 1122 found
found_brca = found_vars %>% filter(Mutation_Gene1 == "BRCA1" | Mutation_Gene1 == "BRCA2")
id1 <- found_brca$Sample.ID
# 37 likely
likely_brca = likely_vars %>% filter(Mutation_Gene1 == "BRCA1" | Mutation_Gene1 == "BRCA2")
id2 <- likely_brca$Sample.ID

ids <- unique(c(id1, id2))
# 1159 unique ppl, none overlapping makes sense

ss2 <- read_excel(here("ch", "CREPTotal_wRegeneronv3_050125_w_mutations_w_PMBBID.xlsx"), sheet = "Clean - V3")
ss2$Sample.ID <- sprintf("UP%04d", ss2$Participant)

# filter by gnomAD




ss2 = ss2 %>% filter(Sample.ID %in% ids)
ss2 <- ss2[match(ids, ss2$Sample.ID), ]

found_lines = found_lines %>% filter(Sample_ID %in% ids)
found_lines$Sample_ID <- sub(".*(UP\\d+).*", "\\1", found_lines$Sample.ID)

duplicates <- found_lines[duplicated(found_lines$Sample.ID) | duplicated(found_lines$Sample.ID, fromLast = TRUE), ]









