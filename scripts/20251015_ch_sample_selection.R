# ========================
# CH SAMPLE SELECTION
# ========================
library(here)
setwd( "/Users/jennyzli/Documents/Nathanson")
source(here("R", "config.R"))
library(VennDiagram)

# ========================
# LOAD DATA
# ========================
### PMBB ###
progeny$SampNum <- as.numeric(progeny$SampNum)
consent <- read_excel(here("PMBB", "pmbb_flag_tables.xlsx")) %>% filter(PMBB_Consent == 1)
# 53062 bruhhh
up <- read.csv(here("simplexo", "data", "simplexo_up_map.csv"))

PMBB_DIR = here("PMBB", "3.0")
OCC <- here(PMBB_DIR, "PMBB-Release-2024-3.0_phenotype_condition_occurrence.txt")
COV <- here(PMBB_DIR, "PMBB-Release-2024-3.0_covariates.txt")
PER <- here(PMBB_DIR, "PMBB-Release-2024-3.0_phenotype_person.txt")
cov <- fread(COV, header = TRUE)
person <- fread(PER, header = TRUE)

### BRCA1/2 ###
progeny <- read_excel(here("ch", "ss", "brca_carriers_ch_freq_w_seen_in_crep_20251020.xlsx"), sheet = "Data_from_master_table")
progeny <- progeny %>% filter(DNA == "D")

brad <- read_excel(here("ch", "ss", "BW-Regeneron_Workbook_20250627.xlsx"), sheet = "CREP_workbook")
old_brca1 <- read_excel(here("ch", "ss", "BRCA1.BRCA2 P.LP_04.11.21.xlsx"), sheet = "BRCA1")
old_brca2 <- read_excel(here("ch", "ss", "BRCA1.BRCA2 P.LP_04.11.21.xlsx"), sheet = "BRCA2")

### VARIANTS
# var <- read.csv(here("ch", "data", "all_genes_filtered_variants.csv"))
# x <- read.csv(here("ch", "data", "all_genes_with_filter_info.csv"))

# ========================
# PROGENY QC
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
table(progeny_pmbb$BRCA1_Presence_of_Mutation)
table(progeny_pmbb$BRCA2_Presence_of_Mutation)

# ========================
# DATA CAPTURE BY LOCATION
# ========================
progeny_tumor <- read_excel(here("ch", "ss", "brca_carriers_ch_freq_w_seen_in_crep_20251020.xlsx"), sheet = "Tumor") %>% filter(globalid %in% progeny_pmbb$globalid)
progeny_treat <- read_excel(here("ch", "ss", "brca_carriers_ch_freq_w_seen_in_crep_20251020.xlsx"), sheet = "Treatment") %>% filter(globalid %in% progeny_pmbb$globalid)
progeny_rec <- read_excel(here("ch", "ss", "brca_carriers_ch_freq_w_seen_in_crep_20251020.xlsx"), sheet = "Recurrence") %>% filter(globalid %in% progeny_pmbb$globalid)
progeny_path <- read_excel(here("ch", "ss", "brca_carriers_ch_freq_w_seen_in_crep_20251020.xlsx"), sheet = "Pathology") %>% filter(globalid %in% progeny_pmbb$globalid)

# length(unique(progeny_treat$globalid))
# length(unique(progeny_rec$globalid))
# length(unique(progeny_path$globalid))
# length(unique(progeny_tumor$globalid))
# length(unique(progeny_pmbb$globalid))

# smoking?
table(progeny_pmbb$SmokingEver)
# NA   Never      No Unknown     Yes
# 52     775      62      22     367
# 52 + 22 = 74 unknowns not terrible

progeny_pmbb <- progeny_pmbb %>%
    mutate(
        seen_at_penn = case_when(
            Seen_In_CREP %in% c("Yes", "HUP CREP", "Pennsylvania Hospital",
                                "Chester County Hospital", "Penn Medicine Cherry Hill",
                                "Penn Medicine Princeton") ~ "Yes",
            Seen_In_CREP == "No" ~ "No",
            TRUE ~ "Unknown/Missing"
        )
    )

data_capture <- progeny_pmbb %>%
    mutate(
        has_tumor = globalid %in% progeny_tumor$globalid,
        has_treatment = globalid %in% progeny_treat$globalid,
        has_recurrence = globalid %in% progeny_rec$globalid,
        has_pathology = globalid %in% progeny_path$globalid
    ) %>%
    group_by(seen_at_penn) %>%
    summarise(
        n = n(),
        n_tumor = sum(has_tumor),
        pct_tumor = round(mean(has_tumor) * 100, 1),
        n_treatment = sum(has_treatment),
        pct_treatment = round(mean(has_treatment) * 100, 1),
        n_recurrence = sum(has_recurrence),
        pct_recurrence = round(mean(has_recurrence) * 100, 1),
        n_pathology = sum(has_pathology),
        pct_pathology = round(mean(has_pathology) * 100, 1)
    )
print(data_capture)
write.csv(data_capture, here("ch", "data", "progeny_data_capture.csv"))

capture_table <- table(
    progeny_pmbb$seen_at_penn,
    progeny_pmbb$globalid %in% progeny_tumor$globalid
)
print(capture_table)
print(chisq.test(capture_table))
# data:  capture_table
# X-squared = 3.1429, df = 2, p-value = 0.2077
write.csv(capture_table, here("ch", "data", "progeny_tumor_table.csv"))

# ========================
# PROGENY - BRCA1/2
# ========================
brca1_progeny_pmbb <- progeny_pmbb %>% filter(BRCA1_Presence_of_Mutation %in% c("Yes", "Obligate Carrier", "VUS Positive"))
dim(brca1_progeny_pmbb) #682
brca2_progeny_pmbb <- progeny_pmbb %>% filter(BRCA2_Presence_of_Mutation %in% c("Yes", "Obligate Carrier", "VUS Positive"))
dim(brca2_progeny_pmbb) #607

length(intersect(brca1_progeny_pmbb$globalid, brca2_progeny_pmbb$globalid))
# 11 with both and in PMBB

# these are the tumors of the ones w/ dna in pmbb
progeny_pmbb_tumor <- progeny_tumor %>% filter(globalid %in% progeny_pmbb$globalid)

# ========================
# BRAD - BRCA1/2 CARRIERS
# ========================
# merge w/ map to PMBB ID...
length(unique(brad$Participant)) #2869

brad_pmbb <- left_join(brad, up, by = "VCFID")
dim(brad_pmbb) # 2869

brca1_brad <- brad_pmbb %>% filter(Mutation_Gene1 %in% c("BRCA1")) %>% filter(!(is.na(Participant)))
brca2_brad <- brad_pmbb %>% filter(Mutation_Gene1 %in% c("BRCA2")) %>% filter(!(is.na(Participant)))
length(unique(brca1_brad$Participant)) #678
length(unique(brca2_brad$Participant)) #592

# ========================
# FREEZE 2.0 - BRCA1/2 CARRIERS
# ========================
dim(old_brca1) # 153
dim(old_brca2) # 263
length(unique(old_brca1$SampleID)) # 152
length(unique(old_brca2$SampleID)) # 263
old_brca1 <- merge_duplicates(old_brca1, "SampleID")
dim(old_brca1) # 152

# ========================
# VAR CALLING - BRCA1/2 CARRIERS
# ========================
brca1 <- read.csv(here("ch", "data", "PMBB-Release-2024-3.0_genetic_exome_BRCA1_NF.consented_pmbb-only.norm.vep.report.csv"))
brca2 <- read.csv(here("ch", "data", "PMBB-Release-2024-3.0_genetic_exome_BRCA2_NF.consented_pmbb-only.norm.vep.report.csv"))

# BW-Regeneron_Workbook_20250627
#
flags <- read.csv(here("PMBB", "3.0", "rgcname_pmbbid_metadata_flags.csv"))
library(UpSetR)
library(ggplot2)

# Convert flags to logical matrix (UpSetR needs TRUE/FALSE or 1/0)
upset_data <- flags %>%
    select(PMBB_Consent, CREP, Melanoma, Testicular_CA, Para_Pheo, Ky) %>%
    mutate(across(everything(), ~ as.numeric(. == 1)))

# Create UpSet plot
upset_plot <- upset(
    upset_data,
    sets = c("PMBB_Consent", "CREP", "Melanoma", "Testicular_CA", "Para_Pheo", "Ky"),
    order.by = "freq",
    keep.order = TRUE,
    sets.bar.color = "#4472C4",
    main.bar.color = "#C55A11",
    point.size = 3.5,
    line.size = 1.5,
    text.scale = c(1.5, 1.3, 1.3, 1.2, 1.5, 1.2),  # Axis labels, numbers, set names
    mb.ratio = c(0.6, 0.4),  # Ratio of main bar to set size bar
    mainbar.y.label = "Intersection Size",
    sets.x.label = "Total per Category"
)

# Save
png(here("ch", "figures", "flags_upset_plot.png"),
    width = 12, height = 7, units = "in", res = 300)
upset_plot
dev.off()


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



# ========================
# VENN DIAGRAMS - 5 SOURCES (including VUS)
# ========================
library(ggVennDiagram)
library(ggplot2)

# Prepare data with 5 sources
brca1_sources <- list(
    Progeny = unique(brca1_progeny_pmbb$PMBB_ID),
    Brad = unique(brca1_brad$PMBB_ID),
    Freeze_2.0 = unique(old_brca1$SampleID),
    VCF_Path = unique(brca1_path$Sample.ID),      # LoF level 1
    VCF_VUS = unique(brca1_vus$Sample.ID)         # LoF level 2
)

brca2_sources <- list(
    Progeny = unique(brca2_progeny_pmbb$PMBB_ID),
    Brad = unique(brca2_brad$PMBB_ID),
    Freeze_2.0 = unique(old_brca2$SampleID),
    VCF_Path = unique(brca2_path$Sample.ID),      # LoF level 1
    VCF_VUS = unique(brca2_vus$Sample.ID)         # LoF level 2
)

# Print sample sizes for each source
cat("\n=== BRCA1 Sample Sizes ===\n")
sapply(brca1_sources, length)

cat("\n=== BRCA2 Sample Sizes ===\n")
sapply(brca2_sources, length)

cat("\n=== Total Unique BRCA1 Carriers (including VUS) ===\n")
length(unique(unlist(brca1_sources)))

cat("\n=== Total Unique BRCA2 Carriers (including VUS) ===\n")
length(unique(unlist(brca2_sources)))

##### BRCA1 VENN (5-way) #####
# Note: 5-way Venn diagrams are complex and may not display perfectly
brca1_venn <- ggVennDiagram(
    brca1_sources,
    label = "count",
    label_alpha = 0,
    label_size = 2.5,      # Smaller text for 5-way
    set_size = 3,
    edge_size = 0.5
) +
    scale_fill_gradient(low = "#F4F4F4", high = "#4472C4") +
    labs(title = "BRCA1 Carriers: Clinical Sources + VCF (P/LP & VUS)",
         subtitle = paste0("Total unique carriers: ",
                           length(unique(unlist(brca1_sources))))) +
    theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        plot.subtitle = element_text(hjust = 0.5, size = 11),
        legend.position = "right",
        legend.title = element_text(size = 9),
        plot.margin = margin(10, 10, 10, 10)
    )

ggsave(here("ch", "figures", "brca1_sources_venn_5way.png"),
       plot = brca1_venn,
       width = 12, height = 9,
       dpi = 300, bg = "white")

##### BRCA2 VENN (5-way) #####
brca2_venn <- ggVennDiagram(
    brca2_sources,
    label = "count",
    label_alpha = 0,
    label_size = 2.5,
    set_size = 3,
    edge_size = 0.5
) +
    scale_fill_gradient(low = "#F4F4F4", high = "#C55A11") +
    labs(title = "BRCA2 Carriers: Clinical Sources + VCF (P/LP & VUS)",
         subtitle = paste0("Total unique carriers: ",
                           length(unique(unlist(brca2_sources))))) +
    theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        plot.subtitle = element_text(hjust = 0.5, size = 11),
        legend.position = "right",
        legend.title = element_text(size = 9),
        plot.margin = margin(10, 10, 10, 10)
    )

ggsave(here("ch", "figures", "brca2_sources_venn_5way.png"),
       plot = brca2_venn,
       width = 12, height = 9,
       dpi = 300, bg = "white")

# ========================
# DETAILED OVERLAP ANALYSIS
# ========================

# Combine clinical sources for comparison
brca1_clinical <- unique(c(brca1_sources$Progeny, brca1_sources$Brad, brca1_sources$Freeze_2.0))
brca2_clinical <- unique(c(brca2_sources$Progeny, brca2_sources$Brad, brca2_sources$Freeze_2.0))

# Check overlap between VCF pathogenic and clinical
brca1_path_in_clinical <- intersect(brca1_sources$VCF_Path, brca1_clinical)
brca1_path_only <- setdiff(brca1_sources$VCF_Path, brca1_clinical)

brca2_path_in_clinical <- intersect(brca2_sources$VCF_Path, brca2_clinical)
brca2_path_only <- setdiff(brca2_sources$VCF_Path, brca2_clinical)

# Check overlap between VCF VUS and clinical
brca1_vus_in_clinical <- intersect(brca1_sources$VCF_VUS, brca1_clinical)
brca1_vus_only <- setdiff(brca1_sources$VCF_VUS, brca1_clinical)

brca2_vus_in_clinical <- intersect(brca2_sources$VCF_VUS, brca2_clinical)
brca2_vus_only <- setdiff(brca2_sources$VCF_VUS, brca2_clinical)

# Summary
vcf_summary <- data.frame(
    Category = c(
        "VCF P/LP: Total",
        "VCF P/LP: In clinical sources",
        "VCF P/LP: NOT in clinical sources",
        "VCF P/LP: % in clinical",
        "",
        "VCF VUS: Total",
        "VCF VUS: In clinical sources",
        "VCF VUS: NOT in clinical sources",
        "VCF VUS: % in clinical",
        "",
        "Clinical sources: Total",
        "Clinical: Have VCF P/LP",
        "Clinical: Have VCF VUS",
        "Clinical: No VCF variant"
    ),
    BRCA1 = c(
        length(brca1_sources$VCF_Path),
        length(brca1_path_in_clinical),
        length(brca1_path_only),
        round(length(brca1_path_in_clinical) / length(brca1_sources$VCF_Path) * 100, 1),
        "",
        length(brca1_sources$VCF_VUS),
        length(brca1_vus_in_clinical),
        length(brca1_vus_only),
        round(length(brca1_vus_in_clinical) / length(brca1_sources$VCF_VUS) * 100, 1),
        "",
        length(brca1_clinical),
        length(intersect(brca1_clinical, brca1_sources$VCF_Path)),
        length(intersect(brca1_clinical, brca1_sources$VCF_VUS)),
        length(setdiff(brca1_clinical, c(brca1_sources$VCF_Path, brca1_sources$VCF_VUS)))
    ),
    BRCA2 = c(
        length(brca2_sources$VCF_Path),
        length(brca2_path_in_clinical),
        length(brca2_path_only),
        round(length(brca2_path_in_clinical) / length(brca2_sources$VCF_Path) * 100, 1),
        "",
        length(brca2_sources$VCF_VUS),
        length(brca2_vus_in_clinical),
        length(brca2_vus_only),
        round(length(brca2_vus_in_clinical) / length(brca2_sources$VCF_VUS) * 100, 1),
        "",
        length(brca2_clinical),
        length(intersect(brca2_clinical, brca2_sources$VCF_Path)),
        length(intersect(brca2_clinical, brca2_sources$VCF_VUS)),
        length(setdiff(brca2_clinical, c(brca2_sources$VCF_Path, brca2_sources$VCF_VUS)))
    )
)

print(vcf_summary)
write.csv(vcf_summary, here("ch", "data", "brca_vcf_clinical_overlap.csv"), row.names = FALSE)

cat("\n=== Key Findings ===\n")
cat("BRCA1 P/LP variants found by VCF but NOT in clinical sources:", length(brca1_path_only), "\n")
cat("BRCA2 P/LP variants found by VCF but NOT in clinical sources:", length(brca2_path_only), "\n")
cat("\nBRCA1 VUS found by VCF but NOT in clinical sources:", length(brca1_vus_only), "\n")
cat("BRCA2 VUS found by VCF but NOT in clinical sources:", length(brca2_vus_only), "\n")



# ========================
# VENN DIAGRAMS
# ========================
library(ggVennDiagram)
library(ggplot2)

brca1_sources <- list(
    Progeny = unique(brca1_progeny_pmbb$PMBB_ID),
    Brad = unique(brca1_brad$PMBB_ID),
    Freeze_2.0 = unique(old_brca1$SampleID)
)

brca2_sources <- list(
    Progeny = unique(brca2_progeny_pmbb$PMBB_ID),
    Brad = unique(brca2_brad$PMBB_ID),
    Freeze_2.0 = unique(old_brca2$SampleID)
)

# Print sample sizes for each source
cat("\n=== BRCA1 Sample Sizes ===\n")
sapply(brca1_sources, length)

cat("\n=== BRCA2 Sample Sizes ===\n")
sapply(brca2_sources, length)

##### BRCA1 VENN #####
brca1_venn <- ggVennDiagram(
    brca1_sources,
    label = "count",           # Show counts in each region
    label_alpha = 0,           # Transparent label background
    label_size = 4,            # Label text size
    set_size = 4,              # Category name text size
    edge_size = 0.5            # Border thickness
) +
    scale_fill_gradient(low = "#F4F4F4", high = "#4472C4") +
    scale_color_manual(values = rep("black", 7)) +  # Black borders
    labs(title = "BRCA1 Pathogenic/Likely Pathogenic Variant Carriers",
         subtitle = paste0("Total unique carriers: ",
                           length(unique(unlist(brca1_sources))))) +
    theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
        plot.subtitle = element_text(hjust = 0.5, size = 12),
        legend.position = "right",
        legend.title = element_text(size = 10),
        plot.margin = margin(10, 10, 10, 10)
    )

ggsave(here("ch", "figures", "brca1_sources_venn.png"),
       plot = brca1_venn,
       width = 10, height = 7,
       dpi = 300, bg = "white")

##### BRCA2 VENN #####
brca2_venn <- ggVennDiagram(
    brca2_sources,
    label = "count",
    label_alpha = 0,
    label_size = 4,
    set_size = 4,
    edge_size = 0.5
) +
    scale_fill_gradient(low = "#F4F4F4", high = "#C55A11") +
    scale_color_manual(values = rep("black", 7)) +
    labs(title = "BRCA2 Pathogenic/Likely Pathogenic Variant Carriers",
         subtitle = paste0("Total unique carriers: ",
                           length(unique(unlist(brca2_sources))))) +
    theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
        plot.subtitle = element_text(hjust = 0.5, size = 12),
        legend.position = "right",
        legend.title = element_text(size = 10),
        plot.margin = margin(10, 10, 10, 10)
    )

ggsave(here("ch", "figures", "brca2_sources_venn.png"),
       plot = brca2_venn,
       width = 10, height = 7,
       dpi = 300, bg = "white")

# Optional: Create a combined summary table
overlap_summary <- data.frame(
    Source = c("Progeny only", "Brad only", "Freeze 2.0 only",
               "Progeny & Brad", "Progeny & Freeze 2.0", "Brad & Freeze 2.0",
               "All three", "Total unique"),
    BRCA1 = c(
        length(setdiff(setdiff(brca1_sources$Progeny, brca1_sources$Brad), brca1_sources$Freeze_2.0)),
        length(setdiff(setdiff(brca1_sources$Brad, brca1_sources$Progeny), brca1_sources$Freeze_2.0)),
        length(setdiff(setdiff(brca1_sources$Freeze_2.0, brca1_sources$Progeny), brca1_sources$Brad)),
        length(setdiff(intersect(brca1_sources$Progeny, brca1_sources$Brad), brca1_sources$Freeze_2.0)),
        length(setdiff(intersect(brca1_sources$Progeny, brca1_sources$Freeze_2.0), brca1_sources$Brad)),
        length(setdiff(intersect(brca1_sources$Brad, brca1_sources$Freeze_2.0), brca1_sources$Progeny)),
        length(Reduce(intersect, brca1_sources)),
        length(unique(unlist(brca1_sources)))
    ),
    BRCA2 = c(
        length(setdiff(setdiff(brca2_sources$Progeny, brca2_sources$Brad), brca2_sources$Freeze_2.0)),
        length(setdiff(setdiff(brca2_sources$Brad, brca2_sources$Progeny), brca2_sources$Freeze_2.0)),
        length(setdiff(setdiff(brca2_sources$Freeze_2.0, brca2_sources$Progeny), brca2_sources$Brad)),
        length(setdiff(intersect(brca2_sources$Progeny, brca2_sources$Brad), brca2_sources$Freeze_2.0)),
        length(setdiff(intersect(brca2_sources$Progeny, brca2_sources$Freeze_2.0), brca2_sources$Brad)),
        length(setdiff(intersect(brca2_sources$Brad, brca2_sources$Freeze_2.0), brca2_sources$Progeny)),
        length(Reduce(intersect, brca2_sources)),
        length(unique(unlist(brca2_sources)))
    )
)

print(overlap_summary)
write.csv(overlap_summary, here("ch", "data", "brca_sources_overlap_summary.csv"), row.names = FALSE)






##### INTERSECTION OF F3 + F2 IDS #####
length(unique(c(brca2_progeny_pmbb$PMBB_ID, brca2_brad$PMBB_ID))) # 608
length(brca2_brad$PMBB_ID) # 592 so mostly the same, so just take the intersection of brad's and the progeny crep sheets

brca2_progeny_brad <- intersect(brca2_progeny_pmbb$PMBB_ID, brca2_brad$PMBB_ID)
brca1_progeny_brad <- intersect(brca1_progen_pmbb$PMBB_ID, brca1_brad$PMBB_ID)

all_brca2_crep_ids <- unique(c(brca2_progeny_brad, old_brca2$SampleID))
length(all_brca2_crep_ids)
# 592

all_brca1_crep_ids <- unique(c(brca1_progeny_brad, old_brca1$SampleID))
length(all_brca1_crep_ids)
# 555

all_brca12_crep_ids <- sort(unique(c(all_brca1_crep_ids, all_brca2_crep_ids)))
length(all_brca12_crep_ids)
# 1412

# ##### INTERSECTION INCLUDING NON CREP - FOR EXCLUSION/CONTROLS #####
#
# all_brca2_ids <- unique(c(
#     brca2_progeny_pmbb$PMBB_ID,
#     brca2_brad$PMBB_ID,
#     old_brca2$SampleID
# ))
# length(all_brca2_ids)
# # 870
#
# all_brca1_ids <- unique(c(
#     brca1_progeny_pmbb$PMBB_ID,
#     brca1_brad$PMBB_ID,
#     old_brca1$SampleID
# ))
# length(all_brca1_ids)
# # 834

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












# progeny_crep <- progeny %>% filter(!(Seen_In_CREP %in% c("No", "Unknown", NA)))
# progeny_nocrep <- progeny %>% filter((Seen_In_CREP %in% c("No", "Unknown", NA)))
# length(unique(progeny_crep$globalid))
# length(progeny_crep$globalid %in% tumor_ids)
# # 1648/1648 have cancer info or tumors...
# length(progeny_nocrep$globalid %in% tumor_ids)
# # 1168
# length(unique(progeny_nocrep$globalid))
# # 1168
# # how
#
# table(progeny_pmbb$Seen_In_CREP)
# dim(progeny_crep)
# # 1648
# progeny_crep_merged <- merge_duplicates(progeny_crep, "SampNum")
#
# progeny_crep_pmbb <- progeny_crep_merged %>%
#     left_join(up, by = "SampNum") %>%
#     filter(!is.na(PMBB_ID))
# dim(progeny_crep_merged)
# # 783 are in the PMBB
#
# brca1_progeny_crep_pmbb <- progeny_crep_pmbb %>% filter(BRCA1_Presence_of_Mutation %in% c("Yes", "Obligate Carrier", "Yes-unconfirmed"))
# dim(brca1_progeny_crep_pmbb) #409
# brca2_progeny_crep_pmbb <- progeny_crep_pmbb %>% filter(BRCA2_Presence_of_Mutation %in% c("Yes", "Obligate Carrier", "Yes-unconfirmed"))
# dim(brca2_progeny_crep_pmbb) #380


### ==== USER INPUTS ====
excel_path <- file.path("ch", "ss", "ch_genes.xlsx")
sheet <- 1
folder_path <- "~/Downloads/All_variants"
output_file <- file.path("ch", "ss", "gene_matches.csv")

genes <- read_excel(excel_path, sheet = sheet, col_names = FALSE)
colnames(genes) <- c("Gene")

file_list <- list.files(folder_path, recursive = TRUE, full.names = TRUE)

genes$Found <- sapply(genes$Gene, function(g) {
    any(grepl(g, file_list, ignore.case = TRUE))
})

genes$Found <- ifelse(genes$Found, 1, 0)

write.csv(genes, output_file, row.names = FALSE)


