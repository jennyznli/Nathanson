library(here)
setwd("/Users/jennyzli/Documents/Nathanson")
source(here("R", "config.R"))

library(gridExtra)
library(grid)

# ========================
# LOAD DATA
# ========================
progeny <- read_excel(here("ch", "ss", "brca_carriers_ch_freq_w_seen_in_crep_20251020.xlsx"), sheet = "Data_from_master_table") %>% filter(DNA == "D")
up <- read.csv(here("simplexo", "data", "simplexo_up_map.csv"))
up$SampNum <- as.numeric(up$SampNum)

PMBB_DIR = here("PMBB", "3.0")
OCC <- here(PMBB_DIR, "PMBB-Release-2024-3.0_phenotype_condition_occurrence.txt")
COV <- here(PMBB_DIR, "PMBB-Release-2024-3.0_covariates.txt")
PER <- here(PMBB_DIR, "PMBB-Release-2024-3.0_phenotype_person.txt")
cov <- fread(COV, header = TRUE)
person <- fread(PER, header = TRUE)

# ========================
# PROGENY
# ========================
dim(progeny)
# 2816

progeny_merged <- merge_duplicates(progeny, "SampNum")
progeny_merged$SampNum <- as.numeric(progeny_merged$SampNum)

progeny_pmbb <- progeny_merged %>%
    left_join(up, by = "SampNum") %>%
    filter(!is.na(PMBB_ID))
dim(progeny_pmbb)
# 1278 are in the PMBB

# more relaxed
brca1_progeny_pmbb <- progeny_pmbb %>% filter(BRCA1_Presence_of_Mutation != "No")
dim(brca1_progeny_pmbb) #885
brca2_progeny_pmbb <- progeny_pmbb %>% filter(BRCA2_Presence_of_Mutation != "No")
dim(brca2_progeny_pmbb) #818
# > dim(brca1_progeny_pmbb) #885
# [1] 885  37
# > brca2_progeny_pmbb <- progeny_pmbb %>% filter(BRCA2_Presence_of_Mutation != "No")
# > dim(brca2_progeny_pmbb) #818
# [1] 818  37

# stringent filtering...
brca1_progeny_pmbb <- progeny_pmbb %>% filter(BRCA1_Presence_of_Mutation == "Yes")
dim(brca1_progeny_pmbb) #885
brca2_progeny_pmbb <- progeny_pmbb %>% filter(BRCA2_Presence_of_Mutation == "Yes")
dim(brca2_progeny_pmbb) #818
# > dim(brca1_progeny_pmbb) #885
# [1] 680  37
# > brca2_progeny_pmbb <- progeny_pmbb %>% filter(BRCA2_Presence_of_Mutation == "Yes")
# > dim(brca2_progeny_pmbb) #818
# [1] 605  37

length(intersect(brca1_progeny_pmbb$globalid, brca2_progeny_pmbb$globalid))
# 9 with both ??

progeny_pmbb_cov <- progeny_pmbb %>% left_join(cov, by = c("PMBB_ID" = "person_id"))
unique(progeny_pmbb_cov$Batch)
# both 2 and 1...
table(progeny_pmbb_cov$Batch)
# 1    2
# 1 1277
# ========================
# FREEZE 3 ONLY
# ========================
f3_brca1 <- read_excel(here("ch", "ss", "BW-Regeneron_Workbook_20250627.xlsx"), sheet = "CREP_workbook") %>% filter(Mutation_Gene1 == "BRCA1")
f3_brca1_df <- f3_brca1 %>% left_join(up, by = "VCFID") %>% filter(!is.na(PMBB_ID))
length(unique(f3_brca1_df$PMBB_ID)) # 674 - must be freeze 3.0 only....?
x <- f3_brca1_df %>% left_join(cov, by = c("PMBB_ID" = "person_id"))
f3_brca1_fnd <- read_excel(here("ch", "ss", "BW-Regeneron_Workbook_20250627.xlsx"), sheet = "CREP.found_lines") %>% filter(Gene == "BRCA1") %>% filter(!is.na(Sample.ID))


f3_brca1_pos <- read_excel(here("ch", "ss", "BW-Regeneron_Workbook_20250627.xlsx"), sheet = "CREP.possible_lines") %>% filter(Gene == "BRCA1") %>% filter(!is.na(Sample.ID))
f3_brca1_out <- read_excel(here("ch", "ss", "BW-Regeneron_Workbook_20250627.xlsx"), sheet = "CREP.filtered_out_variant_lines") %>% filter(Gene == "BRCA1") %>% filter(!is.na(Sample.ID))
f3_brca1_und <- read_excel(here("ch", "ss", "BW-Regeneron_Workbook_20250627.xlsx"), sheet = "CREP.undetectable_variants") %>% filter(Mutation_Gene1 == "BRCA1")
f3_brca1_und <- f3_brca1_und %>% left_join(up, by = c("Sample.ID" = "VCFID")) %>% filter(!is.na(PMBB_ID))
f3_brca1_mis <- read_excel(here("ch", "ss", "BW-Regeneron_Workbook_20250627.xlsx"), sheet = "CREP.missing_variants") %>% filter(Mutation_Gene1 == "BRCA1")
f3_brca1_mis <- f3_brca1_mis %>% left_join(up, by = c("Sample.ID" = "VCFID")) %>% filter(!is.na(PMBB_ID))

f3_brca2_fnd <- read_excel(here("ch", "ss", "BW-Regeneron_Workbook_20250627.xlsx"), sheet = "CREP.found_lines") %>% filter(Gene == "BRCA2") %>% filter(!is.na(Sample.ID))

set_f3_brca1 <- unique(f3_brca1_df$match_col)
set_f3_brca1_fnd <- unique(f3_brca1_fnd$Sample.ID)
set_f3_brca1_pos <- unique(f3_brca1_pos$Sample.ID)
set_f3_brca1_out <- unique(f3_brca1_out$Sample.ID)
set_f3_brca1_und <- unique(f3_brca1_und$match_col)
set_f3_brca1_mis <- unique(f3_brca1_mis$match_col)

print(length(set_f3_brca1))
print(length(set_f3_brca1_fnd))
print(length(set_f3_brca1_pos))
print(length(set_f3_brca1_out))
print(length(set_f3_brca1_und))
print(length(set_f3_brca1_mis))
length(unique(c(set_f3_brca1_mis, set_f3_brca1_und, set_f3_brca1_pos, set_f3_brca1_fnd)))

# > print(length(set_f3_brca1))
# [1] 674
# > print(length(set_f3_brca1_fnd))
# [1] 595
# > print(length(set_f3_brca1_pos))
# [1] 95
# > print(length(set_f3_brca1_out))
# [1] 4
# > print(length(set_f3_brca1_und))
# [1] 66
# > print(length(set_f3_brca1_mis))
# [1] 11

length(intersect(brca1_progeny_pmbb$PMBB_ID, f3_brca1_df$PMBB_ID)) # all are intersecting which makes sense?

unique(f3_brca1_fnd$Bio.type)
unique(f3_brca1_fnd$Variant.LoF_level) # 1 2 3
table(f3_brca1_fnd$Variant.LoF_level)
table(f3_brca2_fnd$Variant.LoF_level)
# 1   2   3
# 582  13   1
# > table(f3_brca2_fnd$Variant.LoF_level)
#
# 1   2   3   4
# 558   3   1   6

f3_brca1_fnd_1 <- f3_brca1_fnd %>%
    filter(Gene == "BRCA1") %>%
    filter(Bio.type == "protein_coding") %>%
    filter(Sample.AltFrac != ".", Sample.Depth != ".", Sample.AltDepth != ".") %>%
    mutate(
        Sample.AltFrac = as.numeric(Sample.AltFrac),
        Sample.Depth = as.numeric(Sample.Depth),
        Sample.AltDepth = as.numeric(Sample.AltDepth)
    ) %>%
    filter(!is.na(Sample.AltFrac), !is.na(Sample.Depth), !is.na(Sample.AltDepth))

f3_brca2_fnd_1 <- f3_brca2_fnd %>%
    filter(Gene == "BRCA2") %>%
    filter(Bio.type == "protein_coding") %>%
    filter(Sample.AltFrac != ".", Sample.Depth != ".", Sample.AltDepth != ".") %>%
    mutate(
        Sample.AltFrac = as.numeric(Sample.AltFrac),
        Sample.Depth = as.numeric(Sample.Depth),
        Sample.AltDepth = as.numeric(Sample.AltDepth)
    ) %>%
    filter(!is.na(Sample.AltFrac), !is.na(Sample.Depth), !is.na(Sample.AltDepth))

remove_pattern <- c("synonymous", "upstream", "3_prime_UTR_variant", "5_prime_UTR_variant", "inframe")
f3_brca1_fnd_2 <- f3_brca1_fnd_1 %>% filter(!grepl(paste(remove_pattern, collapse="|"), Variant.Consequence), Variant.Consequence != "intron_variant")
f3_brca2_fnd_2 <- f3_brca2_fnd_1 %>% filter(!grepl(paste(remove_pattern, collapse="|"), Variant.Consequence), Variant.Consequence != "intron_variant")

f3_brca1_fnd_3 <- f3_brca1_fnd_2 %>%
    filter(Sample.AltFrac > 0.3, Sample.AltFrac < 0.7) %>%
    filter(Sample.Depth > 20, Sample.AltDepth > 5)
f3_brca2_fnd_3 <- f3_brca2_fnd_2 %>%
    filter(Sample.AltFrac > 0.3, Sample.AltFrac < 0.7) %>%
    filter(Sample.Depth > 20, Sample.AltDepth > 5)

table(f3_brca1_fnd_3$Variant.LoF_level)
table(f3_brca2_fnd_3$Variant.LoF_level)
# 1   2   3
# 560  10   1
# > table(f3_brca2_fnd_3$Variant.LoF_level)
#
# 1   2   3   4
# 533   2   1   2

# ========================
# ALL PAIRWISE INTERSECTIONS
# ========================

# cat("\n=== PAIRWISE INTERSECTIONS ===\n\n")
#
# # Main workbook vs all others
# cat("Main ∩ Found:         ", length(intersect(set_f3_brca1, set_f3_brca1_fnd)), "\n") # 594
# cat("Main ∩ Possible:      ", length(intersect(set_f3_brca1, set_f3_brca1_pos)), "\n") # 8
# cat("Main ∩ Filtered Out:  ", length(intersect(set_f3_brca1, set_f3_brca1_out)), "\n") # 0
# cat("Main ∩ Undetectable:  ", length(intersect(set_f3_brca1, set_f3_brca1_und)), "\n") # 65
# cat("Main ∩ Missing:       ", length(intersect(set_f3_brca1, set_f3_brca1_mis)), "\n") # 10
#
# # Found vs all others - none
# cat("Found ∩ Possible:     ", length(intersect(set_f3_brca1_fnd, set_f3_brca1_pos)), "\n")
# cat("Found ∩ Filtered Out: ", length(intersect(set_f3_brca1_fnd, set_f3_brca1_out)), "\n")
# cat("Found ∩ Undetectable: ", length(intersect(set_f3_brca1_fnd, set_f3_brca1_und)), "\n")
# cat("Found ∩ Missing:      ", length(intersect(set_f3_brca1_fnd, set_f3_brca1_mis)), "\n")
#
# # Possible vs remaining
# cat("Possible ∩ Filtered Out: ", length(intersect(set_f3_brca1_pos, set_f3_brca1_out)), "\n")
# cat("Possible ∩ Undetectable: ", length(intersect(set_f3_brca1_pos, set_f3_brca1_und)), "\n") # 1
# cat("Possible ∩ Missing:      ", length(intersect(set_f3_brca1_pos, set_f3_brca1_mis)), "\n") # 3
# length(set_f3_brca1_und)
# length(set_f3_brca1_pos)
#
# # Filtered Out vs remaining - none
# cat("Filtered Out ∩ Undetectable: ", length(intersect(set_f3_brca1_out, set_f3_brca1_und)), "\n")
# cat("Filtered Out ∩ Missing:      ", length(intersect(set_f3_brca1_out, set_f3_brca1_mis)), "\n")
#
# # Undetectable vs Missing - none
# cat("Undetectable ∩ Missing:      ", length(intersect(set_f3_brca1_und, set_f3_brca1_mis)), "\n")
#
# cat("Total unique across all except out sets: ", length(unique(c(set_f3_brca1, set_f3_brca1_fnd, set_f3_brca1_pos,
#                                                       set_f3_brca1_und, set_f3_brca1_mis))), "\n")
# # so i think that we first split found vs possibly
# # and those possible, we looked at which ones matched or didn't
#
# # so the main workbook has all categories except for filtered out
# # so what i think happened:
# # get found vs. possible
# # from the possible we split it into udetectable and missing
# # anyway i'm pretty sure we just use the work book itself...
#
# f3_brca2 <- read_excel(here("ch", "ss", "BW-Regeneron_Workbook_20250627.xlsx"), sheet = "CREP_workbook") %>% filter(Mutation_Gene1 == "BRCA2")
# f3_brca2_df <- f3_brca2 %>% left_join(up, by = "VCFID") %>% filter(!is.na(PMBB_ID))
# f3_brca2_id <- unique(f3_brca2_df$PMBB_ID)
# f3_brca1_id <- unique(f3_brca1_df$PMBB_ID)
# length(f3_brca2_id)
# length(f3_brca1_id)

# ========================
# FREEZE 2 ONLY
# ========================
f2_brca1 <- read_excel(here("ch", "ss", "BRCA1.BRCA2 P.LP_04.11.21.xlsx"), sheet = "BRCA1")
f2_brca2 <- read_excel(here("ch", "ss", "BRCA1.BRCA2 P.LP_04.11.21.xlsx"), sheet = "BRCA2")
# dim(f2_brca1)
# [1] 153  15
# > dim(f2_brca2)
# [1] 263  10

# > range(f2_brca1$ALT_AlleleFrac)
# [1] 0.259 0.684
range(f2_brca2$ALT_AlleleFrac)
# 0.185

f2_brca1_fil <- f2_brca1 %>% filter(Total_Depth > 20, ALT_AlleleDepth > 5, ALT_AlleleFrac > 0.3, ALT_AlleleFrac < 0.7)
f2_brca2_fil <- f2_brca2 %>% filter(Total_Depth > 20, ALT_AlleleDepth > 5, ALT_AlleleFrac > 0.3, ALT_AlleleFrac < 0.7)
length(unique(f2_brca1_fil$SampleID))
length(unique(f2_brca2_fil$SampleID))
# > length(unique(f2_brca1_fil$SampleID))
# [1] 148
# > length(unique(f2_brca2_fil$SampleID))
# [1] 256

# ========================
# FILTERING FUNCTION
# ========================
process_brca_data <- function(data, gene_name) {
    remove_pattern <- c("synonymous", "upstream", "3_prime_UTR_variant",
                        "5_prime_UTR_variant", "inframe")

    step1 <- data %>%
        filter(Gene == gene_name) %>%
        filter(Bio.type == "protein_coding") %>%
        filter(Sample.AltFrac != ".", Sample.Depth != ".", Sample.AltDepth != ".") %>%
        mutate(
            Sample.AltFrac = as.numeric(Sample.AltFrac),
            Sample.Depth = as.numeric(Sample.Depth),
            Sample.AltDepth = as.numeric(Sample.AltDepth)
        ) %>%
        filter(!is.na(Sample.AltFrac), !is.na(Sample.Depth), !is.na(Sample.AltDepth))

    step2 <- step1 %>%
        filter(!grepl(paste(remove_pattern, collapse="|"), Variant.Consequence),
               Variant.Consequence != "intron_variant")

    step3 <- step2 %>%
        filter(Sample.AltFrac > 0.3, Sample.AltFrac < 0.7,
               Sample.Depth > 20, Sample.AltDepth > 5)

    return(list(step1 = step1, step2 = step2, step3 = step3))
}

# ========================
# PLOTTING FUNCTION
# ========================
create_diagnostic_pdf <- function(step1, step2, step3, gene_name, output_file) {

    plots <- list()
    metrics <- list(
        list(col = "Sample.Depth", name = "Total Depth", xlab = "Total Depth (DP)",
             xlim = c(0, 200), vlines = 20, color = "steelblue"),
        list(col = "Sample.AltDepth", name = "Alt Depth", xlab = "Alt Depth (AD)",
             xlim = c(0, 100), vlines = 5, color = "darkgreen"),
        list(col = "Sample.AltFrac", name = "VAF", xlab = "Variant Allele Fraction",
             xlim = c(0, 1), vlines = c(0.3, 0.7), color = "purple")
    )

    steps <- list(
        list(data = step1, label = "Initial"),
        list(data = step2, label = "After VC Filter"),
        list(data = step3, label = "Final")
    )

    # Change xlim for final VAF plot
    metrics_final <- metrics
    metrics_final[[3]]$xlim <- c(0.3, 0.7)

    idx <- 1
    for (s in 1:3) {
        step_data <- steps[[s]]$data
        step_label <- steps[[s]]$label
        n <- length(unique(step_data$Sample.ID))

        # Use special xlim for final VAF
        metrics_to_use <- if (s == 3) metrics_final else metrics

        for (m in metrics_to_use) {
            plots[[idx]] <- ggplot(step_data, aes(x = .data[[m$col]])) +
                geom_histogram(bins = 100, fill = m$color, alpha = 0.7) +
                geom_vline(xintercept = m$vlines, color = "red",
                           linetype = "dashed", size = 1) +
                labs(title = paste0(step_label, ": ", m$name),
                     subtitle = paste("n =", n, "samples"),
                     x = m$xlab, y = "Count") +
                xlim(m$xlim) +
                theme_minimal() +
                theme(plot.title = element_text(size = 10),
                      plot.subtitle = element_text(size = 8))
            idx <- idx + 1
        }
    }

    pdf(output_file, width = 15, height = 12)
    do.call(grid.arrange, c(plots, ncol = 3,
                            top = textGrob(paste(gene_name, "Filtering Pipeline"),
                                           gp = gpar(fontsize = 16, fontface = "bold"))))
    dev.off()
    cat("✓ Saved:", output_file, "\n")
}

# ========================
# READ DATA
# ========================
all_brca1 <- read.csv(file.path("ch", "data", "brca_vep.BRCA1.vep.report.csv"))
all_brca2 <- read.csv(file.path("ch", "data", "brca_vep.BRCA2.vep.report.csv"))
cons_brca1 <- read.csv(here("ch", "data", "PMBB-Release-2024-3.0_genetic_exome_BRCA1_NF.consented_pmbb-only.norm.vep.report.csv"))
cons_brca2 <- read.csv(here("ch", "data", "PMBB-Release-2024-3.0_genetic_exome_BRCA2_NF.consented_pmbb-only.norm.vep.report.csv"))

cat("  all_brca1:", length(unique(all_brca1$Sample.ID)), "\n")
cat("  all_brca2:", length(unique(all_brca2$Sample.ID)), "\n")
cat("  cons_brca1:", length(unique(cons_brca1$Sample.ID)), "\n")
cat("  cons_brca2:", length(unique(cons_brca2$Sample.ID)), "\n\n")

# ========================
# PROCESS DATA
# ========================
all_b1 <- process_brca_data(all_brca1, "BRCA1")
all_b2 <- process_brca_data(all_brca2, "BRCA2")
cons_b1 <- process_brca_data(cons_brca1, "BRCA1")
cons_b2 <- process_brca_data(cons_brca2, "BRCA2")

# ========================
# PRINT SUMMARY STATISTICS
# ========================
cat("ALL BRCA1:\n")
cat("  Step 1 (Initial + QC):", length(unique(all_b1$step1$Sample.ID)), "\n")
cat("  Step 2 (+ Consequence):", length(unique(all_b1$step2$Sample.ID)), "\n")
cat("  Step 3 (+ VAF/Depth):", length(unique(all_b1$step3$Sample.ID)), "\n")
cat("  LoF levels:", paste(names(table(all_b1$step3$Variant.LoF_level)), "=",
                           table(all_b1$step3$Variant.LoF_level), collapse = ", "), "\n\n")

cat("ALL BRCA2:\n")
cat("  Step 1 (Initial + QC):", length(unique(all_b2$step1$Sample.ID)), "\n")
cat("  Step 2 (+ Consequence):", length(unique(all_b2$step2$Sample.ID)), "\n")
cat("  Step 3 (+ VAF/Depth):", length(unique(all_b2$step3$Sample.ID)), "\n")
cat("  LoF levels:", paste(names(table(all_b2$step3$Variant.LoF_level)), "=",
                           table(all_b2$step3$Variant.LoF_level), collapse = ", "), "\n\n")

cat("CONSENTED BRCA1:\n")
cat("  Step 1 (Initial + QC):", length(unique(cons_b1$step1$Sample.ID)), "\n")
cat("  Step 2 (+ Consequence):", length(unique(cons_b1$step2$Sample.ID)), "\n")
cat("  Step 3 (+ VAF/Depth):", length(unique(cons_b1$step3$Sample.ID)), "\n")
cat("  LoF levels:", paste(names(table(cons_b1$step3$Variant.LoF_level)), "=",
                           table(cons_b1$step3$Variant.LoF_level), collapse = ", "), "\n\n")

cat("CONSENTED BRCA2:\n")
cat("  Step 1 (Initial + QC):", length(unique(cons_b2$step1$Sample.ID)), "\n")
cat("  Step 2 (+ Consequence):", length(unique(cons_b2$step2$Sample.ID)), "\n")
cat("  Step 3 (+ VAF/Depth):", length(unique(cons_b2$step3$Sample.ID)), "\n")
cat("  LoF levels:", paste(names(table(cons_b2$step3$Variant.LoF_level)), "=",
                           table(cons_b2$step3$Variant.LoF_level), collapse = ", "), "\n\n")

# ========================
# DATA QUALITY CHECKS
# ========================
cat("ALL BRCA1 - VAF > 1 issues:\n")
cat("  Rows with VAF > 1:", sum(all_b1$step1$Sample.AltFrac > 1, na.rm = TRUE),
    "(", round(sum(all_b1$step1$Sample.AltFrac > 1, na.rm = TRUE) / nrow(all_b1$step1) * 100, 3), "%)\n\n")

cat("ALL BRCA2 - VAF > 1 issues:\n")
cat("  Rows with VAF > 1:", sum(all_b2$step1$Sample.AltFrac > 1, na.rm = TRUE),
    "(", round(sum(all_b2$step1$Sample.AltFrac > 1, na.rm = TRUE) / nrow(all_b2$step1) * 100, 3), "%)\n\n")

cat("CONSENTED BRCA1 - VAF > 1 issues:\n")
cat("  Rows with VAF > 1:", sum(cons_b1$step1$Sample.AltFrac > 1, na.rm = TRUE),
    "(", round(sum(cons_b1$step1$Sample.AltFrac > 1, na.rm = TRUE) / nrow(cons_b1$step1) * 100, 3), "%)\n\n")

cat("CONSENTED BRCA2 - VAF > 1 issues:\n")
cat("  Rows with VAF > 1:", sum(cons_b2$step1$Sample.AltFrac > 1, na.rm = TRUE),
    "(", round(sum(cons_b2$step1$Sample.AltFrac > 1, na.rm = TRUE) / nrow(cons_b2$step1) * 100, 3), "%)\n\n")

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

write.csv(all_b1$step3, file.path("ch", "data", "brca1_all_filtered.csv"), row.names = FALSE)
write.csv(all_b2$step3, file.path("ch", "data", "brca2_all_filtered.csv"), row.names = FALSE)
write.csv(cons_b1$step3, file.path("ch", "data", "brca1_cons_filtered.csv"), row.names = FALSE)
write.csv(cons_b2$step3, file.path("ch", "data", "brca2_cons_filtered.csv"), row.names = FALSE)


