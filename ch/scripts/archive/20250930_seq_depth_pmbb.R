# ========================
# Packages
# ========================
library(here)
setwd("/Users/jennyzli/Documents/Nathanson")
source(here("R", "load_packages.R"))
source(here("R", "sample_selection.R"))
source(here("R", "control_selection.R"))

# ============================================================================
# COMBINE ALL VEP
# ============================================================================

# head -n 1 crep.ASXL1.norm.vep.report.csv > crep.all.norm.vep.report.csv
# tail -n +2 -q crep.*.norm.vep.report.csv >> crep.all.norm.vep.report.csv

gene_names <- str_extract(basename(vep_files), "(?<=crep\\.).*(?=\\.norm)")
gene_names <- gene_names[2:74]

all_data <- read.csv(here("ch", "vep", "crep.all.norm.vep.report.csv"))
# 916450     58

# ============================================================================
# GENE CHECK
# ============================================================================
vep_files <- list.files(
    path = here("ch", "vep"),
    pattern = "^crep\\.[^a][^l][^l].*\\.norm\\.vep\\.report\\.csv$",  # Exclude crep.all.norm...
    full.names = TRUE
)

gene_check <- map_dfr(vep_files, function(file) {
    data <- read.csv(file, stringsAsFactors = FALSE)
    unique_genes <- unique(data$Gene)
    filename <- basename(file)
    data.frame(
        Filename = filename,
        N_variants = nrow(data),
        N_unique_genes = length(unique_genes),
        Genes = paste(sort(unique_genes), collapse = ", "),
        stringsAsFactors = FALSE
    )
})
print(gene_check)

write.csv(gene_check,
          here("ch", "data", "gene_check_per_file.csv"),
          row.names = FALSE)
# some have 2 <= genes..should prob filter the all file to only the genes we're looking at

# ============================================================================
# FILTERING
# ============================================================================
test_filtered <- all_data %>%
    mutate(
        filter_reason = NA_character_,
        pass_all = TRUE
    )

# Filter 1: Total read depth >= 20
test_filtered <- test_filtered %>%
    mutate(
        pass_depth = Sample.Depth >= 20,
        filter_reason = if_else(!pass_depth & is.na(filter_reason),
                                "DP<20", filter_reason),
        pass_all = pass_all & pass_depth
    )

# Filter 2: minAD >= 3
test_filtered <- test_filtered %>%
    mutate(
        pass_minAD = Sample.AltDepth >= 3,
        filter_reason = if_else(!pass_minAD & is.na(filter_reason),
                                "minAD<3", filter_reason),
        pass_all = pass_all & pass_minAD
    )

# Filter 3: Indels with >=10 alt reads
test_filtered <- test_filtered %>%
    mutate(
        is_indel = Variant.Class %in% c("insertion", "deletion"),
        pass_indel = !(is_indel & Sample.AltDepth <= 10),
        filter_reason = if_else(!pass_indel & is.na(filter_reason),
                                "Indel_AD<=10", filter_reason),
        pass_all = pass_all & pass_indel
    )

# Filter 4: VAF >= 0.02
test_filtered <- test_filtered %>%
    mutate(
        pass_vaf_low = Sample.AltFrac >= 0.02,
        filter_reason = if_else(!pass_vaf_low & is.na(filter_reason),
                                "VAF<0.02", filter_reason),
        pass_all = pass_all & pass_vaf_low
    )

# Filter 5: VAF <= 0.35 (likely germline)
test_filtered <- test_filtered %>%
    mutate(
        pass_germline = Sample.AltFrac <= 0.35,
        filter_reason = if_else(!pass_germline & is.na(filter_reason),
                                "Likely_germline", filter_reason),
        pass_all = pass_all & pass_germline
    )

# Filter 6: gnomAD AF > 0.01
test_filtered <- test_filtered %>%
    mutate(
        # First convert to character, then handle special cases
        gnomAD.MAX_AF_char = as.character(gnomAD.MAX_AF),
        gnomAD.MAX_AF.clean = case_when(
            is.na(gnomAD.MAX_AF_char) ~ 0,
            gnomAD.MAX_AF_char == "." ~ 0,
            gnomAD.MAX_AF_char == "" ~ 0,
            gnomAD.MAX_AF_char == "NA" ~ 0,
            TRUE ~ suppressWarnings(as.numeric(gnomAD.MAX_AF_char))
        ),
        # Replace any remaining NAs with 0
        gnomAD.MAX_AF.clean = if_else(is.na(gnomAD.MAX_AF.clean), 0, gnomAD.MAX_AF.clean),
        pass_gnomad = gnomAD.MAX_AF.clean <= 0.01,
        filter_reason = if_else(!pass_gnomad & is.na(filter_reason),
                                "gnomAD_AF>0.01", filter_reason),
        pass_all = pass_all & pass_gnomad
    ) %>%
    select(-gnomAD.MAX_AF_char)


# ============================================================================
# PER-GENE SUMMARY TABLE
# ============================================================================
genes <- read.csv(here("ch", "data", "ch_genes.csv"))$ch_genes

per_gene_summary <- test_filtered %>%
    group_by(Gene) %>%
    summarise(
        N_raw = n(),
        N_filtered = sum(pass_all),
        Pass_rate = round(sum(pass_all)/n()*100, 1),
        Median_depth_raw = median(Sample.Depth),
        Median_depth_filtered = median(Sample.Depth[pass_all]),
        Median_VAF_raw = round(median(Sample.AltFrac), 3),
        Median_VAF_filtered = round(median(Sample.AltFrac[pass_all]), 3),
        N_failed_depth = sum(!pass_depth),
        N_failed_minAD = sum(!pass_minAD),
        N_failed_indel = sum(!pass_indel),
        N_failed_vaf_low = sum(!pass_vaf_low),
        N_failed_germline = sum(!pass_germline),
        N_failed_gnomad = sum(!pass_gnomad)
    ) %>%
    arrange(desc(N_filtered)) %>%
    filter(Gene %in% genes)
dim(per_gene_summary)
# 70 14

print(per_gene_summary)

write.csv(per_gene_summary,
          here("ch", "data", "per_gene_filtering_summary.csv"),
          row.names = FALSE)

# ============================================================================
# KEY VISUALIZATIONS
# ============================================================================

# 1. Per-gene filtering results (bar chart)
p1 <- ggplot(per_gene_summary,
             aes(x = reorder(Gene, N_filtered), y = N_filtered)) +
    geom_col(aes(y = N_raw), fill = "gray80", alpha = 0.5) +
    geom_col(fill = "steelblue") +
    coord_flip() +
    theme_minimal() +
    labs(title = "Variants per Gene (Gray = Raw, Blue = Filtered)",
         x = "Gene", y = "Number of Variants")

# 2. Pass rate by gene (lollipop plot)
p2 <- ggplot(per_gene_summary,
             aes(x = reorder(Gene, Pass_rate), y = Pass_rate)) +
    geom_segment(aes(xend = Gene, y = 0, yend = Pass_rate),
                 color = "gray70") +
    geom_point(size = 3, color = "steelblue") +
    coord_flip() +
    theme_minimal() +
    labs(title = "Filtering Pass Rate by Gene",
         x = "Gene", y = "Pass Rate (%)")

# 3. Depth distribution by gene (boxplot)
p3 <- test_filtered %>%
    filter(pass_all) %>%
    ggplot(aes(x = reorder(Gene, Sample.Depth, median),
               y = Sample.Depth)) +
    geom_boxplot(fill = "lightblue", outlier.size = 0.5) +
    coord_flip() +
    theme_minimal() +
    labs(title = "Depth Distribution by Gene (Filtered)",
         x = "Gene", y = "Sample Depth")

# 4. VAF distribution by gene (violin plot)
p4 <- test_filtered %>%
    filter(pass_all) %>%
    ggplot(aes(x = reorder(Gene, Sample.AltFrac, median),
               y = Sample.AltFrac)) +
    geom_violin(fill = "lightcoral") +
    geom_boxplot(width = 0.1, fill = "white", outlier.size = 0.5) +
    coord_flip() +
    theme_minimal() +
    labs(title = "VAF Distribution by Gene (Filtered)",
         x = "Gene", y = "VAF")

# 6. Overall filtering cascade (all genes combined)
filter_cascade <- data.frame(
    Step = c("Raw", "DP>=20", "minAD>=3", "Indel", "VAF>=0.02",
             "Germline", "gnomAD", "Final"),
    N = c(
        nrow(all_data),
        sum(test_filtered$pass_depth),
        sum(test_filtered$pass_depth & test_filtered$pass_minAD),
        sum(test_filtered$pass_depth & test_filtered$pass_minAD &
                test_filtered$pass_indel),
        sum(test_filtered$pass_depth & test_filtered$pass_minAD &
                test_filtered$pass_indel & test_filtered$pass_vaf_low),
        sum(test_filtered$pass_depth & test_filtered$pass_minAD &
                test_filtered$pass_indel & test_filtered$pass_vaf_low &
                test_filtered$pass_germline),
        sum(test_filtered$pass_depth & test_filtered$pass_minAD &
                test_filtered$pass_indel & test_filtered$pass_vaf_low &
                test_filtered$pass_germline & test_filtered$pass_gnomad),
        sum(test_filtered$pass_all)
    )
) %>%
    mutate(
        Percent = round(N / nrow(all_data) * 100, 1),
        N_removed = lag(N, default = nrow(all_data)) - N,
        Percent_removed = round(N_removed / nrow(all_data) * 100, 1)
    )

write.csv(filter_cascade,
          here("ch", "data", "filtering_cascade_summary.csv"),
          row.names = FALSE)

print(filter_cascade)

### save everything
ggsave(here("ch", "figures", "per_gene_variants.png"), p1,
       width = 10, height = 14, dpi = 300)
ggsave(here("ch", "figures", "per_gene_passrate.png"), p2,
       width = 10, height = 14, dpi = 300)
ggsave(here("ch", "figures", "depth_by_gene.png"), p3,
       width = 10, height = 14, dpi = 300)
ggsave(here("ch", "figures", "vaf_by_gene.png"), p4,
       width = 10, height = 14, dpi = 300)
# ggsave(here("ch", "figures", "filter_heatmap.png"), p5,
#        width = 12, height = 14, dpi = 300)
# ggsave(here("ch", "figures", "filtering_cascade.png"), p6,
#        width = 10, height = 8, dpi = 300)

# ============================================================================
# EXPORT DATA
# ============================================================================
final_variants <- test_filtered %>%
    filter(pass_all) %>%
    select(-starts_with("pass_"), -filter_reason, -is_indel)

write.csv(final_variants,
          here("ch", "data", "all_genes_filtered_variants.csv"),
          row.names = FALSE)

write.csv(test_filtered,
          here("ch", "data", "all_genes_with_filter_info.csv"),
          row.names = FALSE)

# ============================================================================
# OVERALL COHORT QUALITY PLOTS - SEPARATE BEFORE/AFTER
# ============================================================================
data_before <- all_data %>% mutate(Status = "Before Filtering")
data_after <- test_filtered %>% filter(pass_all) %>% mutate(Status = "After Filtering")

# ============================================================================
# 1. TOTAL DEPTH DISTRIBUTION
# ============================================================================

p_depth_before <- ggplot(data_before, aes(x = Sample.Depth)) +
    geom_histogram(bins = 50, fill = "gray70", color = "black") +
    geom_vline(xintercept = 20, color = "red", linetype = "dashed", size = 1) +
    annotate("text", x = 20, y = Inf, label = "DP = 20",
             vjust = 1.5, hjust = -0.1, color = "red", size = 4) +
    theme_minimal() +
    labs(title = "Total Depth Distribution - Before Filtering",
         subtitle = sprintf("N = %s variants", format(nrow(data_before), big.mark = ",")),
         x = "Sample Depth", y = "Count") +
    xlim(0, 200)

p_depth_after <- ggplot(data_after, aes(x = Sample.Depth)) +
    geom_histogram(bins = 50, fill = "steelblue", color = "black") +
    geom_vline(xintercept = 20, color = "red", linetype = "dashed", size = 1) +
    annotate("text", x = 20, y = Inf, label = "DP = 20",
             vjust = 1.5, hjust = -0.1, color = "red", size = 4) +
    theme_minimal() +
    labs(title = "Total Depth Distribution - After Filtering",
         subtitle = sprintf("N = %s variants", format(nrow(data_after), big.mark = ",")),
         x = "Sample Depth", y = "Count") +
    xlim(0, 200)

# ============================================================================
# 2. ALT ALLELE DEPTH DISTRIBUTION
# ============================================================================
p_altdepth_before <- ggplot(data_before, aes(x = Sample.AltDepth)) +
    geom_histogram(bins = 50, fill = "gray70", color = "black") +
    geom_vline(xintercept = 3, color = "red", linetype = "dashed", size = 1) +
    geom_vline(xintercept = 10, color = "blue", linetype = "dashed", size = 1) +
    annotate("text", x = 3, y = Inf, label = "all = 3",
             vjust = 1.5, hjust = -0.1, color = "red", size = 3.5) +
    annotate("text", x = 10, y = Inf, label = "Indel = 10",
             vjust = 3, hjust = -0.1, color = "blue", size = 3.5) +
    theme_minimal() +
    labs(title = "Alt Allele Depth Distribution - Before Filtering",
         subtitle = sprintf("N = %s variants", format(nrow(data_before), big.mark = ",")),
         x = "Alt Allele Depth", y = "Count") +
    xlim(0, 100)

p_altdepth_after <- ggplot(data_after, aes(x = Sample.AltDepth)) +
    geom_histogram(bins = 50, fill = "steelblue", color = "black") +
    geom_vline(xintercept = 3, color = "red", linetype = "dashed", size = 1) +
    geom_vline(xintercept = 10, color = "blue", linetype = "dashed", size = 1) +
    annotate("text", x = 3, y = Inf, label = "all = 3",
             vjust = 1.5, hjust = -0.1, color = "red", size = 3.5) +
    annotate("text", x = 10, y = Inf, label = "Indel = 10",
             vjust = 3, hjust = -0.1, color = "blue", size = 3.5) +
    theme_minimal() +
    labs(title = "Alt Allele Depth Distribution - After Filtering",
         subtitle = sprintf("N = %s variants", format(nrow(data_after), big.mark = ",")),
         x = "Alt Allele Depth", y = "Count") +
    xlim(0, 100)

# ============================================================================
# 3. VAF DISTRIBUTION
# ============================================================================

x <- data_before %>% filter(Sample.AltFrac <= 1)
y <- head(data_before %>% filter(Sample.AltFrac == 1), 50)
# write.csv(x, here("ch", "data", "x.csv"))
# write.csv(y, here("ch", "data", "y.csv"))

p_vaf_before <- ggplot(x, aes(x = Sample.AltFrac)) +
    geom_histogram(bins = 60, fill = "gray70", color = "black") +
    geom_vline(xintercept = 0.02, color = "red", linetype = "dashed", size = 1) +
    geom_vline(xintercept = 0.35, color = "orange", linetype = "dashed", size = 1) +
    annotate("text", x = 0.02, y = Inf, label = "0.02",
             vjust = 1.5, hjust = -0.1, color = "red", size = 3.5) +
    annotate("text", x = 0.35, y = Inf, label = "0.35",
             vjust = 1.5, hjust = 1.1, color = "orange", size = 3.5) +
    theme_minimal() +
    labs(title = "Variant Allele Fraction Distribution - Before Filtering",
         subtitle = sprintf("N = %s variants | CHIP range: 0.02 - 0.35",
                            format(nrow(data_before), big.mark = ",")),
         x = "VAF", y = "Count")

p_vaf_after <- ggplot(data_after, aes(x = Sample.AltFrac)) +
    geom_histogram(bins = 60, fill = "steelblue", color = "black") +
    geom_vline(xintercept = 0.02, color = "red", linetype = "dashed", size = 1) +
    geom_vline(xintercept = 0.35, color = "orange", linetype = "dashed", size = 1) +
    annotate("text", x = 0.02, y = Inf, label = "0.02",
             vjust = 1.5, hjust = -0.1, color = "red", size = 3.5) +
    annotate("text", x = 0.35, y = Inf, label = "0.35",
             vjust = 1.5, hjust = 1.1, color = "orange", size = 3.5) +
    theme_minimal() +
    labs(title = "Variant Allele Fraction Distribution - After Filtering",
         subtitle = sprintf("N = %s variants | CHIP range: 0.02 - 0.35",
                            format(nrow(data_after), big.mark = ",")),
         x = "VAF", y = "Count")

ggsave(here("ch", "figures", "vaf_before.png"),
       p_vaf_before, width = 6, height = 6, dpi = 300)
ggsave(here("ch", "figures", "vaf_after.png"),
       p_vaf_after, width = 6, height = 6, dpi = 300)

# ============================================================================
# 4. VARIANT CLASS DISTRIBUTION
# ============================================================================

# Prepare variant class data
variant_class_before <- data_before %>%
    count(Variant.Class) %>%
    mutate(
        Total = sum(n),
        Percent = round(n / Total * 100, 1)
    )

variant_class_after <- data_after %>%
    count(Variant.Class) %>%
    mutate(
        Total = sum(n),
        Percent = round(n / Total * 100, 1)
    )

# Bar plot - Before
p_varclass_before <- ggplot(variant_class_before,
                            aes(x = Variant.Class, y = n)) +
    geom_col(fill = "gray70") +
    geom_text(aes(label = paste0(format(n, big.mark = ","), "\n(", Percent, "%)")),
              vjust = -0.3, size = 3.5) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Variant Class Distribution - Before Filtering",
         subtitle = sprintf("Total: N = %s", format(sum(variant_class_before$n), big.mark = ",")),
         x = "Variant Class", y = "Count") +
    ylim(0, max(variant_class_before$n) * 1.15)

# Bar plot - After
p_varclass_after <- ggplot(variant_class_after,
                           aes(x = Variant.Class, y = n)) +
    geom_col(fill = "steelblue") +
    geom_text(aes(label = paste0(format(n, big.mark = ","), "\n(", Percent, "%)")),
              vjust = -0.3, size = 3.5) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Variant Class Distribution - After Filtering",
         subtitle = sprintf("Total: N = %s", format(sum(variant_class_after$n), big.mark = ",")),
         x = "Variant Class", y = "Count") +
    ylim(0, max(variant_class_after$n) * 1.15)

# Stacked percentage bar - keep this as comparison
variant_class_data <- bind_rows(
    variant_class_before %>% mutate(Status = "Before Filtering"),
    variant_class_after %>% mutate(Status = "After Filtering")
)

p_varclass_pct <- ggplot(variant_class_data,
                         aes(x = Status, y = Percent, fill = Variant.Class)) +
    geom_col(position = "stack") +
    geom_text(aes(label = paste0(Percent, "%")),
              position = position_stack(vjust = 0.5),
              size = 3, color = "white", fontface = "bold") +
    scale_fill_brewer(palette = "Set2") +
    theme_minimal() +
    theme(legend.position = "right") +
    labs(title = "Variant Class Composition (%)",
         x = "", y = "Percentage", fill = "Variant Class")

# ============================================================================
# SUMMARY STATISTICS TABLE
# ============================================================================

summary_stats <- data.frame(
    Metric = c("Total Variants",
               "Median Depth",
               "Median Alt Depth",
               "Median VAF",
               "SNVs",
               "Insertions",
               "Deletions",
               "Other"),
    Before = c(
        nrow(all_data),
        round(median(all_data$Sample.Depth), 1),
        round(median(all_data$Sample.AltDepth), 1),
        round(median(all_data$Sample.AltFrac), 3),
        sum(all_data$Variant.Class == "substitution", na.rm = TRUE),
        sum(all_data$Variant.Class == "insertion", na.rm = TRUE),
        sum(all_data$Variant.Class == "deletion", na.rm = TRUE),
        sum(!all_data$Variant.Class %in% c("substitution", "insertion", "deletion"), na.rm = TRUE)
    ),
    After = c(
        sum(test_filtered$pass_all),
        round(median(test_filtered$Sample.Depth[test_filtered$pass_all]), 1),
        round(median(test_filtered$Sample.AltDepth[test_filtered$pass_all]), 1),
        round(median(test_filtered$Sample.AltFrac[test_filtered$pass_all]), 3),
        sum(test_filtered$Variant.Class[test_filtered$pass_all] == "substitution", na.rm = TRUE),
        sum(test_filtered$Variant.Class[test_filtered$pass_all] == "insertion", na.rm = TRUE),
        sum(test_filtered$Variant.Class[test_filtered$pass_all] == "deletion", na.rm = TRUE),
        sum(!test_filtered$Variant.Class[test_filtered$pass_all] %in%
                c("substitution", "insertion", "deletion"), na.rm = TRUE)
    )
)

print(summary_stats)

write.csv(summary_stats,
          here("ch", "data", "overall_cohort_summary_stats.csv"),
          row.names = FALSE)

# ============================================================================
# SAVE ALL PLOTS SEPARATELY
# ============================================================================
# Depth plots
ggsave(here("ch", "figures", "depth_before.png"),
       p_depth_before, width = 6, height = 6, dpi = 300)
ggsave(here("ch", "figures", "depth_after.png"),
       p_depth_after, width = 6, height = 6, dpi = 300)

# Alt depth plots
ggsave(here("ch", "figures", "altdepth_before.png"),
       p_altdepth_before, width = 6, height = 6, dpi = 300)
ggsave(here("ch", "figures", "altdepth_after.png"),
       p_altdepth_after, width = 6, height = 6, dpi = 300)

# Variant class plots
ggsave(here("ch", "figures", "varclass_before.png"),
       p_varclass_before, width = 6, height = 6, dpi = 300)
ggsave(here("ch", "figures", "varclass_after.png"),
       p_varclass_after, width = 6, height = 6, dpi = 300)
ggsave(here("ch", "figures", "varclass_comparison_pct.png"),
       p_varclass_pct, width = 6, height = 6, dpi = 300)

# unique
length(unique(data_before$Sample.ID))
# 2840 samples - in CREP?





