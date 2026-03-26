library(here)
setwd("/Users/jennyzli/Documents/Nathanson")
source(here("R", "config.R"))
source(here("ch", "config.R"))
library(data.table, quietly=T)


# Get all unique genes called across all samples and callers
all_genes_called <- map_dfr(sample_ids, function(id) {
    map_dfr(names(vep[[id]]), function(caller) {
        vep[[id]][[caller]] %>%
            select(Gene) %>%
            distinct() %>%
            mutate(sample = id, caller = caller)
    })
})

# Unique genes overall
unique_genes_overall <- all_genes_called %>%
    pull(Gene) %>%
    unique()

# How many are CH genes
ch_genes_found_overall <- intersect(unique_genes_overall, ch_genes)
non_ch_genes <- setdiff(unique_genes_overall, ch_genes)

# Summary
cat("\n=== GENE SUMMARY ===\n")
cat(paste("Total unique genes called across all samples/callers:", length(unique_genes_overall), "\n"))
cat(paste("CH genes detected:", length(ch_genes_found_overall), "\n"))
cat(paste("Non-CH genes:", length(non_ch_genes), "\n"))
cat(paste("CH genes in list:", length(ch_genes), "\n"))
cat(paste("CH genes NOT detected:", length(ch_genes) - length(ch_genes_found_overall), "\n"))
cat(paste("Percentage of called genes that are CH genes:",
          round(length(ch_genes_found_overall) / length(unique_genes_overall) * 100, 1), "%\n"))


# =============================================================================
# FIX VEP DATA TYPES - RUN THIS ONCE AT THE BEGINNING
# =============================================================================

library(tidyverse)

# Function to standardize column types in a dataframe
standardize_vep_columns <- function(df) {
    df %>%
        mutate(
            # Numeric columns
            Start = as.numeric(Start),
            STRAND = as.numeric(STRAND),
            Sample.Depth = as.numeric(Sample.Depth),
            Sample.AltDepth = as.numeric(Sample.AltDepth),
            Sample.AltFrac = as.numeric(Sample.AltFrac),
            Variant.LoF_level = as.numeric(Variant.LoF_level),

            # Character columns (in case they're factors)
            Chr = as.character(Chr),
            REF = as.character(REF),
            ALT = as.character(ALT),
            Gene = as.character(Gene),
            Sample.ID = as.character(Sample.ID)
        )
}

# Apply to all dataframes in vep
cat("Standardizing data types in vep object...\n")

vep_fixed <- map(names(vep), function(sample_id) {
    cat("Processing sample:", sample_id, "\n")

    sample_data <- map(names(vep[[sample_id]]), function(caller) {
        cat("  - Caller:", caller, "\n")
        standardize_vep_columns(vep[[sample_id]][[caller]])
    })

    names(sample_data) <- names(vep[[sample_id]])
    return(sample_data)
})

names(vep_fixed) <- names(vep)

# Replace vep with the fixed version
vep <- vep_fixed

cat("Data type standardization complete!\n")

# Verify it worked
cat("\nColumn types for first sample, first caller:\n")
str(vep[[1]][[1]])

# =============================================================================
# SECTION 1: BASIC VARIANT CALLING STATISTICS
# =============================================================================

library(tidyverse)
library(patchwork)

# Create output directory for this analysis
dir.create(file.path("ch", "figures", "section1_basic_stats"),
           recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------------------------------
# 1.1: Total variants per caller across all samples
# -----------------------------------------------------------------------------

# You already have variant_counts from your code
# Let's expand it with more detail

variant_summary <- variant_counts %>%
    group_by(caller) %>%
    summarise(
        total_variants = sum(n_variants),
        mean_per_sample = mean(n_variants),
        median_per_sample = median(n_variants),
        sd_per_sample = sd(n_variants),
        min_per_sample = min(n_variants),
        max_per_sample = max(n_variants)
    ) %>%
    arrange(desc(total_variants))

print("Variant calling summary by caller:")
print(variant_summary)
# caller  total_variants mean_per_sample median_per_sample sd_per_sample min_per_sample
# <chr>            <int>           <dbl>             <dbl>         <dbl>          <int>
#     1 vardict           5987            748.              706.         323.             396
# 2 mutect2           1815            227.              228.          21.6            198
# 3 pmbb                 0              0                 0            0                0

write_csv(variant_summary,
          file.path("ch", "data", "variant_summary_by_caller.csv"))

# -----------------------------------------------------------------------------
# 1.2: Variants per caller - visualizations
# -----------------------------------------------------------------------------

# Bar plot: Total variants by caller
p1_1 <- ggplot(variant_summary, aes(x = reorder(caller, -total_variants),
                                    y = total_variants, fill = caller)) +
    geom_col() +
    geom_text(aes(label = scales::comma(total_variants)),
              vjust = -0.5, size = 4) +
    labs(title = "Total Variants Called Across All Samples",
         x = "Caller",
         y = "Total Number of Variants") +
    theme_bw() +
    theme(legend.position = "none",
          axis.text.x = element_text(size = 12))

# Box plot: Variants per sample distribution
p1_2 <- ggplot(variant_counts, aes(x = caller, y = n_variants, fill = caller)) +
    geom_boxplot(alpha = 0.7) +
    geom_jitter(width = 0.2, alpha = 0.5, size = 2) +
    labs(title = "Variants per Sample Distribution",
         x = "Caller",
         y = "Number of Variants per Sample") +
    theme_bw() +
    theme(legend.position = "none")

# Faceted bar plot: Variants per sample
p1_3 <- ggplot(variant_counts, aes(x = caller, y = n_variants, fill = caller)) +
    geom_col() +
    facet_wrap(~sample, ncol = 3, scales = "free_y") +
    labs(title = "Variants Called per Sample and Caller",
         x = "Caller",
         y = "Number of Variants") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none")

# Combine plots
p1_combined <- p1_1 / p1_2

ggsave(file.path("ch", "figures", "section1_basic_stats", "variants_overview.pdf"),
       p1_combined, width = 5, height = 10)

ggsave(file.path("ch", "figures", "section1_basic_stats", "variants_per_sample.pdf"),
       p1_3, width = 10, height = 10)

# -----------------------------------------------------------------------------
# 1.3: Variant type analysis (SNVs vs INDELs)
# -----------------------------------------------------------------------------

# Extract variant types from all samples
variant_types <- map_dfr(sample_ids, function(id) {
    map_dfr(names(vep[[id]]), function(caller) {
        vep[[id]][[caller]] %>%
            mutate(
                sample = id,
                caller = caller,
                # Classify variant type based on REF and ALT length
                variant_type = case_when(
                    nchar(as.character(REF)) == 1 & nchar(as.character(ALT)) == 1 ~ "SNV",
                    nchar(as.character(REF)) > nchar(as.character(ALT)) ~ "Deletion",
                    nchar(as.character(REF)) < nchar(as.character(ALT)) ~ "Insertion",
                    TRUE ~ "Complex"
                ),
                indel_length = abs(nchar(as.character(REF)) - nchar(as.character(ALT)))
            ) %>%
            select(sample, caller, variant_type, indel_length)
    })
})

# Summary by variant type
variant_type_summary <- variant_types %>%
    count(caller, variant_type) %>%
    group_by(caller) %>%
    mutate(
        total = sum(n),
        percentage = n / total * 100
    )

print("Variant type distribution:")
print(variant_type_summary)
# Groups:   caller [2]
# caller  variant_type     n total percentage
# <chr>   <chr>        <int> <int>      <dbl>
#     1 mutect2 Complex         19  1815       1.05
# 2 mutect2 Deletion       714  1815      39.3
# 3 mutect2 Insertion      427  1815      23.5
# 4 mutect2 SNV            655  1815      36.1
# 5 vardict Complex         83  5987       1.39
# 6 vardict Deletion      2268  5987      37.9
# 7 vardict Insertion     2302  5987      38.4
# 8 vardict SNV           1334  5987      22.3

# Save
write_csv(variant_type_summary,
          file.path("ch", "data", "variant_type_summary.csv"))

# Stacked bar plot: Variant types by caller
p1_4 <- ggplot(variant_type_summary, aes(x = caller, y = n, fill = variant_type)) +
    geom_col(position = "stack") +
    geom_text(aes(label = paste0(round(percentage, 1), "%")),
              position = position_stack(vjust = 0.5)) +
    labs(title = "Variant Type Distribution by Caller",
         x = "Caller",
         y = "Number of Variants",
         fill = "Variant Type") +
    theme_bw()

# Side-by-side comparison
p1_5 <- ggplot(variant_type_summary, aes(x = variant_type, y = n, fill = caller)) +
    geom_col(position = "dodge") +
    labs(title = "SNVs vs INDELs by Caller",
         x = "Variant Type",
         y = "Number of Variants") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

# INDEL length distribution
p1_6 <- variant_types %>%
    filter(variant_type %in% c("Insertion", "Deletion"), indel_length <= 50) %>%
    ggplot(aes(x = indel_length, fill = caller)) +
    geom_histogram(bins = 50, alpha = 0.6, position = "identity") +
    facet_wrap(~variant_type, scales = "free_y") +
    labs(title = "INDEL Length Distribution by Caller",
         x = "INDEL Length (bp)",
         y = "Count") +
    theme_bw()

ggsave(file.path("ch", "figures", "section1_basic_stats", "variant_types.pdf"),
       p1_4 / p1_5, width = 10, height = 10)

ggsave(file.path("ch", "figures", "section1_basic_stats", "indel_lengths.pdf"),
       p1_6, width = 10, height = 6)

# # -----------------------------------------------------------------------------
# # Print summary
# # -----------------------------------------------------------------------------
# cat("\n=== SECTION 1 COMPLETE ===\n")
# cat("Files created:\n")
# cat("- ch/figures/section1_basic_stats/variants_overview.pdf\n")
# cat("- ch/figures/section1_basic_stats/variants_per_sample.pdf\n")
# cat("- ch/figures/section1_basic_stats/variant_types.pdf\n")
# cat("- ch/figures/section1_basic_stats/indel_lengths.pdf\n")
# cat("- ch/data/variant_summary_by_caller.csv\n")
# cat("- ch/data/variant_type_summary.csv\n")

# =============================================================================
# SECTION 2: CONCORDANCE ANALYSIS
# =============================================================================

library(tidyverse)
library(ggVennDiagram)
library(patchwork)
library(UpSetR)  # For more complex overlap visualization

# Create output directory
dir.create(file.path("ch", "figures", "section2_concordance"),
           recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------------------------------
# 2.1: Create combined dataset with variant IDs
# -----------------------------------------------------------------------------

# You already have create_variant_id function
# Let's create a master dataset of all variants with caller info

combined_variants <- map_dfr(sample_ids, function(id) {
    map_dfr(names(vep[[id]]), function(caller) {
        create_variant_id(vep[[id]][[caller]]) %>%
            mutate(
                sample = id,
                caller = caller
            ) %>%
            select(sample, caller, variant_id, Gene, Chr, Start, REF, ALT)
    })
})

# Add concordance information
combined_variants_with_concordance <- combined_variants %>%
    group_by(sample, variant_id) %>%
    mutate(
        n_callers = n_distinct(caller),
        callers_list = paste(sort(unique(caller)), collapse = ", ")
    ) %>%
    ungroup()

# Save this for later use
write_csv(combined_variants_with_concordance,
          file.path("ch", "data", "combined_variants_all.csv"))

# -----------------------------------------------------------------------------
# 2.2: Overall concordance statistics
# -----------------------------------------------------------------------------

# Count variants by concordance level
concordance_summary <- combined_variants_with_concordance %>%
    distinct(sample, variant_id, .keep_all = TRUE) %>%
    count(callers_list, n_callers) %>%
    arrange(desc(n))

print("Concordance Summary:")
print(concordance_summary)
# callers_list     n_callers     n
# <chr>                <int> <int>
#     1 vardict                  1  5233
# 2 mutect2                  1  1090
# 3 mutect2, vardict         2   725

# Calculate concordance percentages
total_unique_variants <- sum(concordance_summary$n)
concordance_summary <- concordance_summary %>%
    mutate(
        percentage = n / total_unique_variants * 100,
        concordance_category = case_when(
            n_callers == 3 ~ "All 3 callers",
            n_callers == 2 ~ "2 callers",
            n_callers == 1 ~ "1 caller only"
        )
    )

print("Concordance by category:")
concordance_summary %>%
    group_by(concordance_category) %>%
    summarise(
        n_variants = sum(n),
        percentage = sum(percentage)
    ) %>%
    print()
# concordance_category n_variants percentage
# <chr>                     <int>      <dbl>
#     1 1 caller only              6323       89.7
# 2 2 callers                   725       10.3

write_csv(concordance_summary,
          file.path("ch", "data", "concordance_summary.csv"))

# -----------------------------------------------------------------------------
# 2.3: Concordance visualization - pie chart
# -----------------------------------------------------------------------------

concordance_by_category <- concordance_summary %>%
    group_by(concordance_category) %>%
    summarise(n_variants = sum(n)) %>%
    mutate(
        percentage = n_variants / sum(n_variants) * 100,
        label = paste0(concordance_category, "\n",
                       n_variants, " variants\n",
                       round(percentage, 1), "%")
    )

p2_1 <- ggplot(concordance_by_category,
               aes(x = "", y = n_variants, fill = concordance_category)) +
    geom_col(width = 1, color = "white", size = 1) +
    coord_polar(theta = "y") +
    geom_text(aes(label = label),
              position = position_stack(vjust = 0.5),
              size = 4) +
    labs(title = "Variant Concordance Across All Samples",
         fill = "Concordance Level") +
    theme_void() +
    theme(legend.position = "right")

# Bar plot version
p2_2 <- ggplot(concordance_summary,
               aes(x = reorder(callers_list, -n), y = n, fill = callers_list)) +
    geom_col() +
    geom_text(aes(label = paste0(n, "\n(", round(percentage, 1), "%)")),
              vjust = -0.3, size = 3) +
    labs(title = "Variants by Caller Combination",
         x = "Caller Combination",
         y = "Number of Variants") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none")

ggsave(file.path("ch", "figures", "section2_concordance", "concordance_overview.pdf"),
       p2_1 / p2_2, width = 10, height = 12)

# -----------------------------------------------------------------------------
# 2.4: Sample-level concordance
# -----------------------------------------------------------------------------

# Concordance per sample
sample_concordance <- combined_variants_with_concordance %>%
    distinct(sample, variant_id, .keep_all = TRUE) %>%
    count(sample, n_callers) %>%
    group_by(sample) %>%
    mutate(
        total_variants = sum(n),
        percentage = n / total_variants * 100
    ) %>%
    ungroup()

print("Concordance by sample:")
print(sample_concordance)

p2_3 <- ggplot(sample_concordance,
               aes(x = sample, y = n, fill = factor(n_callers))) +
    geom_col(position = "fill") +
    scale_y_continuous(labels = scales::percent) +
    labs(title = "Concordance Rate by Sample",
         x = "Sample",
         y = "Percentage of Variants",
         fill = "Number of\nCallers") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

p2_4 <- ggplot(sample_concordance,
               aes(x = sample, y = n, fill = factor(n_callers))) +
    geom_col(position = "stack") +
    labs(title = "Variant Count by Sample and Concordance",
         x = "Sample",
         y = "Number of Variants",
         fill = "Number of\nCallers") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path("ch", "figures", "section2_concordance", "concordance_by_sample.pdf"),
       p2_3 / p2_4, width = 12, height = 10)

# -----------------------------------------------------------------------------
# 2.5: Pairwise concordance between callers
# -----------------------------------------------------------------------------

# Calculate Jaccard index for each pair of callers
pairwise_concordance <- map_dfr(sample_ids, function(sample_id) {

    callers <- names(vep[[sample_id]])

    # Get variant sets
    variant_sets <- map(callers, function(caller) {
        create_variant_id(vep[[sample_id]][[caller]]) %>%
            pull(variant_id)
    })
    names(variant_sets) <- callers

    # Calculate pairwise overlaps
    results <- list()
    for (i in 1:(length(callers)-1)) {
        for (j in (i+1):length(callers)) {
            caller1 <- callers[i]
            caller2 <- callers[j]

            set1 <- variant_sets[[caller1]]
            set2 <- variant_sets[[caller2]]

            intersection <- length(intersect(set1, set2))
            union <- length(union(set1, set2))

            jaccard <- intersection / union

            results[[length(results) + 1]] <- tibble(
                sample = sample_id,
                caller1 = caller1,
                caller2 = caller2,
                n_caller1 = length(set1),
                n_caller2 = length(set2),
                n_shared = intersection,
                n_union = union,
                jaccard_index = jaccard,
                pct_caller1_in_caller2 = intersection / length(set1) * 100,
                pct_caller2_in_caller1 = intersection / length(set2) * 100
            )
        }
    }

    bind_rows(results)
})

print("Pairwise concordance:")
print(pairwise_concordance)

write_csv(pairwise_concordance,
          file.path("ch", "data", "pairwise_concordance.csv"))

# Visualize pairwise concordance
p2_5 <- pairwise_concordance %>%
    mutate(comparison = paste(caller1, "vs", caller2)) %>%
    ggplot(aes(x = sample, y = jaccard_index, fill = comparison)) +
    geom_col(position = "dodge") +
    labs(title = "Pairwise Concordance (Jaccard Index) by Sample",
         subtitle = "Jaccard Index = Intersection / Union",
         x = "Sample",
         y = "Jaccard Index",
         fill = "Comparison") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Average concordance across samples
avg_concordance <- pairwise_concordance %>%
    group_by(caller1, caller2) %>%
    summarise(
        avg_jaccard = mean(jaccard_index),
        avg_pct_shared = mean((pct_caller1_in_caller2 + pct_caller2_in_caller1) / 2),
        .groups = "drop"
    )

print("Average concordance across samples:")
print(avg_concordance)
# caller1 caller2 avg_jaccard avg_pct_shared
# <chr>   <chr>         <dbl>          <dbl>
#     1 mutect2 pmbb          0              NaN
# 2 mutect2 vardict       0.110           26.6
# 3 pmbb    vardict       0              NaN

p2_6 <- avg_concordance %>%
    mutate(comparison = paste(caller1, "vs", caller2)) %>%
    ggplot(aes(x = comparison, y = avg_jaccard, fill = comparison)) +
    geom_col() +
    geom_text(aes(label = round(avg_jaccard, 3)), vjust = -0.5) +
    labs(title = "Average Pairwise Concordance Across All Samples",
         x = "Caller Comparison",
         y = "Average Jaccard Index") +
    ylim(0, 1) +
    theme_bw() +
    theme(legend.position = "none")

ggsave(file.path("ch", "figures", "section2_concordance", "pairwise_concordance.pdf"),
       p2_5 / p2_6, width = 10, height = 10)

# -----------------------------------------------------------------------------
# 2.6: What's unique to each caller?
# -----------------------------------------------------------------------------

unique_to_caller <- combined_variants_with_concordance %>%
    filter(n_callers == 1) %>%
    count(caller, name = "n_unique") %>%
    arrange(desc(n_unique))

print("Variants unique to each caller:")
print(unique_to_caller)
# caller  n_unique
# <chr>      <int>
#     1 vardict     5262
# 2 mutect2     1090

# Which genes have unique variants?
unique_genes_by_caller <- combined_variants_with_concordance %>%
    filter(n_callers == 1) %>%
    count(caller, Gene, name = "n_unique_variants") %>%
    group_by(caller) %>%
    arrange(desc(n_unique_variants)) %>%
    slice_head(n = 20) %>%
    ungroup()

p2_7 <- ggplot(unique_genes_by_caller,
               aes(x = reorder(Gene, n_unique_variants),
                   y = n_unique_variants, fill = caller)) +
    geom_col() +
    coord_flip() +
    facet_wrap(~caller, scales = "free_y") +
    labs(title = "Top 20 Genes with Caller-Unique Variants",
         x = "Gene",
         y = "Number of Unique Variants") +
    theme_bw() +
    theme(legend.position = "none")

ggsave(file.path("ch", "figures", "section2_concordance", "unique_variants_by_gene.pdf"),
       p2_7, width = 10, height = 8)

# -----------------------------------------------------------------------------
# Print summary
# -----------------------------------------------------------------------------
cat("\n=== SECTION 2 COMPLETE ===\n")
cat("Files created:\n")
cat("- ch/figures/section2_concordance/concordance_overview.pdf\n")
cat("- ch/figures/section2_concordance/concordance_by_sample.pdf\n")
cat("- ch/figures/section2_concordance/pairwise_concordance.pdf\n")
cat("- ch/figures/section2_concordance/unique_variants_by_gene.pdf\n")
cat("- ch/data/combined_variants_all.csv\n")
cat("- ch/data/concordance_summary.csv\n")
cat("- ch/data/pairwise_concordance.csv\n")



# =============================================================================
# SECTION 3: ALLELE FREQUENCY (VAF) ANALYSIS
# =============================================================================

library(tidyverse)
library(patchwork)

# Create output directory
dir.create(file.path("ch", "figures", "section3_vaf"),
           recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------------------------------
# 3.1: VAF distribution by caller
# -----------------------------------------------------------------------------

# You already have af_data from earlier, but let's recreate it with more info
af_data_full <- map_dfr(sample_ids, function(id) {
    map_dfr(names(vep[[id]]), function(caller) {
        create_variant_id(vep[[id]][[caller]]) %>%
            mutate(
                sample = id,
                caller = caller,
                AF = as.numeric(Sample.AltFrac),
                TotalDepth = as.numeric(Sample.Depth),
                AltDepth = as.numeric(Sample.AltDepth)
            ) %>%
            select(variant_id, sample, caller, AF, TotalDepth, AltDepth, Gene, Chr)
    })
})

# Add concordance info
af_data_full <- af_data_full %>%
    group_by(sample, variant_id) %>%
    mutate(
        n_callers = n_distinct(caller),
        concordance_level = case_when(
            n_callers == 3 ~ "High (all 3)",
            n_callers == 2 ~ "Medium (2)",
            n_callers == 1 ~ "Low (1)"
        )
    ) %>%
    ungroup()

# Summary statistics
vaf_summary <- af_data_full %>%
    group_by(caller) %>%
    summarise(
        n_variants = n(),
        mean_vaf = mean(AF, na.rm = TRUE),
        median_vaf = median(AF, na.rm = TRUE),
        sd_vaf = sd(AF, na.rm = TRUE),
        q25_vaf = quantile(AF, 0.25, na.rm = TRUE),
        q75_vaf = quantile(AF, 0.75, na.rm = TRUE),
        min_vaf = min(AF, na.rm = TRUE),
        max_vaf = max(AF, na.rm = TRUE)
    )

print("VAF summary by caller:")
print(vaf_summary)
# caller  n_variants mean_vaf median_vaf sd_vaf q25_vaf q75_vaf min_vaf max_vaf
# <chr>        <int>    <dbl>      <dbl>  <dbl>   <dbl>   <dbl>   <dbl>   <dbl>
#     1 mutect2       1815   0.108       0.048 0.156    0.031  0.0955    0          1
# 2 vardict       5987   0.0494      0.029 0.0862   0.021  0.043     0.01       1

write_csv(vaf_summary, file.path("ch", "data", "vaf_summary_by_caller.csv"))

# -----------------------------------------------------------------------------
# 3.2: VAF distribution plots
# -----------------------------------------------------------------------------

# Density plot
p3_1 <- ggplot(af_data_full, aes(x = AF, fill = caller)) +
    geom_density(alpha = 0.5) +
    geom_vline(xintercept = c(0.5, 1.0), linetype = "dashed", color = "red") +
    labs(title = "VAF Distribution by Caller",
         subtitle = "Dashed lines at 0.5 (heterozygous) and 1.0 (homozygous)",
         x = "Variant Allele Frequency",
         y = "Density") +
    theme_bw()

# Histogram
p3_2 <- ggplot(af_data_full, aes(x = AF, fill = caller)) +
    geom_histogram(bins = 50, alpha = 0.6, position = "identity") +
    facet_wrap(~caller, scales = "free_y") +
    geom_vline(xintercept = c(0.5, 1.0), linetype = "dashed", color = "red") +
    labs(title = "VAF Histogram by Caller",
         x = "Variant Allele Frequency",
         y = "Count") +
    theme_bw()

# Box plot
p3_3 <- ggplot(af_data_full, aes(x = caller, y = AF, fill = caller)) +
    geom_boxplot(alpha = 0.7, outlier.alpha = 0.3) +
    geom_hline(yintercept = c(0.5, 1.0), linetype = "dashed", color = "red") +
    labs(title = "VAF Distribution Comparison",
         x = "Caller",
         y = "Variant Allele Frequency") +
    theme_bw() +
    theme(legend.position = "none")

# Violin plot
p3_4 <- ggplot(af_data_full, aes(x = caller, y = AF, fill = caller)) +
    geom_violin(alpha = 0.7) +
    geom_boxplot(width = 0.2, alpha = 0.5, outlier.alpha = 0.3) +
    geom_hline(yintercept = c(0.5, 1.0), linetype = "dashed", color = "red") +
    labs(title = "VAF Distribution (Violin Plot)",
         x = "Caller",
         y = "Variant Allele Frequency") +
    theme_bw() +
    theme(legend.position = "none")

ggsave(file.path("ch", "figures", "section3_vaf", "vaf_distribution.pdf"),
       (p3_1 / p3_3), width = 10, height = 10)

ggsave(file.path("ch", "figures", "section3_vaf", "vaf_histograms.pdf"),
       p3_2, width = 12, height = 6)

ggsave(file.path("ch", "figures", "section3_vaf", "vaf_violin.pdf"),
       p3_4, width = 8, height = 6)

# -----------------------------------------------------------------------------
# 3.3: Low VAF variants analysis
# -----------------------------------------------------------------------------

# Define VAF bins
af_data_binned <- af_data_full %>%
    mutate(
        vaf_bin = case_when(
            AF < 0.05 ~ "<5%",
            AF < 0.10 ~ "5-10%",
            AF < 0.20 ~ "10-20%",
            AF < 0.30 ~ "20-30%",
            AF < 0.40 ~ "30-40%",
            AF < 0.50 ~ "40-50%",
            AF < 0.60 ~ "50-60%",
            AF < 0.70 ~ "60-70%",
            AF < 0.80 ~ "70-80%",
            AF < 0.90 ~ "80-90%",
            TRUE ~ "≥90%"
        ),
        vaf_bin = factor(vaf_bin, levels = c("<5%", "5-10%", "10-20%", "20-30%",
                                             "30-40%", "40-50%", "50-60%", "60-70%",
                                             "70-80%", "80-90%", "≥90%"))
    )

# Count variants in each bin
vaf_bin_counts <- af_data_binned %>%
    count(caller, vaf_bin) %>%
    group_by(caller) %>%
    mutate(percentage = n / sum(n) * 100)

print("VAF bin distribution:")
print(vaf_bin_counts)

p3_5 <- ggplot(vaf_bin_counts, aes(x = vaf_bin, y = n, fill = caller)) +
    geom_col(position = "dodge") +
    labs(title = "Variant Count by VAF Bin and Caller",
         x = "VAF Bin",
         y = "Number of Variants") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

p3_6 <- ggplot(vaf_bin_counts, aes(x = vaf_bin, y = percentage, fill = caller)) +
    geom_col(position = "dodge") +
    labs(title = "Percentage of Variants by VAF Bin",
         x = "VAF Bin",
         y = "Percentage of Variants (%)") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path("ch", "figures", "section3_vaf", "vaf_bins.pdf"),
       p3_5 / p3_6, width = 12, height = 10)

# Low VAF variants (<10%)
low_vaf_summary <- af_data_full %>%
    filter(AF < 0.10) %>%
    count(caller, name = "n_low_vaf") %>%
    left_join(
        af_data_full %>% count(caller, name = "total"),
        by = "caller"
    ) %>%
    mutate(percentage_low_vaf = n_low_vaf / total * 100)

print("Low VAF variants (<10%):")
print(low_vaf_summary)
# caller  n_low_vaf total percentage_low_vaf
# <chr>       <int> <int>              <dbl>
#     1 mutect2      1373  1815               75.6
# 2 vardict      5635  5987               94.1

# -----------------------------------------------------------------------------
# 3.4: VAF agreement for shared variants
# -----------------------------------------------------------------------------

# For variants called by multiple callers, compare their VAFs
shared_vaf <- af_data_full %>%
    filter(n_callers > 1) %>%
    select(sample, variant_id, caller, AF, n_callers) %>%
    pivot_wider(names_from = caller, values_from = AF)

# Mutect2 vs VarDict scatter
if("mutect2" %in% names(shared_vaf) && "vardict" %in% names(shared_vaf)) {

    p3_7 <- shared_vaf %>%
        filter(!is.na(mutect2) & !is.na(vardict)) %>%
        ggplot(aes(x = mutect2, y = vardict)) +
        geom_point(alpha = 0.5, size = 2) +
        geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
        geom_smooth(method = "lm", se = TRUE, color = "blue") +
        labs(title = "VAF Agreement: Mutect2 vs VarDict",
             subtitle = "For variants called by both callers",
             x = "Mutect2 VAF",
             y = "VarDict VAF") +
        theme_bw()

    # Calculate correlation
    vaf_cor <- shared_vaf %>%
        filter(!is.na(mutect2) & !is.na(vardict)) %>%
        summarise(
            correlation = cor(mutect2, vardict, use = "complete.obs"),
            n_shared = n(),
            mean_diff = mean(abs(mutect2 - vardict), na.rm = TRUE),
            median_diff = median(abs(mutect2 - vardict), na.rm = TRUE)
        )

    print("VAF correlation between Mutect2 and VarDict:")
    print(vaf_cor)

    # Difference plot (Bland-Altman style)
    p3_8 <- shared_vaf %>%
        filter(!is.na(mutect2) & !is.na(vardict)) %>%
        mutate(
            mean_vaf = (mutect2 + vardict) / 2,
            diff_vaf = mutect2 - vardict
        ) %>%
        ggplot(aes(x = mean_vaf, y = diff_vaf)) +
        geom_point(alpha = 0.5) +
        geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
        geom_smooth(method = "loess", se = TRUE) +
        labs(title = "VAF Difference: Mutect2 - VarDict",
             subtitle = "Bland-Altman style plot",
             x = "Mean VAF (average of both callers)",
             y = "Difference in VAF (Mutect2 - VarDict)") +
        theme_bw()

    ggsave(file.path("ch", "figures", "section3_vaf", "vaf_agreement.pdf"),
           p3_7 / p3_8, width = 10, height = 12)
}

# -----------------------------------------------------------------------------
# 3.5: VAF by concordance level
# -----------------------------------------------------------------------------

p3_9 <- ggplot(af_data_full, aes(x = concordance_level, y = AF, fill = concordance_level)) +
    geom_violin(alpha = 0.7) +
    geom_boxplot(width = 0.2, alpha = 0.5, outlier.alpha = 0.3) +
    facet_wrap(~caller) +
    labs(title = "VAF Distribution by Concordance Level",
         subtitle = "Do concordant variants have different VAFs?",
         x = "Concordance Level",
         y = "Variant Allele Frequency") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none")

ggsave(file.path("ch", "figures", "section3_vaf", "vaf_by_concordance.pdf"),
       p3_9, width = 12, height = 6)

# Statistical comparison
vaf_by_concordance_stats <- af_data_full %>%
    group_by(caller, concordance_level) %>%
    summarise(
        n = n(),
        mean_vaf = mean(AF, na.rm = TRUE),
        median_vaf = median(AF, na.rm = TRUE),
        sd_vaf = sd(AF, na.rm = TRUE),
        .groups = "drop"
    )

print("VAF by concordance level:")
print(vaf_by_concordance_stats)

write_csv(vaf_by_concordance_stats,
          file.path("ch", "data", "vaf_by_concordance.csv"))

# -----------------------------------------------------------------------------
# 3.6: VAF stratified by coverage depth
# -----------------------------------------------------------------------------

# Define depth bins
af_data_with_depth_bins <- af_data_full %>%
    mutate(
        depth_bin = case_when(
            TotalDepth < 30 ~ "<30x",
            TotalDepth < 50 ~ "30-50x",
            TotalDepth < 100 ~ "50-100x",
            TotalDepth < 200 ~ "100-200x",
            TRUE ~ "≥200x"
        ),
        depth_bin = factor(depth_bin, levels = c("<30x", "30-50x", "50-100x", "100-200x", "≥200x"))
    )

p3_10 <- ggplot(af_data_with_depth_bins, aes(x = depth_bin, y = AF, fill = caller)) +
    geom_boxplot(alpha = 0.7, outlier.alpha = 0.2) +
    facet_wrap(~caller) +
    labs(title = "VAF by Coverage Depth",
         subtitle = "Does VAF estimation change with depth?",
         x = "Coverage Depth Bin",
         y = "Variant Allele Frequency") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none")

ggsave(file.path("ch", "figures", "section3_vaf", "vaf_by_depth.pdf"),
       p3_10, width = 12, height = 6)

# =============================================================================
# SECTION 4: READ DEPTH ANALYSIS
# =============================================================================

library(tidyverse)
library(patchwork)

# Create output directory
dir.create(file.path("ch", "figures", "section4_depth"),
           recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------------------------------
# 4.1: Depth data with concordance info
# -----------------------------------------------------------------------------

# You already have depth_data from earlier, let's expand it
depth_data_full <- map_dfr(sample_ids, function(id) {
    map_dfr(names(vep[[id]]), function(caller) {
        create_variant_id(vep[[id]][[caller]]) %>%
            mutate(
                sample = id,
                caller = caller,
                TotalDepth = as.numeric(Sample.Depth),
                AltDepth = as.numeric(Sample.AltDepth),
                AF = as.numeric(Sample.AltFrac)
            ) %>%
            select(variant_id, sample, caller, TotalDepth, AltDepth, AF, Gene, Chr)
    })
})

# Add concordance info
depth_data_full <- depth_data_full %>%
    group_by(sample, variant_id) %>%
    mutate(
        n_callers = n_distinct(caller),
        concordance_level = case_when(
            n_callers == 3 ~ "High (all 3)",
            n_callers == 2 ~ "Medium (2)",
            n_callers == 1 ~ "Low (1)"
        )
    ) %>%
    ungroup()

# Summary statistics
depth_summary <- depth_data_full %>%
    group_by(caller) %>%
    summarise(
        n_variants = n(),
        mean_total_depth = mean(TotalDepth, na.rm = TRUE),
        median_total_depth = median(TotalDepth, na.rm = TRUE),
        sd_total_depth = sd(TotalDepth, na.rm = TRUE),
        mean_alt_depth = mean(AltDepth, na.rm = TRUE),
        median_alt_depth = median(AltDepth, na.rm = TRUE),
        min_total_depth = min(TotalDepth, na.rm = TRUE),
        max_total_depth = max(TotalDepth, na.rm = TRUE),
        pct_low_depth = sum(TotalDepth < 30, na.rm = TRUE) / n() * 100,
        pct_high_depth = sum(TotalDepth >= 100, na.rm = TRUE) / n() * 100
    )

print("Depth summary by caller:")
print(depth_summary)
# caller  n_variants mean_total_depth median_total_depth sd_total_depth mean_alt_depth
# <chr>        <int>            <dbl>              <dbl>          <dbl>          <dbl>
#     1 mutect2       1815             76.8                 64           68.1           6.07
# 2 vardict       5987             88.7                 78           54.2           3.71

write_csv(depth_summary, file.path("ch", "data", "depth_summary_by_caller.csv"))

# -----------------------------------------------------------------------------
# 4.2: Total depth distribution
# -----------------------------------------------------------------------------

# Density plot (log scale)
p4_1 <- ggplot(depth_data_full, aes(x = TotalDepth, fill = caller)) +
    geom_density(alpha = 0.5) +
    scale_x_log10(labels = scales::comma) +
    geom_vline(xintercept = c(30, 100), linetype = "dashed", color = "red") +
    labs(title = "Total Read Depth Distribution (log scale)",
         subtitle = "Dashed lines at 30x and 100x coverage",
         x = "Total Depth (log scale)",
         y = "Density") +
    theme_bw()

# Histogram
p4_2 <- ggplot(depth_data_full, aes(x = TotalDepth, fill = caller)) +
    geom_histogram(bins = 50, alpha = 0.6, position = "identity") +
    scale_x_log10(labels = scales::comma) +
    facet_wrap(~caller, scales = "free_y") +
    geom_vline(xintercept = c(30, 100), linetype = "dashed", color = "red") +
    labs(title = "Total Depth Histogram by Caller",
         x = "Total Depth (log scale)",
         y = "Count") +
    theme_bw()

# Box plot
p4_3 <- ggplot(depth_data_full, aes(x = caller, y = TotalDepth, fill = caller)) +
    geom_boxplot(alpha = 0.7, outlier.alpha = 0.2) +
    scale_y_log10(labels = scales::comma) +
    geom_hline(yintercept = c(30, 100), linetype = "dashed", color = "red") +
    labs(title = "Total Depth Distribution Comparison",
         x = "Caller",
         y = "Total Depth (log scale)") +
    theme_bw() +
    theme(legend.position = "none")

ggsave(file.path("ch", "figures", "section4_depth", "total_depth_distribution.pdf"),
       (p4_1 / p4_3), width = 10, height = 10)

ggsave(file.path("ch", "figures", "section4_depth", "total_depth_histograms.pdf"),
       p4_2, width = 12, height = 6)

# -----------------------------------------------------------------------------
# 4.3: Alt depth distribution
# -----------------------------------------------------------------------------

p4_4 <- ggplot(depth_data_full, aes(x = AltDepth, fill = caller)) +
    geom_density(alpha = 0.5) +
    scale_x_log10(labels = scales::comma) +
    geom_vline(xintercept = c(3, 5, 10), linetype = "dashed", color = "red") +
    labs(title = "Alternate Allele Depth Distribution (log scale)",
         subtitle = "Dashed lines at 3, 5, and 10 alt reads",
         x = "Alt Depth (log scale)",
         y = "Density") +
    theme_bw()

p4_5 <- ggplot(depth_data_full, aes(x = caller, y = AltDepth, fill = caller)) +
    geom_boxplot(alpha = 0.7, outlier.alpha = 0.2) +
    scale_y_log10(labels = scales::comma) +
    geom_hline(yintercept = c(3, 5, 10), linetype = "dashed", color = "red") +
    labs(title = "Alt Depth Distribution Comparison",
         x = "Caller",
         y = "Alt Depth (log scale)") +
    theme_bw() +
    theme(legend.position = "none")

ggsave(file.path("ch", "figures", "section4_depth", "alt_depth_distribution.pdf"),
       p4_4 / p4_5, width = 10, height = 10)

# Alt read support analysis
alt_support_summary <- depth_data_full %>%
    mutate(
        alt_support_category = case_when(
            AltDepth < 3 ~ "<3 reads",
            AltDepth < 5 ~ "3-4 reads",
            AltDepth < 10 ~ "5-9 reads",
            AltDepth < 20 ~ "10-19 reads",
            TRUE ~ "≥20 reads"
        ),
        alt_support_category = factor(alt_support_category,
                                      levels = c("<3 reads", "3-4 reads", "5-9 reads",
                                                 "10-19 reads", "≥20 reads"))
    ) %>%
    count(caller, alt_support_category) %>%
    group_by(caller) %>%
    mutate(percentage = n / sum(n) * 100)

print("Alt read support summary:")
print(alt_support_summary)

p4_6 <- ggplot(alt_support_summary, aes(x = alt_support_category, y = percentage, fill = caller)) +
    geom_col(position = "dodge") +
    labs(title = "Distribution of Alternate Read Support",
         x = "Alt Read Support",
         y = "Percentage of Variants (%)") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path("ch", "figures", "section4_depth", "alt_support_categories.pdf"),
       p4_6, width = 10, height = 6)

# -----------------------------------------------------------------------------
# 4.4: Depth vs VAF relationship
# -----------------------------------------------------------------------------

p4_7 <- ggplot(depth_data_full, aes(x = TotalDepth, y = AF, color = caller)) +
    geom_point(alpha = 0.3, size = 1) +
    scale_x_log10(labels = scales::comma) +
    geom_hline(yintercept = c(0.5, 1.0), linetype = "dashed", color = "red") +
    facet_wrap(~caller) +
    labs(title = "VAF vs Total Depth by Caller",
         subtitle = "Does VAF depend on coverage?",
         x = "Total Depth (log scale)",
         y = "Variant Allele Frequency") +
    theme_bw() +
    theme(legend.position = "none")

# Alt depth vs total depth
p4_8 <- ggplot(depth_data_full, aes(x = TotalDepth, y = AltDepth, color = caller)) +
    geom_point(alpha = 0.3, size = 1) +
    scale_x_log10(labels = scales::comma) +
    scale_y_log10(labels = scales::comma) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
    facet_wrap(~caller) +
    labs(title = "Alt Depth vs Total Depth",
         x = "Total Depth (log scale)",
         y = "Alt Depth (log scale)") +
    theme_bw() +
    theme(legend.position = "none")

ggsave(file.path("ch", "figures", "section4_depth", "depth_vs_vaf.pdf"),
       p4_7, width = 12, height = 6)

ggsave(file.path("ch", "figures", "section4_depth", "alt_vs_total_depth.pdf"),
       p4_8, width = 12, height = 6)

# -----------------------------------------------------------------------------
# 4.5: Depth by concordance level
# -----------------------------------------------------------------------------

p4_9 <- ggplot(depth_data_full, aes(x = concordance_level, y = TotalDepth, fill = concordance_level)) +
    geom_violin(alpha = 0.7) +
    geom_boxplot(width = 0.2, alpha = 0.5, outlier.alpha = 0.2) +
    scale_y_log10(labels = scales::comma) +
    facet_wrap(~caller) +
    labs(title = "Total Depth by Concordance Level",
         subtitle = "Do shared variants have higher coverage?",
         x = "Concordance Level",
         y = "Total Depth (log scale)") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none")

p4_10 <- ggplot(depth_data_full, aes(x = concordance_level, y = AltDepth, fill = concordance_level)) +
    geom_violin(alpha = 0.7) +
    geom_boxplot(width = 0.2, alpha = 0.5, outlier.alpha = 0.2) +
    scale_y_log10(labels = scales::comma) +
    facet_wrap(~caller) +
    labs(title = "Alt Depth by Concordance Level",
         subtitle = "Do shared variants have more alt reads?",
         x = "Concordance Level",
         y = "Alt Depth (log scale)") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none")

ggsave(file.path("ch", "figures", "section4_depth", "depth_by_concordance.pdf"),
       p4_9 / p4_10, width = 12, height = 12)

# Statistical summary
depth_by_concordance_stats <- depth_data_full %>%
    group_by(caller, concordance_level) %>%
    summarise(
        n = n(),
        mean_total_depth = mean(TotalDepth, na.rm = TRUE),
        median_total_depth = median(TotalDepth, na.rm = TRUE),
        mean_alt_depth = mean(AltDepth, na.rm = TRUE),
        median_alt_depth = median(AltDepth, na.rm = TRUE),
        .groups = "drop"
    )

print("Depth by concordance level:")
print(depth_by_concordance_stats)
# caller  concordance_level     n mean_total_depth median_total_depth mean_alt_depth
# <chr>   <chr>             <int>            <dbl>              <dbl>          <dbl>
#     1 mutect2 Low (1)            1090             78.4                 64           4.73
# 2 mutect2 Medium (2)          725             74.4                 65           8.09
# 3 vardict Low (1)            5262             90.4                 80           3.11
# 4 vardict Medium (2)          725             76.9                 62           8.08

write_csv(depth_by_concordance_stats,
          file.path("ch", "data", "depth_by_concordance.csv"))

# -----------------------------------------------------------------------------
# 4.6: Depth agreement for shared variants
# -----------------------------------------------------------------------------

# For variants called by multiple callers, compare their depths
shared_depth <- depth_data_full %>%
    filter(n_callers > 1) %>%
    select(sample, variant_id, caller, TotalDepth, AltDepth, n_callers) %>%
    pivot_wider(names_from = caller,
                values_from = c(TotalDepth, AltDepth),
                names_sep = "_")

# Mutect2 vs VarDict depth comparison
if("TotalDepth_mutect2" %in% names(shared_depth) && "TotalDepth_vardict" %in% names(shared_depth)) {

    p4_11 <- shared_depth %>%
        filter(!is.na(TotalDepth_mutect2) & !is.na(TotalDepth_vardict)) %>%
        ggplot(aes(x = TotalDepth_mutect2, y = TotalDepth_vardict)) +
        geom_point(alpha = 0.5, size = 2) +
        geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
        scale_x_log10(labels = scales::comma) +
        scale_y_log10(labels = scales::comma) +
        labs(title = "Total Depth Agreement: Mutect2 vs VarDict",
             subtitle = "For variants called by both",
             x = "Mutect2 Total Depth (log scale)",
             y = "VarDict Total Depth (log scale)") +
        theme_bw()

    p4_12 <- shared_depth %>%
        filter(!is.na(AltDepth_mutect2) & !is.na(AltDepth_vardict)) %>%
        ggplot(aes(x = AltDepth_mutect2, y = AltDepth_vardict)) +
        geom_point(alpha = 0.5, size = 2) +
        geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
        scale_x_log10(labels = scales::comma) +
        scale_y_log10(labels = scales::comma) +
        labs(title = "Alt Depth Agreement: Mutect2 vs VarDict",
             subtitle = "For variants called by both",
             x = "Mutect2 Alt Depth (log scale)",
             y = "VarDict Alt Depth (log scale)") +
        theme_bw()

    ggsave(file.path("ch", "figures", "section4_depth", "depth_agreement.pdf"),
           p4_11 / p4_12, width = 10, height = 12)

    # Calculate correlation
    depth_cor <- shared_depth %>%
        filter(!is.na(TotalDepth_mutect2) & !is.na(TotalDepth_vardict)) %>%
        summarise(
            total_depth_cor = cor(TotalDepth_mutect2, TotalDepth_vardict, use = "complete.obs"),
            alt_depth_cor = cor(AltDepth_mutect2, AltDepth_vardict, use = "complete.obs"),
            n_shared = n()
        )

    print("Depth correlation between Mutect2 and VarDict:")
    print(depth_cor)
}



# =============================================================================
# SECTION 6: CH GENE-SPECIFIC ANALYSIS
# =============================================================================

library(tidyverse)
library(patchwork)
library(readxl)

# Create output directory
dir.create(file.path("ch", "figures", "section6_ch_genes"),
           recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------------------------------
# 6.1: Load CH gene list (you already have this)
# -----------------------------------------------------------------------------

# CH genes already loaded as ch_genes
print(paste("Number of CH genes:", length(ch_genes)))
print("CH genes:")
print(ch_genes)

# -----------------------------------------------------------------------------
# 6.2: Identify CH gene variants in each caller
# -----------------------------------------------------------------------------

# Extract all variants and flag if they're in CH genes
all_variants_ch_flagged <- map_dfr(sample_ids, function(id) {
    map_dfr(names(vep[[id]]), function(caller) {
        vep[[id]][[caller]] %>%
            mutate(
                sample = id,
                caller = caller,
                is_ch_gene = Gene %in% ch_genes,
                TotalDepth = as.numeric(Sample.Depth),
                AltDepth = as.numeric(Sample.AltDepth),
                AF = as.numeric(Sample.AltFrac)
            ) %>%
            select(sample, caller, Gene, is_ch_gene, Chr, Start, REF, ALT,
                   AF, TotalDepth, AltDepth, Variant.Consequence)
    })
})

# Add variant ID and concordance
all_variants_ch_flagged <- all_variants_ch_flagged %>%
    mutate(variant_id = paste(Chr, Start, REF, ALT, sep = ":")) %>%
    group_by(sample, variant_id) %>%
    mutate(
        n_callers = n_distinct(caller),
        concordance_level = case_when(
            n_callers == 3 ~ "High (all 3)",
            n_callers == 2 ~ "Medium (2)",
            n_callers == 1 ~ "Low (1)"
        )
    ) %>%
    ungroup()

# Summary: CH gene detection
ch_gene_detection <- all_variants_ch_flagged %>%
    filter(is_ch_gene) %>%
    count(caller, name = "n_ch_variants") %>%
    left_join(
        all_variants_ch_flagged %>% count(caller, name = "total_variants"),
        by = "caller"
    ) %>%
    mutate(
        pct_ch_variants = n_ch_variants / total_variants * 100
    )

print("CH gene variant detection by caller:")
print(ch_gene_detection)
# caller  n_ch_variants total_variants pct_ch_variants
# <chr>           <int>          <int>           <dbl>
#     1 mutect2            20           1815           1.10
# 2 vardict            46           5987           0.768

write_csv(ch_gene_detection,
          file.path("ch", "data", "ch_gene_detection_summary.csv"))

# -----------------------------------------------------------------------------
# 6.3: Which CH genes are detected?
# -----------------------------------------------------------------------------

ch_genes_found <- all_variants_ch_flagged %>%
    filter(is_ch_gene) %>%
    distinct(Gene, caller) %>%
    count(Gene, name = "n_callers_found") %>%
    arrange(desc(n_callers_found), Gene)

print("CH genes found and by how many callers:")
print(ch_genes_found)

# Which CH genes were NOT found by any caller?
ch_genes_not_found <- setdiff(ch_genes, ch_genes_found$Gene)
print("CH genes NOT detected by any caller:")
print(ch_genes_not_found)

# Create a presence/absence matrix
ch_gene_presence_matrix <- all_variants_ch_flagged %>%
    filter(is_ch_gene) %>%
    distinct(Gene, caller) %>%
    mutate(present = TRUE) %>%
    pivot_wider(names_from = caller, values_from = present, values_fill = FALSE) %>%
    arrange(Gene)

print("CH gene presence by caller:")
print(ch_gene_presence_matrix)

write_csv(ch_gene_presence_matrix,
          file.path("ch", "data", "ch_gene_presence_matrix.csv"))

# -----------------------------------------------------------------------------
# 6.4: Visualizations - CH gene detection
# -----------------------------------------------------------------------------

# Bar plot: CH variants per caller
p6_1 <- ggplot(ch_gene_detection, aes(x = caller, y = n_ch_variants, fill = caller)) +
    geom_col() +
    geom_text(aes(label = n_ch_variants), vjust = -0.5, size = 5) +
    labs(title = "CH Gene Variants Detected by Caller",
         x = "Caller",
         y = "Number of CH Gene Variants") +
    theme_bw() +
    theme(legend.position = "none")

# Percentage
p6_2 <- ggplot(ch_gene_detection, aes(x = caller, y = pct_ch_variants, fill = caller)) +
    geom_col() +
    geom_text(aes(label = paste0(round(pct_ch_variants, 1), "%")),
              vjust = -0.5, size = 5) +
    labs(title = "Percentage of Variants in CH Genes",
         x = "Caller",
         y = "Percentage (%)") +
    theme_bw() +
    theme(legend.position = "none")


ggsave(file.path("ch", "figures", "section6_ch_genes", "ch_detection_overview.pdf"),
       p6_1 / p6_2, width = 10, height = 10)

# Heatmap of CH gene presence
ch_gene_heatmap_data <- ch_gene_presence_matrix %>%
    pivot_longer(-Gene, names_to = "caller", values_to = "detected") %>%
    mutate(detected_label = ifelse(detected, "✓", ""))

p6_3 <- ggplot(ch_gene_heatmap_data, aes(x = caller, y = Gene, fill = detected)) +
    geom_tile(color = "white", size = 0.5) +
    geom_text(aes(label = detected_label), size = 5) +
    scale_fill_manual(values = c("FALSE" = "gray90", "TRUE" = "steelblue")) +
    labs(title = "CH Gene Detection Matrix",
         subtitle = "Which CH genes are detected by which callers?",
         x = "Caller",
         y = "CH Gene") +
    theme_minimal() +
    theme(axis.text.y = element_text(size = 10),
          legend.position = "none")

ggsave(file.path("ch", "figures", "section6_ch_genes", "ch_gene_heatmap.pdf"),
       p6_3, width = 8, height = max(6, length(unique(ch_gene_heatmap_data$Gene)) * 0.3))

# -----------------------------------------------------------------------------
# 6.5: CH gene variant characteristics
# -----------------------------------------------------------------------------

# Compare characteristics of CH gene variants vs non-CH variants
variant_characteristics <- all_variants_ch_flagged %>%
    mutate(gene_type = ifelse(is_ch_gene, "CH Gene", "Non-CH Gene")) %>%
    group_by(caller, gene_type) %>%
    summarise(
        n_variants = n(),
        mean_vaf = mean(AF, na.rm = TRUE),
        median_vaf = median(AF, na.rm = TRUE),
        mean_depth = mean(TotalDepth, na.rm = TRUE),
        median_depth = median(TotalDepth, na.rm = TRUE),
        mean_alt_depth = mean(AltDepth, na.rm = TRUE),
        .groups = "drop"
    )

print("CH gene vs non-CH gene variant characteristics:")
print(variant_characteristics)
# caller  gene_type   n_variants mean_vaf median_vaf mean_depth median_depth mean_alt_depth
# <chr>   <chr>            <int>    <dbl>      <dbl>      <dbl>        <dbl>          <dbl>
#     1 mutect2 CH Gene             20   0.163      0.0985       47.0           44           7.45
# 2 mutect2 Non-CH Gene       1795   0.107      0.048        77.1           65           6.06
# 3 vardict CH Gene             46   0.0795     0.029        74.1           75           4.02
# 4 vardict Non-CH Gene       5941   0.0491     0.029        88.9           78           3.71

write_csv(variant_characteristics,
          file.path("ch", "data", "ch_vs_nonch_characteristics.csv"))

# VAF comparison
p6_4 <- all_variants_ch_flagged %>%
    mutate(gene_type = ifelse(is_ch_gene, "CH Gene", "Non-CH Gene")) %>%
    ggplot(aes(x = gene_type, y = AF, fill = gene_type)) +
    geom_violin(alpha = 0.7) +
    geom_boxplot(width = 0.2, alpha = 0.5, outlier.alpha = 0.2) +
    facet_wrap(~caller) +
    labs(title = "VAF Distribution: CH Genes vs Non-CH Genes",
         x = "Gene Type",
         y = "Variant Allele Frequency") +
    theme_bw() +
    theme(legend.position = "none")

# Depth comparison
p6_5 <- all_variants_ch_flagged %>%
    mutate(gene_type = ifelse(is_ch_gene, "CH Gene", "Non-CH Gene")) %>%
    ggplot(aes(x = gene_type, y = TotalDepth, fill = gene_type)) +
    geom_violin(alpha = 0.7) +
    geom_boxplot(width = 0.2, alpha = 0.5, outlier.alpha = 0.2) +
    scale_y_log10(labels = scales::comma) +
    facet_wrap(~caller) +
    labs(title = "Depth Distribution: CH Genes vs Non-CH Genes",
         x = "Gene Type",
         y = "Total Depth (log scale)") +
    theme_bw() +
    theme(legend.position = "none")

ggsave(file.path("ch", "figures", "section6_ch_genes", "ch_variant_characteristics.pdf"),
       p6_4 / p6_5, width = 12, height = 10)

# -----------------------------------------------------------------------------
# 6.6: CH gene concordance
# -----------------------------------------------------------------------------

ch_concordance <- all_variants_ch_flagged %>%
    filter(is_ch_gene) %>%
    distinct(sample, variant_id, Gene, .keep_all = TRUE) %>%
    count(Gene, n_callers) %>%
    pivot_wider(names_from = n_callers, values_from = n, values_fill = 0,
                names_prefix = "callers_") %>%
    mutate(
        total_variants = rowSums(select(., starts_with("callers_"))),
        pct_concordant = if_else("callers_3" %in% names(.), callers_3 / total_variants * 100, 0)
    ) %>%
    arrange(desc(total_variants))

print("CH gene concordance:")
print(ch_concordance)

write_csv(ch_concordance,
          file.path("ch", "data", "ch_gene_concordance.csv"))

# Plot
p6_6 <- ch_concordance %>%
    pivot_longer(starts_with("callers_"), names_to = "concordance", values_to = "count") %>%
    mutate(concordance = str_replace(concordance, "callers_", "")) %>%
    ggplot(aes(x = reorder(Gene, -total_variants), y = count, fill = concordance)) +
    geom_col(position = "stack") +
    labs(title = "CH Gene Variant Concordance",
         subtitle = "How many variants per gene are detected by 1, 2, or 3 callers?",
         x = "CH Gene",
         y = "Number of Variants",
         fill = "Number of\nCallers") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path("ch", "figures", "section6_ch_genes", "ch_concordance.pdf"),
       p6_6, width = 12, height = 6)

# -----------------------------------------------------------------------------
# 6.7: Detailed CH gene variant table
# -----------------------------------------------------------------------------

ch_variant_details <- all_variants_ch_flagged %>%
    filter(is_ch_gene) %>%
    select(sample, Gene, Chr, Start, REF, ALT, AF, TotalDepth, AltDepth,
           caller, n_callers, concordance_level, Variant.Consequence) %>%
    arrange(Gene, sample, Chr, Start)

write_csv(ch_variant_details,
          file.path("ch", "data", "ch_gene_variants_detailed.csv"))

print(paste("Total CH gene variants:", nrow(ch_variant_details)))
print(paste("Unique CH gene variant positions:",
            n_distinct(ch_variant_details$sample, ch_variant_details$Chr,
                       ch_variant_details$Start, ch_variant_details$REF,
                       ch_variant_details$ALT)))

# -----------------------------------------------------------------------------
# 6.8: Top CH genes by variant count
# -----------------------------------------------------------------------------

top_ch_genes <- all_variants_ch_flagged %>%
    filter(is_ch_gene) %>%
    count(Gene, caller) %>%
    group_by(Gene) %>%
    summarise(total_variants = sum(n)) %>%
    arrange(desc(total_variants)) %>%
    head(15)

p6_7 <- all_variants_ch_flagged %>%
    filter(is_ch_gene, Gene %in% top_ch_genes$Gene) %>%
    count(Gene, caller) %>%
    ggplot(aes(x = reorder(Gene, -n), y = n, fill = caller)) +
    geom_col(position = "dodge") +
    labs(title = "Top 15 CH Genes by Variant Count",
         x = "CH Gene",
         y = "Number of Variants") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path("ch", "figures", "section6_ch_genes", "top_ch_genes.pdf"),
       p6_7, width = 12, height = 6)




