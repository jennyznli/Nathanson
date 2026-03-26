library(here)
setwd("/Users/jennyzli/Documents/Nathanson")
source(here("R", "config.R"))
library(data.table, quietly=T)

library(tidyverse)
library(ggVennDiagram)
library(patchwork)
library(tidyverse)
library(fs)
library(UpSetR)
library(xfun)

# ========================
# LOAD DATA
# ========================
ss <- read_excel(file.path("ch", "ss", "20251110_aml_ngs.xlsx"))
up <- read.csv(here("simplexo", "data", "simplexo_up_map.csv"))
flags <- fread(here("PMBB", "3.0", "rgcname_pmbbid_metadata_flags.csv"))

up_pattern <- "^UPENN-PMBB_UP[0-9]+_UP[0-9]+$"
key <- flags %>% select(PMBB_ID, RGC_sample_name)

# ========================
# EXTRACT ALL FILES FROM FOLDERS
# ========================
base_dir <- "/Users/jennyzli/Documents/Nathanson/ch/data/all"
output_dir <- "/Users/jennyzli/Documents/Nathanson/ch/data/all2"
dir_create(output_dir)

# Find all CSV files recursively
csv_files <- dir_ls(base_dir, recurse = TRUE, glob = "*.csv")

# Process each file
results <- tibble(
    original_path = csv_files,
    filename = path_file(csv_files)
) %>%
    mutate(
        # Extract the parts from filename pattern: ID_caller_caller.norm.vep.report.csv
        # Split by underscore and extract caller (appears twice)
        parts = str_split(filename, "_"),
        # Get caller from the duplicated part (second to last before the dot)
        caller = map_chr(parts, ~{
            # Remove .norm.vep.report.csv and split what remains
            base_name = str_remove(.x[length(.x)], "\\.norm\\.vep\\.report\\.csv$")
            # Split by underscore to get caller (it appears as "caller_caller")
            caller_parts = str_split(base_name, "_")[[1]]
            caller_parts[1]  # Take first occurrence
        }),
        # Extract original ID (everything before the last occurrence of _caller)
        original_id = str_remove(filename, paste0("_", caller, "_", caller, "\\.norm\\.vep\\.report\\.csv$")),
        # Map to new ID
        new_id = key$PMBB_ID[match(original_id, key$RGC_sample_name)]
    )

# Copy files with corrected names
for(i in 1:nrow(results)) {
    original_file <- results$original_path[i]
    caller <- results$caller[i]
    original_id <- results$original_id[i]
    new_id <- results$new_id[i]

    # New filename with caller appearing only once
    new_filename <- glue::glue("{new_id}.{caller}.norm.vep.report.csv")
    print(new_filename)

    file_copy(original_file,
              path(output_dir, new_filename),
              overwrite = TRUE)
    cat("Copied:", original_id, "from", caller, "\n")
}

# ========================
# COMPARISON
# ========================
# Function to read and convert numeric columns
read_and_convert <- function(filepath) {
    df <- read_csv(filepath, show_col_types = FALSE)

    # Define columns that should be numeric
    numeric_cols <- c(
        "Start",
        "SpliceAI.DS_AG", "SpliceAI.DS_AL", "SpliceAI.DS_DG", "SpliceAI.DS_DL",
        "REVEL",
        "gnomAD.AF", "gnomAD.AFR", "gnomAD.AMR", "gnomAD.ASJ",
        "gnomAD.EAS", "gnomAD.FIN", "gnomAD.NFE", "gnomAD.OTH", "gnomAD.SAS",
        "gnomAD.MAX_AF",
        "MaveDB.score",
        "Sample.Depth", "Sample.AltDepth", "Sample.AltFrac"
    )

    # Define columns that should be character (even if they look numeric)
    char_cols <- c(
        "Sample.ID", "Chr", "REF", "ALT", "FILTER", "ID", "Gene",
        "Gene.Accession", "Variant.LoF_level", "Variant.Category",
        "Variant.Class", "Variant.Consequence", "HGVSc", "HGVSp",
        "Feature.Type", "Feature.Accession", "Bio.type", "Existing.variation",
        "EXON", "INTRON", "STRAND", "cDNA.position", "CDS.position",
        "Protein.position", "Amino.acids", "Codons",
        "gnomAD.MAX_POPS", "ClinVar", "ClinVar.SIG", "ClinVar.REVSTAT",
        "ClinVar.DN", "AutoGVP", "AM.class", "AM.pathogenicity",
        "MaveDB.nt", "MaveDB.pro", "MaveDB.urn", "Sample.Zyg"
    )

    # Convert numeric columns
    df <- df %>%
        mutate(suppressWarnings(across(
            all_of(numeric_cols[numeric_cols %in% names(.)]),
            ~as.numeric(.)
        ))) %>%
        # Ensure character columns stay character
        mutate(across(
            all_of(char_cols[char_cols %in% names(.)]),
            ~as.character(.)
        ))

    return(df)
}

# Set your paths
base_dir <- "/Users/jennyzli/Documents/Nathanson/ch/data/all2"  # or wherever your renamed files are

# Create file_info tibble
file_info <- tibble(
    path = dir_ls(base_dir, glob = "*.csv")
) %>%
    mutate(
        filename = path_file(path),
        # Extract ID and caller from filename pattern: ID.caller.norm.vep.report.csv
        id = str_extract(filename, "^[^\\.]+"),
        caller = str_extract(filename, "(?<=\\.)[^\\.]+(?=\\.norm\\.vep\\.report\\.csv)")
    )

# Check it looks right
print(file_info)

# Then continue with your read_and_convert function and vep creation...

vep <- file_info %>%
    mutate(data = map(path, read_and_convert)) %>%
    select(id, caller, data) %>%
    group_by(id) %>%
    summarise(
        callers = list(set_names(data, caller))
    ) %>%
    deframe()

sample_ids <- names(vep)

# -----------------
# 1. BASIC STATS: Variants per caller
# -----------------
variant_counts <- map_dfr(sample_ids, function(id) {
    map_dfr(names(vep[[id]]), function(caller) {
        tibble(
            sample = id,
            caller = caller,
            n_variants = nrow(vep[[id]][[caller]])
        )
    })
})

dir_create(file.path("ch", "figures"))

p1_individual <- ggplot(variant_counts, aes(x = caller, y = n_variants, fill = caller)) +
    geom_col() +
    facet_wrap(~sample, scales = "free_y") +
    labs(title = "Variants Called per Method (Individual Samples)",
         y = "Number of Variants",
         x = "Caller") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

p1_combined <- variant_counts %>%
    group_by(caller) %>%
    summarise(
        total_variants = sum(n_variants),
        mean_variants = mean(n_variants),
        sd_variants = sd(n_variants)
    ) %>%
    ggplot(aes(x = caller, y = total_variants, fill = caller)) +
    geom_col() +
    geom_text(aes(label = total_variants), vjust = -0.5) +
    labs(title = "Total Variants Called per Method (All Samples)",
         y = "Total Number of Variants",
         x = "Caller") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path("ch", "figures", "variant_counts_individual_all.pdf"), p1_individual, width = 10, height = 6)
ggsave(file.path("ch", "figures", "variant_counts_combined_all.pdf"), p1_combined, width = 8, height = 6)

# -----------------
# 2. VARIANT OVERLAP (Venn diagrams)
# -----------------
create_variant_id <- function(df) {
    df %>%
        mutate(variant_id = paste(Chr, Start, REF, ALT, sep = "_"))
}

make_venn <- function(sample_id) {
    callers <- names(vep[[sample_id]])

    variant_lists <- map(callers, function(caller) {
        create_variant_id(vep[[sample_id]][[caller]]) %>%
            pull(variant_id)
    })
    names(variant_lists) <- callers

    ggVennDiagram(variant_lists,
                  label = "count",
                  label_alpha = 0,
                  label_size = 5,
                  label_color = "black",
                  edge_size = 1) +
        scale_fill_distiller(palette = "YlOrRd", direction = 1) +
        labs(title = paste("Sample:", sample_id)) +
        theme(legend.position = "bottom")
}

venn_plots <- map(sample_ids, make_venn)
p2_individual <- wrap_plots(venn_plots, ncol = 3)

all_variants_by_caller <- map(names(vep[[1]]), function(caller) {
    # Get all variants from this caller across all samples
    map_dfr(sample_ids, function(id) {
        if (caller %in% names(vep[[id]])) {
            create_variant_id(vep[[id]][[caller]])
        }
    }) %>%
        distinct(variant_id) %>%
        pull(variant_id)
})
names(all_variants_by_caller) <- names(vep[[1]])

p2_combined <- ggVennDiagram(all_variants_by_caller,
                             label = "count",
                             label_alpha = 0,
                             label_size = 6,
                             label_color = "black",
                             edge_size = 1.5) +
    scale_fill_distiller(palette = "YlOrRd", direction = 1) +
    labs(title = "Variant Overlap Across All Samples Combined") +
    theme(legend.position = "bottom")

ggsave(file.path("ch", "figures", "venn_diagrams_individual_all.pdf"), p2_individual, width = 15, height = 15)
ggsave(file.path("ch", "figures", "venn_diagram_combined_all.pdf"), p2_combined, width = 8, height = 6)

# -----------------
# 2B. UPSET PLOTS FOR VARIANT OVERLAP
# -----------------
# Combined across all samples
all_variants_by_caller_list <- map(names(vep[[1]]), function(caller) {
    map_dfr(sample_ids, function(id) {
        if (caller %in% names(vep[[id]])) {
            create_variant_id(vep[[id]][[caller]])
        }
    }) %>%
        distinct(variant_id) %>%
        pull(variant_id)
})
names(all_variants_by_caller_list) <- names(vep[[1]])

# Create presence/absence matrix
all_unique_variants <- unique(unlist(all_variants_by_caller_list))
upset_data_combined <- data.frame(variant_id = all_unique_variants)

for(caller in names(all_variants_by_caller_list)) {
    upset_data_combined[[caller]] <- as.numeric(
        upset_data_combined$variant_id %in% all_variants_by_caller_list[[caller]]
    )
}

upset_data_combined <- upset_data_combined %>% select(-variant_id)

pdf(file.path("ch", "figures", "upset_combined_all.pdf"), width = 10, height = 6)
upset(upset_data_combined,
      sets = colnames(upset_data_combined),
      order.by = "freq",
      main.bar.color = "steelblue",
      sets.bar.color = "darkblue",
      text.scale = 1.5,
      point.size = 4,
      line.size = 1.5,
      mainbar.y.label = "Number of Unique Variants",
      sets.x.label = "Variants per Caller")
dev.off()

# -----------------
# 3. VARIANTS PER GENE (histogram)
# -----------------
gene_counts <- map_dfr(sample_ids, function(id) {
    map_dfr(names(vep[[id]]), function(caller) {
        vep[[id]][[caller]] %>%
            count(Gene, name = "n_variants") %>%
            mutate(sample = id, caller = caller)
    })
})

top_genes <- map_dfr(sample_ids, function(id) {
    map_dfr(names(vep[[id]]), function(caller) {
        create_variant_id(vep[[id]][[caller]]) %>%
            select(variant_id, Gene)
    })
}) %>%
    distinct(variant_id, Gene) %>%  # Get unique variants per gene
    count(Gene, name = "total_variants") %>%
    arrange(desc(total_variants)) %>%
    head(20)

p4 <- ggplot(top_genes, aes(x = reorder(Gene, total_variants), y = total_variants)) +
    geom_col(fill = "steelblue") +
    geom_text(aes(label = total_variants), hjust = -0.2) +
    coord_flip() +
    labs(title = "Top 20 Genes by Unique Variant Count (All Samples Combined)",
         x = "Gene", y = "Number of Unique Variants") +
    theme_bw()

ggsave(file.path("ch", "figures", "top_genes_combined_all.pdf"), p4, width = 8, height = 10)

### Count unique GENES per caller
genes_per_caller <- map_dfr(sample_ids, function(id) {
    map_dfr(names(vep[[id]]), function(caller) {
        vep[[id]][[caller]] %>%
            select(Gene) %>%
            mutate(caller = caller)  # Add caller column
    })
}) %>%
    distinct(Gene, caller) %>%  # Unique genes per caller
    count(caller, name = "n_genes")

print(genes_per_caller)

# If you want to see WHICH genes are unique to each caller:
genes_by_caller <- map_dfr(sample_ids, function(id) {
    map_dfr(names(vep[[id]]), function(caller) {
        vep[[id]][[caller]] %>%
            select(Gene) %>%
            mutate(caller = caller)
    })
}) %>%
    distinct(Gene, caller)

# See which genes each caller detects
genes_mutect2 <- genes_by_caller %>% filter(caller == "mutect2") %>% pull(Gene) %>% unique()
genes_vardict <- genes_by_caller %>% filter(caller == "vardict") %>% pull(Gene) %>% unique()
genes_pmbb <- genes_by_caller %>% filter(caller == "pmbb") %>% pull(Gene) %>% unique()

# Genes unique to each caller
unique_to_mutect2 <- setdiff(genes_mutect2, c(genes_vardict, genes_pmbb))
unique_to_vardict <- setdiff(genes_vardict, c(genes_mutect2, genes_pmbb))
unique_to_pmbb <- setdiff(genes_pmbb, c(genes_mutect2, genes_vardict))

# Genes detected by all three
genes_all_three <- Reduce(intersect, list(genes_mutect2, genes_vardict, genes_pmbb))

cat("\nGene detection by caller:\n")
cat("Mutect2:", length(genes_mutect2), "genes\n")
cat("VarDict:", length(genes_vardict), "genes\n")
cat("PMBB:", length(genes_pmbb), "genes\n")

cat("\nUnique to Mutect2:", length(unique_to_mutect2), "genes\n")
cat("Unique to VarDict:", length(unique_to_vardict), "genes\n")
cat("Unique to PMBB:", length(unique_to_pmbb), "genes\n")
cat("Detected by all three:", length(genes_all_three), "genes\n")

# -----------------
# 4. ALLELE FREQUENCY COMPARISON
# -----------------
af_data <- map_dfr(sample_ids, function(id) {
    map_dfr(names(vep[[id]]), function(caller) {
        create_variant_id(vep[[id]][[caller]]) %>%
            mutate(
                sample = id,
                caller = caller,
                AF = as.numeric(Sample.AltFrac)
            ) %>%
            select(variant_id, sample, caller, AF)
    })
})

# Compare AF across callers for shared variants
af_comparison <- af_data %>%
    group_by(sample, variant_id, caller) %>%
    summarise(AF = mean(AF, na.rm = TRUE), .groups = "drop") %>%
    group_by(sample, variant_id) %>%
    filter(n() > 1) %>%  # Only variants called by multiple callers
    ungroup() %>%
    pivot_wider(names_from = caller, values_from = AF)

# Combined scatter plot
af_comparison_combined <- af_data %>%
    group_by(variant_id, caller) %>%
    summarise(AF = mean(AF, na.rm = TRUE), .groups = "drop") %>%
    group_by(variant_id) %>%
    filter(n() > 1) %>%
    ungroup() %>%
    pivot_wider(names_from = caller, values_from = AF)

p5_combined <- ggplot(af_comparison_combined, aes(x = mutect2, y = vardict)) +
    geom_point(alpha = 0.5, size = 2) +
    geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed", size = 1) +
    labs(title = "Allele Frequency: Mutect2 vs VarDict (All Samples)",
         x = "Mutect2 AF", y = "VarDict AF") +
    theme_bw()
ggsave(file.path("ch", "figures", "af_scatter_mutect2_vardict_combined_all.pdf"), p5_combined, width = 6, height = 4)

p5_combined <- ggplot(af_comparison_combined, aes(x = mutect2, y = pmbb)) +
    geom_point(alpha = 0.5, size = 2) +
    geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed", size = 1) +
    labs(title = "Allele Frequency: Mutect2 vs PMBB (All Samples)",
         x = "Mutect2 AF", y = "PMBB AF") +
    theme_bw()
ggsave(file.path("ch", "figures", "af_scatter_mutect2_pmbb_combined_all.pdf"), p5_combined, width = 6, height = 4)

p5_combined <- ggplot(af_comparison_combined, aes(x = vardict, y = pmbb)) +
    geom_point(alpha = 0.5, size = 2) +
    geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed", size = 1) +
    labs(title = "Allele Frequency: Vardict vs PMBB (All Samples)",
         x = "Vardict AF", y = "PMBB AF") +
    theme_bw()
ggsave(file.path("ch", "figures", "af_scatter_vardict_pmbb_combined_all.pdf"), p5_combined, width = 6, height = 4)

### EXAMINE OUTLIERS
# ========================
# IDENTIFY DISCORDANT VARIANTS
# ========================

# Get full data with annotations
af_data_full <- map_dfr(sample_ids, function(id) {
    map_dfr(names(vep[[id]]), function(caller) {
        create_variant_id(vep[[id]][[caller]]) %>%
            mutate(sample = id, caller = caller, AF = as.numeric(Sample.AltFrac))
    })
})

# Function to find discordant variants
find_discordant <- function(comparison_df, caller1, caller2, threshold = 0.1) {
    comparison_df %>%
        filter(!is.na(.data[[caller1]]) & !is.na(.data[[caller2]])) %>%
        mutate(AF_diff = abs(.data[[caller1]] - .data[[caller2]])) %>%
        filter(AF_diff > threshold) %>%
        arrange(desc(AF_diff)) %>%
        # Add full annotations
        left_join(
            af_data_full %>%
                filter(caller == caller1) %>%
                select(variant_id, sample, Gene, Variant.Consequence, HGVSp,
                       Sample.Depth, Sample.AltDepth, gnomAD.AF, ClinVar.SIG),
            by = "variant_id"
        ) %>%
        select(variant_id, !!caller1 := .data[[caller1]], !!caller2 := .data[[caller2]],
               AF_diff, Gene, Variant.Consequence, HGVSp, Sample.Depth,
               Sample.AltDepth, gnomAD.AF, ClinVar.SIG)
}

# Get discordant variants (AF difference > 0.1)
discordant_mutect2_vardict <- find_discordant(af_comparison_combined, "mutect2", "vardict")
discordant_mutect2_pmbb <- find_discordant(af_comparison_combined, "mutect2", "pmbb")
discordant_vardict_pmbb <- find_discordant(af_comparison_combined, "vardict", "pmbb")


# Summary
cat("Discordant variants (AF diff > 0.1):\n")
cat("  Mutect2 vs VarDict:", nrow(discordant_mutect2_vardict), "\n")
cat("  Mutect2 vs PMBB:", nrow(discordant_mutect2_pmbb), "\n")
cat("  VarDict vs PMBB:", nrow(discordant_vardict_pmbb), "\n\n")

# View top 10 of each
View(head(discordant_mutect2_vardict, 10))
View(head(discordant_mutect2_pmbb, 10))
View(head(discordant_vardict_pmbb, 10))

# Save to CSV
write.csv(discordant_mutect2_vardict,
          file.path("ch", "output", "discordant_mutect2_vardict_all.csv"),
          row.names = FALSE)
write.csv(discordant_mutect2_pmbb,
          file.path("ch", "output", "discordant_mutect2_pmbb_all.csv"),
          row.names = FALSE)
write.csv(discordant_vardict_pmbb,
          file.path("ch", "output", "discordant_vardict_pmbb_all.csv"),
          row.names = FALSE)

### DENSITY PLOTS ###
# Density plot of AF distribution - per sample
p6_individual <- ggplot(af_data, aes(x = AF, fill = caller)) +
    geom_density(alpha = 0.5) +
    facet_wrap(~sample) +
    labs(title = "Allele Frequency Distribution by Caller (Individual)") +
    theme_bw()
ggsave(file.path("ch", "figures", "af_density_individual_all.pdf"), p6_individual, width = 8, height = 6)

# Combined density plot
p6_combined <- ggplot(af_data, aes(x = AF, fill = caller)) +
    geom_density(alpha = 0.5) +
    labs(title = "Allele Frequency Distribution by Caller (All Samples)",
         x = "Allele Frequency",
         y = "Density") +
    theme_bw() +
    theme(legend.position = "bottom")
ggsave(file.path("ch", "figures", "af_density_combined_all.pdf"), p6_combined, width = 8, height = 6)

# -----------------
# 5. READ DEPTH COMPARISON
# -----------------
depth_data <- map_dfr(sample_ids, function(id) {
    map_dfr(names(vep[[id]]), function(caller) {
        create_variant_id(vep[[id]][[caller]]) %>%
            mutate(
                sample = id,
                caller = caller,
                TotalDepth = as.numeric(Sample.Depth),
                AltDepth = as.numeric(Sample.AltDepth)
            ) %>%
            select(variant_id, sample, caller, TotalDepth, AltDepth)
    })
})

# Per sample depth distribution
p7_individual <- ggplot(depth_data, aes(x = TotalDepth, fill = caller)) +
    geom_histogram(bins = 50, alpha = 0.6, position = "identity") +
    scale_x_log10() +
    facet_grid(sample ~ caller) +
    labs(title = "Read Depth Distribution by Sample and Caller",
         x = "log10(Total Depth)", y = "Count") +
    theme_bw()

# Combined depth distribution
p7_combined <- ggplot(depth_data, aes(x = TotalDepth, fill = caller)) +
    geom_histogram(bins = 50, alpha = 0.6, position = "identity") +
    facet_wrap(~caller) +
    labs(title = "Read Depth Distribution (All Samples)",
         x = "log10(Total Depth)", y = "Count") +
    theme_bw()
ggsave(file.path("ch", "figures", "depth_histogram_individual_all.pdf"), p7_individual, width = 13, height = 10)
ggsave(file.path("ch", "figures", "depth_histogram_combined_all.pdf"), p7_combined, width = 10, height = 6)

p7_box <- ggplot(depth_data, aes(x = caller, y = TotalDepth, fill = caller)) +
    geom_boxplot(alpha = 0.7) +
    scale_y_log10() +
    labs(title = "Read Depth Distribution by Caller (All Samples)",
         x = "Caller", y = "log10(Total Depth)") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none")
ggsave(file.path("ch", "figures", "depth_boxplot_combined_all.pdf"), p7_box, width = 8, height = 6)


# -----------------
# 5. READ DEPTH COMPARISON (Total and Alt Depth)
# -----------------
depth_data <- map_dfr(sample_ids, function(id) {
    map_dfr(names(vep[[id]]), function(caller) {
        create_variant_id(vep[[id]][[caller]]) %>%
            mutate(
                sample = id,
                caller = caller,
                TotalDepth = as.numeric(Sample.Depth),
                AltDepth = as.numeric(Sample.AltDepth)
            ) %>%
            select(variant_id, sample, caller, TotalDepth, AltDepth)
    })
})

# ========================
# TOTAL DEPTH
# ========================
# Per sample depth distribution
p7_individual <- ggplot(depth_data, aes(x = TotalDepth, fill = caller)) +
    geom_histogram(bins = 50, alpha = 0.6, position = "identity") +
    scale_x_log10() +
    facet_grid(sample ~ caller) +
    labs(title = "Total Read Depth Distribution by Sample and Caller",
         x = "Total Depth (log10)", y = "Count") +
    theme_bw()

# Combined depth distribution
p7_combined <- ggplot(depth_data, aes(x = TotalDepth, fill = caller)) +
    geom_histogram(bins = 50, alpha = 0.6, position = "identity") +
    scale_x_log10() +
    facet_wrap(~caller) +
    labs(title = "Total Read Depth Distribution (All Samples)",
         x = "Total Depth (log10)", y = "Count") +
    theme_bw() +
    theme(legend.position = "bottom")

# Box plot
p7_box <- ggplot(depth_data, aes(x = caller, y = TotalDepth, fill = caller)) +
    geom_boxplot(alpha = 0.7) +
    scale_y_log10() +
    labs(title = "Total Read Depth Distribution by Caller (All Samples)",
         x = "Caller", y = "Total Depth (log10)") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none")

ggsave(file.path("ch", "figures", "total_depth_histogram_individual_all.pdf"), p7_individual, width = 13, height = 10)
ggsave(file.path("ch", "figures", "total_depth_histogram_combined_all.pdf"), p7_combined, width = 10, height = 6)
ggsave(file.path("ch", "figures", "total_depth_boxplot_combined_all.pdf"), p7_box, width = 8, height = 6)

# ========================
# ALT DEPTH
# ========================
# Per sample alt depth distribution
p8_individual <- ggplot(depth_data, aes(x = AltDepth, fill = caller)) +
    geom_histogram(bins = 50, alpha = 0.6, position = "identity") +
    scale_x_log10() +
    facet_grid(sample ~ caller) +
    labs(title = "Alt Read Depth Distribution by Sample and Caller",
         x = "Alt Depth (log10)", y = "Count") +
    theme_bw()

# Combined alt depth distribution
p8_combined <- ggplot(depth_data, aes(x = AltDepth, fill = caller)) +
    geom_histogram(bins = 50, alpha = 0.6, position = "identity") +
    scale_x_log10() +
    facet_wrap(~caller) +
    labs(title = "Alt Read Depth Distribution (All Samples)",
         x = "Alt Depth (log10)", y = "Count") +
    theme_bw() +
    theme(legend.position = "bottom")

# Box plot
p8_box <- ggplot(depth_data, aes(x = caller, y = AltDepth, fill = caller)) +
    geom_boxplot(alpha = 0.7) +
    scale_y_log10() +
    labs(title = "Alt Read Depth Distribution by Caller (All Samples)",
         x = "Caller", y = "Alt Depth (log10)") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none")

ggsave(file.path("ch", "figures", "alt_depth_histogram_individual_all.pdf"), p8_individual, width = 13, height = 10)
ggsave(file.path("ch", "figures", "alt_depth_histogram_combined_all.pdf"), p8_combined, width = 10, height = 6)
ggsave(file.path("ch", "figures", "alt_depth_boxplot_combined_all.pdf"), p8_box, width = 8, height = 6)

# -----------------
# FILTER VEP FOR AML GENES AND COMBINE INTO MEGATABLE
# -----------------
aml_genes <- read_excel(here("ch", "ss", "20251110_aml_ngs.xlsx"))
aml_genes <- aml_genes$Gene
aml_genes <- aml_genes[!is.na(aml_genes)]
aml_genes <- as.character(aml_genes)

aml_all <- map_dfr(names(vep), function(sample_id) {
    map_dfr(names(vep[[sample_id]]), function(caller) {
        vep[[sample_id]][[caller]] %>%
            filter(Gene %in% aml_genes) %>%
            mutate(
                vep_sample_id = sample_id,
                caller = caller,
                # Ensure numeric columns are numeric
                across(c(Start, Variant.LoF_level, STRAND, Sample.Depth,
                         Sample.AltDepth, Sample.AltFrac), as.numeric)
            )
    })
})

aml_all2 <- aml_all %>%
    left_join(key, by = c("Sample.ID" = "RGC_sample_name")) %>%
    relocate(vep_sample_id, .before = 1)  # Move vep_sample_id to first column

write_csv(aml_all2, file.path("ch", "data", "pilot_somatic_aml.csv"))

aml_all2_mutect2 <- aml_all2 %>% filter(caller == "mutect2")
aml_all2_vardict <- aml_all2 %>% filter(caller == "vardict")
aml_all2_pmbb <- aml_all2 %>% filter(caller == "pmbb")

write_csv(aml_all2_mutect2, file.path("ch", "data", "pilot_mutect2_aml.csv"))
write_csv(aml_all2_vardict, file.path("ch", "data", "pilot_vardict_aml.csv"))
write_csv(aml_all2_pmbb, file.path("ch", "data", "pilot_pmbb_aml.csv"))

# -----------------
# FILTER VEP FOR CH GENES AND COMBINE INTO MEGATABLE
# -----------------
ch_genes <- read_excel(here("ch", "ss", "ch_genes.xlsx"), col_names = FALSE)
ch_genes <- as.factor(ch_genes$...1)

ch_all <- map_dfr(names(vep), function(sample_id) {
    map_dfr(names(vep[[sample_id]]), function(caller) {
        vep[[sample_id]][[caller]] %>%
            filter(Gene %in% ch_genes) %>%
            mutate(
                vep_sample_id = sample_id,
                caller = caller,
                # Ensure numeric columns are numeric
                across(c(Start, Variant.LoF_level, STRAND, Sample.Depth,
                         Sample.AltDepth, Sample.AltFrac), as.numeric)
            )
    })
})

ch_all2 <- ch_all %>%
    left_join(key, by = c("Sample.ID" = "RGC_sample_name")) %>%
    relocate(vep_sample_id, .before = 1)  # Move vep_sample_id to first column

write_csv(ch_all2, file.path("ch", "data", "pilot_somatic_ch_lof1.csv"))

ch_all2_mutect2 <- ch_all2 %>% filter(caller == "mutect2")
ch_all2_vardict <- ch_all2 %>% filter(caller == "vardict")
ch_all2_pmbb <- ch_all2 %>% filter(caller == "pmbb")

write_csv(ch_all2_mutect2, file.path("ch", "data", "pilot_mutect2_ch_lof1.csv"))
write_csv(ch_all2_vardict, file.path("ch", "data", "pilot_vardict_ch_lof1.csv"))
write_csv(ch_all2_pmbb, file.path("ch", "data", "pilot_pmbb_ch_lof1.csv"))

# -----------------
# SAVE MEGA VERSION - SEPARATED BY CALLER
# -----------------
all <- map_dfr(names(vep), function(sample_id) {
    map_dfr(names(vep[[sample_id]]), function(caller) {
        vep[[sample_id]][[caller]] %>%
            mutate(
                vep_sample_id = sample_id,
                caller = caller,
                # Ensure numeric columns are numeric
                across(c(Start, Variant.LoF_level, STRAND, Sample.Depth,
                         Sample.AltDepth, Sample.AltFrac), as.numeric)
            )
    })
})

x <- all %>%
    left_join(key, by = c("Sample.ID" = "RGC_sample_name")) %>%
    relocate(vep_sample_id, .before = 1)  # Move vep_sample_id to first column

mutect2_variants <- x %>% filter(caller == "mutect2")
vardict_variants <- x %>% filter(caller == "vardict")
pmbb_variants <- x %>% filter(caller == "pmbb")

write_csv(mutect2_variants, file.path("ch", "data", "pilot_mutect2_lof1_all.csv"))
write_csv(vardict_variants, file.path("ch", "data", "pilot_vardict_lof1_all.csv"))
write_csv(pmbb_variants, file.path("ch", "data", "pilot_pmbb_lof1_all.csv"))
write_csv(x, file.path("ch", "data", "pilot_somatic_lof1_all.csv"))

mutect2_vep <- map(vep, ~ list(mutect2 = .x$mutect2))
vardict_vep <- map(vep, ~ list(vardict = .x$vardict))
pmbb_vep <- map(vep, ~ list(pmbb = .x$pmbb))

saveRDS(mutect2_vep, file.path("ch", "data", "pilot_mutect2_lof1_all.rds"))
saveRDS(vardict_vep, file.path("ch", "data", "pilot_vardict_lof1_all.rds"))
saveRDS(pmbb_vep, file.path("ch", "data", "pilot_pmbb_lof1_all.rds"))
saveRDS(vep, file.path("ch", "data", "pilot_somatic_lof1_all.rds"))
