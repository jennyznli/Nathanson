# ========================
# FILTERING FUNCTION WITH GENE-SPECIFIC LOGGING
# ========================
process_brca_data <- function(data, gene_name, dataset_label) {

    # Generate datetime stamp
    datetime_stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")

    # Set up gene-specific log file with datetime
    log_file <- file.path("ch", "log",
                          paste0(datetime_stamp, "_postvep_qc_", gene_name, "_", dataset_label, ".txt"))
    dir.create(file.path("ch", "log"), showWarnings = FALSE, recursive = TRUE)

    # Function to write to both console and log
    log_cat <- function(...) {
        message <- paste0(...)
        cat(message)
        cat(message, file = log_file, append = TRUE)
    }

    # Function to calculate and log stats (excluding -1 values)
    log_stats <- function(data, label) {
        log_cat(label, ":\n")

        # VAF stats (exclude -1)
        vaf_clean <- data$Sample.AltFrac[data$Sample.AltFrac != -1]
        if (length(vaf_clean) > 0) {
            vaf_stats <- summary(vaf_clean)
            log_cat("  VAF - Min: ", round(vaf_stats[1], 3),
                    ", Median: ", round(vaf_stats[3], 3),
                    ", Max: ", round(vaf_stats[6], 3), "\n")
        } else {
            log_cat("  VAF - No valid data\n")
        }

        # Total Depth stats (exclude -1)
        depth_clean <- data$Sample.Depth[data$Sample.Depth != -1]
        if (length(depth_clean) > 0) {
            depth_stats <- summary(depth_clean)
            log_cat("  Total Depth - Min: ", round(depth_stats[1], 1),
                    ", Median: ", round(depth_stats[3], 1),
                    ", Max: ", round(depth_stats[6], 1), "\n")
        } else {
            log_cat("  Total Depth - No valid data\n")
        }

        # Alt Depth stats (exclude -1)
        altdepth_clean <- data$Sample.AltDepth[data$Sample.AltDepth != -1]
        if (length(altdepth_clean) > 0) {
            altdepth_stats <- summary(altdepth_clean)
            log_cat("  Alt Depth - Min: ", round(altdepth_stats[1], 1),
                    ", Median: ", round(altdepth_stats[3], 1),
                    ", Max: ", round(altdepth_stats[6], 1), "\n")
        } else {
            log_cat("  Alt Depth - No valid data\n")
        }
    }

    # Start fresh log
    cat("", file = log_file, append = FALSE)
    log_cat("=== POST-VEP GENE QC LOG ===\n")
    log_cat("Gene: ", gene_name, "\n")
    log_cat("Dataset: ", dataset_label, "\n\n")

    # Initial file info
    log_cat("=== INITIAL DATA ===\n")
    log_cat("Total rows: ", nrow(data), "\n")
    log_cat("Total individuals: ", length(unique(data$Sample.ID)), "\n")
    log_cat("Genomic range: ", min(data$Start, na.rm = TRUE), " - ", max(data$Start, na.rm = TRUE), "\n\n")

    remove_pattern <- c("synonymous", "upstream", "3_prime_UTR_variant",
                        "5_prime_UTR_variant", "inframe")

    # Step 1: Gene and protein coding filter
    log_cat("=== FILTERING STEPS ===\n\n")

    step1 <- data %>%
        filter(Gene == gene_name) %>%
        filter(Bio.type == "protein_coding") %>%
        filter(Sample.AltFrac != ".", Sample.Depth != ".", Sample.AltDepth != ".") %>%
        mutate(
            Sample.AltFrac = as.numeric(Sample.AltFrac),
            Sample.Depth = as.numeric(Sample.Depth),
            Sample.AltDepth = as.numeric(Sample.AltDepth)
        ) %>%
        filter(!is.na(Sample.AltFrac), !is.na(Sample.Depth), !is.na(Sample.AltDepth)) %>%
        filter(Sample.AltFrac != -1, Sample.Depth != -1, Sample.AltDepth != -1)  # Remove -1 values

    log_cat("Step 1 (Gene + Protein Coding): ",
            nrow(step1), " variants, ",
            length(unique(step1$Sample.ID)), " individuals\n")
    log_stats(step1, "  Quality Metrics")
    log_cat("\n")

    # Step 2: Read depth and VAF filter
    step2 <- step1 %>%
        filter(Sample.AltFrac > 0.3, Sample.AltFrac < 0.7,
               Sample.Depth > 20, Sample.AltDepth > 5)

    log_cat("Step 2 (+ Depth/VAF Filter): ",
            nrow(step2), " variants, ",
            length(unique(step2$Sample.ID)), " individuals")
    log_cat(" (removed ", nrow(step1) - nrow(step2), " variants)\n")
    log_stats(step2, "  Quality Metrics")
    log_cat("\n")

    # Step 3: Consequence filter
    step3 <- step2 %>%
        filter(!grepl(paste(remove_pattern, collapse="|"), Variant.Consequence),
               Variant.Consequence != "intron_variant")

    log_cat("Step 3 (+ Consequence Filter): ",
            nrow(step3), " variants, ",
            length(unique(step3$Sample.ID)), " individuals")
    log_cat(" (removed ", nrow(step2) - nrow(step3), " variants)\n")
    log_stats(step3, "  Quality Metrics")

    # Final summary
    log_cat("=== FINAL RESULTS ===\n")
    log_cat("Final genomic range: ", min(step3$Start, na.rm = TRUE), " - ", max(step3$Start, na.rm = TRUE), "\n\n")
    log_cat("All consequences: ", paste(unique(step3$Variant.Consequence), collapse = ", "), "\n\n")

    if (nrow(step3) > 0) {
        # LoF levels with both variant and individual counts, sorted 1-4
        log_cat("LoF levels (variants | unique individuals):\n")
        lof_levels <- sort(unique(step3$Variant.LoF_level))
        for (lof in lof_levels) {
            variant_count <- sum(step3$Variant.LoF_level == lof)
            individuals_in_lof <- step3 %>%
                filter(Variant.LoF_level == lof) %>%
                pull(Sample.ID) %>%
                unique() %>%
                length()
            log_cat("  ", lof, ": ", variant_count, " variants | ", individuals_in_lof, " individuals\n")
        }
    } else {
        log_cat("No variants remaining after filtering\n")
    }

    log_cat("\nLog saved to: ", log_file, "\n")

    return(list(step1 = step1, step2 = step2, step3 = step3))
}

# ========================
# PLOTTING FUNCTION (NO XLIM)
# ========================
create_diagnostic_pdf <- function(step1, step2, step3, gene_name, output_file) {

    plots <- list()
    metrics <- list(
        list(col = "Sample.Depth", name = "Total Depth", xlab = "Total Depth (DP)",
             vlines = 20, color = "steelblue"),
        list(col = "Sample.AltDepth", name = "Alt Depth", xlab = "Alt Depth (AD)",
             vlines = 5, color = "darkgreen"),
        list(col = "Sample.AltFrac", name = "VAF", xlab = "Variant Allele Fraction",
             vlines = c(0.3, 0.7), color = "purple")
    )

    steps <- list(
        list(data = step1, label = "Initial"),
        list(data = step2, label = "After Depth/VAF"),
        list(data = step3, label = "Final")
    )

    idx <- 1
    for (s in 1:3) {
        step_data <- steps[[s]]$data
        step_label <- steps[[s]]$label
        n <- length(unique(step_data$Sample.ID))

        for (m in metrics) {
            plots[[idx]] <- ggplot(step_data, aes(x = .data[[m$col]])) +
                geom_histogram(bins = 100, fill = m$color, alpha = 0.7) +
                geom_vline(xintercept = m$vlines, color = "red",
                           linetype = "dashed", linewidth = 1) +
                labs(title = paste0(step_label, ": ", m$name),
                     subtitle = paste("n =", n, "individuals"),
                     x = m$xlab, y = "Count") +
                theme_minimal() +
                theme(plot.title = element_text(size = 10),
                      plot.subtitle = element_text(size = 8))
            idx <- idx + 1
        }
    }

    pdf(output_file, width = 15, height = 12)

    grid.arrange(
        grobs = plots,
        ncol = 3,
        top = textGrob(paste(gene_name, "Filtering Pipeline"),
                       gp = gpar(fontsize = 16, fontface = "bold"))
    )

    dev.off()
}
