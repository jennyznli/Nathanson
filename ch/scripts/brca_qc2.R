# ========================
# POST VEP FILTERING
# ========================
process_brca_data <- function(data, gene_name, dataset_label, AF_LOWER, AF_UPPER,
                              DP, AD, GQ = NULL) {

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

    # Check if GQ column exists
    has_gq <- "Sample.GQ" %in% names(data)

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

        # GQ (only if column exists)
        if (has_gq) {
            gq_clean <- data$Sample.GQ[data$Sample.GQ != -1]
            if (length(gq_clean) > 0) {
                gq_stats <- summary(gq_clean)
                log_cat("  GQ - Min: ", round(gq_stats[1], 1),
                        ", Median: ", round(gq_stats[3], 1),
                        ", Max: ", round(gq_stats[6], 1), "\n")
            } else {
                log_cat("  Genome Quality - No valid data\n")
            }
        } else {
            log_cat("  GQ - Column not available\n")
        }
    }

    # Start fresh log
    cat("", file = log_file, append = FALSE)
    log_cat("=== POST-VEP GENE QC LOG ===\n")
    log_cat("Gene: ", gene_name, "\n")
    log_cat("Dataset: ", dataset_label, "\n")
    log_cat("Timestamp: ", datetime_stamp, "\n\n")

    # LOG INPUT PARAMETERS
    log_cat("=== FILTERING PARAMETERS ===\n")
    log_cat("VAF Lower Threshold: ", AF_LOWER, "\n")
    log_cat("VAF Upper Threshold: ", AF_UPPER, "\n")
    log_cat("Minimum Total Depth (DP): ", DP, "\n")
    log_cat("Minimum Alternate Depth (AD): ", AD, "\n")
    if (has_gq && !is.null(GQ)) {
        log_cat("Minimum Genotype Quality (GQ): ", GQ, "\n\n")
    } else {
        log_cat("Minimum Genotype Quality (GQ): Not applied (column not available)\n\n")
    }

    # Initial file info
    log_cat("=== INITIAL DATA ===\n")
    log_cat("Total rows: ", nrow(data), "\n")
    log_cat("Total individuals: ", length(unique(data$Sample.ID)), "\n")
    log_cat("Genomic range: ", min(data$Start, na.rm = TRUE), " - ", max(data$Start, na.rm = TRUE), "\n\n")

    # Step 1: Gene and protein coding filter
    log_cat("=== FILTERING STEPS ===\n\n")

    step1 <- data %>%
        filter(Gene == gene_name) %>%
        filter(Bio.type == "protein_coding")

    # Conditional filtering based on GQ availability
    if (has_gq) {
        step1 <- step1 %>%
            filter(Sample.AltFrac != ".", Sample.Depth != ".", Sample.AltDepth != ".", Sample.GQ != ".") %>%
            mutate(
                Sample.AltFrac = as.numeric(Sample.AltFrac),
                Sample.Depth = as.numeric(Sample.Depth),
                Sample.AltDepth = as.numeric(Sample.AltDepth),
                Sample.GQ = as.numeric(gsub("\\[|\\]", "", Sample.GQ))
            ) %>%
            filter(!is.na(Sample.AltFrac), !is.na(Sample.Depth), !is.na(Sample.AltDepth), !is.na(Sample.GQ)) %>%
            filter(Sample.AltFrac != -1, Sample.Depth != -1, Sample.AltDepth != -1, Sample.GQ > 0)
    } else {
        step1 <- step1 %>%
            filter(Sample.AltFrac != ".", Sample.Depth != ".", Sample.AltDepth != ".") %>%
            mutate(
                Sample.AltFrac = as.numeric(Sample.AltFrac),
                Sample.Depth = as.numeric(Sample.Depth),
                Sample.AltDepth = as.numeric(Sample.AltDepth)
            ) %>%
            filter(!is.na(Sample.AltFrac), !is.na(Sample.Depth), !is.na(Sample.AltDepth)) %>%
            filter(Sample.AltFrac != -1, Sample.Depth != -1, Sample.AltDepth != -1)
    }

    log_cat("Step 1 (Gene + Protein Coding): ",
            nrow(step1), " variants, ",
            length(unique(step1$Sample.ID)), " individuals\n")
    log_stats(step1, "  Quality Metrics")
    log_cat("\n")

    # Step 2: Read depth and VAF filter
    if (has_gq && !is.null(GQ)) {
        step2 <- step1 %>%
            filter(Sample.AltFrac >= AF_LOWER, Sample.AltFrac <= AF_UPPER,
                   Sample.Depth >= DP, Sample.AltDepth >= AD,
                   Sample.GQ >= GQ)
    } else {
        step2 <- step1 %>%
            filter(Sample.AltFrac >= AF_LOWER, Sample.AltFrac <= AF_UPPER,
                   Sample.Depth >= DP, Sample.AltDepth >= AD)
    }

    log_cat("Step 2 (+ Depth/VAF Filter): ",
            nrow(step2), " variants, ",
            length(unique(step2$Sample.ID)), " individuals")
    log_cat(" (removed ", nrow(step1) - nrow(step2), " variants)\n")
    log_stats(step2, "  Quality Metrics")
    log_cat("\n")

    # Step 3: Consequence filter
    step3 <- step2 %>%
        filter(!(Variant.Consequence %in% c("intron_variant", "synonymous_variant",
                                            "3_prime_UTR_variant", "5_prime_UTR_variant")))


    # log_cat("Step 3 (+ Consequence Filter): ",
    #         nrow(step3), " variants, ",
    #         length(unique(step3$Sample.ID)), " individuals")
    # log_cat(" (removed ", nrow(step2) - nrow(step3), " variants)\n")
    # log_stats(step3, "  Quality Metrics")
    # log_cat("\n")

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

    # Check if GQ column exists in any of the steps
    has_gq <- "Sample.GQ" %in% names(step1) ||
        "Sample.GQ" %in% names(step2) ||
        "Sample.GQ" %in% names(step3)

    # Define metrics - conditionally include GQ
    metrics <- list(
        list(col = "Sample.Depth", name = "Total Depth", xlab = "Total Depth (DP)",
             color = "steelblue"),
        list(col = "Sample.AltDepth", name = "Alt Depth", xlab = "Alt Depth (AD)",
             color = "darkgreen"),
        list(col = "Sample.AltFrac", name = "VAF", xlab = "Variant Allele Fraction",
             color = "purple")
    )

    # Add GQ metric only if column exists
    if (has_gq) {
        metrics[[4]] <- list(col = "Sample.GQ", name = "Genome Quality",
                             xlab = "Genome Quality",
                             color = "pink")
    }

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
            # Check if this specific metric column exists in this step's data
            if (m$col %in% names(step_data)) {
                plots[[idx]] <- ggplot(step_data, aes(x = .data[[m$col]])) +
                    geom_histogram(bins = 100, fill = m$color, alpha = 0.7) +
                    labs(title = paste0(step_label, ": ", m$name),
                         subtitle = paste("n =", n, "individuals"),
                         x = m$xlab, y = "Count") +
                    theme_minimal() +
                    theme(plot.title = element_text(size = 10),
                          plot.subtitle = element_text(size = 8))
                idx <- idx + 1
            }
        }
    }

    # Determine number of columns based on whether GQ exists
    ncol <- if (has_gq) 4 else 3

    pdf(output_file, width = 15, height = 12)
    grid.arrange(
        grobs = plots,
        ncol = ncol,
        top = textGrob(paste(gene_name, "Filtering Pipeline"),
                       gp = gpar(fontsize = 16, fontface = "bold"))
    )
    dev.off()
}
