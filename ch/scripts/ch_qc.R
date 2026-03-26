ch_qc <- function(data, dataset_label, AF_LOWER, AF_UPPER, DP, AD, GQ) {

    datetime_stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
    log_file <- file.path("ch", "log", paste0(datetime_stamp, "_postvep_qc_", dataset_label, ".txt"))
    dir.create(file.path("ch", "log"), showWarnings = FALSE, recursive = TRUE)

    log_cat <- function(...) {
        message <- paste0(...)
        cat(message)
        cat(message, file = log_file, append = TRUE)
    }

    log_stats <- function(data, label) {
        log_cat(label, ":\n")
        for (col in c("Sample.AltFrac", "Sample.Depth", "Sample.AltDepth", "Sample.GQ")) {
            vals <- data[[col]][data[[col]] != -1]
            if (length(vals) > 0) {
                s <- summary(vals)
                log_cat("  ", col, " - Min: ", round(s[1], 3), ", Median: ", round(s[3], 3), ", Max: ", round(s[6], 3), "\n")
            }
        }
    }

    cat("", file = log_file, append = FALSE)
    log_cat("=== POST-VEP QC LOG ===\n")
    log_cat("Dataset: ", dataset_label, "\n\n")
    log_cat("=== INITIAL DATA ===\n")
    log_cat("Total rows: ", nrow(data), "\n")
    log_cat("Total individuals: ", length(unique(data$Sample.ID)), "\n\n")

    # Step 1: protein coding + convert columns + remove missing
    step1 <- data %>%
        filter(Bio.type == "protein_coding") %>%
        filter(Sample.AltFrac != ".", Sample.Depth != ".", Sample.AltDepth != ".", Sample.GQ != ".") %>%
        mutate(
            Sample.AltFrac   = as.numeric(Sample.AltFrac),
            Sample.Depth     = as.numeric(Sample.Depth),
            Sample.AltDepth  = as.numeric(Sample.AltDepth),
            Sample.GQ        = as.numeric(gsub("\\[|\\]", "", Sample.GQ))
        ) %>%
        filter(!is.na(Sample.AltFrac), !is.na(Sample.Depth), !is.na(Sample.AltDepth), !is.na(Sample.GQ)) %>%
        filter(Sample.AltFrac != -1, Sample.Depth != -1, Sample.AltDepth != -1, Sample.GQ > 0)

    log_cat("Step 1 (Protein Coding + valid values): ", nrow(step1), " variants, ", length(unique(step1$Sample.ID)), " individuals\n")
    log_stats(step1, "  Quality Metrics")
    log_cat("\n")

    # Step 2: GQ filter
    step2 <- step1 %>%
        filter(Sample.GQ >= GQ)

    log_cat("Step 2 (GQ >= ", GQ, "): ", nrow(step2), " variants, ", length(unique(step2$Sample.ID)), " individuals")
    log_cat(" (removed ", nrow(step1) - nrow(step2), " variants)\n")
    log_stats(step2, "  Quality Metrics")
    log_cat("\n")

    # Step 3: VAF filter
    step3 <- step2 %>%
        filter(Sample.AltFrac >= AF_LOWER, Sample.AltFrac <= AF_UPPER)

    log_cat("Step 3 (VAF ", AF_LOWER, "-", AF_UPPER, "): ", nrow(step3), " variants, ", length(unique(step3$Sample.ID)), " individuals")
    log_cat(" (removed ", nrow(step2) - nrow(step3), " variants)\n")
    log_stats(step3, "  Quality Metrics")
    log_cat("\n")

    # Step 4: DP + AD filter
    step4 <- step3 %>%
        filter(Sample.Depth >= DP, Sample.AltDepth >= AD)

    log_cat("Step 4 (DP >= ", DP, ", AD >= ", AD, "): ", nrow(step4), " variants, ", length(unique(step4$Sample.ID)), " individuals")
    log_cat(" (removed ", nrow(step3) - nrow(step4), " variants)\n")
    log_stats(step4, "  Quality Metrics")
    log_cat("\n")

    # Final summary
    log_cat("=== FINAL RESULTS ===\n")
    log_cat("Final rows: ", nrow(step4), "\n")
    log_cat("Final individuals: ", length(unique(step4$Sample.ID)), "\n")
    log_cat("All consequences: ", paste(unique(step4$Variant.Consequence), collapse = ", "), "\n\n")

    if (nrow(step4) > 0) {
        log_cat("LoF levels (variants | unique individuals):\n")
        for (lof in sort(unique(step4$Variant.LoF_level))) {
            variant_count <- sum(step4$Variant.LoF_level == lof)
            ind_count <- step4 %>% filter(Variant.LoF_level == lof) %>% pull(Sample.ID) %>% unique() %>% length()
            log_cat("  ", lof, ": ", variant_count, " variants | ", ind_count, " individuals\n")
        }
    } else {
        log_cat("No variants remaining after filtering\n")
    }

    log_cat("\nLog saved to: ", log_file, "\n")

    return(list(step1 = step1, step2 = step2, step3 = step3, step4 = step4))
}
