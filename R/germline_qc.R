# ============================================================
# GERMLINE GENE QC
#
# Post-VEP filtering of germline variant call reports down to
# high-confidence per-gene calls + diagnostic histograms.
#
# Expected input columns (one row per sample x variant):
#   Gene                 - HGNC gene symbol
#   Bio.type             - Ensembl biotype ("protein_coding" kept)
#   Sample.ID            - sample ID
#   Sample.AltFrac        - variant allele fraction (VAF)
#   Sample.Depth          - total read depth (DP)
#   Sample.AltDepth       - alt-allele read depth (AD)
#   Variant.Consequence   - VEP consequence
#   Variant.LoF_level     - Brad's LoF level (1-4)
#   Start                 - genomic position
#
# Missing/uncalled genotype fields are encoded as "." or -1 in these
# reports and are dropped before filtering.
#
# Accepts
#   qc_germline_gene()   one gene   x one data.frame  -> list(step1/2/3)
#   qc_germline_genes()  N genes    x one data.frame  -> named list of the above
#   run_germline_qc()    N genes    x file-or-folder  -> writes *_qc.csv
# ============================================================

#' Should a (possibly VEP-compound, "&"-joined) Variant.Consequence value be
#' dropped by `exclude_consequence`?
#'
#' VEP joins multiple overlapping consequences for one variant with "&"
#' (e.g. "synonymous_variant&splice_region_variant"). A plain
#' `grepl(paste(patterns, collapse = "|"), consequence)` would match that
#' whole compound string on "synonymous" alone and drop a variant that also
#' carries a splice-relevant (potentially pathogenic) annotation. This
#' splits each value on "&" first and only flags it for exclusion if
#' *every* component term matches one of `patterns` -- so a row is kept as
#' long as at least one of its component consequences isn't excluded.
#'
#' @param consequence character vector, possibly with "&"-joined terms.
#' @param patterns character vector of regex substrings to exclude on.
#' @return logical vector, same length as `consequence`.
is_excluded_consequence <- function(consequence, patterns) {
    terms <- strsplit(as.character(consequence), "&", fixed = TRUE)
    pattern <- paste(patterns, collapse = "|")
    vapply(terms, function(t) all(grepl(pattern, t)), logical(1))
}

#' QC-filter one gene's variant calls out of a VEP report csv
#'
#' Applies 3 filtering stages and writes a log of the
#' variant/individual counts and quality-metric summaries at each stage:
#'   1. Gene + protein_coding filter, drop uncalled ("." / -1)
#'   2. Genotype quality filter: keep calls with `af_lower <= VAF <=
#'      af_upper`, `Sample.Depth >= min_depth`, `Sample.AltDepth >=
#'      min_altdepth`.
#'   3. Consequence + LoF level filter: drop rows whose Variant.Consequence,
#'      split on VEP's "&" compound-term separator, has *every* component
#'      term matching a pattern in `exclude_consequence` (a row with even
#'      one non-excluded co-occurring term, e.g.
#'      "synonymous_variant&splice_region_variant", is kept). Also keep
#'      only Variant.LoF_level 1 or 2, defaults to
#'
#'
#' @param data data.frame of VEP report rows; must contain the columns above
#' @param gene_name character(1). Value to match against `data$Gene`.
#' @param dataset_label character(1). label used in file naming/logging
#' @param af_lower,af_upper numeric(1). Inclusive VAF band a call must fall
#'   in to pass step 2. Defaults (0.20-0.80) for heterozygous call.
#' @param min_depth numeric(1). Minimum total read depth (Sample.Depth)
#'   required at step 2. Defaults 14.
#' @param min_altdepth numeric(1). Minimum alt-allele read depth
#'   (Sample.AltDepth) required at step 2. Defaults 5.
#' @param exclude_consequence character vector of substrings matched
#'   against `Variant.Consequence` and dropped at step 3. Default excludes synonymous,
#'   upstream, 3'/5' UTR, inframe, and intronic consequences.
#' @param log_dir directory the per-run QC log is written into (created if
#'   missing). Defaults into "log".
#'
#' @return list(step1 =, step2 =, step3 =) -- data.frames after each
#'   filtering stage. step3 is the final QC'd call set for this gene.
qc_germline_gene <- function(data, gene_name, dataset_label,
                              af_lower = 0.20, af_upper = 0.80,
                              min_depth = 14, min_altdepth = 5,
                              exclude_consequence = c("synonymous", "upstream",
                                                       "3_prime_UTR_variant",
                                                       "5_prime_UTR_variant",
                                                       "inframe", "intron_variant"),
                              log_dir = file.path("log")) {

    # Generate datetime stamp
    datetime_stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")

    # Set up gene-specific log file with datetime
    log_file <- file.path(log_dir,
                          paste0(datetime_stamp, "_postvep_qc_", gene_name, "_", dataset_label, ".txt"))
    dir.create(log_dir, showWarnings = FALSE, recursive = TRUE)

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

    # LOG INPUT PARAMETERS
    log_cat("=== FILTERING PARAMETERS ===\n")
    log_cat("VAF band: [", af_lower, ", ", af_upper, "]\n")
    log_cat("Minimum Total Depth (DP): ", min_depth, "\n")
    log_cat("Minimum Alternate Depth (AD): ", min_altdepth, "\n")
    log_cat("Excluded consequence patterns: ", paste(exclude_consequence, collapse = ", "), "\n")
    log_cat("LoF level filter: 1 or 2 only\n\n")

    # Initial file info
    log_cat("=== INITIAL DATA ===\n")
    log_cat("Total rows: ", nrow(data), "\n")
    log_cat("Total individuals: ", length(unique(data$Sample.ID)), "\n")
    log_cat("Genomic range: ", min(data$Start, na.rm = TRUE), " - ", max(data$Start, na.rm = TRUE), "\n\n")

    # Step 1: Gene and protein coding filter
    log_cat("=== FILTERING STEPS ===\n\n")

    step1 <- data %>%
        filter(Gene == gene_name) %>%
        filter(FILTER == "PASS") %>%
        filter(Bio.type == "protein_coding") %>%
        filter(Sample.AltFrac != ".", Sample.Depth != ".", Sample.AltDepth != ".") %>%
        mutate(
            Sample.AltFrac = as.numeric(Sample.AltFrac),
            Sample.Depth = as.numeric(Sample.Depth),
            Sample.AltDepth = as.numeric(Sample.AltDepth),
            Variant.LoF_level = suppressWarnings(as.numeric(Variant.LoF_level))
        ) %>%
        filter(!is.na(Sample.AltFrac), !is.na(Sample.Depth), !is.na(Sample.AltDepth)) %>%
        filter(Sample.AltFrac != -1, Sample.Depth != -1, Sample.AltDepth != -1)  # Remove uncalled (-1) values

    log_cat("Step 1 (Gene + Protein Coding): ",
            nrow(step1), " variants, ",
            length(unique(step1$Sample.ID)), " individuals\n")
    log_stats(step1, "  Quality Metrics")
    log_cat("\n")

    # Step 2: Read depth and VAF filter
    step2 <- step1 %>%
        filter(Sample.AltFrac >= af_lower, Sample.AltFrac <= af_upper,
               Sample.Depth >= min_depth, Sample.AltDepth >= min_altdepth)

    log_cat("Step 2 (+ Depth/VAF Filter): ",
            nrow(step2), " variants, ",
            length(unique(step2$Sample.ID)), " individuals")
    log_cat(" (removed ", nrow(step1) - nrow(step2), " variants)\n")
    log_stats(step2, "  Quality Metrics")
    log_cat("\n")

    # Step 3: Consequence + LoF level filter
    step3 <- step2 %>%
        filter(!is_excluded_consequence(Variant.Consequence, exclude_consequence)) %>%
        filter(Variant.LoF_level %in% c(1, 2))

    log_cat("Step 3 (+ Consequence + LoF Level Filter): ",
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

#' QC-filter multiple genes' variant calls out of one VEP report data.frame
#'
#' Thin loop around qc_germline_gene(): runs the same three-stage filter
#' once per gene in `gene_names` and returns one result list per gene.
#' Genes with no matching rows in `data$Gene` are skipped with a warning
#' rather than erroring, so a partial gene list against a single-gene
#' report file still works.
#'
#' @param data data.frame of VEP report rows (may contain multiple genes).
#' @param gene_names character vector of Gene values to filter to.
#' @param dataset_label passed through to qc_germline_gene() (see there).
#' @param ... threshold/log overrides passed through to
#'   qc_germline_gene() (af_lower, af_upper, min_depth, min_altdepth,
#'   exclude_consequence, log_dir).
#'
#' @return named list, one element per gene actually found in `data`, each
#'   the list(step1 =, step2 =, step3 =) from qc_germline_gene().
qc_germline_genes <- function(data, gene_names, dataset_label, ...) {
    results <- list()
    for (gene in gene_names) {
        if (!gene %in% unique(data$Gene)) {
            warning("Gene '", gene, "' not found in data$Gene -- skipping")
            next
        }
        results[[gene]] <- qc_germline_gene(data, gene, dataset_label, ...)
    }
    results
}

#' Run germline QC over a CSV file or a folder of CSVs, one output file per
#' (source file, gene) pair
#'
#' Two input modes, chosen automatically from `input_path`:
#'   - Single .csv file: read once, then filtered for every gene in
#'     `genes` that's actually present in its Gene column (a report can
#'     hold more than one gene). Genes not present are skipped with a
#'     warning.
#'   - Folder: every "*.csv" file directly inside the folder whose file
#'     name contains one of `genes` (case-insensitive substring match,
#'     e.g. matches "*.BRCA1.vep.report.csv" against gene "BRCA1") is
#'     read and filtered for that gene. This matches the one-gene-per-file
#'     layout produced by germline_gene_extract.smk. Files that don't
#'     match any requested gene are left untouched.
#'
#' Output: for every (file, gene) pair with a QC'd result, writes the
#' final (step3) filtered rows to
#'   "<original_basename>_qc.csv"
#' in the SAME directory as the source file.
#'
#' @param input_path path to a single .csv file, or a directory containing
#'   .csv files.
#' @param genes character vector of gene symbols to filter to.
#' @param dataset_label free-text label for QC logs (see
#'   qc_germline_gene()). Defaults to the source file's basename (without
#'   extension) if not supplied, so each output/log is traceable to its
#'   input file.
#' @param ... threshold/log overrides passed through to
#'   qc_germline_gene() (af_lower, af_upper, min_depth, min_altdepth,
#'   exclude_consequence, log_dir).
#'
#' @return invisible named list, keyed "<source_basename>::<gene>", of the
#'   qc_germline_gene() results -- useful for feeding create_diagnostic_pdf()
#'   without re-reading the output CSVs.
run_germline_qc <- function(input_path, genes, dataset_label = NULL, ...) {
    stopifnot(length(genes) >= 1)

    is_dir <- dir.exists(input_path)
    if (!is_dir && !file.exists(input_path)) {
        stop("input_path does not exist: ", input_path)
    }

    if (is_dir) {
        csv_files <- list.files(input_path, pattern = "\\.csv$",
                                 full.names = TRUE, ignore.case = TRUE)
        gene_pattern <- paste(genes, collapse = "|")
        csv_files <- csv_files[grepl(gene_pattern, basename(csv_files), ignore.case = TRUE)]
        if (length(csv_files) == 0) {
            stop("No CSV files in ", input_path, " matched any of: ", paste(genes, collapse = ", "))
        }
    } else {
        if (!grepl("\\.csv$", input_path, ignore.case = TRUE)) {
            stop("input_path must be a .csv file or a directory: ", input_path)
        }
        csv_files <- input_path
    }

    all_results <- list()

    for (f in csv_files) {
        data <- read.csv(f, stringsAsFactors = FALSE)
        out_dir <- dirname(f)
        base <- tools::file_path_sans_ext(basename(f))
        label <- if (!is.null(dataset_label)) dataset_label else base

        # In folder mode, only try genes whose name actually matched this
        # file (avoids e.g. re-filtering a BRCA1 file for BRCA2). In
        # single-file mode the file could hold several genes, so try them
        # all -- qc_germline_genes() skips any that aren't present.
        genes_to_try <- if (is_dir) {
            genes[sapply(genes, function(g) grepl(g, basename(f), ignore.case = TRUE))]
        } else {
            genes
        }

        gene_results <- qc_germline_genes(data, genes_to_try, label, ...)

        for (gene in names(gene_results)) {
            res <- gene_results[[gene]]
            out_file <- file.path(out_dir, paste0(base, "_qc.csv"))
            write.csv(res$step3, out_file, row.names = FALSE)
            message("Wrote ", nrow(res$step3), " QC'd variants for ", gene, " -> ", out_file)
            all_results[[paste0(base, "::", gene)]] <- res
        }
    }

    invisible(all_results)
}

# ========================
# PLOTTING FUNCTION (NO XLIM)
# ========================

#' Build a 3x3 grid of QC diagnostic histograms (Depth/AltDepth/VAF x
#' Initial/After-filter/Final) for one gene and write it to a PDF.
#'
#' @param step1,step2,step3 the three stages returned by
#'   qc_germline_gene() (or one element of qc_germline_genes()/
#'   run_germline_qc()'s result list).
#' @param gene_name used only for the plot title.
#' @param output_file path to write the PDF to.
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
        list(data = step3, label = "After Consequence/LoF Filtering")
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
