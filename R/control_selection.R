# ============================================================
# Non-Cancer Control Selection Script
# ============================================================
# This script selects non-cancer control samples from PMBB
# phenotype files with flexible filtering options.
#
# Usage:
#   select_controls(
#       control_name = "healthy_controls",  # Name for the control group
#       exclude_codes = ...                 # ICD codes to exclude (default: malignant neoplasms)
#       gender_filter = NULL,               # "Male", "Female", or NULL
#       crep_filter = NULL,                 # TRUE - include only CREP, FALSE - exclude CREP, NULL - no filter (so including CREP)
#       age_filter = NULL,                  # int - minimum age, vector - range, or NULL
#       data_dir = NULL,                    # Directory for outputs
#       pmbb_dir = NULL,                    # Directory for PMBB input files
#       log_dir = NULL,                     # Directory for log files
#       output_prefix = NULL                # Filename prefix (defaults to control_name)
#   )
#
# Requires:
#   PMBB-Release-2024-3.0_phenotype_condition_occurrence.txt
#   PMBB-Release-2024-3.0_covariates.txt
#   PMBB-Release-2024-3.0_phenotype_person.txt
#
# Output:
#   - Filtered control patient lists and covariates
#   - Log file with selection steps and stats
# ============================================================

select_controls <- function(
        control_name = "controls",
        exclude_codes = c("^C(?!44)", "^Z85(?!\\.828)", "^(?!173)(?:1[4-9][0-9]|20[0-8]|209\\.[0-3])", "^V10(?!\\.83)"),
        gender_filter = NULL,
        crep_filter = NULL,
        age_filter = NULL,
        data_dir = NULL,
        pmbb_dir = NULL,
        log_dir = NULL,
        output_prefix = NULL
) {
    packages <- c("tidyr", "dplyr", "readr", "stringr", "data.table")
    purrr::walk(packages, ~ require(.x, character.only = TRUE))

    # ========================
    # SETUP
    # ========================
    if (is.null(data_dir)) {
        DATA_DIR <- getwd()
    } else {
        DATA_DIR <- data_dir
    }

    if (is.null(pmbb_dir)) {
        PMBB_DIR <- file.path(DATA_DIR, "PMBB", "3.0")
    } else {
        PMBB_DIR <- pmbb_dir
    }

    if (is.null(log_dir)) {
        LOG_DIR <- file.path(DATA_DIR, "log")
    } else {
        LOG_DIR <- log_dir
    }

    if (!dir.exists(DATA_DIR)) dir.create(DATA_DIR, recursive = TRUE)
    if (!dir.exists(LOG_DIR)) dir.create(LOG_DIR, recursive = TRUE)

    OCC <- file.path(PMBB_DIR, "PMBB-Release-2024-3.0_phenotype_condition_occurrence.txt")
    COV <- file.path(PMBB_DIR, "PMBB-Release-2024-3.0_covariates.txt")
    PER <- file.path(PMBB_DIR, "PMBB-Release-2024-3.0_phenotype_person.txt")

    files_to_check <- c(OCC, COV, PER)
    missing_files <- files_to_check[!file.exists(files_to_check)]

    if (length(missing_files) > 0) {
        stop("ERROR: Missing files:\n", paste(missing_files, collapse = "\n"))
    }

    cat("All required PMBB files found successfully.\n")

    occ <- fread(OCC, header = TRUE)
    cov <- fread(COV, header = TRUE)
    per <- fread(PER, header = TRUE)

    # ========================
    # CONFIGURATION
    # ========================
    TIMESTAMP <- format(Sys.time(), "%Y%m%d_%H%M%S")

    if (is.null(output_prefix)) {
        output_prefix <- tolower(gsub("[^A-Za-z0-9]", "_", control_name))
    }

    # Start logging
    LOG_FILE <- file.path(LOG_DIR, paste0(TIMESTAMP, "_", output_prefix, "_control_selection_log.txt"))
    sink(LOG_FILE)

    cat("Non-Cancer Control Selection Log\n")
    cat("================================\n\n")
    cat("Control Name:", control_name, "\n")
    cat("Exclude Codes:", paste(exclude_codes, collapse = ", "), "\n")
    cat("Gender Filter:", ifelse(is.null(gender_filter), "None", gender_filter), "\n")
    cat("CREP Filter:", ifelse(is.null(crep_filter), "None", crep_filter), "\n")
    cat("Age Filter:", ifelse(is.null(age_filter), "None", paste(age_filter, collapse = "-")), "\n\n")

    # ========================
    # HELPER FUNCTIONS
    # ========================
    apply_gender_filter <- function(data, gender_filter) {
        if (is.null(gender_filter)) {
            cat("=== No gender filter applied ===\n")
            return(data)
        }

        cat("=== Applying gender filter:", gender_filter, "===\n")
        cat("Before:", nrow(data), "patients\n")

        filtered <- data[Sequenced_gender == gender_filter]
        cat("After:", nrow(filtered), "patients\n\n")

        return(filtered)
    }

    apply_crep_filter <- function(data, crep_filter) {
        if (is.null(crep_filter)) {
            cat("=== No CREP filter applied ===\n")
            return(data)
        }

        cat("=== Applying CREP filter ===\n")
        cat("Before:", nrow(data), "patients\n")

        if (crep_filter) {
            filtered <- data[CREP_HighRisk_Flag == 1]
            cat("After (CREP only):", nrow(filtered), "patients\n\n")
        } else {
            filtered <- data[CREP_HighRisk_Flag == 0]
            cat("After (no CREP):", nrow(filtered), "patients\n\n")
        }

        return(filtered)
    }

    apply_age_filter <- function(data, age_filter) {
        if (is.null(age_filter)) {
            cat("=== No age filter applied ===\n")
            return(data)
        }

        cat("=== Applying age filter ===\n")
        cat("Before:", nrow(data), "patients\n")

        if (length(age_filter) == 1) {
            filtered <- data[Sample_age >= age_filter[1]]
            cat("After (>=", age_filter[1], "):", nrow(filtered), "patients\n\n")
        } else {
            filtered <- data[Sample_age >= age_filter[1] & Sample_age <= age_filter[2]]
            cat("After (", age_filter[1], "-", age_filter[2], "):", nrow(filtered), "patients\n\n")
        }

        return(filtered)
    }

    # ========================
    # MAIN ANALYSIS
    # ========================

    # Find patients to exclude based on ICD codes
    cat("=== Finding patients to exclude ===\n")
    exclude_pattern <- paste(exclude_codes, collapse = "|")
    exclude_occ <- occ[grepl(exclude_pattern, condition_source_value, perl = TRUE)]
    exclude_patients <- unique(exclude_occ$person_id)

    cat("Patients with exclusion codes:", length(exclude_patients), "\n")

    # Get all patients and exclude those with unwanted codes
    all_patients <- unique(cov$person_id)
    control_patients <- setdiff(all_patients, exclude_patients)

    cat("Total patients in dataset:", length(all_patients), "\n")
    cat("Control patients (after exclusions):", length(control_patients), "\n\n")

    # Create control dataset
    controls_cov <- cov[person_id %in% control_patients]

    # # Combine unknown classes automatically
    # cat("=== Combining unknown classes ===\n")
    # before_unknown <- sum(controls_cov$Class %in% c("UNKNOWN0", "UNKNOWN1", "UNKNOWN2"))
    # controls_cov[Class %in% c("UNKNOWN0", "UNKNOWN1", "UNKNOWN2"), Class := "UNKNOWN"]
    # after_unknown <- sum(controls_cov$Class == "UNKNOWN")
    # cat("Combined", before_unknown, "unknown entries into", after_unknown, "UNKNOWN entries\n\n")

    # Apply filters
    controls_filtered <- apply_gender_filter(controls_cov, gender_filter)
    controls_filtered <- apply_crep_filter(controls_filtered, crep_filter)
    controls_final <- apply_age_filter(controls_filtered, age_filter)

    # ========================
    # SAVE RESULTS
    # ========================

    # File names
    all_controls_file <- paste0(output_prefix, "_all_controls.txt")
    final_controls_file <- paste0(output_prefix, "_final_controls.txt")
    all_ids_file <- paste0(output_prefix, "_all_control_ids.txt")
    final_ids_file <- paste0(output_prefix, "_final_control_ids.txt")

    # Save datasets
    fwrite(controls_cov, file.path(DATA_DIR, all_controls_file), sep = "\t")
    fwrite(controls_final, file.path(DATA_DIR, final_controls_file), sep = "\t")

    # Save ID lists
    writeLines(as.character(controls_cov$person_id), file.path(DATA_DIR, all_ids_file))
    writeLines(as.character(controls_final$person_id), file.path(DATA_DIR, final_ids_file))

    # Summary statistics
    cat("=== FINAL SUMMARY ===\n")
    cat("Total controls (all filters):", nrow(controls_cov), "\n")
    cat("Final controls (after all filters):", nrow(controls_final), "\n")

    if ("Sample_age" %in% colnames(controls_final)) {
        cat("\nAge statistics:\n")
        cat("- Mean age:", round(mean(controls_final$Sample_age, na.rm = TRUE), 1), "\n")
        cat("- Age range:", min(controls_final$Sample_age, na.rm = TRUE), "-",
            max(controls_final$Sample_age, na.rm = TRUE), "\n")
    }

    sink()

    # Return results
    results <- list(
        all_controls = controls_cov,
        final_controls = controls_final,
        summary = list(
            control_name = control_name,
            total_controls = nrow(controls_cov),
            final_controls = nrow(controls_final),
            filters_applied = list(
                exclude_codes = exclude_codes,
                gender = gender_filter,
                crep = crep_filter,
                age = age_filter
            )
        ),
        files_created = c(all_controls_file, final_controls_file, all_ids_file, final_ids_file)
    )

    cat("Control selection complete!\n")
    cat("Log saved to:", LOG_FILE, "\n")
    cat("Files created:\n")
    for(file in results$files_created) {
        cat("  -", file, "\n")
    }

    return(results)
}

