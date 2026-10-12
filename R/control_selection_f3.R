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
#       save_excluded_detail = TRUE         # Writes detailed summary of the excluded patients
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

select_controls_f3 <- function(
        control_name = "controls",
        exclude_codes = NULL,
        gender_filter = NULL,
        crep_filter = NULL,
        age_filter = NULL,
        data_dir = NULL,
        pmbb_dir = NULL,
        log_dir = NULL,
        output_prefix = NULL,
        save_excluded_details = TRUE
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

    OCC <- file.path(PMBB_DIR, "3.0", "PMBB-Release-2024-3.0_phenotype_condition_occurrence.txt")
    COV <- file.path(PMBB_DIR, "3.0", "PMBB-Release-2024-3.0_covariates.txt")
    PER <- file.path(PMBB_DIR, "3.0", "PMBB-Release-2024-3.0_phenotype_person.txt")

    files_to_check <- c(OCC, COV, PER)
    missing_files <- files_to_check[!file.exists(files_to_check)]

    if (length(missing_files) > 0) {
        stop("ERROR: Missing files:\n", paste(missing_files, collapse = "\n"))
    }

    cat("All required PMBB files found successfully.\n")

    # OPTIMIZATION 1: Load only needed columns from occurrence table
    cat("Loading data files...\n")
    occ <- fread(OCC,
                 select = c("person_id", "condition_source_value"),
                 header = TRUE)
    cov <- fread(COV, header = TRUE)
    per <- fread(PER, header = TRUE)
    cat("Data loaded successfully.\n\n")

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
    cat("Exclude Codes:", ifelse(is.null(exclude_codes) || length(exclude_codes) == 0,
                                 "None", paste(exclude_codes, collapse = ", ")), "\n")
    cat("Gender Filter:", ifelse(is.null(gender_filter), "None", gender_filter), "\n")
    cat("CREP Filter:", ifelse(is.null(crep_filter), "None", crep_filter), "\n")
    cat("Age Filter:", ifelse(is.null(age_filter), "None", paste(age_filter, collapse = "-")), "\n\n")

    # ========================
    # MAIN ANALYSIS - OPTIMIZED
    # ========================

    # OPTIMIZATION 2: Handle case where no exclusion codes provided
    if (is.null(exclude_codes) || length(exclude_codes) == 0) {
        cat("=== No exclusion codes provided ===\n")
        cat("Using all patients as potential controls\n\n")
        exclude_patients <- character(0)
        excluded_patient_details <- NULL
    } else {
        # Find patients to exclude based on ICD codes
        cat("=== Finding patients to exclude ===\n")
        exclude_pattern <- paste(exclude_codes, collapse = "|")

        # OPTIMIZATION 3: Use vectorized grepl then get unique person_ids
        has_exclude_code <- grepl(exclude_pattern, occ$condition_source_value,
                                  perl = TRUE, ignore.case = TRUE)

        # NEW: Get all exclusion occurrences with details
        exclude_occurrences <- occ[has_exclude_code]
        exclude_patients <- unique(exclude_occurrences$person_id)

        cat("Patients with exclusion codes:", length(exclude_patients), "\n")
        cat("Total exclusion occurrences:", nrow(exclude_occurrences), "\n\n")

        # NEW: Create summary of excluded patients
        if (save_excluded_details && nrow(exclude_occurrences) > 0) {
            cat("=== Analyzing excluded patients ===\n")

            # Summarize by patient
            excluded_patient_summary <- exclude_occurrences[, .(
                num_exclusion_codes = .N,
                unique_exclusion_codes = length(unique(condition_source_value)),
                all_exclusion_codes = paste(unique(condition_source_value), collapse = "; ")
            ), by = person_id]

            # Add covariates for excluded patients
            setkey(excluded_patient_summary, person_id)
            setkey(cov, person_id)
            excluded_patient_details <- cov[excluded_patient_summary]

            # Summarize ICD code frequencies
            excluded_icd_summary <- exclude_occurrences[, .(
                patient_count = length(unique(person_id)),
                occurrence_count = .N
            ), by = condition_source_value][order(-patient_count)]
        } else {
            excluded_patient_details <- NULL
            excluded_icd_summary <- NULL
        }
    }

    # Get all patients and exclude those with unwanted codes
    all_patients <- unique(cov$person_id)
    control_patients <- setdiff(all_patients, exclude_patients)

    cat("Total patients in dataset:", length(all_patients), "\n")
    cat("Control patients (after exclusions):", length(control_patients), "\n")
    cat("Exclusion rate:", round((length(exclude_patients) / length(all_patients)) * 100, 1), "%\n\n")

    # Get all patients and exclude those with unwanted codes
    all_patients <- unique(cov$person_id)
    control_patients <- setdiff(all_patients, exclude_patients)

    cat("Total patients in dataset:", length(all_patients), "\n")
    cat("Control patients (after exclusions):", length(control_patients), "\n\n")

    # OPTIMIZATION 4: Filter using data.table instead of creating intermediate dataset
    cat("=== Applying demographic filters ===\n")

    # Start with all controls
    controls_final <- cov[person_id %in% control_patients]
    cat("Starting controls:", nrow(controls_final), "\n")

    # OPTIMIZATION 5: Chain all filters in sequence (avoiding function call overhead)
    # Gender filter
    if (!is.null(gender_filter)) {
        if (!"Sequenced_gender" %in% colnames(controls_final)) {
            warning("Sequenced_gender column not found. Skipping gender filter.")
        } else {
            cat("- Applying gender filter (", gender_filter, ")\n")
            cat("  Before:", nrow(controls_final), "\n")
            controls_final <- controls_final[Sequenced_gender == gender_filter]
            cat("  After:", nrow(controls_final), "\n")
        }
    } else {
        cat("- No gender filter applied\n")
    }

    # CREP filter
    if (!is.null(crep_filter)) {
        if (!"CREP_HighRisk_Flag" %in% colnames(controls_final)) {
            warning("CREP_HighRisk_Flag column not found. Skipping CREP filter.")
        } else {
            cat("- Applying CREP filter (", ifelse(crep_filter, "CREP only", "exclude CREP"), ")\n")
            cat("  Before:", nrow(controls_final), "\n")
            controls_final <- controls_final[CREP_HighRisk_Flag == as.integer(crep_filter)]
            cat("  After:", nrow(controls_final), "\n")
        }
    } else {
        cat("- No CREP filter applied\n")
    }

    # Age filter
    if (!is.null(age_filter)) {
        if (!"Sample_age" %in% colnames(controls_final)) {
            warning("Sample_age column not found. Skipping age filter.")
        } else {
            if (length(age_filter) == 1) {
                cat("- Applying age filter (>=", age_filter[1], ")\n")
                cat("  Before:", nrow(controls_final), "\n")
                controls_final <- controls_final[Sample_age >= age_filter[1] & !is.na(Sample_age)]
                cat("  After:", nrow(controls_final), "\n")
            } else {
                cat("- Applying age filter (", age_filter[1], "-", age_filter[2], ")\n")
                cat("  Before:", nrow(controls_final), "\n")
                controls_final <- controls_final[
                    Sample_age >= age_filter[1] &
                        Sample_age <= age_filter[2] &
                        !is.na(Sample_age)
                ]
                cat("  After:", nrow(controls_final), "\n")
            }
        }
    } else {
        cat("- No age filter applied\n")
    }
    cat("\n")

    # OPTIMIZATION 6: Only create "all controls" dataset if different from final
    # (saves memory and disk I/O)
    controls_cov <- cov[person_id %in% control_patients]

    # ========================
    # SAVE RESULTS
    # ========================

    cat("=== Saving results ===\n")

    # File names
    all_controls_file <- paste0(output_prefix, "_all_controls.txt")
    final_controls_file <- paste0(output_prefix, "_final_controls.txt")
    all_ids_file <- paste0(output_prefix, "_all_control_ids.txt")
    final_ids_file <- paste0(output_prefix, "_final_control_ids.txt")

    # NEW: Files for excluded patients
    excluded_patients_file <- paste0(output_prefix, "_excluded_patients.txt")
    excluded_icd_file <- paste0(output_prefix, "_excluded_icd_summary.txt")

    # Save datasets
    fwrite(controls_cov, file.path(DATA_DIR, all_controls_file), sep = "\t")
    fwrite(controls_final, file.path(DATA_DIR, final_controls_file), sep = "\t")

    # Save ID lists
    writeLines(as.character(controls_cov$person_id), file.path(DATA_DIR, all_ids_file))
    writeLines(as.character(controls_final$person_id), file.path(DATA_DIR, final_ids_file))

    # NEW: Save excluded patient details
    files_created <- c(all_controls_file, final_controls_file, all_ids_file, final_ids_file)

    if (save_excluded_details && !is.null(excluded_patient_details)) {
        fwrite(excluded_patient_details, file.path(DATA_DIR, excluded_patients_file), sep = "\t")
        fwrite(excluded_icd_summary, file.path(DATA_DIR, excluded_icd_file), sep = "\t")
        files_created <- c(files_created, excluded_patients_file, excluded_icd_file)
        cat("Saved excluded patient details\n")
    }

    # Summary statistics
    cat("\n=== FINAL SUMMARY ===\n")
    cat("Total controls (no demographic filters):", nrow(controls_cov), "\n")
    cat("Final controls (after all filters):", nrow(controls_final), "\n")
    cat("Excluded patients:", length(exclude_patients), "\n\n")

    if ("Sequenced_gender" %in% colnames(controls_final)) {
        cat("\nGender distribution:\n")
        print(table(controls_final$Sequenced_gender))
    }

    if ("CREP_HighRisk_Flag" %in% colnames(controls_final)) {
        cat("\nCREP status:\n")
        cat("- CREP high-risk:", sum(controls_final$CREP_HighRisk_Flag == 1, na.rm = TRUE), "\n")
        cat("- Non-CREP:", sum(controls_final$CREP_HighRisk_Flag == 0, na.rm = TRUE), "\n")
    }

    if ("Sample_age" %in% colnames(controls_final)) {
        cat("\nAge statistics:\n")
        cat("- Mean age:", round(mean(controls_final$Sample_age, na.rm = TRUE), 1), "\n")
        cat("- Median age:", round(median(controls_final$Sample_age, na.rm = TRUE), 1), "\n")
        cat("- Age range:", min(controls_final$Sample_age, na.rm = TRUE), "-",
            max(controls_final$Sample_age, na.rm = TRUE), "\n")
    }

    if ("Batch" %in% colnames(controls_final)) {
        cat("\nBatch distribution:\n")
        print(table(controls_final$Batch))
    }

    sink()

    # Return results
    results <- list(
        all_controls = controls_cov,
        final_controls = controls_final,
        excluded_patients = excluded_patient_details,  # NEW
        excluded_icd_summary = excluded_icd_summary,   # NEW
        summary = list(
            control_name = control_name,
            total_controls = nrow(controls_cov),
            final_controls = nrow(controls_final),
            excluded_count = length(exclude_patients),  # NEW
            filters_applied = list(
                exclude_codes = exclude_codes,
                gender = gender_filter,
                crep = crep_filter,
                age = age_filter
            )
        ),
        files_created = files_created
    )

    cat("Control selection complete!\n")
    cat("Log saved to:", LOG_FILE, "\n")
    cat("Files created:\n")
    for(file in results$files_created) {
        cat("  -", file.path(DATA_DIR, file), "\n")
    }

    return(results)
}
