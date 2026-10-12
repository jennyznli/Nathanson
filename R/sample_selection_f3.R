# ============================================================
# Sample Selection Script
# ============================================================
# This script selects samples from PMBB phenotype files
# using ICD codes, filtering by visit frequency, time span
# between first and last occurrence, gender, and age.
#
# Usage:
#   select_samples(
#       sample_name = "breast",              # Name of the type as string
#       icd_codes = c("^C50", "^Z85.3"),     # ICD9/10 code patterns (regex-compatible), or NULL if taking all ICD codes
#       gender_filter = NULL,                # "Male", "Female", or NULL
#       crep_filter = TRUE,                  # TRUE - include only CREP samples, FALSE - exclude CREP samples, NULL - no filter for CREP
#       age_filter = NULL,                   # int - minimum age, vector - range of edges, or NULL for no filter
#       min_instances = 3,                   # Minimum instances of ICD codes
#       min_timespan = NULL,                 # Minimum timespan between first and last diagnosis occurrence (days)
#       exclude = FALSE,                     # Whether to exclude patients with certain ICD codes (e. g. those w/ history of other cancer)
#       exclude_codes = ...,                 # If excluding patients, ICD codes to exclude (like cancer codes)
#       min_exclude_instances = 2,           # If excluding patients, minimum instances for exclusion
#       min_exclude_timespan = NULL,         # If excluding patients, minimum timespan between first and last diagnosis occurrence (days)
#       data_dir = NULL,                     # Directory for outputs
#       pmbb_dir = NULL,                     # Directory for PMBB input w/ covariates, condition_occurrences, demographic tables
#       log_dir = NULL,                      # Directory for log files
#       output_prefix = NULL                 # Filename prefix (defaults to sample_name)
#   )
#
# Requires:
#   PMBB-Release-2024-3.0_phenotype_condition_occurrence.txt
#   PMBB-Release-2024-3.0_covariates.txt
#   PMBB-Release-2024-3.0_phenotype_person.txt
#
# Output:
#   - All matched patients and filtered patient lists
#   - Log file w/ selection steps and stats
# ============================================================

select_samples_f3 <- function(
        sample_name,
        icd_codes,
        gender_filter = NULL,
        crep_filter = NULL,
        min_instances = 2,
        min_timespan = NULL,
        age_filter = NULL,
        exclude = FALSE,
        exclude_codes = NULL,
        min_exclude_instances = 2,
        min_exclude_timespan = NULL,
        data_dir = NULL,
        pmbb_dir = NULL,
        log_dir = NULL,
        output_prefix = NULL
) {
    packages <- c("tidyr", "dplyr", "readr",
                  "stringr",  "data.table")
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
        PMBB_DIR <- file.path(DATA_DIR, "PMBB")
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

    setwd(DATA_DIR)

    # File paths
    OCC <- file.path(PMBB_DIR, "3.0", "PMBB-Release-2024-3.0_phenotype_condition_occurrence.txt")
    COV <- file.path(PMBB_DIR, "3.0", "PMBB-Release-2024-3.0_covariates.txt")
    PER <- file.path(PMBB_DIR, "3.0", "PMBB-Release-2024-3.0_phenotype_person.txt")

    # Check if these files exist
    files_to_check <- list(
        "Condition Occurrence" = OCC,
        "Covariates" = COV,
        "Person" = PER
    )

    missing_files <- c()
    for (file_type in names(files_to_check)) {
        file_path <- files_to_check[[file_type]]
        if (!file.exists(file_path)) {
            missing_files <- c(missing_files, paste0(file_type, ": ", file_path))
        }
    }

    if (length(missing_files) > 0) {
        stop("ERROR: The following required files are missing:\n",
             paste(missing_files, collapse = "\n"),
             "\n\nPlease ensure PMBB_DIR is set correctly and contains these files:\n",
             "- PMBB-Release-2024-3.0_phenotype_condition_occurrence.txt\n",
             "- PMBB-Release-2024-3.0_covariates.txt\n",
             "- PMBB-Release-2024-3.0_phenotype_person.txt")
    }

    cat("All required PMBB files found successfully.\n")

    # OPTIMIZATION 1: Load only needed columns
    cat("Loading data files...\n")
    occ <- fread(OCC,
                 select = c("person_id", "condition_source_value",
                            "condition_start_date", "condition_end_date",
                            "visit_occurrence_id"),
                 header = TRUE)
    cov <- fread(COV, header = TRUE)
    per <- fread(PER, header = TRUE)
    cat("Data loaded successfully.\n\n")

    # ========================
    # CONFIGURATION
    # ========================
    DATE <- format(Sys.Date(), "%Y%m%d")
    TIMESTAMP <- format(Sys.time(), "%Y%m%d_%H%M%S")

    # Set output prefix if not provided
    if (is.null(output_prefix)) {
        output_prefix <- tolower(gsub("[^A-Za-z0-9]", "_", sample_name))
    }

    # Start logging
    LOG_FILE <- file.path(LOG_DIR, paste0(TIMESTAMP, "_", output_prefix, "_case_selection_log.txt"))
    sink(LOG_FILE)

    cat("Select cases for", sample_name, "selection log\n")
    cat("========================================", "\n")

    cat("Selection Name:", sample_name, "\n")
    cat("ICD Codes:", paste(icd_codes, collapse = ", "), "\n")
    cat("Gender Filter:", gender_filter, "\n")
    cat("CREP Filter:", crep_filter, "\n")
    cat("Minimum Instances:", min_instances, "\n")
    cat("Minimum Timespan:", min_timespan, "\n")
    cat("Exclude Other Codes:", exclude, "\n")
    if (exclude) {
        cat("Other ICD Codes to Exclude:", paste(exclude_codes, collapse = ", "), "\n")
    }
    if (!is.null(age_filter)) {
        if (length(age_filter) == 1) {
            cat("Age Filter: >=", age_filter, "\n")
        } else {
            cat("Age Filter:", paste(age_filter, collapse = "-"), "\n")
        }
    }
    cat("\n")

    # ========================
    # OPTIMIZED HELPER FUNCTION
    # ========================
    # OPTIMIZATION 2: Optimized patient summary creation
    create_patient_summary <- function(occurrences_data) {
        patient_summary <- occurrences_data[, {
            # Use na.omit for faster NA removal
            all_dates <- unique(na.omit(c(condition_start_date, condition_end_date)))
            all_visits <- unique(na.omit(visit_occurrence_id))

            list(
                num_visits = length(all_visits),
                num_dates = length(all_dates),
                first_date = if(length(all_dates) > 0) min(all_dates) else as.Date(NA),
                last_date = if(length(all_dates) > 0) max(all_dates) else as.Date(NA),
                total_records = .N,
                num_unique_codes = length(unique(condition_source_value)),
                all_matching_codes = paste(unique(condition_source_value), collapse = ";")
            )
        }, by = person_id]

        # Calculate time span in days and months using := for in-place modification
        patient_summary[, `:=`(
            days_span = as.numeric(last_date - first_date),
            months_span = round(as.numeric(last_date - first_date) / 30.44, 1)
        )]

        return(patient_summary)
    }

    # ========================
    # OPTIMIZED EXCLUSION FUNCTION
    # ========================
    exclude_patients <- function(patient_data, condition_data, target_codes,
                                 exclude_codes, min_exclude_visits, min_exclude_timespan) {

        # SAFETY CHECK 1: Is exclude feature even enabled?
        if (!exclude) {
            cat("=== No exclusion applied ===\n\n")
            return(patient_data)
        }

        # SAFETY CHECK 2: Do we have codes to exclude?
        if (is.null(exclude_codes) || length(exclude_codes) == 0) {
            cat("=== No exclude_codes provided. Skipping exclusion. ===\n\n")
            return(patient_data)
        }

        # SAFETY CHECK 3: Do we have target codes?
        if (is.null(target_codes) || length(target_codes) == 0) {
            warning("No target_codes provided. Cannot determine what to exclude.")
            return(patient_data)
        }

        cat("=== Finding patients with ICD codes to exclude ===\n")
        cat("Exclude patterns:", paste(exclude_codes, collapse = ", "), "\n")
        cat("Target patterns:", paste(target_codes, collapse = ", "), "\n\n")

        # Get all patients in our current dataset
        target_patients <- patient_data$person_id
        cat("Starting with", length(target_patients), "patients\n\n")

        # OPTIMIZATION 3: Run both greps once, use boolean logic
        cat("1. Identifying exclusion occurrences...\n")
        exclude_pattern <- paste(exclude_codes, collapse = "|")
        target_pattern <- paste(target_codes, collapse = "|")

        is_exclude_match <- grepl(
            exclude_pattern,
            condition_data$condition_source_value,
            perl = TRUE,
            ignore.case = TRUE
        )

        is_target_match <- grepl(
            target_pattern,
            condition_data$condition_source_value,
            perl = TRUE,
            ignore.case = TRUE
        )

        # Keep: matches exclusion BUT NOT target
        exclude_occurrences <- condition_data[is_exclude_match & !is_target_match]

        cat("- Total exclusion occurrences found:", sum(is_exclude_match), "\n")
        cat("- Exclusion occurrences after removing target codes:", nrow(exclude_occurrences), "\n\n")

        if (nrow(exclude_occurrences) == 0) {
            cat("No exclusion occurrences found. No patients to exclude.\n\n")
            return(patient_data)
        }

        cat("2. Creating patient summary for exclusion occurrences...\n")
        exclude_summary <- create_patient_summary(exclude_occurrences)

        # OPTIMIZATION 4: Apply filters in one step with data.table
        cat("3. Filtering for substantial exclusion evidence...\n")
        substantial_exclude <- exclude_summary

        if (!is.null(min_exclude_visits)) {
            substantial_exclude <- substantial_exclude[num_visits >= min_exclude_visits]
        }
        if (!is.null(min_exclude_timespan)) {
            substantial_exclude <- substantial_exclude[days_span >= min_exclude_timespan]
        }

        # Save exclusion summaries
        fwrite(exclude_summary, "exclude_summary.txt", sep = "\t")
        fwrite(substantial_exclude, "substantial_exclude_summary.txt", sep = "\t")

        cat("Patients with substantial exclusion evidence (>=", min_exclude_visits, "visit occurrences):",
            nrow(substantial_exclude), "\n")

        # Get patients to exclude
        patients_to_exclude <- substantial_exclude$person_id

        # Find overlap with target patients
        overlap_patients <- intersect(target_patients, patients_to_exclude)
        cat("Target patients who also have substantial exclusion codes:", length(overlap_patients), "\n")

        # Exclude these patients
        clean_patients <- setdiff(target_patients, patients_to_exclude)
        cat("Patients remaining after exclusion:", length(clean_patients), "\n")
        cat("- Exclusion rate:", round((length(overlap_patients) / length(target_patients)) * 100, 1), "%\n\n")

        # Save removed patients info
        if (length(overlap_patients) > 0) {
            removed_patients_id <- overlap_patients
            writeLines(as.character(removed_patients_id), "excluded_patients_ids.txt")

            removed_patients <- exclude_summary[person_id %in% removed_patients_id]
            fwrite(removed_patients, "excluded_patients.txt", sep = "\t")
        }

        # Filter the dataset
        filtered_data <- patient_data[person_id %in% clean_patients]

        return(filtered_data)
    }

    # ========================
    # MAIN ANALYSIS - OPTIMIZED
    # ========================

    cat("=== Searching for patients with target ICD codes ===\n")

    # SAFETY CHECK: Validate ICD codes input
    if (is.null(icd_codes) || length(icd_codes) == 0) {
        stop("ERROR: No ICD codes provided. Please specify icd_codes parameter.")
    }

    # Create combined pattern for all ICD codes
    combined_pattern <- paste(icd_codes, collapse = "|")

    # OPTIMIZATION 5: Single grep operation to find matching occurrences
    cat("1. Filtering condition occurrences table...\n")
    matching_occurrences <- occ[
        grepl(combined_pattern, condition_source_value, perl = TRUE, ignore.case = TRUE)
    ]

    cat("- Total occurrences in dataset:", nrow(occ), "\n")
    cat("- Occurrences matching codes:", nrow(matching_occurrences), "\n")

    if(nrow(matching_occurrences) == 0) {
        cat("No matching occurrences found. Stopping analysis.\n")
        sink()
        return(NULL)
    }

    # Create patient summary
    cat("2. Creating patient summary from matching occurrences...\n")
    all_patients <- create_patient_summary(matching_occurrences)
    cat("- Total patients with matching codes:", nrow(all_patients), "\n\n")

    results <- list()

    # OPTIMIZATION 6: Filter BEFORE merging with covariates
    cat("3. Applying visit and timespan filters...\n")
    cat("- Patients before filters:", nrow(all_patients), "\n")

    all_patients_filtered <- all_patients

    if (!is.null(min_instances)) {
        all_patients_filtered <- all_patients_filtered[num_visits >= min_instances]
        cat("- Patients after visit filter (>= ", min_instances, " visits):", nrow(all_patients_filtered), "\n")
    }

    if (!is.null(min_timespan)) {
        all_patients_filtered <- all_patients_filtered[days_span >= min_timespan]
        cat("- Patients after timespan filter (>= ", min_timespan, " days):", nrow(all_patients_filtered), "\n")
    }
    cat("\n")

    # OPTIMIZATION 7: Merge only filtered patients (smaller dataset)
    cat("4. Adding covariate data...\n")
    setkey(all_patients_filtered, person_id)
    setkey(cov, person_id)
    all_patients_filtered <- cov[all_patients_filtered]  # Right join - keeps all filtered patients
    cat("- Patients after merging covariates:", nrow(all_patients_filtered), "\n\n")

    # Apply exclusion filter
    all_patients_clean <- exclude_patients(
        all_patients_filtered,
        occ,
        target_codes = icd_codes,
        exclude_codes = exclude_codes,
        min_exclude_visits = min_exclude_instances,
        min_exclude_timespan = min_exclude_timespan
    )

    # OPTIMIZATION 8: Apply remaining filters in optimized sequence
    cat("5. Applying demographic filters...\n")
    all_patients_final <- all_patients_clean

    # CREP filter
    if (!is.null(crep_filter)) {
        if (!"CREP_HighRisk_Flag" %in% colnames(all_patients_final)) {
            warning("CREP_HighRisk_Flag column not found. Skipping CREP filter.")
        } else {
            cat("- Patients before CREP filter:", nrow(all_patients_final), "\n")
            all_patients_final <- all_patients_final[CREP_HighRisk_Flag == as.integer(crep_filter)]
            cat("- Patients after CREP filter:", nrow(all_patients_final), "\n")
        }
    }

    # Gender filter
    if (!is.null(gender_filter)) {
        if (!"Sequenced_gender" %in% colnames(all_patients_final)) {
            warning("Sequenced_gender column not found. Skipping gender filter.")
        } else {
            cat("- Patients before gender filter:", nrow(all_patients_final), "\n")
            all_patients_final <- all_patients_final[Sequenced_gender == gender_filter]
            cat("- Patients after gender filter (", gender_filter, "):", nrow(all_patients_final), "\n")
        }
    }

    # Age filter
    if (!is.null(age_filter)) {
        if (!"Sample_age" %in% colnames(all_patients_final)) {
            warning("Sample_age column not found. Skipping age filter.")
        } else {
            cat("- Patients before age filter:", nrow(all_patients_final), "\n")

            if (length(age_filter) == 1) {
                # Minimum age
                all_patients_final <- all_patients_final[
                    Sample_age >= age_filter[1] & !is.na(Sample_age)
                ]
                cat("- Patients after age filter (>= ", age_filter[1], "):", nrow(all_patients_final), "\n")
            } else {
                # Age range
                all_patients_final <- all_patients_final[
                    Sample_age >= age_filter[1] &
                        Sample_age <= age_filter[2] &
                        !is.na(Sample_age)
                ]
                cat("- Patients after age filter (", paste(age_filter, collapse = "-"), "):",
                    nrow(all_patients_final), "\n")
            }
        }
    }
    cat("\n")

    # ========================
    # SAVE RESULTS
    # ========================

    cat("6. Saving results...\n")
    all_patients_file <- paste0(output_prefix, "_cancer_all_patients.txt")
    filtered_patients_file <- paste0(output_prefix, "_cancer_filtered_patients.txt")
    all_patients_id_file <- paste0(output_prefix, "_cancer_all_patients_ids.txt")
    filtered_patients_id_file <- paste0(output_prefix, "_cancer_filtered_patients_ids.txt")

    # Save full datasets
    fwrite(all_patients_clean, all_patients_file, sep = "\t")
    fwrite(all_patients_final, filtered_patients_file, sep = "\t")

    # Save ID lists
    writeLines(as.character(all_patients_clean$person_id), all_patients_id_file)
    writeLines(as.character(all_patients_final$person_id), filtered_patients_id_file)

    # Store results
    results$all_patients <- all_patients_clean
    results$filtered_patients <- all_patients_final
    results$files_created <- c(all_patients_file, filtered_patients_file,
                               all_patients_id_file, filtered_patients_id_file)

    results$summary <- list(
        cancer_type = sample_name,
        icd_codes = icd_codes,
        total_patients = nrow(all_patients_clean),
        filtered_patients = nrow(all_patients_final),
        filters_applied = list(
            gender = gender_filter,
            age = age_filter,
            crep = crep_filter,
            min_instances = min_instances,
            min_timespan = min_timespan,
            exclude = exclude,
            exclude_codes = if(exclude) exclude_codes else NULL
        )
    )

    sink()

    cat("Analysis complete\n")
    cat("Log saved to:", LOG_FILE, "\n")

    if (!is.null(results$files_created)) {
        cat("Files created:\n")
        for(file in results$files_created) {
            cat("  -", file, "\n")
        }
    }

    return(results)
}
