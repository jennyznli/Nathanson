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

select_samples <- function(
        sample_name,
        icd_codes,
        gender_filter = NULL,
        crep_filter = NULL,
        min_instances = 3,
        min_timespan = NULL,
        age_filter = NULL,
        exclude = TRUE,
        exclude_codes = c("^C(?!44)", "^Z85(?!\\.828)", "^(?!173)(?:1[4-9][0-9]|20[0-8]|209\\.[0-3])", "^V10(?!\\.83)"),
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

    occ <- fread(OCC, header = TRUE)
    cov <- fread(COV, header = TRUE)
    per <- fread(PER, header = TRUE)

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
    # FILTER HELPER FUNCTIONS
    # ========================

    # Function to apply gender filter
    apply_gender_filter <- function(patient_data, gender_filter) {
        if (is.null(gender_filter)) {
            cat("=== No gender filter applied === \n\n")
            return(patient_data)
        }

        if (!"Sequenced_gender" %in% colnames(patient_data)) {
            warning("Sequenced_gender column not found in covariates. Skipping gender filter.")
            return(patient_data)
        }

        cat("=== Applying gender filter:", gender_filter, "=== \n")
        cat("- Patients before gender filter:", nrow(patient_data), "\n")

        filtered_data <- patient_data[Sequenced_gender == gender_filter]
        cat("- Patients after gender filter:", nrow(filtered_data), "\n\n")

        return(filtered_data)
    }

    # Function to apply age filter
    apply_age_filter <- function(patient_data, age_filter) {
        if (is.null(age_filter)) {
            cat("=== No age filter applied === \n\n")
            return(patient_data)
        }

        if (!"Sample_age" %in% colnames(patient_data)) {
            warning("Sample_age column not found. Skipping age filter.")
            return(patient_data)
        }

        if (length(age_filter) == 1) {
            # Single number - minimum age
            cat("=== Applying minimum age filter: >=", age_filter, "=== \n")
            cat("- Patients before age filter:", nrow(patient_data), "\n")
            age_condition <- patient_data$Sample_age >= age_filter[1]
        } else {
            # Range - minimum and maximum age
            cat("=== Applying age range filter:", paste(age_filter, collapse = "-"), "=== \n")
            cat("- Patients before age filter:", nrow(patient_data), "\n")
            age_condition <- patient_data$Sample_age >= age_filter[1] &
                patient_data$Sample_age <= age_filter[2]
        }

        filtered_data <- patient_data[age_condition & !is.na(age_condition)]
        cat("- Patients after age filter:", nrow(filtered_data), "\n\n")

        return(filtered_data)
    }

    # Function to apply CREP filter
    apply_crep_filter <- function(patient_data, crep_filter = NULL) {
        if (is.null(crep_filter)) {
            cat("=== No CREP filter applied === \n\n")
            return(patient_data)
        }
        if (!"CREP_HighRisk_Flag" %in% colnames(patient_data)) {
            warning("CREP_HighRisk_Flag column not found. Skipping CREP filter.")
            return(patient_data)
        }

        cat("=== Applying CREP filter ===\n")
        cat("- Patients before CREP filter:", nrow(patient_data), "\n")

        if (crep_filter) {
            # keep only CREP high-risk patients
            filtered_data <- patient_data[CREP_HighRisk_Flag == 1]
            cat("- Patients after CREP filter (only high risk):", nrow(filtered_data), "\n\n")
        } else {
            # Exclude CREP high-risk patients
            filtered_data <- patient_data[CREP_HighRisk_Flag == 0]
            cat("- Patients after CREP filter (excluding high risk):", nrow(filtered_data), "\n\n")
        }
        return(filtered_data)
    }

    # Function to filter by minimum visits
    apply_visits_filter <- function(patient_data, min_visits = 3) {
        if (is.null(min_visits)) {
            cat("=== Skipping visit filter ===\n\n")
            return(patient_data)
        }

        cat("=== Applying visit filter === \n")
        cat("Minimum visits:", min_visits, "\n")

        cat("- Patients before visits filter:", nrow(patient_data), "\n")
        filtered_data <- patient_data[num_visits >= min_visits]

        cat("- Patients after visits filter:", nrow(filtered_data), "\n\n")
        return(filtered_data)
    }

    # Function to filter by minimum timespan between first and last diagnosis
    apply_timespan_filter <- function(patient_data, min_days = NULL) {
        if (is.null(min_days)) {
            cat("=== Skipping timespan filter === \n\n")
            return(patient_data)
        }

        cat("=== Applying timespan filter ===\n")
        cat("Minimum time span:", min_days, "days\n")

        cat("- Patients before timespan filter:", nrow(patient_data), "\n")
        filtered_data <- patient_data[days_span >= min_days]

        cat("- Patients after timespan filter:", nrow(filtered_data), "\n\n")
        return(filtered_data)
    }

    # ========================
    # HELPER FUNCTION
    # ========================
    # aggregates condition_occurrences table down to one individual per row
    create_patient_summary <- function(occurrences_data) {
        patient_summary <- occurrences_data[, {
            # All unique dates from both start and end dates
            all_dates <- unique(c(condition_start_date, condition_end_date))
            all_dates <- all_dates[!is.na(all_dates)]
            # All visit occurrences
            all_visits <- unique(visit_occurrence_id)
            all_visits <- all_visits[!is.na(all_visits)]
            list(
                num_visits = length(all_visits),
                num_dates = length(all_dates),
                first_date = min(all_dates, na.rm = TRUE),
                last_date = max(all_dates, na.rm = TRUE),
                total_records = .N,
                num_unique_codes = length(unique(condition_source_value)),
                all_matching_codes = paste(unique(condition_source_value), collapse = ";")
            )
        }, by = person_id]

        # Calculate time span in days and months
        patient_summary[, days_span := as.numeric(last_date - first_date)]
        patient_summary[, months_span := round(days_span / 30.44, 1)]

        return(patient_summary)
    }

    # ========================
    # MAIN FUNCTIONS
    # ========================
    # Function to filter and analyze patients by specific ICD pattern
    analyze_condition_patients <- function(condition_data, icd_pattern, pattern_name) {
        cat("=== Searching for patients w/ ICD codes ===\n")
        cat("Pattern:", icd_pattern, "\n\n")

        cat("1. Filtering condition occurrences table...\n")
        matching_occurrences <- condition_data[
            grepl(icd_pattern, condition_source_value, perl = TRUE, ignore.case = TRUE)
        ]

        cat("- Total occurrences in dataset:", nrow(condition_data), "\n")
        cat("- Occurrences matching codes: ", nrow(matching_occurrences), "\n")
        cat("- Patients with matching codes:", length(unique(matching_occurrences$person_id)), "\n\n")

        if(nrow(matching_occurrences) == 0) {
            cat("No matching occurrences found. Stopping analysis.\n")
            return(NULL)
        }

        cat("2. Creating patient summary from matching occurrences only...\n")
        patient_summary <- create_patient_summary(matching_occurrences)

        cat("3. Adding covariate data...\n\n")
        setkey(patient_summary, person_id)
        setkey(cov, person_id)
        patient_summary <- merge(patient_summary, cov, by = "person_id", all.x = TRUE)

        return(patient_summary)
    }

    # Function to exclude patients with other ICD codes
    exclude_patients <- function(patient_data, condition_data, target_codes, exclude_codes, min_exclude_visits, min_exclude_timespan) {
        if (!exclude) {
            cat("No other exclusion applied\n\n")
            return(patient_data)
        }

        cat("=== Finding patients with ICD codes to exclude ===\n")
        cat("Pattern:", exclude_codes, "\n\n")

        # Get all patients in our current dataset
        target_patients <- patient_data$person_id
        cat("Starting with", length(target_patients), "patients\n\n")

        # Create pattern for other codes (excluding our target codes)
        exclude_icd_pattern <- paste(exclude_codes, collapse = "|")
        target_pattern <- paste(target_codes, collapse = "|")

        # Find all occurrences of other codes
        exclude_occurrences <- condition_data[
            grepl(exclude_icd_pattern, condition_source_value, perl = TRUE, ignore.case = TRUE)
        ]
        cat("- Total other ICD code occurrences found:", nrow(exclude_occurrences), "\n")

        # Remove occurrences that match our target
        exclude_occurrences_filtered <- exclude_occurrences[
            !grepl(target_pattern, condition_source_value, perl = TRUE, ignore.case = TRUE)
        ]
        cat("- Other exclusion occurrences excluding target ICD codes:", nrow(exclude_occurrences_filtered), "\n\n")

        if (nrow(exclude_occurrences_filtered) == 0) {
            cat("No other exclude occurrences found. No patients to exclude.\n\n")
            return(patient_data)
        }

        # cat("5. Creating patient summary for exclusion occurrences...\n")
        exclude_summary <- create_patient_summary(exclude_occurrences_filtered, "other codes")

        ## CHECK THIS ??
        substantial_exclude <- apply_visits_filter(exclude_summary, min_visits = min_exclude_visits)
        substantial_exclude <- apply_timespan_filter(substantial_exclude, min_days = min_exclude_timespan)

        ## add covariates??

        # Save exclusion summaries
        fwrite(exclude_summary, "exclude_summary.txt", sep = "\t")
        fwrite(substantial_exclude, "substantial_exclude_summary.txt", sep = "\t")
        cat("Patients with substantial exclusion evidence (>=", min_exclude_visits, "visit occurrences):", nrow(substantial_exclude), "\n")

        # Get patients to exclude (those with substantial other ICD code evidence)
        patients_to_exclude <- substantial_exclude$person_id

        # Find overlap with our target patients
        overlap_patients <- intersect(target_patients, patients_to_exclude)
        cat("Target patients who also have substantial other ICD codes:", length(overlap_patients), "\n")

        # Exclude these patients
        clean_patients <- setdiff(target_patients, patients_to_exclude)
        cat("Patients remaining after substantial other ICD code exclusion:", length(clean_patients), "\n")
        cat("- Exclusion rate:", round((length(overlap_patients) / length(target_patients)) * 100, 1), "%\n\n")

        # Save removed patients info
        removed_patients_id <- setdiff(target_patients, clean_patients)
        writeLines(as.character(removed_patients_id), "excluded_patients_ids.txt")

        removed_patients <- exclude_summary[person_id %in% removed_patients_id]
        fwrite(removed_patients, "excluded_patients.txt", sep = "\t")

        # Filter the dataset
        filtered_data <- patient_data[person_id %in% clean_patients]

        return(filtered_data)
    }

    # ========================
    # MAIN ANALYSIS
    # ========================
    # Create combined pattern for all ICD codes
    combined_pattern <- paste(icd_codes, collapse = "|")

    # Run main analysis function
    all_patients <- analyze_condition_patients(occ, combined_pattern, paste(sample_name, "cancer"))
    results <- list()

    if(is.null(all_patients)) {
        cat("No patients found matching criteria. Stopping analysis.\n")
        sink()
        return(NULL)
    }

    # Apply ICD codes exclusion filter if it exists
    all_patients_clean <- exclude_patients(all_patients, occ,
                                           c, exclude_codes,
                                           min_exclude_visits = min_exclude_instances,
                                           min_exclude_timespan = min_exclude_timespan)

    # Apply other filters
    all_patients_filtered <- apply_visits_filter(all_patients_clean, min_visits = min_instances)
    all_patients_filtered <- apply_timespan_filter(all_patients_filtered, min_days = min_timespan)
    all_patients_filtered <- apply_crep_filter(all_patients_filtered, crep_filter)
    all_patients_filtered <- apply_gender_filter(all_patients_filtered, gender_filter)
    all_patients_final <- apply_age_filter(all_patients_filtered, age_filter)

    # ========================
    # SAVE RESULTS
    # ========================

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
    results$files_created <- c(all_patients_file, filtered_patients_file, all_patients_id_file, filtered_patients_id_file)

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




