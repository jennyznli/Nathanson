# ============================================================
# Sample Selection Script - Observation Occurrences Version
# ============================================================
# This script selects samples from PMBB phenotype observation files
# using observation_source_value patterns, filtering by visit frequency,
# time span between first and last observation date, gender, CREP status,
# and age.
#
# Main difference from condition_occurrence version:
#   - observation table has observation_date / observation_datetime
#   - no condition_start_date / condition_end_date
#   - patient-level dates are based only on observation_date
#
# Usage:
#   select_samples(
#       sample_name = "breast",
#       icd_codes = c("^C50", "^Z85.3"),     # Regex-compatible observation_source_value patterns
#       gender_filter = NULL,                # "Male", "Female", or NULL
#       crep_filter = TRUE,                  # TRUE = include only CREP samples, FALSE = exclude CREP samples, NULL = no CREP filter
#       age_filter = NULL,                   # int = minimum age, vector = age range, or NULL
#       min_instances = 3,                   # Minimum number of distinct visits
#       min_timespan = NULL,                 # Minimum days between first and last observation_date
#       exclude = FALSE,
#       exclude_codes = NULL,
#       min_exclude_instances = 2,           # Minimum distinct visits for exclusion evidence
#       min_exclude_timespan = NULL,
#       data_dir = NULL,
#       pmbb_dir = NULL,
#       log_dir = NULL,
#       output_prefix = NULL
#   )
#
# Output:
#   - All matched patients and filtered patient lists
#   - Log file with selection steps and stats
# ============================================================

select_samples_obs <- function(
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
    packages <- c(
        "tidyr",
        "dplyr",
        "readr",
        "stringr",
        "data.table",
        "purrr"
    )

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
    OBS <- file.path(PMBB_DIR, "4.0", "PMBB-Release-2026-4.0_phenotype_observation.txt")
    COV <- file.path(PMBB_DIR, "4.0", "PMBB-Release-2026-4.0_phenotype_covariates.txt")
    PER <- file.path(PMBB_DIR, "4.0", "PMBB-Release-2026-4.0_phenotype_person.txt")

    # Check if these files exist
    files_to_check <- list(
        "Observation" = OBS,
        "Covariates" = COV,
        "Person" = PER
    )

    missing_files <- c()

    for (file_type in names(files_to_check)) {
        file_path <- files_to_check[[file_type]]

        if (!file.exists(file_path)) {
            missing_files <- c(
                missing_files,
                paste0(file_type, ": ", file_path)
            )
        }
    }

    if (length(missing_files) > 0) {
        stop(
            "ERROR: The following required files are missing:\n",
            paste(missing_files, collapse = "\n"),
            "\n\nPlease ensure PMBB_DIR is set correctly."
        )
    }

    cat("All required PMBB files found successfully.\n")

    # ========================
    # LOAD DATA
    # ========================

    cat("Loading data files...\n")

    occ <- fread(
        OBS,
        select = c(
            "person_id",
            "observation_source_value",
            "observation_date",
            "observation_datetime",
            "visit_occurrence_id"
        ),
        header = TRUE
    )

    cov <- fread(COV, header = TRUE)
    per <- fread(PER, header = TRUE)

    # Convert observation_date for date math
    occ[, observation_date := as.Date(observation_date)]

    cat("Data loaded successfully.\n\n")

    # ========================
    # CONFIGURATION
    # ========================

    DATE <- format(Sys.Date(), "%Y%m%d")
    TIMESTAMP <- format(Sys.time(), "%Y%m%d_%H%M%S")

    if (is.null(output_prefix)) {
        output_prefix <- tolower(gsub("[^A-Za-z0-9]", "_", sample_name))
    }

    LOG_FILE <- file.path(
        LOG_DIR,
        paste0(TIMESTAMP, "_", output_prefix, "_observation_selection_log.txt")
    )

    sink(LOG_FILE)

    cat("Select cases for", sample_name, "selection log\n")
    cat("========================================\n")

    cat("Selection Name:", sample_name, "\n")
    cat("Observation Code Patterns:", paste(icd_codes, collapse = ", "), "\n")
    cat("Gender Filter:", gender_filter, "\n")
    cat("CREP Filter:", crep_filter, "\n")
    cat("Minimum Distinct Visits:", min_instances, "\n")
    cat("Minimum Timespan:", min_timespan, "\n")
    cat("Exclude Other Codes:", exclude, "\n")

    if (exclude) {
        cat("Other Observation Codes to Exclude:", paste(exclude_codes, collapse = ", "), "\n")
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
    # HELPER FUNCTION
    # ========================

    create_patient_summary <- function(occurrences_data) {
        patient_summary <- occurrences_data[, {
            all_dates <- unique(na.omit(observation_date))
            all_visits <- unique(na.omit(visit_occurrence_id))

            list(
                num_visits = length(all_visits),
                num_dates = length(all_dates),
                first_date = if (length(all_dates) > 0) min(all_dates) else as.Date(NA),
                last_date = if (length(all_dates) > 0) max(all_dates) else as.Date(NA),
                total_records = .N,
                num_unique_codes = length(unique(observation_source_value)),
                all_matching_codes = paste(unique(observation_source_value), collapse = ";")
            )
        }, by = person_id]

        patient_summary[, `:=`(
            days_span = as.numeric(last_date - first_date),
            months_span = round(as.numeric(last_date - first_date) / 30.44, 1)
        )]

        return(patient_summary)
    }

    # ========================
    # EXCLUSION FUNCTION
    # ========================

    exclude_patients <- function(
        patient_data,
        observation_data,
        target_codes,
        exclude_codes,
        min_exclude_visits,
        min_exclude_timespan
    ) {

        if (!exclude) {
            cat("=== No exclusion applied ===\n\n")
            return(patient_data)
        }

        if (is.null(exclude_codes) || length(exclude_codes) == 0) {
            cat("=== No exclude_codes provided. Skipping exclusion. ===\n\n")
            return(patient_data)
        }

        if (is.null(target_codes) || length(target_codes) == 0) {
            warning("No target_codes provided. Cannot determine what to exclude.")
            return(patient_data)
        }

        cat("=== Finding patients with observation codes to exclude ===\n")
        cat("Exclude patterns:", paste(exclude_codes, collapse = ", "), "\n")
        cat("Target patterns:", paste(target_codes, collapse = ", "), "\n\n")

        target_patients <- patient_data$person_id

        cat("Starting with", length(target_patients), "patients\n\n")

        cat("1. Identifying exclusion observations...\n")

        exclude_pattern <- paste(exclude_codes, collapse = "|")
        target_pattern <- paste(target_codes, collapse = "|")

        is_exclude_match <- grepl(
            exclude_pattern,
            observation_data$observation_source_value,
            perl = TRUE,
            ignore.case = TRUE
        )

        is_target_match <- grepl(
            target_pattern,
            observation_data$observation_source_value,
            perl = TRUE,
            ignore.case = TRUE
        )

        exclude_occurrences <- observation_data[
            is_exclude_match & !is_target_match
        ]

        cat("- Total exclusion observations found:", sum(is_exclude_match), "\n")
        cat("- Exclusion observations after removing target codes:", nrow(exclude_occurrences), "\n\n")

        if (nrow(exclude_occurrences) == 0) {
            cat("No exclusion observations found. No patients to exclude.\n\n")
            return(patient_data)
        }

        cat("2. Creating patient summary for exclusion observations...\n")
        exclude_summary <- create_patient_summary(exclude_occurrences)

        cat("3. Filtering for substantial exclusion evidence...\n")

        substantial_exclude <- exclude_summary

        if (!is.null(min_exclude_visits)) {
            substantial_exclude <- substantial_exclude[
                num_visits >= min_exclude_visits
            ]
        }

        if (!is.null(min_exclude_timespan)) {
            substantial_exclude <- substantial_exclude[
                days_span >= min_exclude_timespan
            ]
        }

        fwrite(exclude_summary, "exclude_summary.txt", sep = "\t")
        fwrite(substantial_exclude, "substantial_exclude_summary.txt", sep = "\t")

        cat(
            "Patients with substantial exclusion evidence (>=",
            min_exclude_visits,
            "distinct visits):",
            nrow(substantial_exclude),
            "\n"
        )

        patients_to_exclude <- substantial_exclude$person_id

        overlap_patients <- intersect(target_patients, patients_to_exclude)

        cat(
            "Target patients who also have substantial exclusion codes:",
            length(overlap_patients),
            "\n"
        )

        clean_patients <- setdiff(target_patients, patients_to_exclude)

        cat("Patients remaining after exclusion:", length(clean_patients), "\n")

        if (length(target_patients) > 0) {
            cat(
                "- Exclusion rate:",
                round((length(overlap_patients) / length(target_patients)) * 100, 1),
                "%\n\n"
            )
        } else {
            cat("- Exclusion rate: NA\n\n")
        }

        if (length(overlap_patients) > 0) {
            removed_patients_id <- overlap_patients

            writeLines(
                as.character(removed_patients_id),
                "excluded_patients_ids.txt"
            )

            removed_patients <- exclude_summary[
                person_id %in% removed_patients_id
            ]

            fwrite(removed_patients, "excluded_patients.txt", sep = "\t")
        }

        filtered_data <- patient_data[
            person_id %in% clean_patients
        ]

        return(filtered_data)
    }

    # ========================
    # MAIN ANALYSIS
    # ========================

    cat("=== Searching for patients with target observation codes ===\n")

    if (is.null(icd_codes) || length(icd_codes) == 0) {
        cat("ERROR: No observation code patterns provided.\n")
        sink()
        stop("ERROR: No observation code patterns provided. Please specify icd_codes parameter.")
    }

    combined_pattern <- paste(icd_codes, collapse = "|")

    cat("1. Filtering observation table...\n")

    matching_occurrences <- occ[
        grepl(
            combined_pattern,
            observation_source_value,
            perl = TRUE,
            ignore.case = TRUE
        )
    ]

    cat("- Total observations in dataset:", nrow(occ), "\n")
    cat("- Observations matching codes:", nrow(matching_occurrences), "\n")

    if (nrow(matching_occurrences) == 0) {
        cat("No matching observations found. Stopping analysis.\n")
        sink()
        return(NULL)
    }

    cat("2. Creating patient summary from matching observations...\n")

    all_patients <- create_patient_summary(matching_occurrences)

    cat("- Total patients with matching codes:", nrow(all_patients), "\n\n")

    results <- list()

    # ========================
    # VISIT AND TIMESPAN FILTERS
    # ========================

    cat("3. Applying visit and timespan filters...\n")
    cat("- Patients before filters:", nrow(all_patients), "\n")

    all_patients_filtered <- all_patients

    # Option A: keep old behavior.
    # min_instances means minimum number of distinct visits,
    # not minimum number of observation rows.
    if (!is.null(min_instances)) {
        all_patients_filtered <- all_patients_filtered[
            num_visits >= min_instances
        ]

        cat(
            "- Patients after visit filter (>= ",
            min_instances,
            " distinct visits):",
            nrow(all_patients_filtered),
            "\n"
        )
    }

    if (!is.null(min_timespan)) {
        all_patients_filtered <- all_patients_filtered[
            days_span >= min_timespan
        ]

        cat(
            "- Patients after timespan filter (>= ",
            min_timespan,
            " days):",
            nrow(all_patients_filtered),
            "\n"
        )
    }

    cat("\n")

    # ========================
    # MERGE COVARIATES
    # ========================

    cat("4. Adding covariate data...\n")

    setkey(all_patients_filtered, person_id)
    setkey(cov, person_id)

    all_patients_filtered <- cov[all_patients_filtered]

    cat(
        "- Patients after merging covariates:",
        nrow(all_patients_filtered),
        "\n\n"
    )

    # ========================
    # EXCLUSION FILTER
    # ========================

    all_patients_clean <- exclude_patients(
        patient_data = all_patients_filtered,
        observation_data = occ,
        target_codes = icd_codes,
        exclude_codes = exclude_codes,
        min_exclude_visits = min_exclude_instances,
        min_exclude_timespan = min_exclude_timespan
    )

    # ========================
    # DEMOGRAPHIC FILTERS
    # ========================

    cat("5. Applying demographic filters...\n")

    all_patients_final <- all_patients_clean

    # CREP filter
    if (!is.null(crep_filter)) {
        if (!"CREP_HighRisk_Flag" %in% colnames(all_patients_final)) {
            warning("CREP_HighRisk_Flag column not found. Skipping CREP filter.")
        } else {
            cat("- Patients before CREP filter:", nrow(all_patients_final), "\n")

            all_patients_final <- all_patients_final[
                CREP_HighRisk_Flag == as.integer(crep_filter)
            ]

            cat("- Patients after CREP filter:", nrow(all_patients_final), "\n")
        }
    }

    # Gender filter
    if (!is.null(gender_filter)) {
        if (!"sequenced_gender" %in% colnames(all_patients_final)) {
            warning("Sequenced_gender column not found. Skipping gender filter.")
        } else {
            cat("- Patients before gender filter:", nrow(all_patients_final), "\n")

            all_patients_final <- all_patients_final[
                sequenced_gender == gender_filter
            ]

            cat(
                "- Patients after gender filter (",
                gender_filter,
                "):",
                nrow(all_patients_final),
                "\n"
            )
        }
    }

    # Age filter
    if (!is.null(age_filter)) {
        if (!"Sample_age" %in% colnames(all_patients_final)) {
            warning("Sample_age column not found. Skipping age filter.")
        } else {
            cat("- Patients before age filter:", nrow(all_patients_final), "\n")

            if (length(age_filter) == 1) {
                all_patients_final <- all_patients_final[
                    Sample_age >= age_filter[1] &
                        !is.na(Sample_age)
                ]

                cat(
                    "- Patients after age filter (>= ",
                    age_filter[1],
                    "):",
                    nrow(all_patients_final),
                    "\n"
                )
            } else {
                all_patients_final <- all_patients_final[
                    Sample_age >= age_filter[1] &
                        Sample_age <= age_filter[2] &
                        !is.na(Sample_age)
                ]

                cat(
                    "- Patients after age filter (",
                    paste(age_filter, collapse = "-"),
                    "):",
                    nrow(all_patients_final),
                    "\n"
                )
            }
        }
    }

    cat("\n")

    # ========================
    # SAVE RESULTS
    # ========================

    cat("6. Saving results...\n")

    all_patients_file <- paste0(output_prefix, "_observation_all_patients.txt")
    filtered_patients_file <- paste0(output_prefix, "_observation_filtered_patients.txt")
    all_patients_id_file <- paste0(output_prefix, "_observation_all_patients_ids.txt")
    filtered_patients_id_file <- paste0(output_prefix, "_observation_filtered_patients_ids.txt")

    fwrite(all_patients_clean, all_patients_file, sep = "\t")
    fwrite(all_patients_final, filtered_patients_file, sep = "\t")

    writeLines(
        as.character(all_patients_clean$person_id),
        all_patients_id_file
    )

    writeLines(
        as.character(all_patients_final$person_id),
        filtered_patients_id_file
    )

    results$all_patients <- all_patients_clean
    results$filtered_patients <- all_patients_final

    results$files_created <- c(
        all_patients_file,
        filtered_patients_file,
        all_patients_id_file,
        filtered_patients_id_file
    )

    results$summary <- list(
        sample_name = sample_name,
        observation_code_patterns = icd_codes,
        total_patients = nrow(all_patients_clean),
        filtered_patients = nrow(all_patients_final),
        filters_applied = list(
            gender = gender_filter,
            age = age_filter,
            crep = crep_filter,
            min_distinct_visits = min_instances,
            min_timespan = min_timespan,
            exclude = exclude,
            exclude_codes = if (exclude) exclude_codes else NULL,
            min_exclude_distinct_visits = if (exclude) min_exclude_instances else NULL,
            min_exclude_timespan = if (exclude) min_exclude_timespan else NULL
        )
    )

    sink()

    cat("Analysis complete.\n")
    cat("Log saved to:", LOG_FILE, "\n")

    if (!is.null(results$files_created)) {
        cat("Files created:\n")

        for (file in results$files_created) {
            cat("  -", file, "\n")
        }
    }

    return(results)
}
