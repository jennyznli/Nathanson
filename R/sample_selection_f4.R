# ============================================================
# ICD Sample Selection Script
# PMBB Freeze 4.0
# ============================================================
# Purpose:
#   Select samples from PMBB 4.0 phenotype files using ICD/source code
#   patterns that may appear in either:
#       1. condition_occurrence
#       2. observation
#
#   Some ICD/source codes previously found only in condition_occurrence
#   may now appear in observation. This script integrates both tables
#   into one evidence table and applies the selection logic across both.
#
#   min_instances = 3 means the patient must have the target code pattern
#   in at least 3 distinct visit_occurrence_id values across condition +
#   observation combined.

# Usage:
#   select_samples(
#       sample_name = "breast",
#       icd_codes = c("^C50", "^Z85.3"),
#       gender_filter = "Female",
#       crep_filter = NULL,
#       age_filter = NULL,
#       min_instances = 3,
#       min_timespan = NULL,
#       exclude = TRUE,
#       exclude_codes = cancer_codes,
#       min_exclude_instances = 2,
#       min_exclude_timespan = NULL,
#       data_dir = NULL,
#       pmbb_dir = NULL,
#       log_dir = NULL,
#       output_prefix = NULL
#   )
#
# Inputs:
#   - PMBB-Release-2026-4.0_phenotype_condition_occurrence.txt
#   - PMBB-Release-2026-4.0_phenotype_observation.txt
#   - PMBB-Release-2026-4.0_phenotype_covariates.txt
#   - PMBB-Release-2026-4.0_phenotype_person.txt
# ============================================================

select_samples_f4 <- function(
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
    options(datatable.showProgress = FALSE)

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

    # Allows either:
    #   pmbb_dir = ".../PMBB"
    # or:
    #   pmbb_dir = ".../PMBB/4.0"
    if (dir.exists(file.path(PMBB_DIR, "4.0"))) {
        PMBB_4_DIR <- file.path(PMBB_DIR, "4.0")
    } else {
        PMBB_4_DIR <- PMBB_DIR
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
    COND <- file.path(
        PMBB_4_DIR,
        "PMBB-Release-2026-4.0_phenotype_condition_occurrence.txt"
    )

    OBS <- file.path(
        PMBB_4_DIR,
        "PMBB-Release-2026-4.0_phenotype_observation.txt"
    )

    COV <- file.path(
        PMBB_4_DIR,
        "PMBB-Release-2026-4.0_phenotype_covariates.txt"
    )

    PER <- file.path(
        PMBB_4_DIR,
        "PMBB-Release-2026-4.0_phenotype_person.txt"
    )

    # Check if files exist
    files_to_check <- list(
        "Condition Occurrence" = COND,
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

    cond <- fread(
        COND,
        select = c(
            "person_id",
            "condition_source_value",
            "condition_start_date",
            "condition_end_date",
            "visit_occurrence_id"
        ),
        header = TRUE
    )

    obs <- fread(
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

    cat("Data loaded successfully.\n\n")

    # ========================
    # HARMONIZE CONDITION + OBSERVATION
    # ========================

    cat("Combining condition and observation evidence...\n")

    cond_harmonized <- cond[, .(
        person_id = person_id,
        source_value = as.character(condition_source_value),
        event_date = as.Date(condition_start_date),
        end_date = as.Date(condition_end_date),
        event_datetime = as.POSIXct(NA),
        visit_occurrence_id = visit_occurrence_id,
        source_table = "condition_occurrence"
    )]

    obs_harmonized <- obs[, .(
        person_id = person_id,
        source_value = as.character(observation_source_value),
        event_date = as.Date(observation_date),
        end_date = as.Date(NA),
        event_datetime = as.POSIXct(observation_datetime),
        visit_occurrence_id = visit_occurrence_id,
        source_table = "observation"
    )]

    occ <- rbindlist(
        list(cond_harmonized, obs_harmonized),
        use.names = TRUE,
        fill = TRUE
    )

    occ <- occ[!is.na(source_value) & source_value != ""]

    cat("Combined evidence table created.\n")
    cat("- Condition rows:", nrow(cond_harmonized), "\n")
    cat("- Observation rows:", nrow(obs_harmonized), "\n")
    cat("- Combined rows:", nrow(occ), "\n\n")

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
        paste0(TIMESTAMP, "_", output_prefix, "_case_selection_log.txt")
    )

    sink(LOG_FILE)

    on.exit({
        if (sink.number(type = "output") > 0) {
            sink()
        }
    }, add = TRUE)

    cat("Case selection for ", sample_name, "\n")
    cat("========================================\n\n")

    cat("Selection Name:", sample_name, "\n")
    cat("Code Patterns:", paste(icd_codes, collapse = ", "), "\n")
    cat("Gender Filter:", gender_filter, "\n")
    cat("CREP Filter:", crep_filter, "\n")
    cat("Minimum Distinct Visits:", min_instances, "\n")
    cat("Minimum Timespan:", min_timespan, "\n")
    cat("Exclude Other Codes:", exclude, "\n")

    if (exclude) {
        cat("Other Codes to Exclude:", paste(exclude_codes, collapse = ", "), "\n")
        cat("Minimum Exclusion Distinct Visits:", min_exclude_instances, "\n")
        cat("Minimum Exclusion Timespan:", min_exclude_timespan, "\n")
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
            all_dates <- unique(na.omit(event_date))
            all_visits <- unique(na.omit(visit_occurrence_id))

            list(
                num_visits = length(all_visits),
                num_dates = length(all_dates),
                first_date = if (length(all_dates) > 0) min(all_dates) else as.Date(NA),
                last_date = if (length(all_dates) > 0) max(all_dates) else as.Date(NA),
                total_records = .N,
                num_unique_codes = length(unique(source_value)),
                all_matching_codes = paste(sort(unique(source_value)), collapse = ";"),
                evidence_sources = paste(sort(unique(source_table)), collapse = ";")
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
        occurrence_data,
        target_codes,
        exclude_codes,
        min_exclude_visits,
        min_exclude_timespan
    ) {
        # SAFETY CHECK 1: Is exclude feature enabled?
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

        cat("=== Finding patients with codes to exclude ===\n")
        cat("Exclude patterns:", paste(exclude_codes, collapse = ", "), "\n")
        cat("Target patterns:", paste(target_codes, collapse = ", "), "\n\n")

        target_patients <- patient_data$person_id

        cat("Starting with", length(target_patients), "patients\n\n")

        cat("1. Identifying exclusion occurrences...\n")

        exclude_pattern <- paste(exclude_codes, collapse = "|")
        target_pattern <- paste(target_codes, collapse = "|")

        is_exclude_match <- grepl(
            exclude_pattern,
            occurrence_data$source_value,
            perl = TRUE,
            ignore.case = TRUE
        )

        is_target_match <- grepl(
            target_pattern,
            occurrence_data$source_value,
            perl = TRUE,
            ignore.case = TRUE
        )

        # Keep exclusion codes but do not count target codes as exclusions
        exclude_occurrences <- occurrence_data[
            is_exclude_match & !is_target_match
        ]

        cat("- Total exclusion occurrences found:", sum(is_exclude_match), "\n")
        cat("- Exclusion occurrences after removing target codes:", nrow(exclude_occurrences), "\n\n")

        if (nrow(exclude_occurrences) == 0) {
            cat("No exclusion occurrences found. No patients to exclude.\n\n")
            return(patient_data)
        }

        cat("2. Creating patient summary for exclusion occurrences...\n")
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
    cat("=== Searching for patients with target codes ===\n")

    # SAFETY CHECK: Validate code input
    if (is.null(icd_codes) || length(icd_codes) == 0) {
        cat("ERROR: No code patterns provided.\n")
        sink()
        stop("ERROR: No code patterns provided. Please specify icd_codes parameter.")
    }

    combined_pattern <- paste(icd_codes, collapse = "|")

    cat("1. Filtering combined condition + observation evidence table...\n")
    is_target_match <- grepl(
        combined_pattern,
        occ$source_value,
        perl = TRUE,
        ignore.case = TRUE
    )

    matching_occurrences <- occ[is_target_match]

    cat("- Total combined occurrences in dataset:", nrow(occ), "\n")
    cat("- Occurrences matching target codes:", nrow(matching_occurrences), "\n")

    if (nrow(matching_occurrences) == 0) {
        cat("No matching occurrences found. Stopping analysis.\n")
        sink()
        return(NULL)
    }

    cat("\nBreakdown of matching occurrences by source table:\n")
    print(matching_occurrences[, .N, by = source_table][order(-N)])
    cat("\n")

    cat("2. Creating patient summary from matching occurrences...\n")
    all_patients <- create_patient_summary(matching_occurrences)

    cat("- Total patients with matching codes:", nrow(all_patients), "\n\n")
    results <- list()

    # ========================
    # VISIT AND TIMESPAN FILTERS
    # ========================
    cat("3. Applying visit and timespan filters...\n")
    cat("- Patients before filters:", nrow(all_patients), "\n")

    all_patients_filtered <- all_patients

    # min_instances means minimum number of distinct visits,
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

    # timespan filter
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

    # Right join style from original script:
    # keeps all filtered patients and adds covariate columns.
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
        occurrence_data = occ,
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
        if (!"crep_highrisk_flag" %in% colnames(all_patients_final)) {
            warning("crep_highrisk_flag column not found. Skipping CREP filter.")
        } else {
            cat("- Patients before CREP filter:", nrow(all_patients_final), "\n")

            all_patients_final <- all_patients_final[
                crep_highrisk_flag == as.integer(crep_filter)
            ]

            cat("- Patients after CREP filter:", nrow(all_patients_final), "\n")
        }
    } else {
        cat("- No CREP filter applied\n")
    }

    # Gender filter
    if (!is.null(gender_filter)) {
        if (!"sequenced_gender" %in% colnames(all_patients_final)) {
            warning("sequenced_gender column not found. Skipping gender filter.")
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
    } else {
        cat("- No gender filter applied\n")
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
    } else {
        cat("- No age filter applied\n")
    }

    cat("\n")

    # ========================
    # SAVE RESULTS
    # ========================
    cat("6. Saving results...\n")

    all_patients_file <- paste0(output_prefix, "_all_patients.txt")
    filtered_patients_file <- paste0(output_prefix, "_filtered_patients.txt")
    all_patients_id_file <- paste0(output_prefix, "_all_patients_ids.txt")
    filtered_patients_id_file <- paste0(output_prefix, "_filtered_patients_ids.txt")
    matching_occurrences_file <- paste0(output_prefix, "_matching_occurrences.txt")

    fwrite(
        all_patients_clean,
        file.path(DATA_DIR, all_patients_file),
        sep = "\t"
    )

    fwrite(
        all_patients_final,
        file.path(DATA_DIR, filtered_patients_file),
        sep = "\t"
    )

    fwrite(
        matching_occurrences,
        file.path(DATA_DIR, matching_occurrences_file),
        sep = "\t"
    )

    writeLines(
        as.character(all_patients_clean$person_id),
        file.path(DATA_DIR, all_patients_id_file)
    )

    writeLines(
        as.character(all_patients_final$person_id),
        file.path(DATA_DIR, filtered_patients_id_file)
    )

    # ========================
    # FINAL SUMMARY
    # ========================

    cat("\n=== FINAL SUMMARY ===\n")
    cat("Selection name:", sample_name, "\n")
    cat("Patients after visit/timespan/exclusion filters:", nrow(all_patients_clean), "\n")
    cat("Final patients after demographic filters:", nrow(all_patients_final), "\n")

    if ("sequenced_gender" %in% colnames(all_patients_final)) {
        cat("\nGender distribution:\n")
        print(table(all_patients_final$sequenced_gender))
    }

    if ("crep_highrisk_flag" %in% colnames(all_patients_final)) {
        cat("\nCREP status:\n")
        cat("- CREP high-risk:", sum(all_patients_final$crep_highrisk_flag == 1, na.rm = TRUE), "\n")
        cat("- Non-CREP:", sum(all_patients_final$crep_highrisk_flag == 0, na.rm = TRUE), "\n")
    }

    if ("Sample_age" %in% colnames(all_patients_final)) {
        cat("\nAge statistics:\n")
        cat("- Mean age:", round(mean(all_patients_final$Sample_age, na.rm = TRUE), 1), "\n")
        cat("- Median age:", round(median(all_patients_final$Sample_age, na.rm = TRUE), 1), "\n")
        cat(
            "- Age range:",
            min(all_patients_final$Sample_age, na.rm = TRUE),
            "-",
            max(all_patients_final$Sample_age, na.rm = TRUE),
            "\n"
        )
    }

    if ("Batch" %in% colnames(all_patients_final)) {
        cat("\nBatch distribution:\n")
        print(table(all_patients_final$Batch))
    }

    files_created <- c(
        all_patients_file,
        filtered_patients_file,
        all_patients_id_file,
        filtered_patients_id_file,
        matching_occurrences_file
    )

    results$all_patients <- all_patients_clean
    results$filtered_patients <- all_patients_final
    results$matching_occurrences <- matching_occurrences
    results$person <- per

    results$files_created <- files_created

    results$summary <- list(
        sample_name = sample_name,
        code_patterns = icd_codes,
        total_matching_occurrences = nrow(matching_occurrences),
        total_patients_with_matching_codes = nrow(all_patients),
        patients_after_visit_timespan_filters = nrow(all_patients_filtered),
        patients_after_exclusion = nrow(all_patients_clean),
        final_patients = nrow(all_patients_final),
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

    results$log_file <- LOG_FILE

    sink()

    cat("Analysis complete.\n")
    cat("Log saved to:", LOG_FILE, "\n")

    if (!is.null(results$files_created)) {
        cat("Files created:\n")

        for (file in results$files_created) {
            cat("  -", file.path(DATA_DIR, file), "\n")
        }
    }

    return(results)
}


