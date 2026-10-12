# ============================================================
# ICD Non-Cancer Control Selection Script
# PMBB Freeze 4.0
# ============================================================
# Purpose:
#   Select non-cancer control samples from PMBB phenotype files.
#
# Key PMBB 4.0 adaptation:
#   Exclusion codes may appear in either:
#       1. condition_occurrence
#       2. observation
#   This script combines both tables into one evidence table and
#   excludes patients with matching exclusion code patterns in either file.
#
#   A patient is excluded if they have exclusion-code evidence in at least
#   min_exclude_instances distinct visits (visit_occurrence_id).

# Usage:
#   select_controls(
#       control_name = "healthy_controls",
#       exclude_codes = cancer_codes,
#       gender_filter = NULL,
#       crep_filter = NULL,
#       age_filter = NULL,
#       min_exclude_instances = 1,
#       min_exclude_timespan = NULL,
#       data_dir = NULL,
#       pmbb_dir = NULL,
#       log_dir = NULL,
#       output_prefix = NULL,
#       save_excluded_details = TRUE
#   )
# ============================================================
select_controls_f4 <- function(
        control_name = "controls",
        exclude_codes = NULL,
        gender_filter = NULL,
        crep_filter = NULL,
        age_filter = NULL,
        min_exclude_instances = 1,
        min_exclude_timespan = NULL,
        data_dir = NULL,
        pmbb_dir = NULL,
        log_dir = NULL,
        output_prefix = NULL,
        save_excluded_details = TRUE
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
    options(datatable.verbose = FALSE)

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

    files_to_check <- list(
        "Condition Occurrence" = COND,
        "Observation" = OBS,
        "Covariates" = COV
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
            "visit_occurrence_id"
        ),
        header = TRUE
    )

    cov <- fread(COV, header = TRUE)

    cat("Data loaded successfully.\n\n")

    # ========================
    # HARMONIZE CONDITION + OBSERVATION
    # ========================

    cat("Combining condition and observation evidence...\n")

    cond_harmonized <- cond[, .(
        person_id = person_id,
        source_value = as.character(condition_source_value),
        event_date = as.Date(condition_start_date),
        visit_occurrence_id = visit_occurrence_id,
        source_table = "condition_occurrence"
    )]

    obs_harmonized <- obs[, .(
        person_id = person_id,
        source_value = as.character(observation_source_value),
        event_date = as.Date(observation_date),
        visit_occurrence_id = visit_occurrence_id,
        source_table = "observation"
    )]

    occ <- rbindlist(
        list(cond_harmonized, obs_harmonized),
        use.names = TRUE,
        fill = TRUE
    )

    occ <- occ[!is.na(source_value) & source_value != ""]

    cat("Combined evidence table created.\n\n")

    # ========================
    # CONFIGURATION
    # ========================

    TIMESTAMP <- format(Sys.time(), "%Y%m%d_%H%M%S")

    if (is.null(output_prefix)) {
        output_prefix <- tolower(gsub("[^A-Za-z0-9]", "_", control_name))
    }

    LOG_FILE <- file.path(
        LOG_DIR,
        paste0(TIMESTAMP, "_", output_prefix, "_control_selection_log.txt")
    )

    sink(LOG_FILE)

    on.exit({
        if (sink.number(type = "output") > 0) {
            sink()
        }
    }, add = TRUE)

    cat("Non-Cancer Control Selection Log\n")
    cat("================================\n\n")

    cat("Control Name:", control_name, "\n")
    cat(
        "Exclude Codes:",
        ifelse(
            is.null(exclude_codes) || length(exclude_codes) == 0,
            "None",
            paste(exclude_codes, collapse = ", ")
        ),
        "\n"
    )
    cat("Gender Filter:", ifelse(is.null(gender_filter), "None", gender_filter), "\n")
    cat("CREP Filter:", ifelse(is.null(crep_filter), "None", crep_filter), "\n")
    cat("Age Filter:", ifelse(is.null(age_filter), "None", paste(age_filter, collapse = "-")), "\n")
    cat("Minimum Exclusion Distinct Visits:", min_exclude_instances, "\n")
    cat("Minimum Exclusion Timespan:", ifelse(is.null(min_exclude_timespan), "None", min_exclude_timespan), "\n\n")

    cat("Combined evidence table:\n")
    cat("- Condition rows:", nrow(cond_harmonized), "\n")
    cat("- Observation rows:", nrow(obs_harmonized), "\n")
    cat("- Combined rows:", nrow(occ), "\n\n")

    # ========================
    # HELPER FUNCTION
    # ========================

    summarize_exclusion_patients <- function(exclusion_data) {
        patient_summary <- exclusion_data[, {
            all_dates <- unique(na.omit(event_date))
            all_visits <- unique(na.omit(visit_occurrence_id))

            list(
                num_exclusion_records = .N,
                num_exclusion_visits = length(all_visits),
                num_exclusion_dates = length(all_dates),
                first_exclusion_date = if (length(all_dates) > 0) min(all_dates) else as.Date(NA),
                last_exclusion_date = if (length(all_dates) > 0) max(all_dates) else as.Date(NA),
                unique_exclusion_codes = length(unique(source_value)),
                all_exclusion_codes = paste(sort(unique(source_value)), collapse = "; "),
                evidence_sources = paste(sort(unique(source_table)), collapse = "; ")
            )
        }, by = person_id]

        patient_summary[, `:=`(
            exclusion_days_span = as.numeric(last_exclusion_date - first_exclusion_date),
            exclusion_months_span = round(
                as.numeric(last_exclusion_date - first_exclusion_date) / 30.44,
                1
            )
        )]

        return(patient_summary)
    }

    # ========================
    # MAIN ANALYSIS
    # ========================
    #
    # Exclusion is computed ONCE here, independent of save_excluded_details.
    # The threshold uses DISTINCT VISITS (num_exclusion_visits) to match the
    # sample-selection script. save_excluded_details only controls whether the
    # detail tables are built and written; it never changes who is excluded.

    excluded_patient_details <- NULL
    excluded_code_summary <- NULL
    exclude_patients <- character(0)

    if (is.null(exclude_codes) || length(exclude_codes) == 0) {
        cat("=== No exclusion codes provided ===\n")
        cat("Using all patients as potential controls.\n\n")
    } else {
        cat("=== Finding patients to exclude ===\n")

        exclude_pattern <- paste(exclude_codes, collapse = "|")

        has_exclude_code <- grepl(
            exclude_pattern,
            occ$source_value,
            perl = TRUE,
            ignore.case = TRUE
        )

        exclude_occurrences <- occ[has_exclude_code]

        cat("Total exclusion records:", nrow(exclude_occurrences), "\n")

        if (nrow(exclude_occurrences) == 0) {
            cat("No exclusion records found.\n\n")
            exclude_patients <- character(0)
        } else {
            cat("\nExclusion records by source table:\n")
            print(exclude_occurrences[, .N, by = source_table][order(-N)])
            cat("\n")

            # --- Compute exclusion once (distinct-visit thresholding) ---
            cat("=== Summarizing exclusion evidence ===\n")

            excluded_patient_summary <- summarize_exclusion_patients(
                exclude_occurrences
            )

            substantial_excluded_patient_summary <- excluded_patient_summary

            if (!is.null(min_exclude_instances)) {
                substantial_excluded_patient_summary <-
                    substantial_excluded_patient_summary[
                        num_exclusion_visits >= min_exclude_instances
                    ]
            }

            if (!is.null(min_exclude_timespan)) {
                substantial_excluded_patient_summary <-
                    substantial_excluded_patient_summary[
                        exclusion_days_span >= min_exclude_timespan
                    ]
            }

            exclude_patients <- unique(
                substantial_excluded_patient_summary$person_id
            )

            cat(
                "Patients with exclusion evidence meeting threshold:",
                length(exclude_patients),
                "\n"
            )
            cat(
                "Minimum exclusion distinct visits required:",
                min_exclude_instances,
                "\n"
            )

            if (!is.null(min_exclude_timespan)) {
                cat(
                    "Minimum exclusion timespan required:",
                    min_exclude_timespan,
                    "days\n"
                )
            }

            cat("\n")

            # --- Detail tables: the ONLY thing gated by save_excluded_details ---
            if (save_excluded_details) {
                cat("=== Building excluded patient details ===\n")

                setkey(substantial_excluded_patient_summary, person_id)
                setkey(cov, person_id)

                excluded_patient_details <- cov[
                    substantial_excluded_patient_summary
                ]

                excluded_code_summary <- exclude_occurrences[
                    person_id %in% exclude_patients,
                    .(
                        patient_count = length(unique(person_id)),
                        occurrence_count = .N
                    ),
                    by = .(source_table, source_value)
                ][order(-patient_count, source_table, source_value)]

                cat("Excluded patient detail tables prepared.\n\n")
            }
        }
    }

    # ========================
    # DEFINE CONTROL POOL
    # ========================

    all_patients <- unique(cov$person_id)
    control_patients <- setdiff(all_patients, exclude_patients)

    cat("Total patients in covariate dataset:", length(all_patients), "\n")
    cat("Excluded patients:", length(exclude_patients), "\n")
    cat("Control patients after exclusions:", length(control_patients), "\n")

    if (length(all_patients) > 0) {
        cat(
            "Exclusion rate:",
            round((length(exclude_patients) / length(all_patients)) * 100, 1),
            "%\n\n"
        )
    } else {
        cat("Exclusion rate: NA\n\n")
    }

    # ========================
    # DEMOGRAPHIC FILTERS
    # ========================

    cat("=== Applying demographic filters ===\n")

    controls_cov <- cov[person_id %in% control_patients]
    controls_final <- controls_cov

    cat("Starting controls:", nrow(controls_final), "\n")

    # Gender filter
    if (!is.null(gender_filter)) {
        if (!"sequenced_gender" %in% colnames(controls_final)) {
            warning("sequenced_gender column not found. Skipping gender filter.")
        } else {
            cat("- Applying gender filter (", gender_filter, ")\n")
            cat("  Before:", nrow(controls_final), "\n")

            controls_final <- controls_final[
                sequenced_gender == gender_filter
            ]

            cat("  After:", nrow(controls_final), "\n")
        }
    } else {
        cat("- No gender filter applied\n")
    }

    # CREP filter
    if (!is.null(crep_filter)) {
        if (!"crep_highrisk_flag" %in% colnames(controls_final)) {
            warning("crep_highrisk_flag column not found. Skipping CREP filter.")
        } else {
            cat(
                "- Applying CREP filter (",
                ifelse(crep_filter, "CREP only", "exclude CREP"),
                ")\n"
            )
            cat("  Before:", nrow(controls_final), "\n")

            controls_final <- controls_final[
                crep_highrisk_flag == as.integer(crep_filter)
            ]

            cat("  After:", nrow(controls_final), "\n")
        }
    } else {
        cat("- No CREP filter applied\n")
    }

    # Age filter
    if (!is.null(age_filter)) {
        if (!"sample_age" %in% colnames(controls_final)) {
            warning("sample_age column not found. Skipping age filter.")
        } else {
            if (length(age_filter) == 1) {
                cat("- Applying age filter (>=", age_filter[1], ")\n")
                cat("  Before:", nrow(controls_final), "\n")

                controls_final <- controls_final[
                    sample_age >= age_filter[1] &
                        !is.na(sample_age)
                ]

                cat("  After:", nrow(controls_final), "\n")
            } else {
                cat(
                    "- Applying age filter (",
                    age_filter[1],
                    "-",
                    age_filter[2],
                    ")\n"
                )
                cat("  Before:", nrow(controls_final), "\n")

                controls_final <- controls_final[
                    sample_age >= age_filter[1] &
                        sample_age <= age_filter[2] &
                        !is.na(sample_age)
                ]

                cat("  After:", nrow(controls_final), "\n")
            }
        }
    } else {
        cat("- No age filter applied\n")
    }

    cat("\n")

    # ========================
    # SAVE RESULTS
    # ========================

    cat("=== Saving results ===\n")

    all_controls_file <- paste0(output_prefix, "_all_controls.txt")
    final_controls_file <- paste0(output_prefix, "_final_controls.txt")
    all_ids_file <- paste0(output_prefix, "_all_control_ids.txt")
    final_ids_file <- paste0(output_prefix, "_final_control_ids.txt")

    excluded_patients_file <- paste0(output_prefix, "_excluded_patients.txt")
    excluded_code_file <- paste0(output_prefix, "_excluded_icd_summary.txt")

    fwrite(
        controls_cov,
        file.path(DATA_DIR, all_controls_file),
        sep = "\t"
    )

    fwrite(
        controls_final,
        file.path(DATA_DIR, final_controls_file),
        sep = "\t"
    )

    writeLines(
        as.character(controls_cov$person_id),
        file.path(DATA_DIR, all_ids_file)
    )

    writeLines(
        as.character(controls_final$person_id),
        file.path(DATA_DIR, final_ids_file)
    )

    files_created <- c(
        all_controls_file,
        final_controls_file,
        all_ids_file,
        final_ids_file
    )

    if (
        save_excluded_details &&
        !is.null(excluded_patient_details) &&
        nrow(excluded_patient_details) > 0
    ) {
        fwrite(
            excluded_patient_details,
            file.path(DATA_DIR, excluded_patients_file),
            sep = "\t"
        )

        fwrite(
            excluded_code_summary,
            file.path(DATA_DIR, excluded_code_file),
            sep = "\t"
        )

        files_created <- c(
            files_created,
            excluded_patients_file,
            excluded_code_file
        )

        cat("Saved excluded patient details.\n")
    }

    # ========================
    # FINAL SUMMARY
    # ========================

    cat("\n=== FINAL SUMMARY ===\n")
    cat("Total controls before demographic filters:", nrow(controls_cov), "\n")
    cat("Final controls after all filters:", nrow(controls_final), "\n")
    cat("Excluded patients:", length(exclude_patients), "\n\n")

    if ("sequenced_gender" %in% colnames(controls_final)) {
        cat("\nGender distribution:\n")
        print(table(controls_final$sequenced_gender))
    }

    if ("crep_highrisk_flag" %in% colnames(controls_final)) {
        cat("\nCREP status:\n")
        cat(
            "- CREP high-risk:",
            sum(controls_final$crep_highrisk_flag == 1, na.rm = TRUE),
            "\n"
        )
        cat(
            "- Non-CREP:",
            sum(controls_final$crep_highrisk_flag == 0, na.rm = TRUE),
            "\n"
        )
    }

    if ("sample_age" %in% colnames(controls_final)) {
        cat("\nAge statistics:\n")
        cat(
            "- Mean age:",
            round(mean(controls_final$sample_age, na.rm = TRUE), 1),
            "\n"
        )
        cat(
            "- Median age:",
            round(median(controls_final$sample_age, na.rm = TRUE), 1),
            "\n"
        )
        cat(
            "- Age range:",
            min(controls_final$sample_age, na.rm = TRUE),
            "-",
            max(controls_final$sample_age, na.rm = TRUE),
            "\n"
        )
    }

    if ("Batch" %in% colnames(controls_final)) {
        cat("\nBatch distribution:\n")
        print(table(controls_final$Batch))
    }

    # ========================
    # RETURN RESULTS
    # ========================

    results <- list(
        all_controls = controls_cov,
        final_controls = controls_final,
        excluded_patients = excluded_patient_details,
        excluded_icd_summary = excluded_code_summary,
        summary = list(
            control_name = control_name,
            total_controls = nrow(controls_cov),
            final_controls = nrow(controls_final),
            excluded_count = length(exclude_patients),
            filters_applied = list(
                exclude_codes = exclude_codes,
                min_exclude_distinct_visits = min_exclude_instances,
                min_exclude_timespan = min_exclude_timespan,
                gender = gender_filter,
                crep = crep_filter,
                age = age_filter
            )
        ),
        files_created = files_created,
        log_file = LOG_FILE
    )

    sink()

    cat("Control selection complete.\n")
    cat("Log saved.\n")

    return(results)
}
