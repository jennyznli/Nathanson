#' Process ICD codes from Excel format to regex patterns
#'
#' Takes a dataframe with ICD9 and ICD10 columns containing comma-separated codes
#' and converts them to regex-compatible patterns with ^ prefix
#'
#' @param icd_df Dataframe with ICD9 and ICD10 columns
#' @return Named list with 'icd9' and 'icd10' vectors of regex patterns
process_icd_codes <- function(icd_df) {
    process_column <- function(col_vector) {
        col_vector <- col_vector[!is.na(col_vector)]
        all_codes <- paste(col_vector, collapse = ", ")
        # Split by comma
        code_list <- unlist(strsplit(all_codes, ","))
        # Trim whitespace
        code_list <- trimws(code_list)
        # Remove empty strings
        code_list <- code_list[code_list != ""]
        # Add ^ prefix for regex matching (start of string)
        code_list <- paste0("^", code_list)
        return(code_list)
    }
    result <- list()
    if ("ICD9" %in% colnames(icd_df)) {
        result$icd9 <- process_column(icd_df$ICD9)
    }
    if ("ICD10" %in% colnames(icd_df)) {
        result$icd10 <- process_column(icd_df$ICD10)
    }
    result$combined <- c(result$icd9, result$icd10)
    return(result)
}

merge_duplicates <- function(df, group_col) {
    df %>%
        group_by(across(all_of(group_col))) %>%
        summarise(
            across(-any_of(group_col), ~ paste(unique(as.character(.x)), collapse = ";")),
            .groups = "drop"
        )
}

# Function to get earliest age (your existing function)
get_earliest_age <- function(age_str) {
    if (is.na(age_str) || is.null(age_str) || age_str == "" || age_str == "NA") {
        return(NA)
    }
    age_str <- as.character(age_str)
    age_parts <- strsplit(age_str, ";")[[1]]
    ages <- suppressWarnings(as.numeric(age_parts))
    valid_ages <- ages[!is.na(ages)]
    if (length(valid_ages) == 0) {
        warning(paste("No valid ages found in:", age_str))
        return(NA)
    }
    if (any(valid_ages < 0)) {
        warning(paste("Negative age found in:", age_str))
    }
    return(min(valid_ages))
}

get_earliest_date <- function(date_str) {
    if (is.na(date_str) || is.null(date_str) || date_str == "" || date_str == "NA") {
        return(NA)
    }
    date_str <- as.character(date_str)
    date_parts <- strsplit(date_str, ";")[[1]]
    # Trim whitespace
    date_parts <- trimws(date_parts)
    # Parse dates in YYYY-MM-DD format
    dates <- suppressWarnings(as.Date(date_parts, format = "%Y-%m-%d"))
    valid_dates <- dates[!is.na(dates)]

    if (length(valid_dates) == 0) {
        warning(paste("No valid dates found in:", date_str))
        return(NA)
    }

    return(min(valid_dates))
}
