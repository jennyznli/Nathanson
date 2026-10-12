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
    # group_by() lumps all NA group_col values into a single group, which would
    # merge unrelated rows together -- so rows with a missing id are left as-is
    # instead of being deduplicated
    # merged's non-group columns all come out as character (see across() below),
    # so na_rows needs the same coercion or bind_rows() errors on type mismatch
    na_rows <- df %>%
        filter(if_any(all_of(group_col), is.na)) %>%
        mutate(across(-any_of(group_col), as.character))
    id_rows <- df %>% filter(!if_any(all_of(group_col), is.na))

    merged <- id_rows %>%
        group_by(across(all_of(group_col))) %>%
        summarise(
            across(-any_of(group_col), ~ {
                vals <- unique(trimws(as.character(.x)))
                vals <- vals[!is.na(vals) & vals != ""]
                if (length(vals) == 0) NA_character_ else paste(vals, collapse = ";")
            }),
            .groups = "drop"
        )

    bind_rows(merged, na_rows)
}

#' Clean up an already semicolon-joined field: trim whitespace, drop blank
#' segments, and dedupe repeated values (e.g. "Category 1: Negative;" or
#' ";Yes" from a prior manual merge) without touching real multi-value data.
clean_semicolon_field <- function(x) {
    if (is.na(x)) return(NA_character_)
    parts <- trimws(strsplit(x, ";")[[1]])
    parts <- unique(parts[parts != ""])
    if (length(parts) == 0) NA_character_ else paste(parts, collapse = ";")
}

#' Rename columns to Capital_Snake_Case: periods -> underscores, collapse
#' repeated underscores, capitalize the first letter of each segment (the
#' rest of each segment's casing is left alone so acronyms like BRCA1 or
#' IHC survive). If two source columns tidy to the same name (e.g.
#' "Collection_Type" and "Collection.Type" both already existed), make.unique
#' suffixes the later one with "_1" rather than silently colliding -- that
#' suffix is a signal the two source columns need to be reconciled by hand.
tidy_colnames <- function(df) {
    nm <- names(df)
    nm <- gsub("\\.", "_", nm)
    nm <- gsub("_+", "_", nm)
    nm <- gsub("^_|_$", "", nm)
    nm <- sapply(strsplit(nm, "_"), function(parts) {
        parts <- ifelse(nchar(parts) > 0,
                         paste0(toupper(substring(parts, 1, 1)), substring(parts, 2)),
                         parts)
        paste(parts, collapse = "_")
    })
    names(df) <- make.unique(nm, sep = "_")
    df
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
