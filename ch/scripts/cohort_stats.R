# ============================================================
# COVARIATE BALANCE - UNMATCHED + MATCHED
# ============================================================
library(here)
setwd("/Users/jennyzli/Documents/Nathanson")
source(here("R", "config.R"))
library(writexl)

# ============================================================
# READ DATA
# ============================================================
cov_unmatched <- read.csv(file.path("ch", "data", "pmbb_brca12_cov_df.csv"), row.names = 1) %>%
    filter(Strata %in% c(1, 2), !is.na(Sample_age), !is.na(Smoke_History))
dim(cov_unmatched)
# 28822
table(cov_unmatched$BRCA12_Case)
# 0     1
# 28125   697

cov_matched <- read_excel(file.path("ch", "data", "pmbb_brca12_cov_chip_df.xlsx"))
dim(cov_matched)
# 3004

# ============================================================
# HELPERS
# ============================================================
format_p <- function(p) {
    if (is.na(p)) return(NA_character_)
    if (p == 0)   return("< 2.2e-16")
    if (p < 0.001) return(formatC(p, format = "e", digits = 2))
    as.character(round(p, 3))
}

test_continuous <- function(var, data) {
    test <- wilcox.test(data[[var]] ~ data$BRCA12_Case)

    stats <- data %>%
        group_by(BRCA12_Case) %>%
        summarise(
            median = median(.data[[var]], na.rm = TRUE),
            q1     = quantile(.data[[var]], 0.25, na.rm = TRUE),
            q3     = quantile(.data[[var]], 0.75, na.rm = TRUE),
            .groups = "drop"
        ) %>%
        mutate(summary = paste0(median, " [", q1, "–", q3, "]"))

    data.frame(
        Variable  = var,
        Type      = "Continuous",
        Level     = NA_character_,
        Cases     = stats$summary[stats$BRCA12_Case == 1],
        Controls  = stats$summary[stats$BRCA12_Case == 0],
        Method    = "Wilcoxon",
        P_value = format_p(test$p.value),
        Sig       = ifelse(test$p.value < 0.001, "***",
                           ifelse(test$p.value < 0.01,  "**",
                                  ifelse(test$p.value < 0.05,  "*", "ns")))
    )
}

test_categorical <- function(var, data) {
    df_sub <- data[!is.na(data[[var]]), ]
    tab    <- table(df_sub[[var]], df_sub$BRCA12_Case)

    # counts and percentages per level
    counts <- df_sub %>%
        group_by(BRCA12_Case, .data[[var]]) %>%
        summarise(n = n(), .groups = "drop") %>%
        group_by(BRCA12_Case) %>%
        mutate(pct   = round(100 * n / sum(n), 1),
               label = paste0(n, " (", pct, "%)")) %>%
        ungroup() %>%
        select(BRCA12_Case, level = all_of(var), label) %>%
        pivot_wider(names_from = BRCA12_Case, values_from = label,
                    names_prefix = "group_") %>%
        rename(Cases = group_1, Controls = group_0)

    # test
    expected <- suppressWarnings(chisq.test(tab)$expected)
    if (any(expected < 5) || nrow(tab) > 2) {
        test   <- fisher.test(tab, simulate.p.value = TRUE)
        method <- "Fisher"
    } else {
        test   <- chisq.test(tab)
        method <- "Chi-squared"
    }

    p   <- test$p.value
    sig <- ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "ns")))

    # header row with p-value, level rows with counts
    header <- data.frame(Variable = var, Type = "Categorical", Level = NA_character_,
                         Cases = "", Controls = "",
                         Method = method,
                         P_value = format_p(p), Sig = sig)

    level_rows <- data.frame(
        Variable  = var,
        Type      = "",
        Level     = as.character(counts$level),  # <- add this
        Cases     = counts$Cases,
        Controls  = counts$Controls,
        Method    = "",
        P_value   = NA_character_,  # <- change from NA_real_
        Sig       = ""
    )

    bind_rows(header, level_rows)
}


# ============================================================
# RUN FOR BOTH COHORTS
# ============================================================
continuous_vars  <- c("Sample_age")
categorical_vars <- c("Sequenced_gender", "Smoke_History", "Batch", "Strata", "Class")

run_table <- function(data) {
    cont_results <- lapply(continuous_vars,  test_continuous,  data = data)
    cat_results  <- lapply(categorical_vars, test_categorical, data = data)
    bind_rows(c(cont_results, cat_results))
}

table_unmatched <- run_table(cov_unmatched)
table_matched   <- run_table(cov_matched)

write_xlsx(
    list(
        Unmatched = table_unmatched,
        Matched   = table_matched
    ),
    file.path("ch", "data", "covariate_balance_tables.xlsx")
)




