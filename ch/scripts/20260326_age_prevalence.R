# ========================
# PREVALENCE + AGE ANALYSIS
# ========================
library(here)
setwd("/Users/jennyzli/Documents/Nathanson")
source(here("R", "config.R"))

# ========================
# READ IN DATA
# ========================
all_ids <- read.csv(file.path("ch", "data", "ch_psm_matched4_case_control_ids.csv"))$x

### CHIP variants
vars <- read_excel(file.path("ch", "data", "ch_seq_wl_art_minad4_vars.xlsx"))
dim(vars)
length(unique(vars$Sample.ID))

### covariates
cov_all <- read.csv(file.path("ch", "data", "pmbb_brca12_cov_df.csv"), row.names = 1) %>% filter(person_id %in% all_ids)

cov_all <- cov_all %>%
    mutate(
        CHIP_Binary = person_id %in% vars$Sample.ID,
        CHIP_Count  = sapply(person_id, function(id) sum(vars$Sample.ID == id))
    )

# split into various sections..
cov_f <- cov_all %>% filter(Sequenced_gender == "Female")
cov_m <- cov_all %>% filter(Sequenced_gender == "Male")
cov_s1 <- cov_all %>% filter(Strata == 1)
cov_s1f <- cov_all %>% filter(Strata == 1, Sequenced_gender == "Female")
cov_s1m <- cov_all %>% filter(Strata == 1, Sequenced_gender == "Male")

# ========================
# AGE DISTRIBUTION HISTOGRAM
# ========================
plot_age_hist <- function(df, title_str, fill_var = NULL, binwidth = 2) {

    if (!is.null(fill_var)) {
        p <- ggplot(df, aes(x = Sample_age, fill = as.factor(.data[[fill_var]]))) +
            geom_histogram(binwidth = binwidth, color = "white", linewidth = 0.3,
                           alpha = 0.75, position = "identity") +
            scale_fill_manual(
                values = c("0" = "gray", "1" = "#C2185B",
                           "Female" = "purple", "Male" = "orange"),
                name = fill_var
            )
    } else {
        p <- ggplot(df, aes(x = Sample_age)) +
            geom_histogram(binwidth = binwidth, fill = "#546E7A",
                           color = "white", linewidth = 0.3, alpha = 0.85)
    }

    p <- p +
        geom_density(
            data        = df,
            mapping     = aes(x = Sample_age, y = after_stat(count) * binwidth),
            color       = "darkgray",
            linewidth   = 0.4,
            inherit.aes = FALSE
        ) +
        scale_x_continuous(breaks = seq(20, 90, by = 10), limits = c(18, 92)) +
        scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
        labs(
            title    = title_str,
            subtitle = sprintf("N = %d | median age = %.1f | IQR = %.1f\u2013%.1f",
                               nrow(df),
                               median(df$Sample_age, na.rm = TRUE),
                               quantile(df$Sample_age, 0.25, na.rm = TRUE),
                               quantile(df$Sample_age, 0.75, na.rm = TRUE)),
            x = "Age at Sample Collection (years)",
            y = "Count"
        ) +
        theme_minimal(base_size = 12) +
        theme(
            plot.title         = element_text(face = "bold", hjust = 0.5),
            plot.subtitle      = element_text(hjust = 0.5, color = "grey40", size = 10),
            panel.grid.minor   = element_blank(),
            panel.grid.major.x = element_blank(),
            legend.position    = "bottom"
        )

    p
}

# --- Run and save ---
datasets <- list(
    list(df = cov_all,  tag = "all",  lab = "Full"),
    list(df = cov_f,  tag = "female",  lab = "Full, Female"),
    list(df = cov_m,  tag = "male",  lab = "Full, Male"),
    list(df = cov_s1,  tag = "s1",  lab = "Strata 1"),
    list(df = cov_s1f,  tag = "s1_female",  lab = "Strata 1, Female"),
    list(df = cov_s1m,  tag = "s1_male",  lab = "Strata 1, Male")
)

for (d in datasets) {

    # Overall
    p <- plot_age_hist(d$df, paste("Age Distribution —", d$lab))
    ggsave(file.path("ch", "figures", paste0("age_hist_", d$tag, ".pdf")),
           p, width = 7, height = 4.5, dpi = 300)

    # Stratified by carrier status
    p_carrier <- plot_age_hist(d$df, paste("Age Distribution by Carrier Status —", d$lab),
                               fill_var = "BRCA12_Case")
    ggsave(file.path("ch", "figures", paste0("age_hist_by_carrier_", d$tag, ".pdf")),
           p_carrier, width = 7, height = 4.5, dpi = 300)

    # Stratified by sex
    p_sex <- plot_age_hist(d$df, paste("Age Distribution by Sex —", d$lab),
                           fill_var = "Sequenced_gender")
    ggsave(file.path("ch", "figures", paste0("age_hist_by_sex_", d$tag, ".pdf")),
           p_sex, width = 7, height = 4.5, dpi = 300)
}


# ========================
# QUICK STATS
# ========================
sink(file.path("ch", "data", "ch_summary_stats.log"), split = TRUE)
cat("Run date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")

# --- helper: simple count + % table for any variable ---
basic_crosstab <- function(df, var) {
    tab  <- table(df[[var]], useNA = "ifany")
    prop <- prop.table(tab) * 100
    out  <- data.frame(
        Level = names(tab),
        N     = as.integer(tab),
        Pct   = sprintf("%.1f%%", prop)
    )
    cat("\n", var, ":\n")
    print(out, row.names = FALSE)
    invisible(out)
}

chip_crosstab <- function(df, strat_var) {
    tab_n    <- table(df[[strat_var]], df$CHIP_Binary)
    tab_prop <- prop.table(tab_n, margin = 1) * 100
    combined <- matrix(
        sprintf("%d (%.1f%%)", tab_n, tab_prop),
        nrow = nrow(tab_n),
        dimnames = list(rownames(tab_n), c("No CHIP", "CHIP"))
    )
    total_n    <- colSums(tab_n)
    total_prop <- total_n / sum(total_n) * 100
    total_row  <- matrix(
        sprintf("%d (%.1f%%)", total_n, total_prop),
        nrow = 1,
        dimnames = list("Total", c("No CHIP", "CHIP"))
    )
    out <- rbind(combined, total_row)
    cat("\nCHIP by", strat_var, ":\n")
    print(out, quote = FALSE, right = FALSE)
    invisible(out)
}

for (tag in c("all", "s1", "s1f")) {
    cov <- get(paste0("cov_", tag))
    cat("\n", strrep("=", 40), "\n")
    cat(" Dataset:", tag, "| N =", nrow(cov), "\n")
    cat(strrep("=", 40), "\n")

    # --- basic demographic counts ---
    cat("\n--- Basic counts ---\n")
    for (sv in c("Sequenced_gender", "Carrier", "Batch")) {
        basic_crosstab(cov, sv)
    }

    # --- CHIP prevalence ---
    cat("\nOverall CHIP prevalence:\n")
    ov_n    <- table(cov$CHIP_Binary)
    ov_prop <- prop.table(ov_n) * 100
    cat(sprintf("  No CHIP : %d (%.1f%%)\n", ov_n["FALSE"], ov_prop["FALSE"]))
    cat(sprintf("  CHIP    : %d (%.1f%%)\n", ov_n["TRUE"],  ov_prop["TRUE"]))

    for (sv in c("BRCA12_Case", "Carrier", "Sequenced_gender", "Batch")) {
        chip_crosstab(cov, sv)
    }
}
sink()

# ========================
# HELPER: AGE-DECADE PREVALENCE PLOT
# ========================
plot_prev_by_decade <- function(df, title_str, strat_var = NULL, max_age = Inf,
                                strat_colors = NULL) {

    breaks <- if (is.finite(max_age)) c(-Inf, 30, 40, 50, 60, max_age)
    else                    c(-Inf, 30, 40, 50, 60, 70, Inf)
    labels <- if (is.finite(max_age)) c("≤30", "30-40", "40-50", "50-60", paste0("60-", max_age))
    else                    c("≤30", "30-40", "40-50", "50-60", "60-70", "≥70")

    # --- overall summary (always computed) ---
    overall <- df %>%
        mutate(age_group = cut(Sample_age, breaks = breaks, labels = labels, right = TRUE)) %>%
        group_by(age_group) %>%
        summarise(n_total  = n(),
                  n_chip   = sum(CHIP_Binary),
                  prev_pct = 100 * mean(CHIP_Binary),
                  .groups  = "drop") %>%
        mutate(x_num = as.numeric(age_group))

    p <- ggplot() +
        scale_x_continuous(breaks = overall$x_num, labels = overall$age_group) +
        scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.1))) +
        labs(title = title_str, x = "Age group (years)", y = "CHIP prevalence (%)") +
        theme_minimal(base_size = 13) +
        theme(panel.grid.minor   = element_blank(),
              panel.grid.major.x = element_blank(),
              plot.title         = element_text(face = "bold", hjust = 0.5),
              legend.position    = "bottom")

    # --- stratified layer (optional) ---
    if (!is.null(strat_var)) {
        strat <- df %>%
            mutate(age_group = cut(Sample_age, breaks = breaks, labels = labels, right = TRUE),
                   strat_fac = as.factor(.data[[strat_var]])) %>%
            group_by(age_group, strat_fac) %>%
            summarise(n_total  = n(),
                      n_chip   = sum(CHIP_Binary),
                      prev_pct = 100 * mean(CHIP_Binary),
                      .groups  = "drop") %>%
            mutate(x_num = as.numeric(age_group))

        # default colors if not supplied
        if (is.null(strat_colors)) {
            n_levels     <- nlevels(strat$strat_fac)
            strat_colors <- setNames(
                scales::hue_pal()(n_levels),
                levels(strat$strat_fac)
            )
        }

        p <- p +
            geom_smooth(data = strat,
                        aes(x = x_num, y = prev_pct,
                            color = strat_fac, group = strat_fac),
                        method = "loess", se = FALSE, linewidth = 0.8, linetype = "solid") +
            geom_point(data = strat,
                       aes(x = x_num, y = prev_pct, color = strat_fac),
                       size = 2.5) +
            scale_color_manual(values = strat_colors, name = strat_var)
    }

    # --- overall layer on top (always, in dark grey/black, labeled) ---
    p <- p +
        geom_smooth(data = overall,
                    aes(x = x_num, y = prev_pct, group = 1),
                    method = "loess", se = FALSE,
                    color = "darkgray", linewidth = 0.8) +
        geom_point(data = overall,
                   aes(x = x_num, y = prev_pct),
                   color = "darkgray", size = 3.5, shape = 18) +   # diamond shape distinguishes overall
        annotate("text",
                 x     = max(overall$x_num),
                 y     = overall$prev_pct[which.max(overall$x_num)] + 1.5,
                 label = "Overall",
                 color = "darkgray", size = 3.5, hjust = 1, fontface = "italic")

    p
}

# ========================
# CHIP PREVALENCE BY DECADE — ALL / NO-BRCA12 / BRCA12
# ========================
for (tag in c("all", "s1", "s1f")) {
    cov     <- get(paste0("cov_", tag))
    tag_lab <- if (tag == "all") "All" else "Strata 1"
    sex_lab <- if (tag == "s1f") "Female" else "All sexes"

    # overall only
    p <- plot_prev_by_decade(cov, sprintf("CHIP prevalence by decade (%s | %s)", tag_lab, sex_lab))
    ggsave(file.path("ch", "figures", paste0("chip_prevalence_by_decade_", tag, ".pdf")),
           p, width = 7, height = 5)

    # stratified by carrier, with overall on top
    p <- plot_prev_by_decade(
        cov,
        sprintf("CHIP prevalence by decade — carriers vs controls (%s | %s)", tag_lab, sex_lab),
        strat_var    = "BRCA12_Case",
        strat_colors = c("0" = "#546E7A", "1" = "#C2185B")
    )
    ggsave(file.path("ch", "figures", paste0("chip_prevalence_by_decade_by_case_", tag, ".pdf")),
           p, width = 7, height = 5)

    # stratified by sex — skip for s1f since it's already female-only
    if (tag != "s1f") {
        p <- plot_prev_by_decade(
            cov,
            sprintf("CHIP prevalence by decade — by sex (%s | %s)", tag_lab, sex_lab),
            strat_var    = "Sequenced_gender",
            strat_colors = c("Female" = "#7B2D8B", "Male" = "#0277BD")
        )
        ggsave(file.path("ch", "figures", paste0("chip_prevalence_by_decade_by_sex_", tag, ".pdf")),
               p, width = 7, height = 5)
    }

    # stratified by batch, with overall on top
    p <- plot_prev_by_decade(
        cov,
        sprintf("CHIP prevalence by decade — by batch (%s | %s)", tag_lab, sex_lab),
        strat_var = "Batch"
    )
    ggsave(file.path("ch", "figures", paste0("chip_prevalence_by_decade_by_batch_", tag, ".pdf")),
           p, width = 7, height = 5)
}

