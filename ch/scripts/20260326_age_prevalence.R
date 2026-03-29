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
# write_xlsx(all_ch3, file.path("ch", "data", "ch3.xlsx"))

### CHIP variants
vars <- read_excel(file.path("ch", "data", "ch3.xlsx"))
dim(vars)
# 256
length(unique(vars$Sample.ID))
# 220

### covariates
cov_all <- read.csv(file.path("ch", "data", "pmbb_brca12_cov_df.csv"), row.names = 1) %>% filter(person_id %in% all_ids)

cov_all <- cov_all %>%
    mutate(
        CHIP_Binary = person_id %in% vars$Sample.ID,
        CHIP_Count  = sapply(person_id, function(id) sum(vars$Sample.ID == id))
    )


cov_s1 <- cov_all %>% filter(Strata == 1)
dim(cov_all)
# 3004
dim(cov_s1)
# 2605

# write_xlsx(cov_s1, file.path("ch", "data", "pmbb_brca12_s1_cov_chip_df.xlsx"))
# write_xlsx(cov_all, file.path("ch", "data", "pmbb_brca12_cov_chip_df.xlsx"))

# ========================
# QUICK STATS
# ========================
sink(file.path("ch", "log", "ch_summary_stats.log"), split = TRUE)
cat("Run date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")

chip_crosstab <- function(df, strat_var) {
    tab_n    <- table(df[[strat_var]], df$CHIP_Binary)
    tab_prop <- prop.table(tab_n, margin = 1) * 100

    # Build combined "n (%)" strings
    combined <- matrix(
        sprintf("%d (%.1f%%)", tab_n, tab_prop),
        nrow = nrow(tab_n),
        dimnames = list(rownames(tab_n), c("No CHIP", "CHIP"))
    )

    # Total row
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

for (tag in c("all", "s1")) {
    cov <- get(paste0("cov_", tag))
    cat("\n", strrep("=", 40), "\n")
    cat(" Dataset:", tag, "| N =", nrow(cov), "\n")
    cat(strrep("=", 40), "\n")

    cat("\nOverall CHIP prevalence:\n")
    ov_n    <- table(cov$CHIP_Binary)
    ov_prop <- prop.table(ov_n) * 100
    cat(sprintf("  No CHIP : %d (%.1f%%)\n",  ov_n["FALSE"], ov_prop["FALSE"]))
    cat(sprintf("  CHIP    : %d (%.1f%%)\n",  ov_n["TRUE"],  ov_prop["TRUE"]))

    for (sv in c("BRCA12_Case", "Carrier", "Sequenced_gender", "Batch")) {
        chip_crosstab(cov, sv)
    }
}
sink()

# ========================
# HELPER: AGE-DECADE PREVALENCE PLOT
# ========================
plot_prev_by_decade <- function(df, title_str, max_age = Inf) {
    breaks <- if (is.finite(max_age)) c(-Inf, 30, 40, 50, 60, max_age)
    else                     c(-Inf, 30, 40, 50, 60, 70, Inf)
    labels <- if (is.finite(max_age)) c("≤30", "30-40", "40-50", "50-60", paste0("60-", max_age))
    else                     c("≤30", "30-40", "40-50", "50-60", "60-70", "≥70")

    prev_by_decade <- df %>%
        mutate(age_group = cut(Sample_age, breaks = breaks,
                               labels = labels, right = TRUE)) %>%
        group_by(age_group) %>%
        summarise(n_total  = n(),
                  n_chip   = sum(CHIP_Binary),
                  prev_pct = 100 * mean(CHIP_Binary),
                  .groups  = "drop") %>%
        mutate(x_num = as.numeric(age_group))

    ggplot(prev_by_decade, aes(x = x_num, y = prev_pct)) +
        geom_smooth(method = "loess", se = FALSE, color = "#C0392B", linewidth = 1) +
        geom_point(color = "#C0392B", size = 3) +
        scale_x_continuous(breaks = prev_by_decade$x_num,
                           labels = prev_by_decade$age_group) +
        scale_y_continuous(limits = c(0, NA),
                           expand = expansion(mult = c(0, 0.1))) +
        labs(title = title_str, x = "Age group (years)", y = "CHIP prevalence (%)") +
        theme_minimal(base_size = 13) +
        theme(panel.grid.minor   = element_blank(),
              panel.grid.major.x = element_blank(),
              plot.title         = element_text(face = "bold", hjust = 0.5))
}

# ========================
# HELPER: AGE-DECADE PREVALENCE PLOT STRATIFIED
# ========================
plot_prev_stratified <- function(df, strat_var, title_str, max_age = Inf) {
    breaks <- if (is.finite(max_age)) c(-Inf, 30, 40, 50, 60, max_age)
    else                     c(-Inf, 30, 40, 50, 60, 70, Inf)
    labels <- if (is.finite(max_age)) c("≤30", "30-40", "40-50", "50-60", paste0("60-", max_age))
    else                     c("≤30", "30-40", "40-50", "50-60", "60-70", "≥70")

    prev_by_decade <- df %>%
        mutate(age_group  = cut(Sample_age, breaks = breaks,
                                labels = labels, right = TRUE),
               strat_fac  = as.factor(.data[[strat_var]])) %>%
        group_by(age_group, strat_fac) %>%
        summarise(n_total  = n(),
                  n_chip   = sum(CHIP_Binary),
                  prev_pct = 100 * mean(CHIP_Binary),
                  .groups  = "drop") %>%
        mutate(x_num = as.numeric(age_group))

    ggplot(prev_by_decade, aes(x = x_num, y = prev_pct,
                               color = strat_fac, group = strat_fac)) +
        geom_smooth(method = "loess", se = FALSE, linewidth = 1) +
        geom_point(size = 3) +
        scale_x_continuous(breaks = unique(prev_by_decade$x_num),
                           labels = unique(prev_by_decade$age_group)) +
        scale_y_continuous(limits = c(0, NA),
                           expand = expansion(mult = c(0, 0.1))) +
        labs(title  = title_str,
             x      = "Age group (years)",
             y      = "CHIP prevalence (%)",
             color  = strat_var) +
        theme_minimal(base_size = 13) +
        theme(panel.grid.minor   = element_blank(),
              panel.grid.major.x = element_blank(),
              plot.title         = element_text(face = "bold", hjust = 0.5),
              legend.position    = "bottom")
}

# ========================
# CHIP PREVALENCE BY DECADE — ALL / NO-BRCA12 / BRCA12
# (loop over cov_all and cov_s1)
# ========================
for (tag in c("all", "s1")) {
    cov <- get(paste0("cov_", tag))

    # Overall
    p <- plot_prev_by_decade(cov, "CHIP prevalence by decade")
    ggsave(file.path("ch", "figures",
                     paste0("chip_prevalence_by_decade_", tag, ".pdf")),
           p, width = 7, height = 5)

    # Non-carriers only
    p <- plot_prev_by_decade(cov %>% filter(BRCA12_Case == 0),
                             "CHIP prevalence by decade (non-carriers)")
    ggsave(file.path("ch", "figures",
                     paste0("chip_prevalence_by_decade_no_brca12_", tag, ".pdf")),
           p, width = 7, height = 5)

    # BRCA1/2 carriers only (age capped at 70)
    p <- plot_prev_by_decade(cov %>% filter(BRCA12_Case == 1),
                             "CHIP prevalence by decade (BRCA1/2 carriers)",
                             max_age = 70)
    ggsave(file.path("ch", "figures",
                     paste0("chip_prevalence_by_decade_brca12_", tag, ".pdf")),
           p, width = 7, height = 5)

    # Stratified by Sequenced_gender
    p <- plot_prev_stratified(cov, "Sequenced_gender",
                              "CHIP prevalence by decade (by sex)")
    ggsave(file.path("ch", "figures",
                     paste0("chip_prevalence_by_decade_by_sex_", tag, ".pdf")),
           p, width = 7, height = 5)

    # Stratified by BRCA12_Case
    p <- plot_prev_stratified(cov, "BRCA12_Case",
                              "CHIP prevalence by decade (carriers vs controls)")
    ggsave(file.path("ch", "figures",
                     paste0("chip_prevalence_by_decade_by_case_", tag, ".pdf")),
           p, width = 7, height = 5)

    # Stratified by Batch
    p <- plot_prev_stratified(cov, "Batch",
                              "CHIP prevalence by decade (by batch)")
    ggsave(file.path("ch", "figures",
                     paste0("chip_prevalence_by_decade_by_batch_", tag, ".pdf")),
           p, width = 7, height = 5)
}

# ========================
# COVERAGE DISTRIBUTIONS (Total Depth)
# ========================
# Join variant-level depth info onto covariates
# (vars is expected to have a Total.Depth / Total_Depth column)
depth_col <- "Sample.Depth"

# summarise vars to one row per person first
vars_depth <- vars %>%
    group_by(Sample.ID) %>%
    summarise(Total_Depth = mean(Sample.Depth, na.rm = TRUE), .groups = "drop") %>%
    rename(person_id = Sample.ID)

cov_depth_all <- cov_all %>% left_join(vars_depth, by = "person_id")
cov_depth_s1  <- cov_s1  %>% left_join(vars_depth, by = "person_id")

plot_coverage <- function(df, strat_var, title_str) {
    ggplot(df %>% filter(!is.na(Total_Depth)),
           aes(x = Total_Depth, fill = as.factor(.data[[strat_var]]))) +
        geom_histogram(binwidth = 10, color = "white", alpha = 0.7,
                       position = "identity") +
        scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
        scale_fill_brewer(palette = "Set2") +
        labs(title = title_str,
             x     = "Total Depth",
             y     = "Frequency",
             fill  = strat_var) +
        theme_minimal(base_size = 13) +
        theme(panel.grid.minor = element_blank(),
              plot.title        = element_text(face = "bold", hjust = 0.5),
              legend.position   = "bottom")
}

for (tag in c("all", "s1")) {
    cov_depth <- get(paste0("cov_depth_", tag))

    for (sv in c("Carrier", "Batch", "BRCA12_Case")) {
        p <- plot_coverage(cov_depth, sv,
                           paste0("Coverage distribution by ", sv))
        ggsave(file.path("ch", "figures",
                         paste0("coverage_by_", tolower(sv), "_", tag, ".pdf")),
               p, width = 7, height = 5)
    }
}

# ========================
# COVERAGE DISTRIBUTIONS (Alt Depth)
# ========================
depth_col <- "Sample.AltDepth"

# summarise vars to one row per person first
vars_depth <- vars %>%
    group_by(Sample.ID) %>%
    summarise(Sample.AltDepth = round(mean(Sample.AltDepth, na.rm = TRUE)),  # <- round
              .groups = "drop") %>%
    rename(person_id = Sample.ID)

cov_depth_all <- cov_all %>% left_join(vars_depth, by = "person_id")
cov_depth_s1  <- cov_s1  %>% left_join(vars_depth, by = "person_id")

plot_coverage <- function(df, strat_var, title_str) {
    ggplot(df %>% filter(!is.na(Sample.AltDepth)),
           aes(x = Sample.AltDepth, fill = as.factor(.data[[strat_var]]))) +
        geom_bar(color = "white", alpha = 0.7, position = "identity") +
        scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
        scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
        scale_fill_brewer(palette = "Set2") +
        labs(title = title_str,
             x     = "Alt Depth",
             y     = "Number of individuals",
             fill  = strat_var) +
        theme_minimal(base_size = 13) +
        theme(panel.grid.minor = element_blank(),
              plot.title        = element_text(face = "bold", hjust = 0.5),
              legend.position   = "bottom")
}

for (tag in c("all", "s1")) {
    cov_depth <- get(paste0("cov_depth_", tag))

    for (sv in c("Carrier", "Batch", "BRCA12_Case")) {
        p <- plot_coverage(cov_depth, sv,
                           paste0("Coverage distribution by ", sv))
        ggsave(file.path("ch", "figures",
                         paste0("altdepth_by_", tolower(sv), "_", tag, ".pdf")),
               p, width = 7, height = 5)
    }
}



# ========================
# HISTOGRAM OF AGE
# ========================
for (tag in c("all", "s1")) {
    cov <- get(paste0("cov_", tag))
    p <- ggplot(cov, aes(x = Sample_age)) +
        geom_histogram(binwidth = 5, fill = "steelblue", color = "white") +
        scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
        labs(title = "Age distribution", x = "Age (years)", y = "Frequency") +
        theme_minimal(base_size = 13) +
        theme(panel.grid.minor = element_blank(),
              plot.title        = element_text(face = "bold", hjust = 0.5))
    ggsave(file.path("ch", "figures",
                     paste0("age_hist_", tag, ".png")),
           p, width = 7, height = 5)
}

# ========================
# CHIP COUNTS PER INDIVIDUAL
# ========================
for (tag in c("all", "s1")) {
    cov <- get(paste0("cov_", tag))

    chip_counts <- cov %>%
        filter(CHIP_Binary) %>%
        count(CHIP_Count, name = "n_individuals") %>%
        mutate(CHIP_Count = factor(CHIP_Count))   # <- treat as discrete

    p <- ggplot(chip_counts, aes(x = CHIP_Count, y = n_individuals)) +
        geom_col(width = 0.65, fill = "gray50") +
        geom_text(aes(label = n_individuals), vjust = -0.4,
                  size = 3.5, color = "gray20") +
        scale_y_log10(expand = expansion(mult = c(0, 0.15))) +
        # no scale_x_continuous — factor handles axis labels automatically
        labs(title = "Number of CHIP variants per individual",
             x     = "Number of variants",
             y     = "Number of individuals (log scale)") +
        theme_minimal(base_size = 13) +
        theme(panel.grid.minor   = element_blank(),
              panel.grid.major.x = element_blank(),
              plot.title         = element_text(face = "bold", hjust = 0.5))

    ggsave(file.path("ch", "figures",
                     paste0("chip_variant_count_per_person_", tag, ".pdf")),
           p, width = 7, height = 5)
}

high_freq <- cov_all %>% filter(CHIP_Count > 3)
# the one with 9 prob just have 7 or 8
# PMBB7688103445285
# make a manual review sheet for this...


# ========================
# MEAN VAF SECTION
# ========================
# Assumes vars has a VAF column (adjust name if needed: "VAF", "vaf", "AF")
vaf_col <- "Sample.AltFrac"
cov <- cov_all
vars_vaf <- vars %>%
    rename(person_id = Sample.ID, VAF = !!sym(vaf_col))

# --- 1. Overall VAF distribution ---
p_vaf_dist <- ggplot(vars_vaf, aes(x = VAF)) +
    geom_histogram(binwidth = 0.01, fill = "#2980B9", color = "white") +
    scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    labs(title = "Distribution of CHIP VAF",
         x     = "VAF",
         y     = "Number of variants") +
    theme_minimal(base_size = 13) +
    theme(panel.grid.minor = element_blank(),
          plot.title        = element_text(face = "bold", hjust = 0.5))

ggsave(file.path("ch", "figures", "chip_vaf_distribution.pdf"),
       p_vaf_dist, width = 7, height = 5)

# --- 2. Mean VAF by age decade (both cov_all and cov_s1) ---
plot_vaf_by_decade <- function(cov_df, vars_vaf, title_str, max_age = Inf) {
    breaks <- if (is.finite(max_age)) c(-Inf, 30, 40, 50, 60, max_age)
    else                     c(-Inf, 30, 40, 50, 60, 70, Inf)
    labels <- if (is.finite(max_age)) c("≤30", "30-40", "40-50", "50-60", paste0("60-", max_age))
    else                     c("≤30", "30-40", "40-50", "50-60", "60-70", "≥70")

    df <- cov_df %>%
        dplyr::select(person_id, Sample_age) %>%
        inner_join(vars_vaf %>% dplyr::select(person_id, VAF), by = "person_id") %>%
        mutate(age_group = cut(Sample_age, breaks = breaks,
                               labels = labels, right = TRUE)) %>%
        group_by(age_group) %>%
        summarise(mean_vaf = mean(VAF, na.rm = TRUE),
                  sd_vaf   = sd(VAF,   na.rm = TRUE),
                  n        = n(),
                  se_vaf   = sd_vaf / sqrt(n),
                  .groups  = "drop") %>%
        mutate(x_num = as.numeric(age_group))

    ggplot(df, aes(x = x_num, y = mean_vaf)) +
        geom_smooth(method = "loess", se = FALSE, color = "#8E44AD", linewidth = 1) +
        geom_errorbar(aes(ymin = mean_vaf - se_vaf,
                          ymax = mean_vaf + se_vaf),
                      width = 0.2, color = "#8E44AD") +
        geom_point(color = "#8E44AD", size = 3) +
        scale_x_continuous(breaks = df$x_num, labels = df$age_group) +
        scale_y_continuous(labels = scales::percent_format(accuracy = 0.1),
                           expand = expansion(mult = c(0, 0.1))) +
        labs(title = title_str,
             x     = "Age group (years)",
             y     = "Mean VAF (± SE)") +
        theme_minimal(base_size = 13) +
        theme(panel.grid.minor   = element_blank(),
              panel.grid.major.x = element_blank(),
              plot.title         = element_text(face = "bold", hjust = 0.5))
}

for (tag in c("all", "s1")) {
    cov <- get(paste0("cov_", tag))

    # Overall
    p <- plot_vaf_by_decade(cov, vars_vaf, "Mean CHIP VAF by decade")
    ggsave(file.path("ch", "figures",
                     paste0("chip_mean_vaf_by_decade_", tag, ".pdf")),
           p, width = 7, height = 5)

    # Non-carriers
    p <- plot_vaf_by_decade(cov %>% filter(BRCA12_Case == 0), vars_vaf,
                            "Mean CHIP VAF by decade (non-carriers)")
    ggsave(file.path("ch", "figures",
                     paste0("chip_mean_vaf_by_decade_no_brca12_", tag, ".pdf")),
           p, width = 7, height = 5)

    # Carriers (age capped at 70)
    p <- plot_vaf_by_decade(cov %>% filter(BRCA12_Case == 1), vars_vaf,
                            "Mean CHIP VAF by decade (BRCA1/2 carriers)",
                            max_age = 70)
    ggsave(file.path("ch", "figures",
                     paste0("chip_mean_vaf_by_decade_brca12_", tag, ".pdf")),
           p, width = 7, height = 5)
}

# --- 3. Mean VAF by age — stratified (sex, BRCA12_Case, Batch) ---
plot_vaf_stratified <- function(cov_df, vars_vaf, strat_var, title_str, max_age = Inf) {
    breaks <- if (is.finite(max_age)) c(-Inf, 30, 40, 50, 60, max_age)
    else                     c(-Inf, 30, 40, 50, 60, 70, Inf)
    labels <- if (is.finite(max_age)) c("≤30", "30-40", "40-50", "50-60", paste0("60-", max_age))
    else                     c("≤30", "30-40", "40-50", "50-60", "60-70", "≥70")

    df <- cov_df %>%
        dplyr::select(person_id, Sample_age, !!sym(strat_var)) %>%
        inner_join(vars_vaf %>% dplyr::select(person_id, VAF), by = "person_id") %>%
        mutate(age_group = cut(Sample_age, breaks = breaks,
                               labels = labels, right = TRUE),
               strat_fac = as.factor(.data[[strat_var]])) %>%
        group_by(age_group, strat_fac) %>%
        summarise(mean_vaf = mean(VAF, na.rm = TRUE),
                  se_vaf   = sd(VAF, na.rm = TRUE) / sqrt(n()),
                  .groups  = "drop") %>%
        mutate(x_num = as.numeric(age_group))

    ggplot(df, aes(x = x_num, y = mean_vaf,
                   color = strat_fac, group = strat_fac)) +
        geom_smooth(method = "loess", se = FALSE, linewidth = 1) +
        geom_errorbar(aes(ymin = mean_vaf - se_vaf,
                          ymax = mean_vaf + se_vaf),
                      width = 0.2) +
        geom_point(size = 3) +
        scale_x_continuous(breaks = unique(df$x_num),
                           labels = unique(df$age_group)) +
        scale_y_continuous(labels = scales::percent_format(accuracy = 0.1),
                           expand = expansion(mult = c(0, 0.1))) +
        labs(title = title_str,
             x     = "Age group (years)",
             y     = "Mean VAF (± SE)",
             color = strat_var) +
        theme_minimal(base_size = 13) +
        theme(panel.grid.minor   = element_blank(),
              panel.grid.major.x = element_blank(),
              plot.title         = element_text(face = "bold", hjust = 0.5),
              legend.position    = "bottom")
}

for (tag in c("all", "s1")) {
    cov <- get(paste0("cov_", tag))
    for (sv in c("Sequenced_gender", "BRCA12_Case", "Batch")) {
        p <- plot_vaf_stratified(cov, vars_vaf, sv,
                                 paste0("Mean CHIP VAF by decade (by ", sv, ")"))
        ggsave(file.path("ch", "figures",
                         paste0("chip_mean_vaf_by_decade_by_",
                                tolower(sv), "_", tag, ".pdf")),
               p, width = 7, height = 5)
    }
}

# ========================
# BATCH ANALYSIS
# ========================
# Are age distributions different between batches?
cov_all %>% group_by(Batch) %>%
    summarise(mean_age = mean(Sample_age), sd_age = sd(Sample_age), n = n())
# no bc i mathed..

# Quick test
kruskal.test(Sample_age ~ Batch, data = cov_all)
table(cov_all$Batch, cov_all$Sequenced_gender)

m <- glm(CHIP_Binary ~ Batch + Sample_age + Sequenced_gender,
         data = cov_all, family = binomial)
summary(m)
exp(cbind(OR = coef(m), confint(m)))  # odds ratios

gene_batch <- vars %>%
    left_join(cov_all %>% dplyr::select(person_id),
              by = c("Sample.ID" = "person_id")) %>%
    count(Batch, Gene) %>%
    group_by(Batch) %>%
    mutate(pct = 100 * n / sum(n)) %>%
    ungroup() %>%
    pivot_wider(
        names_from  = Batch,
        values_from = c(n, pct),
        values_fill = 0
    ) %>%
    arrange(desc(pct_2))
write_xlsx(gene_batch, file.path("ch", "data", "ch_gene_batch_counts.xlsx"))

# vars %>%
#     left_join(cov_all %>% dplyr::select(person_id),
#               by = c("Sample.ID" = "person_id")) %>%
#     ggplot(aes(x = Sample.AltFrac, fill = Batch)) +
#     geom_histogram(binwidth = 0.01, position = "identity", alpha = 0.6) +
#     facet_wrap(~Batch, ncol = 1) +
#     scale_x_continuous(labels = scales::percent_format()) +
#     theme_minimal()



# ========================
# ALT DEPTH ANALYSIS
# ========================







