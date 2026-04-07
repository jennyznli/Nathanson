# ============================================================
# CHIP prevalence + age analysis
# ============================================================
library(here)
setwd("/Users/jennyzli/Documents/Nathanson")
source(here("R", "config.R"))
library(binom)

# ============================================================
# READ DATA
# ============================================================
cov <- read_excel(file.path("ch", "data", "pmbb_brca12_cov_chip_df.xlsx"))
dim(cov)

vars <- read_excel(file.path("ch", "data", "ch_seq_wl_art_minad4_vars.xlsx"))
vars <- vars %>% filter(Sample.ID %in% cov$person_id)

cov <- cov %>%
    mutate(
        CHIP_Binary = person_id %in% vars$Sample.ID,
        CHIP_Count  = sapply(person_id, function(id) sum(vars$Sample.ID == id))
    )

cat("Overall cohort: N =", nrow(cov),
    "| CHIP+ =", sum(cov$CHIP_Binary),
    sprintf("(%.1f%%)\n", 100 * mean(cov$CHIP_Binary)))
# Overall cohort: N = 3004 | CHIP+ = 193 (6.4%)

# ============================================================
# CHIP BINARY - CHI-SQUARED LOOP
# ============================================================
covariates <- c("BRCA12_Case", "Batch", "Sequenced_gender",
                "Smoke_History", "Strata", "Class")

chisq_results <- lapply(covariates, function(var) {

    df_sub <- cov[!is.na(cov[[var]]) & !is.na(cov$CHIP_Binary), ]

    tryCatch({
        tab    <- table(df_sub[[var]], df_sub$CHIP_Binary)
        n_levels <- nrow(tab)

        # use Fisher if any expected cell < 5
        expected <- chisq.test(tab)$expected

        if (any(expected < 5) || n_levels > 2) {
            # Fisher with simulation for tables larger than 2x2
            test   <- fisher.test(tab, simulate.p.value = (n_levels > 2))
            method <- ifelse(n_levels > 2, "Fisher (simulated)", "Fisher")
            xsq    <- NA
            df     <- NA
        } else {
            test   <- chisq.test(tab)
            method <- "Chi-squared"
            xsq    <- round(test$statistic, 4)
            df     <- test$parameter
        }

        data.frame(
            Covariate = var,
            Method    = method,
            N_levels  = n_levels,
            X_squared = xsq,
            df        = df,
            P_value   = round(test$p.value, 4),
            Sig       = ifelse(test$p.value < 0.001, "***",
                               ifelse(test$p.value < 0.01,  "**",
                                      ifelse(test$p.value < 0.05,  "*", "ns")))
        )
    }, error = function(e) {
        data.frame(Covariate = var, Method = "ERROR", N_levels = NA,
                   X_squared = NA, df = NA, P_value = NA, Sig = NA)
    })
})

chisq_df <- do.call(rbind, chisq_results)

# --- Print ---
cat("\n--- CHIP Binary: Chi-squared/Fisher Results ---\n")
print(chisq_df, row.names = FALSE)
# Covariate             Method N_levels X_squared df P_value Sig
# BRCA12_Case        Chi-squared        2    3.8408  1  0.0500  ns
# Batch        Chi-squared        2    0.8318  1  0.3618  ns
# Sequenced_gender        Chi-squared        2    2.3128  1  0.1283  ns
# Smoke_History        Chi-squared        2    1.9259  1  0.1652  ns
# Strata        Chi-squared        2    1.6537  1  0.1985  ns
# Class Fisher (simulated)        6        NA NA  0.2284  ns

# --- Export ---
write.csv(chisq_df,
          file.path("ch", "data", "chip_binary_chisq_results.csv"),
          row.names = FALSE)

# ============================================================
# WILCOXON TESTS - CHIP COUNT BY COVARIATES
# ============================================================
covariates <- c("BRCA12_Case", "Sequenced_gender", "Smoke_History",
                "Batch", "Strata", "Class")

# --- Run all Wilcoxon tests ---
wilcox_results <- lapply(covariates, function(var) {

    # drop NAs for this variable
    df_sub <- cov[!is.na(cov[[var]]) & !is.na(cov$CHIP_Count), ]

    n_levels <- length(unique(df_sub[[var]]))

    tryCatch({
        if (n_levels == 2) {
            # standard Wilcoxon for binary vars
            test <- wilcox.test(CHIP_Count ~ get(var), data = df_sub)
            pval <- test$p.value
            method <- "Wilcoxon"
        } else {
            # Kruskal-Wallis for 3+ categories
            test <- kruskal.test(CHIP_Count ~ get(var), data = df_sub)
            pval <- test$p.value
            method <- "Kruskal-Wallis"
        }

        data.frame(
            Covariate = var,
            Method    = method,
            N_levels  = n_levels,
            P_value   = round(pval, 4),
            Sig       = ifelse(pval < 0.001, "***",
                               ifelse(pval < 0.01,  "**",
                                      ifelse(pval < 0.05,  "*", "ns")))
        )
    }, error = function(e) {
        data.frame(Covariate = var, Method = "ERROR", N_levels = n_levels,
                   P_value = NA, Sig = NA)
    })
})

results_df <- do.call(rbind, wilcox_results)

# --- Print summary ---
cat("\n--- CHIP Count: Wilcoxon/Kruskal-Wallis Results ---\n")
print(results_df, row.names = FALSE)
# Covariate         Method N_levels P_value Sig
# BRCA12_Case       Wilcoxon        2  0.0362   *
#     Sequenced_gender       Wilcoxon        2  0.1051  ns
# Smoke_History       Wilcoxon        2  0.1417  ns
# Batch       Wilcoxon        2  0.3233  ns
# Strata       Wilcoxon        2  0.1440  ns
# Class Kruskal-Wallis        6  0.2535  ns

# --- Export ---
write.csv(results_df,
          file.path("ch", "data", "chip_count_wilcoxon_results.csv"),
          row.names = FALSE)

# ============================================================
# SHARED THEME
# ============================================================
theme_ch <- theme_minimal(base_size = 12) +
    theme(
        plot.title         = element_text(face = "bold", hjust = 0.5),
        plot.subtitle      = element_text(hjust = 0.5, color = "grey40", size = 10),
        panel.grid.minor   = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.position    = "bottom"
    )
# ============================================================
# TABLE 1 - PVALs
# ============================================================
# assumes your table1 data is in `cov` with BRCA12_Case as 0/1
# and columns: Sample_age, Sequenced_gender, Smoke_History, Batch, Strata, Ancestry

# --- Continuous: Age ---
age_test <- wilcox.test(Sample_age ~ BRCA12_Case, data = cov)

cat("Age Wilcoxon p =", round(age_test$p.value, 4), "\n")
# ge Wilcoxon p = 0.8702

# --- Categorical variables ---
cat_vars <- c("Sequenced_gender", "Smoke_History", "Batch", "Strata", "Class")

table1_tests <- lapply(cat_vars, function(var) {
    df_sub <- cov[!is.na(cov[[var]]), ]
    tab <- table(df_sub[[var]], df_sub$BRCA12_Case)
    expected <- chisq.test(tab)$expected

    if (any(expected < 5) || nrow(tab) > 2) {
        test   <- fisher.test(tab, simulate.p.value = TRUE)
        method <- "Fisher"
    } else {
        test   <- chisq.test(tab)
        method <- "Chi-squared"
    }

    data.frame(
        Variable = var,
        Method   = method,
        P_value  = round(test$p.value, 4),
        Sig      = ifelse(test$p.value < 0.001, "***",
                          ifelse(test$p.value < 0.01,  "**",
                                 ifelse(test$p.value < 0.05,  "*", "ns")))
    )
}) %>% do.call(rbind, .)

print(table1_tests)
write_xlsx(table1_tests,
           file.path("ch", "data", "table1.xlsx"))


### POST HOC FISHER
class_levels <- unique(cov$Class[!is.na(cov$Class)])

pairwise_fisher <- combn(class_levels, 2, simplify = FALSE) %>%
    lapply(function(pair) {
        df_sub <- cov %>% filter(Class %in% pair)
        tab <- table(df_sub$Class, df_sub$BRCA12_Case)
        test <- fisher.test(tab)
        data.frame(
            Class1  = pair[1],
            Class2  = pair[2],
            P_value = round(test$p.value, 4),
            OR      = round(test$estimate, 3)
        )
    }) %>%
    do.call(rbind, .) %>%
    mutate(P_adj = round(p.adjust(P_value, method = "BH"), 4),
           Sig   = ifelse(P_adj < 0.001, "***",
                          ifelse(P_adj < 0.01,  "**",
                                 ifelse(P_adj < 0.05,  "*", "ns"))))

print(pairwise_fisher)


# ============================================================
# TABLE 1
# ============================================================
library(gt)
# --- Age summary---
age_summary <- cov %>%
    group_by(BRCA12_Case) %>%
    summarise(
        median_age = median(Sample_age, na.rm = TRUE),
        q1         = quantile(Sample_age, 0.25, na.rm = TRUE),
        q3         = quantile(Sample_age, 0.75, na.rm = TRUE),
        .groups = "drop"
    ) %>%
    mutate(label = paste0(median_age, " [", q1, "–", q3, "]"))

age_wilcox_p <- wilcox.test(Sample_age ~ BRCA12_Case, data = cov)$p.value
cat("Age: median [IQR]\n")
cat("  Cases:    ", age_summary$label[age_summary$BRCA12_Case == 1], "\n")
cat("  Controls: ", age_summary$label[age_summary$BRCA12_Case == 0], "\n")
cat("  Wilcoxon p =", round(age_wilcox_p, 4), "\n")

# --- Helper to get n (%) for a categorical variable ---
get_counts <- function(var, data = cov) {
    data %>%
        filter(!is.na(.data[[var]])) %>%
        group_by(BRCA12_Case, .data[[var]]) %>%
        summarise(n = n(), .groups = "drop") %>%
        group_by(BRCA12_Case) %>%
        mutate(pct = round(100 * n / sum(n), 1),
               label = paste0(n, " (", pct, "%)")) %>%
        ungroup() %>%
        select(BRCA12_Case, level = all_of(var), label) %>%
        pivot_wider(names_from = BRCA12_Case,
                    values_from = label,
                    names_prefix = "group_") %>%
        rename(Cases = group_1, Controls = group_0)
}

# --- Helper to get p-value ---
get_pval <- function(var, data = cov) {
    df_sub <- data[!is.na(data[[var]]), ]
    tab <- table(df_sub[[var]], df_sub$BRCA12_Case)
    expected <- suppressWarnings(chisq.test(tab)$expected)
    if (any(expected < 5) || nrow(tab) > 2) {
        p <- fisher.test(tab, simulate.p.value = TRUE)$p.value
    } else {
        p <- chisq.test(tab)$p.value
    }
    sig <- ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "")))
    paste0(ifelse(p < 0.001, "<0.001", round(p, 3)), sig)
}

# --- Build each section ---
make_section <- function(var, label, data = cov) {
    rows <- get_counts(var, data)
    pval <- get_pval(var, data)
    header <- tibble(level = paste0(label, "  ", pval), Cases = "", Controls = "")
    bind_rows(header, rows %>% mutate(level = paste0("  ", level)))
}

# --- Assemble table ---
table1 <- bind_rows(
    tibble(level = paste0("Age—median [IQR], year  p=", round(age_wilcox_p, 3)),
           Cases    = age_summary$label[age_summary$BRCA12_Case == 1],
           Controls = age_summary$label[age_summary$BRCA12_Case == 0]),
    make_section("Sequenced_gender", "Sex"),
    make_section("Smoke_History",    "Smoking history"),
    make_section("Batch",            "Sequencing freeze"),
    make_section("Strata",           "Strata"),
    make_section("Class",            "Ancestry†")
)

# --- Render with gt ---
table1 %>%
    gt() %>%
    cols_label(level = "Characteristic",
               Cases = "Cases (N = 675)",
               Controls = "Controls (N = 2,327)") %>%
    tab_footnote("*P<0.05; **P<0.01; ***P<0.001. Chi-squared or Fisher's exact test for categorical variables; Wilcoxon rank-sum for age.") %>%
    tab_footnote("†Ancestry determined by genetic principal components; p-value shown for completeness but ancestry was accounted for via PCs in matching and regression.") %>%
    tab_style(style = cell_text(weight = "bold"),
              locations = cells_column_labels()) %>%
    gtsave(file.path("ch", "figures", "table1.docx"))


# ============================================================
# CHIP PREVALENCE BY AGE
# ============================================================
age_breaks <- c(-Inf, 30, 40, 50, 60, Inf)
age_labels <- c("≤30", "30–40", "40–50", "50–60", "≥60")

prev_overall <- cov %>%
    mutate(age_group = cut(Sample_age, breaks = age_breaks,
                           labels = age_labels, right = TRUE)) %>%
    group_by(age_group) %>%
    summarise(n_total  = n(),
              n_chip   = sum(CHIP_Binary),
              prev_pct = 100 * mean(CHIP_Binary),
              .groups  = "drop") %>%
    mutate(x_num = as.numeric(age_group))

prev_by_carrier <- cov %>%
    mutate(age_group = cut(Sample_age, breaks = age_breaks,
                           labels = age_labels, right = TRUE),
           Carrier_label = ifelse(BRCA12_Case == 1,
                                  "gBRCA1/2 Carrier", "Non-carrier")) %>%
    group_by(age_group, Carrier_label) %>%
    summarise(n_total  = n(),
              n_chip   = sum(CHIP_Binary),
              prev_pct = 100 * mean(CHIP_Binary),
              .groups  = "drop") %>%
    mutate(x_num = as.numeric(age_group))

# fig_prev <- ggplot(prev_overall, aes(x = x_num, y = prev_pct)) +
#     geom_smooth(method = "loess", span = 1.2, se = FALSE,
#                 color = "#546E7A", linewidth = 0.9) +
#     geom_point(size = 3, color = "#546E7A") +
#     scale_x_continuous(breaks = prev_overall$x_num,
#                        labels = prev_overall$age_group) +
#     scale_y_continuous(limits = c(0, NA),
#                        labels = function(x) paste0(x, "%"),
#                        expand = expansion(mult = c(0.02, 0.1))) +
#     labs(
#         title   = "CHIP prevalence by age decade",
#         x       = "Age group (years)",
#         y       = "CHIP prevalence (%)"
#     ) +
#     theme_ch +
#     theme(plot.caption = element_text(size = 8, color = "grey50", hjust = 0))
#
# ggsave(file.path("ch", "figures", "fig_chip_prevalence_by_decade.pdf"),
#        fig_prev, width = 5, height = 4)

# shared y-axis and x-axis scale for both figures
scale_y_prev <- scale_y_continuous(
    limits = c(0, NA),
    labels = function(x) paste0(x, "%"),
    expand = expansion(mult = c(0.02, 0.1))
)
scale_x_prev <- scale_x_continuous(
    breaks = prev_overall$x_num,
    labels = prev_overall$age_group
)

# ── Overall ──
fig_prev_overall <- ggplot(prev_overall, aes(x = x_num, y = prev_pct)) +
    geom_smooth(method = "loess", span = 1.2, se = FALSE,
                color = "#546E7A", linewidth = 0.9) +
    geom_point(size = 3, color = "#546E7A") +
    scale_x_prev + scale_y_prev +
    labs(title = "CHIP prevalence by age decade",
         x = "Age group (years)", y = "CHIP prevalence (%)") +
    theme_ch

ggsave(file.path("ch", "figures", "fig_chip_prevalence_by_decade.pdf"),
       fig_prev_overall, width = 6, height = 5)

# ── Stratified by carrier ──
fig_prev_carrier <- ggplot() +
    geom_smooth(data = prev_by_carrier,
                aes(x = x_num, y = prev_pct, color = Carrier_label),
                method = "loess", span = 1.2, se = FALSE, linewidth = 0.9) +
    geom_point(data = prev_by_carrier,
               aes(x = x_num, y = prev_pct, color = Carrier_label),
               size = 3) +
    geom_smooth(data = prev_overall,
                aes(x = x_num, y = prev_pct, group = 1),
                method = "loess", span = 1.2, se = FALSE,
                color = "darkgray", linewidth = 0.8) +
    geom_point(data = prev_overall,
               aes(x = x_num, y = prev_pct),
               color = "darkgray", size = 3, shape = 18) +
    scale_color_manual(
        values = c("gBRCA1/2 Carrier" = "#C2185B", "Non-carrier" = "#3e77c1"),
        name = NULL
    ) +
    scale_x_prev + scale_y_prev +
    labs(title = "CHIP prevalence by age decade",
         x = "Age group (years)", y = "CHIP prevalence (%)") +
    theme_ch

ggsave(file.path("ch", "figures", "fig_chip_prevalence_by_decade_carrier.pdf"),
       fig_prev_carrier, width = 6, height = 5)

fig_prev_carrier_nokey <- fig_prev_carrier + theme(legend.position = "none")

ggsave(file.path("ch", "figures", "fig_chip_prevalence_by_decade_carrier_nokey.pdf"),
       fig_prev_carrier_nokey, width = 6, height = 5)

# ============================================================
# INTERACTION
# ============================================================
library(lme4)
cov <- cov %>%
    mutate(age_group = cut(Sample_age,
                           breaks = c(-Inf, 30, 40, 50, 60, Inf),
                           labels = c("≤30", "30–40", "40–50", "50–60", "≥60"),
                           right  = TRUE))

m_interaction <- glm(CHIP_Binary ~ BRCA12_Case * Sample_age,
                     data = cov, family = binomial)
m_main        <- glm(CHIP_Binary ~ BRCA12_Case + Sample_age,
                     data = cov, family = binomial)

anova(m_main, m_interaction, test = "Chisq")
# Analysis of Deviance Table
#
# Model 1: CHIP_Binary ~ BRCA12_Case + age_group
# Model 2: CHIP_Binary ~ BRCA12_Case * age_group
# Resid. Df Resid. Dev Df Deviance Pr(>Chi)
# 1      2998     1383.6
# 2      2994     1383.1  4  0.52051   0.9715

# ============================================================
# BY GENE - AGE PREVALENCE
# ============================================================
# after creating prev_by_gene, add CI bounds manually
prev_by_gene <- prev_by_gene %>%
    mutate(
        se       = sqrt((prev_pct/100 * (1 - prev_pct/100)) / n_total) * 100,
        ci_lower = pmax(prev_pct - 1.96 * se, 0),    # floor at 0
        ci_upper = prev_pct + 1.96 * se
    ) %>% filter(Carrier_label != "Non-carrier")

# fig_supp_gene <- ggplot() +
#     # stratified lines with CI
#     geom_ribbon(data = prev_by_gene,
#                 aes(x = x_num, ymin = ci_lower, ymax = ci_upper,
#                     fill = Carrier_label),
#                 alpha = 0.15) +
#     geom_line(data = prev_by_gene,
#               aes(x = x_num, y = prev_pct, color = Carrier_label),
#               linewidth = 0.9) +
#     geom_point(data = prev_by_gene,
#                aes(x = x_num, y = prev_pct, color = Carrier_label),
#                size = 3) +
#     scale_color_manual(
#         values = c(
#             "gBRCA1 Carrier" = "#C2185B",
#             "gBRCA2 Carrier" = "#7B1FA2",
#             "Non-carrier"    = "#546E7A"
#         ),
#         name = NULL
#     ) +
#     scale_fill_manual(
#         values = c(
#             "gBRCA1 Carrier" = "#C2185B",
#             "gBRCA2 Carrier" = "#7B1FA2",
#             "Non-carrier"    = "#546E7A"
#         ),
#         name = NULL
#     ) +
#     scale_x_prev + scale_y_prev +
#     labs(
#         title    = "CHIP prevalence by age decade by germline variant",
#         x        = "Age group (years)",
#         y        = "CHIP prevalence (%)",
#     ) +
#     theme_ch +
#     theme(
#         legend.position  = "bottom",
#         plot.caption     = element_text(size = 8, color = "grey50", hjust = 0),
#         plot.subtitle    = element_text(size = 10, color = "grey30")
#     )
#
# ggsave(file.path("ch", "figures", "fig_supp_chip_prevalence_by_gene.pdf"),
#        fig_supp_gene, width = 6, height = 5)


fig_supp_gene <- ggplot() +
    geom_smooth(data = prev_by_gene,
                aes(x = x_num, y = prev_pct, color = Carrier_label),
                method = "loess", span = 1.0, se = FALSE, linewidth = 0.9) +
    geom_point(data = prev_by_gene,
               aes(x = x_num, y = prev_pct, color = Carrier_label),
               size = 3) +
    scale_color_manual(
        values = c(
            "gBRCA1 Carrier" = "#C2185B",
            "gBRCA2 Carrier" = "#7B1FA2",
            "Non-carrier"    = "#546E7A"
        ),
        name = NULL
    ) +
    scale_x_prev + scale_y_prev +
    labs(
        title = "CHIP prevalence by age decade by germline variant",
        x     = "Age group (years)",
        y     = "CHIP prevalence (%)"
    ) +
    theme_ch +
    theme(
        legend.position = "bottom",
        plot.caption    = element_text(size = 8, color = "grey50", hjust = 0),
        plot.subtitle   = element_text(size = 10, color = "grey30")
    )

ggsave(file.path("ch", "figures", "fig_supp_chip_prevalence_by_gene.pdf"),
       fig_supp_gene, width = 6, height = 5)

# ============================================================
# GENE FREQUENCY ANALYSIS
# ============================================================
vars_cov <- vars %>%
    left_join(
        cov %>% dplyr::select(person_id, BRCA12_Case, BRCA1_Case,
                              BRCA2_Case, Carrier, Sample_age,
                              Batch, Smoke_History, starts_with("PC")),
        by = c("Sample.ID" = "person_id")
    ) %>%
    filter(!is.na(BRCA12_Case))

gene_person <- vars_cov %>%
    dplyr::select(Sample.ID, Gene, BRCA12_Case) %>%
    distinct()

n_carrier_total    <- sum(cov$BRCA12_Case == 1)
n_noncarrier_total <- sum(cov$BRCA12_Case == 0)

gene_freq <- gene_person %>%
    group_by(Gene, BRCA12_Case) %>%
    summarise(n_mutated = n(), .groups = "drop") %>%
    pivot_wider(names_from  = BRCA12_Case,
                values_from = n_mutated,
                values_fill = 0,
                names_prefix = "n_") %>%
    rename(n_noncarrier = n_0, n_carrier = n_1) %>%
    mutate(
        pct_noncarrier = n_noncarrier / n_noncarrier_total * 100,
        pct_carrier    = n_carrier    / n_carrier_total    * 100,
        n_total        = n_carrier + n_noncarrier
    ) %>%
    arrange(desc(n_total))

top_genes <- gene_freq %>%
    slice_max(order_by = n_total, n = 15, with_ties = FALSE) %>%
    arrange(desc(n_total)) %>%
    pull(Gene)

dim(gene_freq)
# 46 6
cat("\nGene frequency table (top 15):\n")
print(head(gene_freq, 15))

# ============================================================
# FISHER'S EXACT TEST PER GENE
# ============================================================
gene_fisher <- gene_freq %>%
    rowwise() %>%
    mutate(
        fisher_p = fisher.test(matrix(
            c(n_carrier,
              n_carrier_total    - n_carrier,
              n_noncarrier,
              n_noncarrier_total - n_noncarrier),
            nrow = 2
        ))$p.value,
        OR_crude = (n_carrier / (n_carrier_total - n_carrier)) /
            (n_noncarrier / (n_noncarrier_total - n_noncarrier))
    ) %>%
    ungroup() %>%
    arrange(fisher_p)

# Separate reportable genes BEFORE computing FDR
gene_fisher_reportable <- gene_fisher %>%
    filter(is.finite(OR_crude), n_noncarrier >= 1) %>%
    mutate(fdr = p.adjust(fisher_p, method = "BH"))

# Keep full table with FDR for saving, but compute separately
gene_fisher <- gene_fisher %>%
    mutate(fdr = p.adjust(fisher_p, method = "BH"))

cat("\nFisher results (top 10):\n")
print(head(gene_fisher_reportable %>% dplyr::select(Gene, n_carrier, n_noncarrier,
                                                    pct_carrier, pct_noncarrier,
                                                    OR_crude, fisher_p, fdr), 10))

# ============================================================
# GENE PREVALENCE PLOT - STRATIFIED
# ============================================================
# shared y scale for both gene figures
scale_y_gene <- scale_y_continuous(
    labels = function(x) paste0(x, "%"),
    expand = expansion(mult = c(0, 0.1))
)

# shared theme additions for both
theme_gene <- list(
    theme_ch,
    theme(
        axis.text.x        = element_text(angle = 90, hjust = 1),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(color = "grey92", linewidth = 0.3)
    )
)

# ── Overall ──
bar_data_overall <- gene_freq %>%
    filter(Gene %in% top_genes) %>%
    mutate(
        pct_total = n_total / (n_carrier_total + n_noncarrier_total) * 100,
        Gene      = factor(Gene, levels = top_genes)
    )

fig_gene_overall <- ggplot(bar_data_overall, aes(x = Gene, y = pct_total)) +
    geom_col(fill = "#546E7A", width = 0.7, alpha = 0.9) +
    scale_y_gene +
    labs(title = "CHIP gene prevalence",
         x = NULL, y = "Prevalence (%)") +
    theme_gene

ggsave(file.path("ch", "figures", "fig_gene_prevalence_bar.pdf"),
       fig_gene_overall, width = 6, height = 5)

# ── Stratified ──
bar_data_strat <- gene_freq %>%
    filter(Gene %in% top_genes) %>%
    pivot_longer(cols = c(pct_carrier, pct_noncarrier),
                 names_to  = "group",
                 values_to = "pct") %>%
    mutate(
        group = recode(group,
                       pct_carrier    = "gBRCA1/2 Carrier",
                       pct_noncarrier = "Non-carrier"),
        Gene  = factor(Gene, levels = top_genes)
    )

fig_gene_strat <- ggplot(bar_data_strat,
                         aes(x = Gene, y = pct, fill = group)) +
    geom_col(position = "dodge", width = 0.7, alpha = 0.9) +
    scale_fill_manual(
        values = c("gBRCA1/2 Carrier" = "#C2185B",
                   "Non-carrier"       = "#546E7A"),
        name = NULL
    ) +
    scale_y_gene +
    labs(title = "CHIP gene prevalence by carrier status",
         x = NULL, y = "Prevalence (%)") +
    theme_gene +
    theme(legend.position = "none")

ggsave(file.path("ch", "figures", "fig_gene_prevalence_bar_stratified_noleg.pdf"),
       fig_gene_strat, width = 6, height = 5)

fig_gene_strat <- ggplot(bar_data_strat,
                         aes(x = Gene, y = pct, fill = group)) +
    geom_col(position = "dodge", width = 0.7, alpha = 0.9) +
    scale_fill_manual(
        values = c("gBRCA1/2 Carrier" = "#C2185B",
                   "Non-carrier"       = "#546E7A"),
        name = NULL
    ) +
    scale_y_gene +
    labs(title = "CHIP gene prevalence by carrier status",
         x = NULL, y = "Prevalence (%)") +
    theme_gene +
    theme(legend.position = "bottom")

ggsave(file.path("ch", "figures", "fig_gene_prevalence_bar_stratified.pdf"),
       fig_gene_strat, width = 6, height = 5)

# ============================================================
# SUPP TABLE 5
# ============================================================
variant_table <- vars %>%
    group_by(variant_id) %>%
    summarise(
        Count   = n_distinct(Sample.ID),
        Gene            = first(Gene),
        Chr             = first(Chr),
        Start           = first(Start),
        REF             = first(REF),
        ALT             = first(ALT),
        HGVSc           = first(HGVSc),
        HGVSp           = first(HGVSp),
        Class   = first(Variant.Class),
        LoF_level = first(Variant.LoF_level),
        Mean_VAF        = round(mean(Sample.AltFrac,   na.rm = TRUE), 3),

        .groups = "drop"
    ) %>%
    arrange(desc(Count))
write.csv(variant_table,
          file.path("ch", "data", "supp5_final_variant_table.csv"),
          row.names = FALSE)

# ============================================================
# COUNT BY INDIVIDIAULS
# ============================================================
n_total <- nrow(cov)

chip_count_dist <- cov %>%
    filter(CHIP_Binary) %>%
    count(CHIP_Count) %>%
    mutate(
        CHIP_Count = factor(CHIP_Count),
        pct        = round(100 * n / n_total, 1),
        label      = paste0(n, "\n(", pct, "%)")
    )

write_xlsx(chip_count_dist,
           file.path("ch", "data", "supp_ind_variant_count_table.xlsx"))

fig_chip_count <- ggplot(chip_count_dist, aes(x = CHIP_Count, y = n)) +
    geom_col(fill = "#546E7A", width = 0.7, alpha = 0.9) +
    geom_text(aes(label = label), vjust = -0.4, size = 3.5, color = "grey30") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
    labs(
        title = "Distribution of CHIP mutation count per individual",
        x     = "Number of CHIP mutations",
        y     = "Number of individuals"
    ) +
    theme_ch


ggsave(file.path("ch", "figures", "fig_chip_count_dist.pdf"),
       fig_chip_count, width = 6, height = 5)
