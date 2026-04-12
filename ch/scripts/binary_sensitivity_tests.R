# ============================================================
# CHIP BINARY ASSOCIATION - MINAD SENSITIVITY ANALYSIS
# ALL SEXES
# ============================================================
library(here)
setwd("/Users/jennyzli/Documents/Nathanson")
source(here("R", "config.R"))
source(here("ch", "scripts", "binary_sensitivity.R"))

library(data.table)
library(MatchIt)
library(survival)
library(splines)
library(sandwich)
library(lmtest)
library(pROC)
library(logistf)
library(lmtest)
library(sandwich)

# ============================================================
# READ DATA
# ============================================================
m.out4 <- readRDS(file.path("ch", "data", "ch_psm_matched4.rds"))
m.data <- match.data(m.out4)
cov_base <- m.data %>%
    mutate(
        Sample_age   = as.numeric(Sample_age),
        Sample_age2  = Sample_age^2,
        subclass_num = as.numeric(as.factor(subclass))
    )
dim(cov_base)
# 3004

vars <- read_excel(file.path("ch", "data", "ch_seq_wl_art_minad4_vars.xlsx")) %>%
    filter(Sample.ID %in% cov_base$person_id)
dim(vars)
# 217

cov_base <- cov_base %>%
    mutate(
        CHIP_Binary = person_id %in% vars$Sample.ID,
        CHIP_Count  = sapply(person_id, function(id) sum(vars$Sample.ID == id))
    )

write_xlsx(cov_base, file.path("ch", "data", "chip_cov_minad4.xlsx"))

# ============================================================
# RUN LOOPS
# ============================================================
cov_all_s12 <- cov_base
cov_f_s12 <- cov_base %>% filter(Sequenced_gender == "Female")

res_all_s12 <- run_chip_analysis(
    cov_base      = cov_all_s12,
    cohort_label  = "all_s12",
    females_only = FALSE
)

res_f_s12 <- run_chip_analysis(
    cov_base      = cov_f_s12,
    cohort_label  = "female_s12",
    females_only = TRUE
)

res_all_age_s12 <- run_chip_age_analysis(
    cov_base      = cov_all_s12,
    cohort_label  = "all_s12",
    females_only = FALSE
)

write_xlsx(res_all_s12,    file.path("ch", "data", "chip_binary_all_s12_minad_sensitivity.xlsx"))
write_xlsx(res_f_s12,    file.path("ch", "data", "chip_binary_female_s12_minad_sensitivity.xlsx"))
write_xlsx(res_all_age_s12, file.path("ch", "data", "chip_binary_all_s12_age_sensitivity.xlsx"))

# ============================================================
# FIGURE 2 - FINAL
# ============================================================
fit_b12 <- glm(CHIP_Binary ~ BRCA12_Case + Sample_age + I(Sample_age^2) + Sequenced_gender +
                   PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + Smoke_History + Batch,
               data = cov_all_s12, family = quasibinomial())

fit_b1 <- glm(CHIP_Binary ~ BRCA1_Case + Sample_age + I(Sample_age^2) + Sequenced_gender +
                  PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + Smoke_History + Batch,
              data = cov_all_s12, family = quasibinomial())

fit_b2 <- glm(CHIP_Binary ~ BRCA2_Case + Sample_age + I(Sample_age^2) + Sequenced_gender +
                  PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + Smoke_History + Batch,
              data = cov_all_s12, family = quasibinomial())

# Helper to extract OR row
extract_or <- function(fit, case_col, label) {
    coef_table <- summary(fit)$coefficients
    data.frame(
        model = label,
        OR    = exp(coef_table[case_col, "Estimate"]),
        CI_lo = exp(coef_table[case_col, "Estimate"] - 1.96 * coef_table[case_col, "Std. Error"]),
        CI_hi = exp(coef_table[case_col, "Estimate"] + 1.96 * coef_table[case_col, "Std. Error"]),
        p     = coef_table[case_col, "Pr(>|t|)"]
    )
}

# Bind all three
final_df <- bind_rows(
    extract_or(fit_b12, "BRCA12_Case", "BRCA1/2"),
    extract_or(fit_b1,  "BRCA1_Case",  "BRCA1"),
    extract_or(fit_b2,  "BRCA2_Case",  "BRCA2")
) %>%
    mutate(
        label      = factor(model, levels = rev(c("BRCA1/2", "BRCA1", "BRCA2"))),
        p_fmt      = ifelse(p < 0.001, "p<0.001", sprintf("p=%.3f", p)),
        annotation = sprintf("%.2f (%.2f, %.2f), %s", OR, CI_lo, CI_hi, p_fmt)
    )

fig_final <- ggplot(final_df,
                    aes(x = OR, y = label, xmin = CI_lo, xmax = CI_hi)) +
    geom_vline(xintercept = 1, linetype = "dashed",
               color = "grey50", linewidth = 0.5) +
    geom_errorbarh(height = 0.15, linewidth = 0.7, color = "#C2185B") +
    geom_point(size = 4, color = "#C2185B") +
    geom_text(aes(x = 4, label = annotation),
              hjust = 0, size = 4, color = "grey30") +
    scale_x_log10(
        breaks = c(0.5, 1, 2, 4),
        labels = c("0.5", "1", "2", "4"),
        limits = c(0.4, 9)
    ) +
    coord_cartesian(clip = "off") +
    labs(
        x        = "Log Odds Ratio",
        y        = NULL
    ) +
    theme_minimal(base_size = 13) +
    theme(
        plot.title         = element_text(face = "bold", hjust = 0, size = 13),
        plot.subtitle      = element_text(color = "grey40", hjust = 0, size = 12),
        plot.margin        = margin(5, 200, 5, 5),
        panel.grid.minor   = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(color = "grey92", linewidth = 0.3),
        axis.text.y        = element_text(size = 12),
        axis.text.x        = element_text(size = 12)
    )

ggsave(
    file.path("ch", "figures", "fig_final_brca_chip_agesq_minad4.pdf"),
    fig_final, width = 8, height = 1.5
)

# ============================================================
# FREEZE 2.0 vs. 3.0
# ============================================================
cov_b1 <- cov_all_s12 %>% filter(Batch == "1")
cov_b2 <- cov_all_s12 %>% filter(Batch == "2")
dim(cov_b1)
# 760
dim(cov_b2)
# 2244

run_three_genes <- function(data, include_sex) {
    sex_cov <- if (include_sex) "+ Sequenced_gender" else ""

    fits <- list(
        list(formula = paste("CHIP_Binary ~ BRCA12_Case + Sample_age + I(Sample_age^2)", sex_cov,
                             "+ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + Smoke_History"),
             case_col = "BRCA12_Case", label = "BRCA1/2"),
        list(formula = paste("CHIP_Binary ~ BRCA1_Case + Sample_age + I(Sample_age^2)", sex_cov,
                             "+ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + Smoke_History"),
             case_col = "BRCA1_Case",  label = "BRCA1"),
        list(formula = paste("CHIP_Binary ~ BRCA2_Case + Sample_age + I(Sample_age^2)", sex_cov,
                             "+ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + Smoke_History"),
             case_col = "BRCA2_Case",  label = "BRCA2")
    )

    bind_rows(lapply(fits, function(f) {
        fit <- glm(as.formula(f$formula), data = data, family = quasibinomial())
        coef_table <- summary(fit)$coefficients
        data.frame(
            model = f$label,
            OR    = exp(coef_table[f$case_col, "Estimate"]),
            CI_lo = exp(coef_table[f$case_col, "Estimate"] - 1.96 * coef_table[f$case_col, "Std. Error"]),
            CI_hi = exp(coef_table[f$case_col, "Estimate"] + 1.96 * coef_table[f$case_col, "Std. Error"]),
            p     = coef_table[f$case_col, "Pr(>|t|)"]
        )
    }))
}

# Run both cohorts
final_df <- bind_rows(
    run_three_genes(cov_b1, include_sex = TRUE) %>% mutate(cohort = "Freeze 2.0"),
    run_three_genes(cov_b2,  include_sex = FALSE) %>% mutate(cohort = "Freeze 3.0"),
    run_three_genes(cov_all_s12,  include_sex = FALSE) %>% mutate(cohort = "All")
) %>%
    mutate(
        model      = factor(model, levels = c("BRCA1/2", "BRCA1", "BRCA2")),
        cohort     = factor(cohort, levels = c("Freeze 2.0", "Freeze 3.0", "All")),
        p_fmt      = ifelse(p < 0.001, "p<0.001", sprintf("p=%.3f", p)),
        annotation = sprintf("%.2f (%.2f, %.2f), %s", OR, CI_lo, CI_hi, p_fmt),
        # y axis: gene + cohort for each row
        label      = factor(
            paste0(model, " (", cohort, ")"),
            levels = rev(paste0(
                rep(c("BRCA1/2", "BRCA1", "BRCA2"), each = 3),
                " (",
                rep(c("Freeze 2.0", "Freeze 3.0", "All"), 3),
                ")"
            ))
        )
    )

fig_freeze <- ggplot(final_df,
                  aes(x = OR, y = label, xmin = CI_lo, xmax = CI_hi,
                      color = cohort)) +
    geom_vline(xintercept = 1, linetype = "dashed",
               color = "grey50", linewidth = 0.5) +
    geom_errorbarh(height = 0.15, linewidth = 0.7) +
    geom_point(size = 4) +
    geom_text(aes(x = 6, label = annotation),
              hjust = 0, size = 3.0, color = "grey30") +
    scale_x_log10(
        breaks = c(0.5, 1, 2, 4),
        labels = c("0.5", "1", "2", "4"),
        limits = c(0.4, 10)
    ) +
    coord_cartesian(clip = "off") +
    labs(
        x        = "Odds ratio (log scale)",
        y        = NULL
    ) +
    theme_minimal(base_size = 13) +
    theme(
        plot.title         = element_text(face = "bold", hjust = 0, size = 13),
        plot.subtitle      = element_text(color = "grey40", hjust = 0, size = 10),
        plot.margin        = margin(5, 220, 5, 5),
        panel.grid.minor   = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(color = "grey92", linewidth = 0.3),
        axis.text.y        = element_text(size = 10),
        axis.text.x        = element_text(size = 10),
        legend.position    = "bottom"
    )

ggsave(
    file.path("ch", "figures", "fig_freeze_comparison_minad4.pdf"),
    fig_freeze, width = 9, height = 4.5
)


# ============================================================
# MINAD
# ============================================================
vars3 <- read_excel(file.path("ch", "data", "ch_seq_wl_art_minad3_vars.xlsx")) %>%
    filter(Sample.ID %in% cov_base$person_id)
vars4 <- read_excel(file.path("ch", "data", "ch_seq_wl_art_minad4_vars.xlsx")) %>%
    filter(Sample.ID %in% cov_base$person_id)
vars5 <- read_excel(file.path("ch", "data", "ch_seq_wl_art_minad5_vars.xlsx")) %>%
    filter(Sample.ID %in% cov_base$person_id)

cov3 <- cov_base %>%
    mutate(
        CHIP_Binary = person_id %in% vars3$Sample.ID,
        CHIP_Count  = sapply(person_id, function(id) sum(vars3$Sample.ID == id))
    )

cov4 <- cov_base %>%
    mutate(
        CHIP_Binary = person_id %in% vars4$Sample.ID,
        CHIP_Count  = sapply(person_id, function(id) sum(vars4$Sample.ID == id))
    )
cov5 <- cov_base %>%
    mutate(
        CHIP_Binary = person_id %in% vars5$Sample.ID,
        CHIP_Count  = sapply(person_id, function(id) sum(vars5$Sample.ID == id))
    )

run_three_genes <- function(data, include_sex) {
    sex_cov <- if (include_sex) "+ Sequenced_gender" else ""

    fits <- list(
        list(formula = paste("CHIP_Binary ~ BRCA12_Case + Sample_age + I(Sample_age^2)", sex_cov,
                             "+ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + Smoke_History  + Batch"),
             case_col = "BRCA12_Case", label = "BRCA1/2"),
        list(formula = paste("CHIP_Binary ~ BRCA1_Case + Sample_age + I(Sample_age^2)", sex_cov,
                             "+ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + Smoke_History  + Batch"),
             case_col = "BRCA1_Case",  label = "BRCA1"),
        list(formula = paste("CHIP_Binary ~ BRCA2_Case + Sample_age + I(Sample_age^2)", sex_cov,
                             "+ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + Smoke_History  + Batch"),
             case_col = "BRCA2_Case",  label = "BRCA2")
    )

    bind_rows(lapply(fits, function(f) {
        fit <- glm(as.formula(f$formula), data = data, family = quasibinomial())
        coef_table <- summary(fit)$coefficients
        data.frame(
            model = f$label,
            OR    = exp(coef_table[f$case_col, "Estimate"]),
            CI_lo = exp(coef_table[f$case_col, "Estimate"] - 1.96 * coef_table[f$case_col, "Std. Error"]),
            CI_hi = exp(coef_table[f$case_col, "Estimate"] + 1.96 * coef_table[f$case_col, "Std. Error"]),
            p     = coef_table[f$case_col, "Pr(>|t|)"]
        )
    }))
}

final_df <- bind_rows(
    run_three_genes(cov3, include_sex = TRUE) %>% mutate(cohort = "minAD3"),
    run_three_genes(cov4,  include_sex = FALSE) %>% mutate(cohort = "minAD4"),
    run_three_genes(cov5,  include_sex = FALSE) %>% mutate(cohort = "minAD5")
) %>%
    mutate(
        model      = factor(model, levels = c("BRCA1/2", "BRCA1", "BRCA2")),
        cohort     = factor(cohort, levels = c("minAD3", "minAD4", "minAD5")),
        p_fmt      = ifelse(p < 0.001, "p<0.001", sprintf("p=%.3f", p)),
        annotation = sprintf("%.2f (%.2f, %.2f), %s", OR, CI_lo, CI_hi, p_fmt),
        # y axis: gene + cohort for each row
        label      = factor(
            paste0(model, " (", cohort, ")"),
            levels = rev(paste0(
                rep(c("BRCA1/2", "BRCA1", "BRCA2"), each = 3),
                " (",
                rep(c("minAD3", "minAD4", "minAD5"), 3),
                ")"
            ))
        )
    )

fig_minad <- ggplot(final_df,
                  aes(x = OR, y = label, xmin = CI_lo, xmax = CI_hi,
                      color = cohort)) +
    geom_vline(xintercept = 1, linetype = "dashed",
               color = "grey50", linewidth = 0.5) +
    geom_errorbarh(height = 0.15, linewidth = 0.7) +
    geom_point(size = 4) +
    geom_text(aes(x = 6, label = annotation),
              hjust = 0, size = 3.0, color = "grey30") +
    # scale_color_manual(
    #     values = c("All" = "darkgray", "Female" = "#C2185B"),
    #     name   = "Cohort"
    # ) +
    scale_x_log10(
        breaks = c(0.5, 1, 2, 4),
        labels = c("0.5", "1", "2", "4"),
        limits = c(0.4, 10)
    ) +
    coord_cartesian(clip = "off") +
    labs(
        x        = "Odds ratio (log scale)",
        y        = NULL
    ) +
    theme_minimal(base_size = 13) +
    theme(
        plot.title         = element_text(face = "bold", hjust = 0, size = 13),
        plot.subtitle      = element_text(color = "grey40", hjust = 0, size = 10),
        plot.margin        = margin(5, 220, 5, 5),
        panel.grid.minor   = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(color = "grey92", linewidth = 0.3),
        axis.text.y        = element_text(size = 10),
        axis.text.x        = element_text(size = 10),
        legend.position    = "bottom"
    )

ggsave(
    file.path("ch", "figures", "fig_minad_comparison.pdf"),
    fig_minad, width = 9, height = 4.5
)



# ============================================================
# FEMALE VS ALL SEXES
# ============================================================
# Female-only cohort — filter and drop Sequenced_gender
cov_female <- cov_all_s12 %>% filter(Sequenced_gender == "Female")  # adjust 0/1 to your coding

# Helper to run all three genes for a given dataset
run_three_genes <- function(data, include_sex) {
    sex_cov <- if (include_sex) "+ Sequenced_gender" else ""

    fits <- list(
        list(formula = paste("CHIP_Binary ~ BRCA12_Case + Sample_age + I(Sample_age^2)", sex_cov,
                             "+ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + Smoke_History + Batch"),
             case_col = "BRCA12_Case", label = "BRCA1/2"),
        list(formula = paste("CHIP_Binary ~ BRCA1_Case + Sample_age + I(Sample_age^2)", sex_cov,
                             "+ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + Smoke_History + Batch"),
             case_col = "BRCA1_Case",  label = "BRCA1"),
        list(formula = paste("CHIP_Binary ~ BRCA2_Case + Sample_age + I(Sample_age^2)", sex_cov,
                             "+ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + Smoke_History + Batch"),
             case_col = "BRCA2_Case",  label = "BRCA2")
    )

    bind_rows(lapply(fits, function(f) {
        fit <- glm(as.formula(f$formula), data = data, family = quasibinomial())
        coef_table <- summary(fit)$coefficients
        data.frame(
            model = f$label,
            OR    = exp(coef_table[f$case_col, "Estimate"]),
            CI_lo = exp(coef_table[f$case_col, "Estimate"] - 1.96 * coef_table[f$case_col, "Std. Error"]),
            CI_hi = exp(coef_table[f$case_col, "Estimate"] + 1.96 * coef_table[f$case_col, "Std. Error"]),
            p     = coef_table[f$case_col, "Pr(>|t|)"]
        )
    }))
}

# Run both cohorts
final_df <- bind_rows(
    run_three_genes(cov_all_s12, include_sex = TRUE) %>% mutate(cohort = "All"),
    run_three_genes(cov_female,  include_sex = FALSE) %>% mutate(cohort = "Female")
) %>%
    mutate(
        model      = factor(model, levels = c("BRCA1/2", "BRCA1", "BRCA2")),
        cohort     = factor(cohort, levels = c("All", "Female")),
        p_fmt      = ifelse(p < 0.001, "p<0.001", sprintf("p=%.3f", p)),
        annotation = sprintf("%.2f (%.2f, %.2f), %s", OR, CI_lo, CI_hi, p_fmt),
        # y axis: gene + cohort for each row
        label      = factor(
            paste0(model, " (", cohort, ")"),
            levels = rev(paste0(
                rep(c("BRCA1/2", "BRCA1", "BRCA2"), each = 2),
                " (",
                rep(c("All", "Female"), 3),
                ")"
            ))
        )
    )

fig_sex <- ggplot(final_df,
                  aes(x = OR, y = label, xmin = CI_lo, xmax = CI_hi,
                      color = cohort)) +
    geom_vline(xintercept = 1, linetype = "dashed",
               color = "grey50", linewidth = 0.5) +
    geom_errorbarh(height = 0.15, linewidth = 0.7) +
    geom_point(size = 4) +
    geom_text(aes(x = 6, label = annotation),
              hjust = 0, size = 3.0, color = "grey30") +
    scale_color_manual(
        values = c("All" = "darkgray", "Female" = "#C2185B"),
        name   = "Cohort"
    ) +
    scale_x_log10(
        breaks = c(0.5, 1, 2, 4),
        labels = c("0.5", "1", "2", "4"),
        limits = c(0.4, 10)
    ) +
    coord_cartesian(clip = "off") +
    labs(
        x        = "Odds ratio (log scale)",
        y        = NULL
    ) +
    theme_minimal(base_size = 13) +
    theme(
        plot.title         = element_text(face = "bold", hjust = 0, size = 13),
        plot.subtitle      = element_text(color = "grey40", hjust = 0, size = 10),
        plot.margin        = margin(5, 220, 5, 5),
        panel.grid.minor   = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(color = "grey92", linewidth = 0.3),
        axis.text.y        = element_text(size = 10),
        axis.text.x        = element_text(size = 10),
        legend.position    = "bottom"
    )

ggsave(
    file.path("ch", "figures", "fig_sex_comparison_minad4.pdf"),
    fig_sex, width = 9, height = 3.5
)

# ============================================================
# INTERACTION
# ============================================================
fit_age_b12 <- glm(CHIP_Binary ~ BRCA12_Case * Sample_age + I(Sample_age^2) + Sequenced_gender +
                   PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + Smoke_History + Batch,
               data = cov_all_s12, family = quasibinomial())

fit_age_b1 <- glm(CHIP_Binary ~ BRCA1_Case * Sample_age + I(Sample_age^2) + Sequenced_gender +
                  PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + Smoke_History + Batch,
              data = cov_all_s12, family = quasibinomial())

fit_age_b2 <- glm(CHIP_Binary ~ BRCA2_Case * Sample_age + I(Sample_age^2) + Sequenced_gender +
                  PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + Smoke_History + Batch,
              data = cov_all_s12, family = quasibinomial())

