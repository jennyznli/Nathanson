# ========================
# Age association for artifacts
# ========================
library(here)
setwd("/Users/jennyzli/Documents/Nathanson")
source(here("R", "config.R"))
library(data.table, quietly = T)
library(ggVennDiagram)
library(UpSetR)

# ========================
# READ IN DATA
# ========================
all_ch <- read_excel(file.path("ch", "data", "ch_seq_wl_art_manual_vars.xlsx"))
dim(all_ch)
# 477

cov <- read.csv(file.path("ch", "data", "pmbb_brca12_cov_df.csv"), row.names = 1)
all_ids <- read.csv(file.path("ch", "data", "ch_psm_matched4_case_control_ids.csv"))$x
cov <- cov %>% filter(person_id %in% all_ids)
table(cov$Batch)
dim(cov)
# 3004

# ========================
# MINAD THRESHOLDS
# ========================
for (thresh in MIN_AD_THRESHOLDS) {
    all_ch[[paste0("minad", thresh)]] <- all_ch$Sample.AltDepth < thresh
}

write_xlsx(all_ch, file.path("ch", "data", "ch_seq_wl_art_manual_minad_vars.xlsx"))

# ========================
# DEFINE FILTER CONFIGURATIONS
# ========================
extra_covs <- "Smoke_History + Sequenced_gender + PC1 + PC2 + PC3 + PC4 + PC5 + PC6"

minad_configs <- list(
    list(name = "minad3_full_covs", flags = c("minad3"), formula_extra = extra_covs),
    list(name = "minad4_full_covs", flags = c("minad4"), formula_extra = extra_covs),
    list(name = "minad5_full_covs", flags = c("minad5"), formula_extra = extra_covs),
    list(name = "minad6_full_covs", flags = c("minad6"), formula_extra = extra_covs),
    list(name = "minad7_full_covs", flags = c("minad7"), formula_extra = extra_covs),
    list(name = "minad8_full_covs", flags = c("minad8"), formula_extra = extra_covs)
)

# ========================
# AGE ASSOCIATION TESTER
# ========================
test_age_association <- function(config, all_ch, cov,
                                 default_formula_extra = NULL,
                                 stratify_by = NULL) {
    if (!is.null(stratify_by)) {
        strata_levels <- sort(unique(na.omit(all_ch[[stratify_by]])))
        strata_list   <- setNames(as.list(strata_levels), strata_levels)
    } else {
        strata_list <- list("all" = NULL)
    }

    results <- lapply(config, function(cfg) {
        nm        <- cfg$name
        flags     <- cfg$flags
        fml_extra <- if (!is.null(cfg$formula_extra)) cfg$formula_extra else default_formula_extra

        lapply(names(strata_list), function(stratum) {
            stratum_val <- strata_list[[stratum]]

            ch_filtered <- all_ch
            if (!is.null(stratum_val)) {
                ch_filtered <- ch_filtered %>% filter(.data[[stratify_by]] == stratum_val)
            }
            if (length(flags) > 0) {
                remove_mask <- Reduce(`|`, lapply(flags, function(f) {
                    col <- ch_filtered[[f]]
                    if (is.logical(col)) col else as.logical(col)
                }))
                remove_mask[is.na(remove_mask)] <- FALSE
                ch_filtered <- ch_filtered[!remove_mask, ]
            }

            carriers <- ch_filtered %>%
                distinct(Sample.ID) %>%
                mutate(carrier = 1L)

            df <- cov %>%
                left_join(carriers, by = c("person_id" = "Sample.ID")) %>%
                mutate(carrier = as.integer(ifelse(is.na(carrier), 0L, carrier)))

            n_carriers <- sum(df$carrier)
            n_total    <- nrow(df)

            cat(sprintf("[%s | stratum: %s] %d carriers / %d total\n",
                        nm, stratum, n_carriers, n_total))

            base_rhs <- if (!is.null(stratify_by) && stratify_by == "Batch") {
                "Sample_age"
            } else {
                "Sample_age + as.factor(Batch)"
            }
            rhs     <- if (!is.null(fml_extra)) paste(base_rhs, "+", fml_extra) else base_rhs
            formula <- as.formula(paste("carrier ~", rhs))

            tryCatch({
                m <- withCallingHandlers(
                    glm(formula, data = df, family = binomial()),
                    warning = function(w) {
                        if (grepl("fitted probabilities numerically 0 or 1", conditionMessage(w)))
                            invokeRestart("muffleWarning")
                    }
                )
                separation <- any(fitted(m) %in% c(0, 1))  | !m$converged
                s  <- coef(summary(m))
                ci <- confint.default(m)
                data.frame(
                    name       = nm, stratum = stratum,
                    flags      = paste(flags, collapse = ", "),
                    n_carriers = n_carriers, n_total = n_total,
                    OR         = exp(coef(m)["Sample_age"]),
                    CI_lo      = exp(ci["Sample_age", 1]),
                    CI_hi      = exp(ci["Sample_age", 2]),
                    p_age      = s["Sample_age", "Pr(>|z|)"],
                    separation = separation,
                    stringsAsFactors = FALSE
                )
            }, error = function(e) {
                warning(sprintf("Model failed for [%s | %s]: %s", nm, stratum, conditionMessage(e)))
                data.frame(name = nm, stratum = stratum,
                           flags = paste(flags, collapse = ", "),
                           n_carriers = n_carriers, n_total = n_total,
                           OR = NA_real_, CI_lo = NA_real_, CI_hi = NA_real_,
                           p_age = NA_real_, separation = NA,
                           stringsAsFactors = FALSE)
            })
        }) %>% bind_rows()
    }) %>% bind_rows()

    results
}

# ========================
# PLOT FUNCTION
# ========================
plot_age_association <- function(results, title = "Age association under different QC filters",
                                 n_total = NULL) {
    stratified <- length(unique(results$stratum)) > 1

    plot_df <- results %>%
        filter(!is.na(OR)) %>%
        mutate(
            name    = factor(name, levels = rev(unique(name))),
            sig     = !is.na(p_age) & p_age < 0.05,
            stratum = factor(stratum)
        )

    subtitle <- if (!is.null(n_total)) sprintf("n = %d total individuals", n_total) else NULL

    p <- ggplot(plot_df, aes(x = OR, y = name, color = if (stratified) stratum else sig)) +
        geom_vline(xintercept = 1, linetype = "dashed", color = "grey60", linewidth = 0.4) +
        geom_errorbar(aes(xmin = CI_lo, xmax = CI_hi, group = if (stratified) stratum else NULL),
                      width = 0.25, linewidth = 0.6, orientation = "y",
                      position = if (stratified) position_dodge(width = 0.6) else "identity") +
        geom_point(size = 3,
                   position = if (stratified) position_dodge(width = 0.6) else "identity") +
        scale_x_continuous(name = "Odds ratio per year of age (95% CI)",
                           breaks = scales::pretty_breaks(n = 6)) +
        labs(title = title, subtitle = subtitle, y = NULL) +
        theme_minimal(base_size = 11) +
        theme(
            plot.title         = element_text(face = "bold", hjust = 0.5),
            plot.subtitle      = element_text(hjust = 0.5, size = 9, color = "grey40"),
            panel.grid.minor   = element_blank(),
            panel.grid.major.y = element_blank(),
            legend.position    = "bottom"
        )

    if (stratified) {
        p <- p +
            scale_color_brewer(palette = "Set2") +
            guides(color = guide_legend(title = "stratum"))
    } else {
        p <- p +
            scale_color_manual(
                values = c("FALSE" = "grey50", "TRUE" = "#D85A30"),
                labels = c("FALSE" = "p \u2265 0.05", "TRUE" = "p < 0.05"),
                name   = NULL
            )
    }
    p
}

# ========================
# RUN
# ========================
minad_res <- test_age_association(minad_configs,  all_ch, cov)

minad_vaf_res <- test_age_association(minad_configs,  all_ch, cov, stratify_by = "VAF_Strata")

minad_freeze_res  <- test_age_association(minad_configs,  all_ch, cov, stratify_by = "Batch")

# ========================
# SAVE RESULTS
# ========================
write_xlsx(minad_res,         file.path("ch", "data", "ch_minad_age_sensitivity.xlsx"))

write_xlsx(minad_vaf_res,       file.path("ch", "data", "ch_minad_age_sensitivity_vaf.xlsx"))

write_xlsx(minad_freeze_res,    file.path("ch", "data", "ch_minad_age_sensitivity_batch.xlsx"))

# ========================
# PLOTS
# ========================
save_assoc_plot <- function(results, filename, title, n_total) {
    p <- plot_age_association(results, title = title, n_total = n_total)
    n_names   <- n_distinct(results$name)
    n_strata  <- n_distinct(results$stratum)
    width  <- if (n_strata > 1) 5 + 1.5 * n_strata else 7
    height <- 0.4 * n_names + 2
    ggsave(filename, p, width = width, height = height, limitsize = FALSE)
}

save_assoc_plot(minad_res,         file.path("ch", "figures", "qc_minad_age_sensitivity2.pdf"),
                "Age association — minAD thresholds", nrow(cov))

save_assoc_plot(minad_vaf_res,     file.path("ch", "figures", "qc_minad_age_sensitivity_vaf2.pdf"),
                "Age association — minAD by VAF stratum", nrow(cov))

save_assoc_plot(minad_freeze_res,  file.path("ch", "figures", "qc_minad_age_sensitivity_batch2.pdf"),
                "Age association — minAD by batch", nrow(cov))

# ========================
# COMPARE VAF 10< vs. ALL VAF
# ========================
vaf_compare_res <- bind_rows(
    test_age_association(minad_configs, all_ch, cov) %>% mutate(vaf_group = "All VAF"),
    test_age_association(minad_configs, all_ch %>% filter(Sample.AltFrac >= 0.10), cov) %>%
        mutate(vaf_group = "VAF ≥ 10%")
)

plot_df <- vaf_compare_res %>%
    filter(!is.na(OR)) %>%
    mutate(
        name      = factor(name, levels = rev(unique(name))),
        sig       = !is.na(p_age) & p_age < 0.05,
        vaf_group = factor(vaf_group, levels = c("All VAF", "VAF ≥ 10%"))
    )

p_vaf_compare <- ggplot(plot_df, aes(x = OR, y = name, color = vaf_group)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey60", linewidth = 0.4) +
    geom_errorbar(aes(xmin = CI_lo, xmax = CI_hi, group = vaf_group),
                  width = 0.25, linewidth = 0.6, orientation = "y",
                  position = position_dodge(width = 0.6)) +
    geom_point(aes(shape = sig), size = 3, position = position_dodge(width = 0.6)) +
    scale_color_manual(
        values = c("All VAF" = "#378ADD", "VAF ≥ 10%" = "#D85A30"),
        name   = NULL
    ) +
    scale_shape_manual(
        values = c("FALSE" = 16, "TRUE" = 8),
        labels = c("FALSE" = "p \u2265 0.05", "TRUE" = "p < 0.05"),
        name   = NULL
    ) +
    scale_x_continuous(
        name   = "Odds ratio per year of age (95% CI)",
        breaks = scales::pretty_breaks(n = 6)
    ) +
    labs(
        title    = "Age association by VAF threshold",
        subtitle = sprintf("n = %d total individuals", nrow(cov)),
        y        = NULL
    ) +
    theme_minimal(base_size = 11) +
    theme(
        plot.title         = element_text(face = "bold", hjust = 0.5),
        plot.subtitle      = element_text(hjust = 0.5, size = 9, color = "grey40"),
        panel.grid.minor   = element_blank(),
        panel.grid.major.y = element_blank(),
        legend.position    = "bottom"
    )

ggsave(file.path("ch", "figures", "qc_filter_age_cov_compare2.pdf"),
       p_vaf_compare,
       width  = 7,
       height = 0.4 * n_distinct(plot_df$name) + 2,
       limitsize = FALSE)


# ========================
# FINAL FILTER
# ========================
all_ch3 <- all_ch
dim(all_ch3)
# 477
write_xlsx(all_ch3, file.path("ch", "data", "ch_seq_wl_art_minad3_vars.xlsx"))

all_ch4 <- all_ch %>% filter(!minad4)
dim(all_ch4)
# 217
write_xlsx(all_ch4, file.path("ch", "data", "ch_seq_wl_art_minad4_vars.xlsx"))

all_ch5 <- all_ch %>% filter(!minad5)
dim(all_ch5)
# 128
write_xlsx(all_ch5, file.path("ch", "data", "ch_seq_wl_art_minad5_vars.xlsx"))

all_ch6 <- all_ch %>% filter(!minad6)
dim(all_ch6)
# 93
write_xlsx(all_ch6, file.path("ch", "data", "ch_seq_wl_art_minad6_vars.xlsx"))
