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
all_ch <- read_excel(file.path("ch", "data", "ch_seq_wl_flags_vars.xlsx"))
dim(all_ch)
# 481

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

write_xlsx(all_ch, file.path("ch", "data", "ch_seq_wl_flags_minad_vars.xlsx"))

# ========================
# DEFINE FILTER CONFIGURATIONS
# ========================
full_qc <- c("flag_freeze_enriched", "flag_not_age_associated",
             "flag_germline_ind", "flag_cluster")
extra_covs <- "Smoke_History + Sequenced_gender + PC1 + PC2 + PC3 + PC4 + PC5 + PC6"

filter_configs <- list(
    list(name = "no_filter",               flags = character(0)),
    list(name = "covs",                    flags = character(0), formula_extra = extra_covs),
    list(name = "flag_freeze_enriched",    flags = c("flag_freeze_enriched")),
    list(name = "flag_not_age_associated", flags = c("flag_not_age_associated")),
    list(name = "flag_germline_ind",       flags = c("flag_germline_ind")),
    list(name = "flag_cluster",            flags = c("flag_cluster")),
    list(name = "flag_freeze_age", flags = c("flag_freeze_enriched", "flag_not_age_associated")),
    list(name = "flag_freeze_age_germ", flags = c("flag_freeze_enriched", "flag_not_age_associated", "flag_germline_ind")),
    list(name = "flag_freeze_age_germ_cluster", flags = c("flag_freeze_enriched", "flag_not_age_associated",
                                                          "flag_germline_ind", "flag_cluster"))
)

minad_configs <- list(
    list(name = "minad3",      flags = c("minad3")),
    list(name = "minad4",      flags = c("minad4")),
    list(name = "minad5",      flags = c("minad5")),
    list(name = "minad6",      flags = c("minad6")),
    list(name = "minad7",      flags = c("minad7")),
    list(name = "minad8",      flags = c("minad8")),
    list(name = "minad3_full", flags = c("minad3", full_qc)),
    list(name = "minad4_full", flags = c("minad4", full_qc)),
    list(name = "minad5_full", flags = c("minad5", full_qc)),
    list(name = "minad6_full", flags = c("minad6", full_qc)),
    list(name = "minad7_full", flags = c("minad7", full_qc)),
    list(name = "minad8_full", flags = c("minad8", full_qc)),
    list(name = "full_covs",       flags = full_qc,                formula_extra = extra_covs),
    list(name = "minad3_full_covs",  flags = c("minad3",  full_qc), formula_extra = extra_covs),
    list(name = "minad4_full_covs",  flags = c("minad4",  full_qc), formula_extra = extra_covs),
    list(name = "minad5_full_covs",  flags = c("minad5",  full_qc), formula_extra = extra_covs),
    list(name = "minad6_full_covs",  flags = c("minad6",  full_qc), formula_extra = extra_covs),
    list(name = "minad7_full_covs",  flags = c("minad7",  full_qc), formula_extra = extra_covs),
    list(name = "minad8_full_covs",  flags = c("minad8",  full_qc), formula_extra = extra_covs)
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
filter_res <- test_age_association(filter_configs, all_ch, cov)
minad_res <- test_age_association(minad_configs,  all_ch, cov)

filter_vaf_res <- test_age_association(filter_configs, all_ch, cov, stratify_by = "VAF_Strata")
minad_vaf_res <- test_age_association(minad_configs,  all_ch, cov, stratify_by = "VAF_Strata")

filter_freeze_res <- test_age_association(filter_configs, all_ch, cov, stratify_by = "Batch")
minad_freeze_res  <- test_age_association(minad_configs,  all_ch, cov, stratify_by = "Batch")

# ========================
# SAVE RESULTS
# ========================
write_xlsx(filter_res,        file.path("ch", "data", "ch_filter_age_sensitivity.xlsx"))
write_xlsx(minad_res,         file.path("ch", "data", "ch_minad_age_sensitivity.xlsx"))

write_xlsx(filter_vaf_res,     file.path("ch", "data", "ch_filter_age_sensitivity_vaf.xlsx"))
write_xlsx(minad_vaf_res,       file.path("ch", "data", "ch_minad_age_sensitivity_vaf.xlsx"))

write_xlsx(filter_freeze_res,  file.path("ch", "data", "ch_filter_age_sensitivity_batch.xlsx"))
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

save_assoc_plot(filter_res,        file.path("ch", "figures", "qc_filter_age_sensitivity.pdf"),
                "Age association — filters", nrow(cov))
save_assoc_plot(minad_res,         file.path("ch", "figures", "qc_minad_age_sensitivity.pdf"),
                "Age association — minAD thresholds", nrow(cov))

save_assoc_plot(filter_vaf_res,    file.path("ch", "figures", "qc_filter_age_sensitivity_vaf.pdf"),
                "Age association — filters by VAF stratum", nrow(cov))
save_assoc_plot(minad_vaf_res,     file.path("ch", "figures", "qc_minad_age_sensitivity_vaf.pdf"),
                "Age association — minAD by VAF stratum", nrow(cov))

save_assoc_plot(filter_freeze_res, file.path("ch", "figures", "qc_filter_age_sensitivity_batch.pdf"),
                "Age association — filters by batch", nrow(cov))
save_assoc_plot(minad_freeze_res,  file.path("ch", "figures", "qc_minad_age_sensitivity_batch.pdf"),
                "Age association — minAD by batch", nrow(cov))

# ========================
# COMPARE VAF 10< vs. ALL VAF
# ========================
vaf_compare_res <- bind_rows(
    test_age_association(filter_configs, all_ch, cov) %>% mutate(vaf_group = "All VAF"),
    test_age_association(filter_configs, all_ch %>% filter(Sample.AltFrac >= 0.10), cov) %>%
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

ggsave(file.path("ch", "figures", "qc_filter_age_cov_compare.pdf"),
       p_vaf_compare,
       width  = 7,
       height = 0.4 * n_distinct(plot_df$name) + 2,
       limitsize = FALSE)


# ========================
# FINAL FILTER
# ========================
all_ch3 <- all_ch %>%
    filter(!flag_freeze_enriched, !flag_not_age_associated,
           !flag_germline_ind, !flag_cluster)
dim(all_ch3)
# 399
write_xlsx(all_ch3, file.path("ch", "data", "ch_seq_wl_art_minad3_vars.xlsx"))

all_ch4 <- all_ch %>%
    filter(!flag_freeze_enriched, !flag_not_age_associated,
           !flag_germline_ind, !flag_cluster,
           !minad4)
dim(all_ch4)
# 175
write_xlsx(all_ch4, file.path("ch", "data", "ch_seq_wl_art_minad4_vars.xlsx"))

all_ch5 <- all_ch %>%
    filter(!flag_freeze_enriched, !flag_not_age_associated,
           !flag_germline_ind, !flag_cluster,
           !minad5)
dim(all_ch5)
# 98
write_xlsx(all_ch5, file.path("ch", "data", "ch_seq_wl_art_minad5_vars.xlsx"))

all_ch6 <- all_ch %>%
    filter(!flag_freeze_enriched, !flag_not_age_associated,
           !flag_germline_ind, !flag_cluster,
           !minad6)
dim(all_ch6)
# 68
write_xlsx(all_ch6, file.path("ch", "data", "ch_seq_wl_art_minad6_vars.xlsx"))

# ========================
# TOP GENES BY FREEZE (POST-QC)
# ========================
top_genes <- all_ch_clean %>%
    count(Gene) %>%
    slice_max(n, n = 20) %>%
    pull(Gene)

gene_counts <- all_ch_clean %>%
    filter(Gene %in% top_genes) %>%
    mutate(
        Freeze = factor(Batch, levels = c(1, 2), labels = c("Freeze 2", "Freeze 3")),
        Gene   = factor(Gene, levels = top_genes)  # ordered by overall count
    ) %>%
    count(Gene, Freeze)

p_genes <- ggplot(gene_counts, aes(x = Gene, y = n, fill = Freeze)) +
    geom_col(position = "dodge", width = 0.7) +
    scale_fill_manual(values = c("Freeze 2" = "#378ADD", "Freeze 3" = "#1D9E75")) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    labs(
        title    = "Top genes after QC filtering",
        subtitle = sprintf("minAD ≥ 4 + full QC | n = %d variants", nrow(all_ch_clean)),
        x        = NULL, y = "Variant count", fill = NULL
    ) +
    theme_minimal(base_size = 11) +
    theme(
        plot.title       = element_text(face = "bold", hjust = 0.5),
        plot.subtitle    = element_text(hjust = 0.5, size = 9, color = "grey40"),
        axis.text.x      = element_text(angle = 40, hjust = 1),
        panel.grid.major.x = element_blank(),
        panel.grid.minor   = element_blank(),
        legend.position    = "top"
    )

ggsave(file.path("ch", "figures", "ch_top_genes_by_freeze.pdf"),
       p_genes, width = 9, height = 5)

# ========================
# DIFF MINAD CUTOFFS for freeze
# ========================
# asymmetric: different cutoffs per batch
# format: minad_f2_{x}_f3_{y} = TRUE if below threshold for that sample's batch
for (t2 in MIN_AD_THRESHOLDS) {
    for (t3 in MIN_AD_THRESHOLDS) {
        col <- paste0("minad_f2_", t2, "_f3_", t3)
        all_ch[[col]] <- ifelse(
            all_ch$Batch == 1,
            all_ch$Sample.AltDepth < t2,   # freeze 2
            all_ch$Sample.AltDepth < t3    # freeze 3
        )
    }
}

write_xlsx(all_ch, file.path("ch", "data", "ch_seq_wl_flags_minad_vars.xlsx"))

# asymmetric minad configs
minad_asym_configs <- lapply(MIN_AD_THRESHOLDS, function(t2) {
    lapply(MIN_AD_THRESHOLDS, function(t3) {
        list(
            name  = sprintf("minad_f2_%d_f3_%d", t2, t3),
            flags = c(paste0("minad_f2_", t2, "_f3_", t3))
        )
    })
}) %>% unlist(recursive = FALSE)

# asymmetric + full qc
minad_asym_full_configs <- lapply(MIN_AD_THRESHOLDS, function(t2) {
    lapply(MIN_AD_THRESHOLDS, function(t3) {
        list(
            name  = sprintf("minad_f2_%d_f3_%d_full", t2, t3),
            flags = c(paste0("minad_f2_", t2, "_f3_", t3), full_qc)
        )
    })
}) %>% unlist(recursive = FALSE)

# run
minad_asym_res      <- test_age_association(minad_asym_configs,      all_ch, cov)
minad_asym_full_res <- test_age_association(minad_asym_full_configs,  all_ch, cov)

write_xlsx(minad_asym_res,      file.path("ch", "data", "ch_minad_asym_age_sensitivity.xlsx"))
write_xlsx(minad_asym_full_res, file.path("ch", "data", "ch_minad_asym_full_age_sensitivity.xlsx"))

# plot — heatmap works better than forest plot for a grid of combinations
asym_heat_df <- minad_asym_full_res %>%
    filter(stratum == "all") %>%
    mutate(
        t2 = as.integer(str_extract(name, "(?<=f2_)\\d")),
        t3 = as.integer(str_extract(name, "(?<=f3_)\\d")),
        sig = !is.na(p_age) & p_age < 0.05
    )

p_asym_heat <- ggplot(asym_heat_df, aes(x = factor(t2), y = factor(t3), fill = -log10(p_age))) +
    geom_tile(color = "white") +
    geom_text(aes(label = ifelse(sig, sprintf("%.3f*", p_age), sprintf("%.3f", p_age))),
              size = 3) +
    scale_fill_gradient(
        low      = "white",
        high     = "#D85A30",
        na.value = "grey90",
        name     = expression(-log[10](p))
    ) +
    geom_tile(data = filter(asym_heat_df, sig),   # outline significant cells
              color = "#D85A30", fill = NA, linewidth = 1) +
    labs(
        title    = "Age association by asymmetric minAD cutoffs (+ full QC)",
        subtitle = "* = p < 0.05 | fill = -log10(p_age)",
        x        = "minAD threshold — Freeze 2",
        y        = "minAD threshold — Freeze 3"
    ) +
    theme_minimal(base_size = 11) +
    theme(
        plot.title    = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5, size = 9, color = "grey40"),
        panel.grid    = element_blank()
    )

ggsave(file.path("ch", "figures", "qc_minad_asym_heatmap.pdf"),
       p_asym_heat, width = 8, height = 6)


# ========================
# COMPARE SPECIFIC ASYMMETRIC VS SYMMETRIC MINAD
# ========================

# define comparison configs
minad_compare_configs <- list(
    list(name = "no_filter",                    flags = character(0)),
    list(name = "minad4 (symmetric)",           flags = c("minad4")),
    list(name = "minad6 (symmetric)",           flags = c("minad6")),
    list(name = "minad f2=4, f3=6 (asymmetric)", flags = c("minad_f2_4_f3_6")),
    list(name = "minad4 + full QC",             flags = c("minad4", full_qc)),
    list(name = "minad6 + full QC",             flags = c("minad6", full_qc)),
    list(name = "minad f2=4, f3=6 + full QC",  flags = c("minad_f2_4_f3_6", full_qc))
)

minad_compare_res <- test_age_association(minad_compare_configs, all_ch, cov)

# optionally also stratified by VAF
minad_compare_vaf_res <- test_age_association(minad_compare_configs, all_ch, cov,
                                              stratify_by = "VAF_Strata")

# plot
save_assoc_plot(minad_compare_res,
                file.path("ch", "figures", "qc_minad_asym_compare.pdf"),
                "Age association — symmetric vs asymmetric minAD", nrow(cov))

save_assoc_plot(minad_compare_vaf_res,
                file.path("ch", "figures", "qc_minad_asym_compare_vaf.pdf"),
                "Age association — symmetric vs asymmetric minAD by VAF", nrow(cov))

# save results
write_xlsx(minad_compare_res,
           file.path("ch", "data", "ch_minad_asym_compare.xlsx"))



