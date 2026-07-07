# ============================================================
# GENE-LEVEL SENSITIVITY - FIRTH ONLy
# ============================================================
make_cov_specs <- function(females_only = FALSE) {
    sex  <- if (females_only) "" else "+ Sequenced_gender"
    base <- glue::glue("Sample_age + Sample_age2 + Batch + Smoke_History {sex}")
    list(
        base = base,
        full = glue::glue("{base} + PC1 + PC2 + PC3 + PC4 + PC5 + PC6")
    ) %>% lapply(\(f) trimws(gsub("\\+\\s*\\+", "+", f)))
}

or_firth <- function(fit, term) {
    i <- which(names(coef(fit)) == term)
    data.frame(OR = exp(coef(fit)[i]), CI_lo = exp(fit$ci.lower[i]),
               CI_hi = exp(fit$ci.upper[i]), p = fit$prob[i])
}

run_one_gene <- function(gene, term, cov, vars_df, cov_specs, min_co = 2) {
    ids    <- vars_df %>% filter(Gene == gene) %>% distinct(Sample.ID) %>% pull()
    dat    <- cov %>% mutate(gene_chip = as.integer(person_id %in% ids))
    n_chip <- sum(dat$gene_chip)
    n_co   <- sum(dat$gene_chip[dat[[term]] == 1])

    if (n_chip < 3 || n_co < min_co || (n_chip - n_co) < min_co) return(NULL)

    purrr::map_dfr(names(cov_specs), \(spec) {
        fml <- as.formula(paste("gene_chip ~", term, "+", cov_specs[[spec]]))
        res <- tryCatch(
            or_firth(logistf(fml, data = dat), term),
            error = \(e) { message(gene, "/", spec, ": ", e$message); NULL }
        )
        if (is.null(res)) return(NULL)
        res %>% mutate(Gene = gene, exposure = term, cov_spec = spec,
                       n_chip = n_chip, n_co = n_co)
    })
}

run_gene_final <- function(
        cov, vars,
        min_carriers = 3,
        females_only = FALSE,
        case_vars    = list("BRCA1/2" = "BRCA12_Case",
                            "BRCA1"   = "BRCA1_Case",
                            "BRCA2"   = "BRCA2_Case"),
        out_dir      = file.path("ch", "data"),
        fig_dir      = file.path("ch", "figures"),
        cohort_label = ""
) {
    if (females_only) {
        cov          <- cov %>% filter(Sequenced_gender == "Female")
        cohort_label <- paste0(cohort_label, "_female")
        cat("females_only: N =", nrow(cov), "\n")
    }

    cov_specs     <- make_cov_specs(females_only)
    genes_to_test <- vars %>%
        left_join(cov %>% dplyr::select(person_id, BRCA12_Case),
                  by = c("Sample.ID" = "person_id")) %>%
        filter(!is.na(BRCA12_Case)) %>%
        distinct(Sample.ID, Gene, BRCA12_Case) %>%
        group_by(Gene) %>%
        summarise(n_carrier = sum(BRCA12_Case == 1), .groups = "drop") %>%
        filter(n_carrier >= min_carriers) %>%
        pull(Gene)

    cat("Genes:", length(genes_to_test), "| Specs:", paste(names(cov_specs), collapse = ", "), "\n")

    results <- purrr::map_dfr(names(case_vars), \(label) {
        term <- case_vars[[label]]
        cat("\nExposure:", label, "\n")
        purrr::map_dfr(genes_to_test, run_one_gene,
                       term = term, cov = cov, vars_df = vars,
                       cov_specs = cov_specs) %>%
            group_by(cov_spec) %>%
            mutate(p_fdr        = p.adjust(p, "BH"),
                   p_bonferroni = p.adjust(p, "bonferroni")) %>%
            ungroup() %>%
            mutate(exposure_label = label,
                   OR_fmt = sprintf("%.2f (%.2f–%.2f)", OR, CI_lo, CI_hi),
                   cohort = cohort_label)
    }) %>% arrange(exposure_label, cov_spec, p_fdr)

    write_xlsx(results, file.path(out_dir,
                                  sprintf("ch_gene_sensitivity_%s.xlsx", cohort_label)))

    plot_sensitivity_forest(results, fig_dir, cohort_label)
    invisible(results)
}

plot_sensitivity_forest <- function(results, fig_dir, cohort_label) {
    for (label in unique(results$exposure_label)) {
        for (spec in unique(results$cov_spec)) {
            df <- results %>%
                filter(exposure_label == label, cov_spec == spec) %>%
                mutate(
                    sig   = p_fdr < 0.05,
                    shape = if_else(sig, "Significant (FDR<0.05)", "Not significant"),
                    annot = sprintf("%.2f (%.2f–%.2f) %s",
                                    OR, CI_lo, CI_hi,
                                    ifelse(p < 0.001, "p<0.001", sprintf("p=%.3f", p)))
                )
            if (nrow(df) == 0) next

            gene_order <- df %>%
                group_by(Gene) %>%
                summarise(med_OR = median(OR, na.rm = TRUE), .groups = "drop") %>%
                arrange(med_OR) %>%
                pull(Gene)

            df <- df %>%
                mutate(Gene = factor(Gene, levels = gene_order))

            n_genes     <- length(gene_order)
            plot_height <- max(3, n_genes * 0.4 + 1.5)   # ~0.4 in per gene, min 3 in

            p <- ggplot(df, aes(OR, Gene, xmin = CI_lo, xmax = CI_hi, shape = shape)) +
                geom_vline(xintercept = 1, linetype = "dashed", color = "grey60") +
                geom_errorbarh(height = 0.3, linewidth = 0.6, color = "#2166ac") +
                geom_point(size = 3, color = "#2166ac") +
                geom_text(aes(x = max(df$CI_hi, na.rm = TRUE) * 1.15, label = annot),
                          hjust = 0, size = 2.5, color = "grey30") +
                scale_x_log10(breaks = c(0.25, 0.5, 1, 2, 4, 8)) +
                scale_shape_manual(
                    values = c("Significant (FDR<0.05)" = 17, "Not significant" = 16),
                    name   = NULL
                ) +
                coord_cartesian(clip = "off") +
                # labs(title = sprintf("%s | covs: %s | %s", label, spec, cohort_label),
                #      x     = "OR (log scale)", y = NULL) +
                theme_classic(base_size = 11) +
                theme(plot.margin     = margin(5, 220, 5, 5),
                      legend.position = "bottom",
                      axis.text.y     = element_text(size = 9))

            ggsave(file.path(fig_dir, sprintf("fig_gene_%s_%s_%s.pdf",
                                              cohort_label, gsub("/", "", label), spec)),
                   p, width = 8, height = plot_height)
        }
    }
}


# # ── covariate specs ───────────────────────────────────────────
# make_cov_specs <- function(females_only = FALSE) {
#     sex <- if (females_only) "" else "+ Sequenced_gender"
#     base <- glue::glue("Sample_age + Sample_age2 + Batch + Smoke_History {sex}")
#     list(
#         base   = base,
#         full   = glue::glue("{base} + PC1 + PC2 + PC3 + PC4 + PC5 + PC6")
#     ) %>% lapply(\(f) trimws(gsub("\\+\\s*\\+", "+", f)))
# }
#
# # ── extract helpers ───────────────────────────────────────────
# or_firth <- function(fit, term) {
#     i <- which(names(coef(fit)) == term)
#     data.frame(OR = exp(coef(fit)[i]), CI_lo = exp(fit$ci.lower[i]),
#                CI_hi = exp(fit$ci.upper[i]), p = fit$prob[i])
# }
#
# or_glm <- function(fit, term) {
#     ct  <- coeftest(fit, vcov = vcovHC(fit, type = "HC3"))
#     i   <- which(rownames(ct) == term)
#     est <- ct[i, "Estimate"]; se <- ct[i, "Std. Error"]
#     data.frame(OR = exp(est), CI_lo = exp(est - 1.96*se),
#                CI_hi = exp(est + 1.96*se), p = ct[i, "Pr(>|z|)"])
# }
#
# or_clogit <- function(fit, term) {
#     i  <- which(rownames(summary(fit)$coefficients) == term)
#     ci <- confint(fit)[i, ]
#     data.frame(OR = exp(coef(fit)[i]), CI_lo = exp(ci[1]),
#                CI_hi = exp(ci[2]), p = summary(fit)$coefficients[i, "Pr(>|z|)"])
# }
#
# # ── per-gene runner ───────────────────────────────────────────
# run_one_gene <- function(gene, term, cov, vars_df, cov_specs,
#                          model_types, min_co = 2) {
#     ids      <- vars_df %>% filter(Gene == gene) %>% distinct(Sample.ID) %>% pull()
#     dat      <- cov %>% mutate(gene_chip = as.integer(person_id %in% ids))
#     n_chip   <- sum(dat$gene_chip)
#     n_co     <- sum(dat$gene_chip[dat[[term]] == 1])
#
#     if (n_chip < 3 || n_co < min_co || (n_chip - n_co) < min_co) return(NULL)
#
#     purrr::map_dfr(names(cov_specs), \(spec) {
#         covs <- cov_specs[[spec]]
#         fml  <- as.formula(paste("gene_chip ~", term, "+", covs))
#
#         purrr::map_dfr(model_types, \(mtype) {
#             res <- tryCatch({
#                 if (mtype == "firth") {
#                     or_firth(logistf(fml, data = dat), term)
#
#                 } else if (mtype == "quasibinom_unweighted") {
#                     or_glm(glm(fml, data = dat, family = quasibinomial()), term)
#
#                 } else if (mtype == "quasibinom_weighted") {
#                     or_glm(glm(fml, data = dat, family = quasibinomial(),
#                                weights = dat$weights), term)
#
#                 } else if (mtype == "clogit") {
#                     cfml <- as.formula(paste("gene_chip ~", term, "+", covs,
#                                              "+ strata(subclass_num)"))
#                     # drop strata with no variation in outcome
#                     dat_c <- dat %>%
#                         group_by(subclass_num) %>%
#                         filter(n_distinct(gene_chip) > 1) %>%
#                         ungroup()
#                     if (nrow(dat_c) == 0) return(NULL)
#                     or_clogit(clogit(cfml, data = dat_c), term)
#                 }
#             }, error = \(e) { message(gene, "/", spec, "/", mtype, ": ", e$message); NULL })
#
#             if (is.null(res)) return(NULL)
#             res %>% mutate(Gene = gene, exposure = term, cov_spec = spec,
#                            model_type = mtype, n_chip = n_chip, n_co = n_co)
#         })
#     })
# }
#
# # ── master wrapper ────────────────────────────────────────────
# run_gene_sensitivity <- function(
#         cov, vars,
#         min_carriers = 3,
#         females_only = FALSE,
#         model_types  = c("firth", "quasibinom_unweighted",
#                          "quasibinom_weighted", "clogit"),
#         case_vars    = list("BRCA1/2" = "BRCA12_Case",
#                             "BRCA1"   = "BRCA1_Case",
#                             "BRCA2"   = "BRCA2_Case"),
#         out_dir      = file.path("ch", "data"),
#         fig_dir      = file.path("ch", "figures"),
#         cohort_label = ""
# ) {
#     if (females_only) {
#         cov          <- cov %>% filter(Sequenced_gender == "Female")
#         cohort_label <- paste0(cohort_label, "_female")
#         cat("females_only: N =", nrow(cov), "\n")
#     }
#
#     cov_specs     <- make_cov_specs(females_only)
#     genes_to_test <- vars %>%
#         left_join(cov %>% dplyr::select(person_id, BRCA12_Case),
#                   by = c("Sample.ID" = "person_id")) %>%
#         filter(!is.na(BRCA12_Case)) %>%
#         distinct(Sample.ID, Gene, BRCA12_Case) %>%
#         group_by(Gene) %>%
#         summarise(n_carrier = sum(BRCA12_Case == 1), .groups = "drop") %>%
#         filter(n_carrier >= min_carriers) %>%
#         pull(Gene)
#
#     cat("Genes:", length(genes_to_test), "| Specs:", paste(names(cov_specs), collapse = ", "),
#         "| Models:", paste(model_types, collapse = ", "), "\n")
#
#     results <- purrr::map_dfr(names(case_vars), \(label) {
#         term <- case_vars[[label]]
#         cat("\nExposure:", label, "\n")
#         purrr::map_dfr(genes_to_test, run_one_gene,
#                        term = term, cov = cov, vars_df = vars,
#                        cov_specs = cov_specs, model_types = model_types) %>%
#             group_by(cov_spec, model_type) %>%
#             mutate(p_fdr        = p.adjust(p, "BH"),
#                    p_bonferroni = p.adjust(p, "bonferroni")) %>%
#             ungroup() %>%
#             mutate(exposure_label = label,
#                    OR_fmt = sprintf("%.2f (%.2f–%.2f)", OR, CI_lo, CI_hi),
#                    cohort = cohort_label)
#     }) %>% arrange(exposure_label, cov_spec, model_type, p_fdr)
#
#     write_xlsx(results, file.path(out_dir,
#                                   sprintf("ch_gene_sensitivity_%s.xlsx", cohort_label)))
#
#     # exclude the clogit when plotting bc skews too much ...
#     plot_sensitivity_forest(results %>% filter(model_type != "clogit"), fig_dir, cohort_label)
#     invisible(results)
# }
#
# # ── forest plot ───────────────────────────────────────────────
# plot_sensitivity_forest <- function(results, fig_dir, cohort_label) {
#     model_colors <- c(firth                = "#2166ac",
#                       quasibinom_unweighted = "#d6604d",
#                       quasibinom_weighted   = "#f4a582",
#                       clogit                = "#4dac26")
#
#     for (label in unique(results$exposure_label)) {
#         for (spec in unique(results$cov_spec)) {
#             df <- results %>%
#                 filter(exposure_label == label, cov_spec == spec) %>%
#                 mutate(sig   = case_when(p_fdr < 0.05 ~ "Padj < 0.05",
#                                          p     < 0.05 ~ "P < 0.05",
#                                          TRUE         ~ "ns"),
#                        annot = sprintf("%.2f (%.2f–%.2f) %s", OR, CI_lo, CI_hi,
#                                        ifelse(p < 0.001, "p<0.001", sprintf("p=%.3f", p))))
#             if (nrow(df) == 0) next
#
#             # order genes by median OR across models, then build stacked y labels
#             gene_order <- df %>%
#                 group_by(Gene) %>%
#                 summarise(med_OR = median(OR, na.rm = TRUE), .groups = "drop") %>%
#                 arrange(med_OR) %>%
#                 pull(Gene)
#
#             df <- df %>%
#                 mutate(Gene = factor(Gene, levels = gene_order)) %>%
#                 arrange(Gene, model_type) %>%
#                 mutate(y_label = paste0(Gene, " — ", model_type),
#                        y_label = factor(y_label, levels = unique(y_label)))
#
#             # gene separator positions (draw lines between gene groups)
#             gene_breaks <- df %>%
#                 group_by(Gene) %>%
#                 summarise(y_max = max(as.integer(y_label)), .groups = "drop") %>%
#                 filter(Gene != last(levels(df$Gene))) %>%
#                 pull(y_max) %>%
#                 `+`(0.5)
#
#             p <- ggplot(df, aes(OR, y_label, xmin = CI_lo, xmax = CI_hi,
#                                 color = model_type, shape = sig)) +
#                 geom_hline(yintercept = gene_breaks, color = "grey88",
#                            linewidth = 0.4) +
#                 geom_vline(xintercept = 1, linetype = "dashed", color = "grey60") +
#                 geom_errorbarh(height = 0.3, linewidth = 0.6) +
#                 geom_point(size = 2.8) +
#                 geom_text(aes(x = max(df$CI_hi, na.rm = TRUE) * 1.15, label = annot),
#                           hjust = 0, size = 2.3, color = "grey30") +
#                 scale_x_log10(breaks = c(0.25, 0.5, 1, 2, 4, 8)) +
#                 scale_color_manual(values = model_colors, name = "Model") +
#                 scale_shape_manual(values = c("FDR<0.05" = 8, "p<0.05" = 17, "ns" = 16),
#                                    name = "Sig") +
#                 coord_cartesian(clip = "off") +
#                 labs(title = sprintf("%s | covs: %s | %s", label, spec, cohort_label),
#                      x = "OR (log scale)", y = NULL) +
#                 theme_classic(base_size = 11) +
#                 theme(strip.background = element_blank(),
#                       plot.margin      = margin(5, 220, 5, 5),
#                       legend.position  = "bottom",
#                       axis.text.y      = element_text(size = 8))
#
#             ggsave(file.path(fig_dir, sprintf("fig_gene_%s_%s_%s.pdf",
#                                               cohort_label, gsub("/", "", label), spec)),
#                    p, width = 8,
#                    height = max(5, 0.28 * nrow(df) + 2))
#         }
#     }
# }
