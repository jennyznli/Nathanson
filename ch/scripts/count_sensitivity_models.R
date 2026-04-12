run_chip_count <- function(
        cov_base,
        minad_thresholds = c(3, 4, 5),
        vars_dir         = file.path("ch", "data"),
        figures_dir      = file.path("ch", "figures"),
        cohort_label     = "all",
        females_only     = FALSE,
        case_vars        = list("BRCA1/2" = "BRCA12_Case",
                                "BRCA1"   = "BRCA1_Case",
                                "BRCA2"   = "BRCA2_Case"),
        base_covs_linear = "Sample_age + Sample_age2 + Batch + Smoke_History +
                            Sequenced_gender + PC1 + PC2 + PC3 + PC4 + PC5 + PC6",
        save_forest      = TRUE
) {
    if (females_only) {
        base_covs_linear <- trimws(gsub("\\+\\s*\\+", "+",
                                        gsub("\\+?\\s*Sequenced_gender\\s*\\+?", "+", base_covs_linear)))
        cat("females_only = TRUE: Sequenced_gender removed\n")
    }

    extract_irr <- function(fit, term, data) {
        vcov_cl <- vcovCL(fit, cluster = ~subclass, data = data)
        ct      <- coeftest(fit, vcov = vcov_cl)
        ci      <- coefci(fit,  vcov = vcov_cl)
        data.frame(
            IRR   = exp(ct[term, "Estimate"]),
            CI_lo = exp(ci[term, 1]),
            CI_hi = exp(ci[term, 2]),
            p     = ct[term, "Pr(>|z|)"]
        ) %>% mutate(IRR_fmt = sprintf("%.2f (%.2f\u2013%.2f)", IRR, CI_lo, CI_hi))
    }

    all_results <- list()

    for (minad in minad_thresholds) {
        cat("\n── minAD =", minad, "| cohort:", cohort_label, "──\n")

        vars        <- read_excel(file.path(vars_dir,
                                            sprintf("ch_seq_wl_art_minad%d_vars.xlsx", minad)))
        chip_counts <- vars %>% count(Sample.ID, name = "CHIP_Count")

        cov <- cov_base %>%
            mutate(CHIP_Count = {
                idx <- match(person_id, chip_counts$Sample.ID)
                ifelse(is.na(idx), 0L, chip_counts$CHIP_Count[idx])
            })

        cat("N:", nrow(cov), "| CHIP+:", sum(cov$CHIP_Count > 0),
            "| mean:", round(mean(cov$CHIP_Count), 3), "\n")

        results_list <- list()

        for (i in seq_along(case_vars)) {
            model_label <- names(case_vars)[i]
            case_col    <- case_vars[[i]]
            fml         <- as.formula(paste("CHIP_Count ~", case_col, "+", base_covs_linear))

            # M1: quasi-Poisson + clustered SE
            fit_qp <- glm(fml, data = cov, family = quasipoisson())

            # M2: negative binomial + clustered SE (with fallback)
            fit_nb <- tryCatch(
                glm.nb(fml, data = cov),
                error   = function(e) { message("NB failed: ", e$message); NULL },
                warning = function(w) { message("NB warning: ", w$message); NULL }
            )

            rows <- list()
            rows[["QP"]] <- extract_irr(fit_qp, case_col, cov) %>%
                mutate(model = model_label, spec = "QuasiPoisson", theta = NA_real_)

            if (!is.null(fit_nb)) {
                rows[["NB"]] <- extract_irr(fit_nb, case_col, cov) %>%
                    mutate(model = model_label, spec = "NegBinom", theta = fit_nb$theta)
            }

            results_list[[model_label]] <- bind_rows(rows)
        }

        results <- bind_rows(results_list) %>%
            mutate(minAD = minad, cohort = cohort_label) %>%
            dplyr::select(cohort, minAD, spec, model, IRR_fmt, IRR, CI_lo, CI_hi, p, theta)

        all_results[[paste0("minAD", minad)]] <- results

        if (save_forest) {
            forest_df <- results %>%
                mutate(
                    p_fmt = ifelse(p < 0.001, "p<0.001", sprintf("p=%.3f", p)),
                    annot = sprintf("%.2f (%.2f, %.2f) %s", IRR, CI_lo, CI_hi, p_fmt),
                    label = factor(paste0(model, " (", spec, ")"),
                                   levels = rev(unique(paste0(model, " (", spec, ")"))))
                )

            fig <- ggplot(forest_df,
                          aes(x = IRR, y = label, xmin = CI_lo, xmax = CI_hi,
                              color = spec)) +
                geom_vline(xintercept = 1, linetype = "dashed",
                           color = "grey50", linewidth = 0.5) +
                geom_errorbarh(height = 0.2, linewidth = 0.6) +
                geom_point(size = 3) +
                geom_text(aes(x = max(CI_hi, na.rm = TRUE) * 1.15,
                              label = annot),
                          hjust = 0, size = 2.8, color = "grey30") +
                scale_x_log10(breaks = c(0.5, 1, 2, 4)) +
                scale_color_manual(values = c(QuasiPoisson = "#2166ac",
                                              NegBinom     = "#d6604d"),
                                   name = "Model") +
                coord_cartesian(clip = "off") +
                labs(title = sprintf("CHIP count | cohort: %s | minAD = %d",
                                     cohort_label, minad),
                     x = "IRR (log scale)", y = NULL) +
                theme_minimal(base_size = 11) +
                theme(plot.title         = element_text(face = "bold", hjust = 0.5),
                      plot.margin        = margin(5, 200, 5, 5),
                      panel.grid.minor   = element_blank(),
                      panel.grid.major.y = element_blank(),
                      panel.grid.major.x = element_line(color = "grey92", linewidth = 0.3),
                      legend.position    = "bottom")

            ggsave(file.path(figures_dir,
                             sprintf("fig_forest_count_%s_minad%d.pdf", cohort_label, minad)),
                   fig, width = 8,
                   height = max(4, 0.5 * nrow(forest_df) + 2))
        }
    }
    return(all_results)
}
