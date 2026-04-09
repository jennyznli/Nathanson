run_chip_age_analysis <- function(
        cov_base,
        minad_thresholds = c(4),
        vars_dir         = file.path("ch", "data"),
        figures_dir      = file.path("ch", "figures"),
        cohort_label     = "all",
        females_only     = FALSE,
        case_vars        = list(
            "BRCA1/2" = "BRCA12_Case",
            "BRCA1"   = "BRCA1_Case",
            "BRCA2"   = "BRCA2_Case"
        ),
        base_covs_linear = "Sample_age + Batch + Smoke_History + Sequenced_gender +
                         PC1 + PC2 + PC3 + PC4 + PC5 + PC6",
        age_specs        = c("linear", "squared", "spline"),   # <-- new
        spline_knots     = 3,                                   # <-- new
        save_forest      = TRUE
) {
    library(splines)

    if (females_only) {
        base_covs_linear <- gsub("\\+?\\s*Sequenced_gender\\s*\\+?", "+", base_covs_linear)
        base_covs_linear <- trimws(gsub("\\+\\s*\\+", "+", base_covs_linear))
        cat("females_only = TRUE: Sequenced_gender removed from formula\n")
    }

    # Helper: swap age term based on spec
    make_covs <- function(base_covs, age_spec) {
        switch(age_spec,
               "linear"  = base_covs,
               "squared" = {
                   # Add age^2 term alongside linear age
                   gsub("Sample_age",
                        "Sample_age + I(Sample_age^2)",
                        base_covs)
               },
               "spline"  = {
                   # Replace linear age with restricted cubic spline
                   gsub("Sample_age",
                        sprintf("ns(Sample_age, df = %d)", spline_knots),
                        base_covs)
               }
        )
    }

    clogit_covs_base <- if (females_only) {
        "Sample_age + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + Smoke_History + Batch"
    } else {
        "Sample_age + Sequenced_gender + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + Smoke_History + Batch"
    }

    all_results <- list()

    for (minad in minad_thresholds) {

        cat("\n============================================================\n")
        cat("Cohort:", cohort_label, "| minAD =", minad, "\n")
        cat("============================================================\n")

        vars_file <- file.path(vars_dir,
                               sprintf("ch_seq_wl_art_minad%d_vars.xlsx", minad))
        vars <- read_excel(vars_file)

        cov <- cov_base %>%
            mutate(
                CHIP_Binary = person_id %in% vars$Sample.ID,
                CHIP_Count  = {
                    chip_counts <- vars %>% count(Sample.ID, name = "CHIP_Count")
                    idx <- match(person_id, chip_counts$Sample.ID)
                    ifelse(is.na(idx), 0L, chip_counts$CHIP_Count[idx])
                }
            )

        cat("N:", nrow(cov), "| CHIP+:", sum(cov$CHIP_Binary), "\n")

        results_list <- list()

        for (i in seq_along(case_vars)) {
            model_label <- names(case_vars)[i]
            case_col    <- case_vars[[i]]

            m1_results <- list()
            for (age_spec in age_specs) {
                covs_this <- make_covs(base_covs_linear, age_spec)
                fit_m1 <- glm(
                    as.formula(paste("CHIP_Binary ~", case_col, "+", covs_this)),
                    data = cov, family = quasibinomial()
                )
                m1_results[[age_spec]] <- extract_or_robust(fit_m1, case_col, cov) %>%
                    mutate(model = model_label, spec = paste0("M1_", age_spec), age_spec = age_spec)  # <-- add this
            }

            results_list[[model_label]] <- bind_rows(m1_results)
        }

        results <- bind_rows(results_list) %>%
            mutate(minAD = minad, cohort = cohort_label) %>%
            dplyr::select(cohort, minAD, age_spec, spec, model,
                          OR_fmt, OR, CI_lo, CI_hi, p, AUC)

        all_results[[paste0("minAD", minad)]] <- results

        # --- Forest plot: one per age_spec ---
        # Remove the for (age_spec in age_specs) loop entirely
        # Just plot all results together
        if (save_forest) {
            forest_df <- results %>%
                mutate(
                    spec_label = case_when(
                        spec == "M1_linear"  ~ "Linear",
                        spec == "M1_squared" ~ "Age²",
                        spec == "M1_spline"  ~ "Spline"
                    ),
                    label = factor(
                        paste0(model, " (", spec_label, ")"),
                        levels = rev(paste0(
                            rep(names(case_vars), each = 3),
                            " (",
                            rep(c("Linear", "Age²", "Spline"), length(case_vars)),
                            ")"
                        ))
                    ),
                    p_fmt      = ifelse(p < 0.001, "p<0.001", sprintf("p=%.3f", p)),
                    annotation = sprintf("%.2f (%.2f, %.2f), %s", OR, CI_lo, CI_hi, p_fmt),
                    color_grp  = factor(spec_label, levels = c("Linear", "Age²", "Spline"))
                )

            # ADD THIS — was missing entirely
            fig_forest <- ggplot(forest_df,
                                 aes(x = OR, y = label, xmin = CI_lo, xmax = CI_hi,
                                     color = color_grp)) +
                geom_vline(xintercept = 1, linetype = "dashed",
                           color = "grey50", linewidth = 0.5) +
                geom_errorbarh(height = 0.15, linewidth = 0.7) +
                geom_point(size = 3.5) +
                geom_text(aes(x = 5, label = annotation),
                          hjust = 0, size = 3.0, color = "grey30") +
                scale_x_log10(breaks = c(0.5, 1, 2, 4),
                              labels = c("0.5", "1", "2", "4"),
                              limits = c(0.4, 5)) +
                coord_cartesian(clip = "off") +
                labs(
                    title = sprintf("CHIP ~ gBRCA1/2 | cohort: %s | minAD = %d | age sensitivity",
                                    cohort_label, minad),
                    x = "Odds ratio (log scale)", y = NULL
                ) +
                theme_minimal(base_size = 12) +
                theme(
                    plot.title         = element_text(face = "bold", hjust = 0.5, size = 12),
                    plot.margin        = margin(5, 180, 5, 5),
                    panel.grid.minor   = element_blank(),
                    panel.grid.major.y = element_blank(),
                    panel.grid.major.x = element_line(color = "grey92", linewidth = 0.3),
                    axis.text.y        = element_text(size = 9),
                    axis.text.x        = element_text(size = 10),
                    legend.position    = "bottom"
                )

            ggsave(
                file.path(figures_dir,
                          sprintf("fig_forest_%s_minad%d_agespec.pdf", cohort_label, minad)),
                fig_forest, width = 8, height = 1.7 * length(case_vars)
            )
            cat("Saved forest plot:", cohort_label, "minAD =", minad, "\n")

        }
    }

    return(all_results)
}


extract_or_robust <- function(fit, term, data) {
    vcov_cl <- vcovCL(fit, cluster = ~subclass, data = data)
    ct      <- coeftest(fit, vcov = vcov_cl)
    ci      <- coefci(fit, vcov = vcov_cl)
    auc     <- tryCatch(as.numeric(pROC::auc(fit$y, fitted(fit))),
                        error = function(e) NA_real_)
    data.frame(
        term  = term,
        OR    = exp(ct[term, "Estimate"]),
        CI_lo = exp(ci[term, 1]),
        CI_hi = exp(ci[term, 2]),
        p     = ct[term, "Pr(>|z|)"],
        AUC   = auc
    ) %>% mutate(OR_fmt = sprintf("%.2f (%.2f\u2013%.2f)", OR, CI_lo, CI_hi))
}

extract_or_clogit <- function(fit, term) {
    s  <- summary(fit)$coefficients
    ci <- confint(fit)
    data.frame(
        term  = term,
        OR    = exp(s[term, "coef"]),
        CI_lo = exp(ci[term, 1]),
        CI_hi = exp(ci[term, 2]),
        p     = s[term, "Pr(>|z|)"],
        AUC   = NA_real_
    ) %>% mutate(OR_fmt = sprintf("%.2f (%.2f\u2013%.2f)", OR, CI_lo, CI_hi))
}

extract_or_firth <- function(fit, term) {
    ci <- confint(fit)
    idx <- which(names(fit$coefficients) == term)
    data.frame(
        term  = term,
        OR    = exp(fit$coefficients[idx]),
        CI_lo = exp(ci[idx, 1]),
        CI_hi = exp(ci[idx, 2]),
        p     = fit$prob[idx],
        AUC   = NA_real_
    ) %>% mutate(OR_fmt = sprintf("%.2f (%.2f\u2013%.2f)", OR, CI_lo, CI_hi))
}
