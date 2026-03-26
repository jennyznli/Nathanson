# ========================
# Age association for artifacts
# ========================
library(here)
setwd("/Users/jennyzli/Documents/Nathanson")
source(here("R", "config.R"))
library(data.table, quietly = T)

# ========================
# READ IN DATA
# ========================
all_ch <- read_excel(file.path("ch", "data", "ch_wl_final.xlsx"))
dim(all_ch)
# 315

cov <- read.csv(file.path("ch", "data", "pmbb_brca12_cov_df.csv"))
all_ids <- read.csv(file.path("ch", "data", "ch_psm_matched4_case_control_ids.csv"))$x
cov <- cov %>% filter(person_id %in% all_ids)
dim(cov)
# 3004

rs <- read_tsv(file.path("ch", "data", "rs7705526_output.raw"))

# ========================
# ASSIGN CATEGORY
# ========================
all_ch <- all_ch %>%
    dplyr::mutate(
        variant_category = dplyr::case_when(
            manualreview == TRUE                                ~ "manual_review",
            wl.exception == TRUE                               ~ "exception",
            wl.mis  == TRUE & !wl.exception & !manualreview   ~ "missense",
            wl.lof  == TRUE & !wl.exception & !manualreview   ~ "lof",
            wl.splice == TRUE & !wl.exception & !manualreview ~ "splice"
            # TRUE                                               ~ "non_whitelist"
        )
    )
print(table(all_ch$variant_category))
# exception           lof manual_review      missense non_whitelist        splice
# 40           407           175            57         22116            95

plot_category_breakdown(all_ch, title = "CHIP variant whitelist breakdown",
                        filename = "ch_wl_categories.pdf")
plot_gene_exceptions(all_ch, filename = "ch_wl_gene_exceptions.pdf")
plot_genes(all_ch, filename = "ch_wl_all_genes.pdf", width = 10)
plot_genes(all_ch, top_n = 10, title = "Top 10 CHIP genes",
           filename = "ch_wl_genes_top10.pdf", width = 7)

# ========================
# BUILD VARIANT_ID
# ========================
trim_splice_hgvsc <- function(hgvsc) {
    if (is.na(hgvsc) || hgvsc == "" || hgvsc == ".") return(NA_character_)
    hgvsc_clean <- sub("^[A-Z]{2}_[0-9]+(?:\\.[0-9]+)?:", "", hgvsc)
    sub("(c\\.\\d+)[+\\-].*", "\\1", hgvsc_clean)
}

clean_hgvsc <- function(hgvsc) {
    if (is.na(hgvsc) || hgvsc == "" || hgvsc == ".") return(NA_character_)
    sub("^[A-Z]{2}_[0-9]+(?:\\.[0-9]+)?:", "", hgvsc)
}

all_ch <- all_ch %>%
    dplyr::mutate(
        HGVSp_clean = dplyr::case_when(
            is.na(HGVSp) | HGVSp == "" | HGVSp == "." ~ NA_character_,
            TRUE ~ sub("^[A-Z]{2}_[0-9]+(?:\\.[0-9]+)?:", "", HGVSp)
        ),
        exon_label = dplyr::case_when(
            !is.na(ExonNumber) ~ paste0("exon", ExonNumber),
            TRUE               ~ "exonUnk"
        ),
        variant_id = dplyr::case_when(
            variant_category == "non_whitelist" ~ NA_character_,
            # if it's a splice, just HGVSc and exon
            variant_category == "splice" ~ paste(
                Gene, MANE_Select_Str, exon_label,
                sapply(HGVSc, trim_splice_hgvsc), sep = ":"
            ),
            # if its anything else, put HGVSp as well
            TRUE ~ paste(
                Gene, MANE_Select_Str, exon_label,
                sapply(HGVSc, clean_hgvsc),
                ifelse(!is.na(HGVSp_clean), HGVSp_clean, ""),
                sep = ":"
            )
        ),
        variant_id = sub(":$", "", variant_id)
    ) %>%
    dplyr::select(-exon_label, -HGVSp_clean)

length(unique((all_ch$variant_id)))
# 213

length(unique((all_ch$Sample.ID)))
# 273

# ========================
# ARTIFACT FILTERING
# Retained if p < 0.10 (associated with age = likely true CHIP).
# ========================
cov <- cov %>% filter(Strata == 1)
N_TOTAL <- nrow(cov)
cat("Total samples:", N_TOTAL, "\n")

FREQ_THRESHOLD <- 4

all_samples_cov <- cov %>%
    rename(Sample.ID = person_id) %>%
    select(Sample.ID, Age = all_of("Sample_age"))

# Count carriers per variant
variant_counts <- all_ch %>%
    filter(!is.na(variant_id)) %>%
    group_by(Chr, Start) %>%
    summarise(n_carriers = n_distinct(Sample.ID), Gene = unique(Gene), .groups = "drop") %>%
    arrange(desc(n_carriers))

# remove known artifact first
all_ch <- all_ch %>% filter(variant_id != "PDS5B:NM_015032:exonUnk:c.2407")

# count carriers per variant_id
variant_counts <- all_ch %>%
    filter(!is.na(variant_id)) %>%
    group_by(variant_id) %>%
    summarise(
        n_carriers = n_distinct(Sample.ID),
        Gene       = first(Gene),
        Chr        = first(Chr),
        Start      = first(Start),
        .groups    = "drop"
    ) %>%
    arrange(desc(n_carriers))

high_freq_variants <- variant_counts %>% filter(n_carriers >= FREQ_THRESHOLD)

vkey <- results$variant_id[1]
carriers <- all_ch %>%
    filter(variant_id == vkey) %>%
    distinct(Sample.ID) %>%
    mutate(carrier = 1)
df <- all_samples_cov %>%
    left_join(carriers, by = "Sample.ID") %>%
    mutate(carrier = ifelse(is.na(carrier), 0, carrier))
m  <- glm(carrier ~ Age, data = df, family = binomial())
s  <- coef(summary(m))
ci <- confint.default(m)
data.frame(
    OR    = exp(coef(m)["Age"]),
    CI_lo = exp(ci["Age", 1]),
    CI_hi = exp(ci["Age", 2]),
    p_age = s["Age", "Pr(>|z|)"]
)

artifact_test_age <- function(vkey, all_samples_cov) {
    carriers <- all_ch %>%
        filter(variant_id == vkey) %>%
        distinct(Sample.ID) %>%
        mutate(carrier = 1)

    df <- all_samples_cov %>%
        left_join(carriers, by = "Sample.ID") %>%
        mutate(carrier = ifelse(is.na(carrier), 0, carrier))

    tryCatch({
        m  <- glm(carrier ~ Age, data = df, family = binomial())
        s  <- coef(summary(m))
        ci <- confint.default(m)
        data.frame(
            OR    = exp(coef(m)["Age"]),
            CI_lo = exp(ci["Age", 1]),
            CI_hi = exp(ci["Age", 2]),
            p_age = s["Age", "Pr(>|z|)"]
        )
    }, error = function(e) {
        data.frame(OR = NA_real_, CI_lo = NA_real_, CI_hi = NA_real_, p_age = NA_real_)
    })
}

results <- high_freq_variants %>%
    rowwise() %>%
    mutate(res = list(artifact_test_age(variant_id, all_samples_cov))) %>%
    unnest(res) %>%
    ungroup() %>%
    mutate(
        status = case_when(
            is.na(p_age) ~ "artifact",
            p_age < 0.10 ~ "retained",
            TRUE         ~ "artifact"
        ),
        OR_fmt = sprintf("%.2f (%.2f\u2013%.2f)", OR, CI_lo, CI_hi)
    )

results %>% arrange(desc(n_carriers)) %>%
    select(variant_id, Gene, n_carriers, OR_fmt, p_age, status)

# ========================
# RS7705526 ASSOCIATION TEST
# ========================
rs <- read_tsv(file.path("ch", "data", "rs7705526_output.raw")) %>%
    select(Sample.ID = IID, rs7705526 = rs7705526_C) %>%
    filter(!is.na(rs7705526))

# merge with age
all_samples_cov_rs <- cov %>%
    filter(Batch == 1) %>%
    rename(Sample.ID = person_id) %>%
    select(Sample.ID, Age = all_of("Sample_age")) %>%
    left_join(rs, by = "Sample.ID") %>%
    filter(!is.na(rs7705526))

cat("Batch 1 samples with rs7705526:", nrow(all_samples_cov_rs), "\n")
cat("Samples with rs7705526:", sum(!is.na(all_samples_cov_rs$rs7705526)), "\n")

# association test — same structure as age test
artifact_test_rs <- function(vkey, all_samples_cov_rs) {
    carriers <- all_ch_clean %>%
        filter(variant_id == vkey) %>%
        distinct(Sample.ID) %>%
        mutate(carrier = 1)

    df <- all_samples_cov_rs %>%
        left_join(carriers, by = "Sample.ID") %>%
        mutate(carrier = ifelse(is.na(carrier), 0, carrier)) %>%
        filter(!is.na(rs7705526))

    n_carriers <- sum(df$carrier)

    # skip if too few carriers to fit reliably
    if (n_carriers < 3) {
        return(data.frame(OR_rs = NA_real_, CI_lo_rs = NA_real_,
                          CI_hi_rs = NA_real_, p_rs = NA_real_))
    }

    tryCatch({
        m  <- suppressWarnings(glm(carrier ~ rs7705526, data = df,
                                   family = binomial(),
                                   control = glm.control(maxit = 100)))  # more iterations

        # check convergence explicitly
        if (!m$converged) {
            return(data.frame(OR_rs = NA_real_, CI_lo_rs = NA_real_,
                              CI_hi_rs = NA_real_, p_rs = NA_real_))
        }

        s  <- coef(summary(m))
        ci <- confint.default(m)
        data.frame(
            OR_rs    = exp(coef(m)["rs7705526"]),
            CI_lo_rs = exp(ci["rs7705526", 1]),
            CI_hi_rs = exp(ci["rs7705526", 2]),
            p_rs     = s["rs7705526", "Pr(>|z|)"]
        )
    }, error = function(e) {
        data.frame(OR_rs = NA_real_, CI_lo_rs = NA_real_,
                   CI_hi_rs = NA_real_, p_rs = NA_real_)
    })
}

# run on same high_freq_variants
results_rs <- high_freq_variants %>%
    rowwise() %>%
    mutate(res = list(artifact_test_rs(variant_id, all_samples_cov_rs))) %>%
    unnest(res) %>%
    ungroup()

# combine with age results and apply joint retention rule:
# retain if EITHER age OR rs association p < 0.10
results_combined <- results %>%
    left_join(results_rs %>% select(variant_id, OR_rs, CI_lo_rs, CI_hi_rs, p_rs),
              by = "variant_id") %>%
    mutate(
        OR_rs_fmt = sprintf("%.2f (%.2f\u2013%.2f)", OR_rs, CI_lo_rs, CI_hi_rs),
        status = case_when(
            is.na(p_age) & is.na(p_rs)          ~ "artifact",
            !is.na(p_age) & p_age < 0.10        ~ "retained",
            !is.na(p_rs)  & p_rs  < 0.10        ~ "retained",
            TRUE                                 ~ "artifact"
        )
    )

results_combined %>%
    arrange(desc(n_carriers)) %>%
    select(variant_id, Gene, n_carriers, OR_fmt, p_age, OR_rs_fmt, p_rs, status) %>%
    print(n = 30)


# ========================
# FOREST PLOT
# ========================
results_plot <- results %>%
    left_join(
        all_ch %>% select(variant_id, HGVSp) %>% distinct(),
        by = "variant_id"
    ) %>%
    mutate(
        label = ifelse(!is.na(HGVSp) & HGVSp != "" & HGVSp != ".",
                       paste0(Gene, ":", HGVSp),
                       variant_id),
        label = factor(label, levels = label[order(OR)])
    )

p_artifact <- ggplot(results_plot, aes(x = OR, y = label, color = status)) +
    geom_vline(xintercept = 1, linetype = "dotted", color = "gray40") +
    geom_errorbar(aes(xmin = CI_lo, xmax = CI_hi), height = 0.3, linewidth = 0.6, orientation = "y") +
    geom_point(size = 2.5, shape = 15) +
    scale_color_manual(
        values = c("retained" = "steelblue", "artifact" = "firebrick"),
        labels = c("retained" = "Associated with age (retained)",
                   "artifact" = "Not associated (flagged)")
    ) +
    labs(
        title    = "Association of common variants with age",
        subtitle = paste0("Variants in ≥", FREQ_THRESHOLD, " individuals | P < 0.10 threshold"),
        x = "Odds ratio for association with age",
        y = NULL, color = NULL
    ) +
    theme_minimal() +
    theme(
        legend.position    = "bottom",
        axis.text.y        = element_text(size = 7),
        panel.grid.major.y = element_blank()
    )

ggsave(file.path("ch", "figures", "artifact_age_association.pdf"),
       plot = p_artifact, width = 8,
       height = max(5, nrow(results_plot) * 0.3), dpi = 300)

# ========================
# REMOVE ARTIFACTS AND SAVE
# ========================
artifacts <- results %>% filter(status == "artifact")

all_ch_clean <- all_ch %>% anti_join(artifacts, by = "variant_id")

cat("Variants before artifact removal:", nrow(all_ch), "\n")
cat("Variants after artifact removal: ", nrow(all_ch_clean), "\n")
# 268
# 230

write_xlsx(all_ch_clean, file.path("ch", "data", "ch_wl_art.xlsx"))
write.csv(results, file.path("ch", "data", "ch_artifact_test_results.csv"), row.names = FALSE)

# ========================
# SUMMARY FIGUREts
# ========================
# cat_summary <- all_ch_clean %>%
#     filter(variant_category != "non_whitelist") %>%
#     count(variant_category, name = "N") %>%
#     mutate(
#         variant_category = factor(variant_category,
#                                   levels = c("splice", "lof", "missense",
#                                              "exception", "manual_review")),
#         fill_color = case_when(
#             variant_category == "manual_review" ~ "#C0392B",
#             variant_category == "exception"     ~ "#BA7517",
#             TRUE                                ~ "steelblue"
#         )
#     )

plot_category_breakdown(all_ch_clean, title = "CHIP variant categories (artifact-filtered)",
                        filename = "ch_wl_art_categories.pdf")
plot_genes(all_ch_clean, title = "Variants per gene (artifact-filtered)",
           filename = "ch_wl_art_genes.pdf", width = 12)
plot_genes(all_ch_clean, top_n = 10, title = "Top 10 genes (artifact-filtered)",
           filename = "ch_wl_art_genes_top10.pdf", width = 7)

# ========================
# STANDARD PLOTTING FUNCTIONS
# ========================

plot_category_breakdown <- function(df, title = "CHIP variant categories", filename = NULL, width = 6, height = 5) {
    cat_summary <- df %>%
        filter(!is.na(variant_category), variant_category != "non_whitelist") %>%
        count(variant_category, name = "N") %>%
        mutate(
            variant_category = factor(variant_category,
                                      levels = c("splice", "lof", "missense", "exception", "manual_review")),
            fill_color = case_when(
                variant_category == "manual_review" ~ "#C0392B",
                variant_category == "exception"     ~ "#BA7517",
                TRUE                                ~ "steelblue"
            )
        )

    p <- ggplot(cat_summary, aes(x = variant_category, y = N, fill = variant_category)) +
        geom_col(width = 0.65) +
        geom_text(aes(label = paste0(N, "\n(", round(100 * N / sum(N), 1), "%)")),
                  vjust = -0.4, size = 3.2, color = "gray30") +
        scale_fill_manual(values = setNames(cat_summary$fill_color, cat_summary$variant_category)) +
        scale_y_continuous(expand = expansion(mult = c(0, 0.25))) +
        scale_x_discrete(labels = c(
            "splice"        = "Splice",
            "lof"           = "LoF / Frameshift",
            "missense"      = "Missense",
            "exception"     = "Gene-specific\nexception",
            "manual_review" = "Manual review"
        )) +
        labs(title = title, x = NULL, y = "Number of variants",
             subtitle = paste0("Total: ", nrow(df))) +
        theme_minimal(base_size = 12) +
        theme(legend.position = "none", panel.grid.major.x = element_blank(),
              panel.grid.minor = element_blank(),
              plot.title    = element_text(face = "bold", size = 13),
              plot.subtitle = element_text(size = 10, color = "gray50"),
              axis.text.x   = element_text(size = 11))

    if (!is.null(filename))
        ggsave(file.path("ch", "figures", filename), p, width = width, height = height)
    p
}

plot_genes <- function(df, title = "CHIP variants by gene", top_n = NULL,
                       fill = "lightblue", filename = NULL, width = 10, height = 5) {
    gene_counts <- df %>%
        filter(!is.na(Gene)) %>%
        count(Gene, name = "N") %>%
        arrange(desc(N))

    if (!is.null(top_n)) gene_counts <- slice_max(gene_counts, N, n = top_n)

    gene_counts <- gene_counts %>%
        mutate(Gene = factor(Gene, levels = (Gene)))  # descending left to right

    p <- ggplot(gene_counts, aes(x = Gene, y = N)) +
        geom_col(width = 0.65, fill = fill) +
        geom_text(aes(label = N), vjust = -0.4, size = 3.2, color = "gray30") +
        scale_y_continuous(expand = expansion(mult = c(0, 0.25))) +
        labs(title = title, x = NULL, y = "Number of variants",
             subtitle = paste0("Total: ", sum(gene_counts$N),
                               " | Genes: ", nrow(gene_counts))) +
        theme_minimal(base_size = 12) +
        theme(panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
              plot.title    = element_text(face = "bold", size = 13),
              plot.subtitle = element_text(size = 10, color = "gray50"),
              axis.text.x   = element_text(size = 9, angle = 45, hjust = 1))

    if (!is.null(filename))
        ggsave(file.path("ch", "figures", filename), p, width = width, height = height)
    p
}

plot_gene_exceptions <- function(df, title = "Gene-specific exceptions", filename = NULL, width = 4, height = 4) {
    gene_exc <- df %>%
        filter(wl.exception) %>%
        count(Gene, name = "N") %>%
        # arrange(desc(N)) %>%
        mutate(Gene = factor(Gene, levels = Gene[order(-N)]))

    p <- ggplot(gene_exc, aes(x = Gene, y = N)) +
        geom_col(width = 0.6, fill = "#BA7517") +
        geom_text(aes(label = N), vjust = -0.4, size = 3.2, color = "gray30") +
        scale_y_continuous(expand = expansion(mult = c(0, 0.25))) +
        labs(title = title, x = NULL, y = "Number of variants") +
        theme_minimal(base_size = 12) +
        theme(panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
              plot.title = element_text(face = "bold", size = 12),
              axis.text.x = element_text(size = 11, angle = 45, hjust = 1))

    if (!is.null(filename))
        ggsave(file.path("ch", "figures", filename), p, width = width, height = height)
    p
}

# ========================
# EXAMINE ASXL1
# ========================
x <- read.csv(file.path("ch", "data", "ch_all_vars.csv"), row.names = 1)
y <- read_excel(file.path("ch", "data", "ch_wl_final.xlsx"))
z <- read_excel(file.path("ch", "data", "ch_wl_art.xlsx"))

View(x %>% filter(Gene == "ASXL1"))
View(y %>% filter(Gene == "ASXL1"))
View(z %>% filter(Gene == "ASXL1"))

