#!/usr/bin/env Rscript
library(data.table)
library(ggplot2)
library(dplyr)
library(biomaRt)
library(data.table)
library(here)

# ========================
# GET REFERENCE FROM EMBL
# ========================

ref <- fread(here("breast", "data","gwas-association-downloaded_2025-07-22-MONDO_0007254-withChildTraits.tsv"))
#  3132   38
ref <- as.data.frame(ref)
cols <- c("SNPS","CHR_ID","CHR_POS", "MAPPED_GENE", "CONTEXT", "REGION", "P-VALUE","DISEASE/TRAIT", "INITIAL SAMPLE SIZE" )
ref1 <- ref[,cols]

ref_p <- ref1 %>% arrange(`P-VALUE`, )
ref_pos <-ref1 %>% arrange(as.numeric(`CHR_ID`), as.numeric(`CHR_POS`))

write.csv(ref_p, file.path("breast", "data", "ref_breast_hits_pval.csv"))
write.csv(ref_pos, file.path("breast", "data", "ref_breast_hits_pos.csv"))

ref_a_pos <- ref_p %>%
    filter(grepl("Africa", `INITIAL SAMPLE SIZE`, ignore.case = TRUE)) %>%
    arrange(as.numeric(`CHR_ID`), as.numeric(`CHR_POS`))
ref_a_p <- ref_p %>%
    filter(grepl("Africa", `INITIAL SAMPLE SIZE`, ignore.case = TRUE)) %>%
    arrange(`P-VALUE`, )
# write.csv(ref_p, file.path("breast", "data", "afr_breast_cancer_hits.csv"))


# Sort data by LOG10P in descending order and select the columns in the required order
# eur_sorted <- eur[order(-LOG10P), .(CHROM, GENPOS, ID, LOG10P)]
# afr_sorted <- afr[order(-LOG10P), .(CHROM, GENPOS, ID, LOG10P)]
#
# # Save sorted data to CSV
# write.csv(eur_sorted, here("breast", "data", "EUR_genome_wide_1e5_sorted.csv"), row.names = FALSE)
# write.csv(afr_sorted, here("breast", "data", "AFR_genome_wide_1e5_sorted.csv"), row.names = FALSE)
#

# ========================
# EUR HITS V. EMBL
# ========================

setwd("~/Documents/Nathanson")
afr <- read_csv(here("breast", "data", "2", "AFR", "AFR_genome_maf05_p1e4.hg38_multianno.csv"),
                col_types = cols(.default = "c"))
eur <- read_csv(here("breast", "data", "2", "EUR", "EUR_genome_maf05_p1e4_annovar.hg38_multianno.csv"),
                      col_types = cols(.default = "c"))
ref <- read_csv(here("breast", "data", "reference_breast_cancer_hits.csv"),
                col_types = cols(.default = "c"))

# clean_gene_names <- function(gene_string) {
#     if (is.na(gene_string) || gene_string == "" || gene_string == ".") {
#         return(character(0))
#     }
#
#     # Split by comma and clean whitespace
#     genes <- str_split(gene_string, ",")[[1]]
#     genes <- str_trim(genes)
#
#     # Remove empty strings
#     genes <- genes[genes != "" & genes != "."]
#
#     return(genes)
# }
#
# colnames(afr) <- make.names(colnames(afr))
# colnames(eur) <- make.names(colnames(eur))
# colnames(ref) <- make.names(colnames(ref))
#
# # Display file information
# # cat("Annotated file dimensions:", nrow(eur), "x", ncol(eur), "\n")
# # cat("ref file dimensions:", nrow(ref), "x", ncol(ref), "\n")
#
# # Preview column names
# # cat("Annotated file columns:", paste(head(colnames(eur), 10), collapse = ", "), "...\n")
# # cat("Reference file columns:", paste(colnames(ref), collapse = ", "), "\n")

# ========================
# EUR HITS V. EMBL
# ========================

setwd("~/Documents/Nathanson")
afr <- read_tsv(here("breast", "data", "2", "AFR", "AFR_genome_maf05_p1e4_anno.txt"),
                col_types = cols(.default = "c"))
eur <- read_tsv(here("breast", "data", "2", "EUR", "EUR_genome_maf05_p1e4_anno.txt"),
                col_types = cols(.default = "c"))
all <- read_tsv(here("breast", "data", "2", "ALL", "all_genome_maf05_p1e4_anno.txt"),
                col_types = cols(.default = "c"))
ref <- read_csv(here("breast", "data", "reference_breast_cancer_hits.csv"),
                col_types = cols(.default = "c"))

col <- c("CHROM", "GENPOS", "avsnp150", "ID", "A1FREQ", "LOG10P", "Func.refGene", "Gene.refGene")
afr <- afr[,col]
eur <- eur[,col]
all <- all[,col]

# sort by pvals
eur_p <- eur %>% arrange(desc(LOG10P))
afr_p <- afr %>% arrange(desc(LOG10P))
all_p <- afr %>% arrange(desc(LOG10P))

eur_pos <- eur %>% arrange(as.numeric(CHROM), as.numeric(GENPOS))
afr_pos <- afr %>% arrange(as.numeric(CHROM), as.numeric(GENPOS))
all_pos <- all %>% arrange(as.numeric(CHROM), as.numeric(GENPOS))


write_csv(eur_p, here("breast", "data", "2", "EUR", "EUR_genome_maf05_p1e4_anno_pval.csv"))
write_csv(afr_p, here("breast", "data", "2", "AFR", "AFR_genome_maf05_p1e4_anno_pval.csv"))
write_csv(all_p, here("breast", "data", "2", "ALL", "ALL_genome_maf05_p1e4_anno_pval.csv"))

write_csv(eur_pos, here("breast", "data", "2", "EUR", "EUR_genome_maf05_p1e4_anno_pos.csv"))
write_csv(afr_pos, here("breast", "data", "2", "AFR", "AFR_genome_maf05_p1e4_anno_pos.csv"))
write_csv(all_pos, here("breast", "data", "2", "ALL", "ALL_genome_maf05_p1e4_anno_pos.csv"))


# ==============================================================================
# GENE MATCHING PIPELINE FOR EUR AND AFR DATASETS
# ==============================================================================

clean_gene_names <- function(gene_string) {
    if (is.na(gene_string) || gene_string == "" || gene_string == ".") {
        return(character(0))
    }
    genes <- str_split(gene_string, ",")[[1]]
    genes <- str_trim(genes)
    genes <- genes[genes != "" & genes != "."]
    return(genes)
}

check_gene_match <- function(gene_string, ref_genes) {
    variant_genes <- clean_gene_names(gene_string)
    if (length(variant_genes) == 0) return(FALSE)
    return(any(variant_genes %in% ref_genes))
}

# ==============================================================================
# STEP 1: PREPARE REFERENCE GENES
# ==============================================================================

# Extract unique genes from reference
ref_genes <- unique(ref$MAPPED_GENE[!is.na(ref$MAPPED_GENE)])
ref_genes <- ref_genes[ref_genes != "" & ref_genes != "."]

# Expand genes (handle commas and dashes)
ref_genes_expanded <- c()
for (gene in ref_genes) {
    if (str_detect(gene, ",")) {
        comma_split <- str_split(gene, ",")[[1]]
        comma_split <- str_trim(comma_split)

        for (comma_gene in comma_split) {
            if (str_detect(comma_gene, " - ")) {
                dash_split <- str_split(comma_gene, " - ")[[1]]
                dash_split <- str_trim(dash_split)
                ref_genes_expanded <- c(ref_genes_expanded, dash_split)
            } else {
                ref_genes_expanded <- c(ref_genes_expanded, comma_gene)
            }
        }
    } else if (str_detect(gene, " - ")) {
        split_genes <- str_split(gene, " - ")[[1]]
        split_genes <- str_trim(split_genes)
        ref_genes_expanded <- c(ref_genes_expanded, split_genes)
    } else {
        ref_genes_expanded <- c(ref_genes_expanded, gene)
    }
}

ref_genes_expanded <- ref_genes_expanded[ref_genes_expanded != "" &
                                             ref_genes_expanded != "." &
                                             !is.na(ref_genes_expanded)]
ref_genes_expanded <- unique(ref_genes_expanded)

cat("Original reference genes:", length(ref_genes), "\n")
cat("Expanded reference genes:", length(ref_genes_expanded), "\n")

# ==============================================================================
# STEP 2: FIND MATCHING VARIANTS
# ==============================================================================
write.csv(eur$avsnp150, file.path("breast", "data", "2", "EUR", "eur_probes.csv"))

# Add match flags
eur$has_match <- sapply(eur$Gene.refGene, function(x) check_gene_match(x, ref_genes_expanded))
afr$has_match <- sapply(afr$Gene.refGene, function(x) check_gene_match(x, ref_genes_expanded))
all$has_match <- sapply(all$Gene.refGene, function(x) check_gene_match(x, ref_genes_expanded))

# Get matching variants
eur_matching <- eur[eur$has_match, ]
afr_matching <- afr[afr$has_match, ]
all_matching <- all[all$has_match, ]

write.csv(eur_matching, file.path("breast", "data", "2", "EUR", "eur_matching.csv"))
write.csv(afr_matching, file.path("breast", "data", "2", "AFR", "afr_matching.csv"))
write.csv(all_matching, file.path("breast", "data", "2", "ALL", "all_matching.csv"))


# ==============================================================================
# EUR
# ==============================================================================
eur_results <- data.frame()

if (nrow(eur_matching) > 0) {
    for (i in 1:nrow(eur_matching)) {
        variant <- eur_matching[i, ]
        variant_genes <- clean_gene_names(variant$Gene.refGene)
        matching_genes <- variant_genes[variant_genes %in% ref_genes_expanded]

        for (gene in matching_genes) {
            ref_matches <- ref[ref$MAPPED_GENE == gene |
                                   str_detect(ref$MAPPED_GENE, paste0("\\b", gene, "\\b")), ]

            for (j in 1:nrow(ref_matches)) {
                result_row <- data.frame(
                    Chr = variant$Chr,
                    Start = variant$Start,
                    End = variant$End,
                    Ref = variant$Ref,
                    Alt = variant$Alt,
                    Gene.refGene = gene,
                    SNPS = ref_matches$SNPS[j],
                    CHR_ID = ref_matches$CHR_ID[j],
                    CHR_POS = ref_matches$CHR_POS[j],
                    MAPPED_GENE = ref_matches$MAPPED_GENE[j],
                    Dataset = "EUR",
                    stringsAsFactors = FALSE
                )
            eur_results <- rbind(eur_results, result_row)
            }
        }
    }
}

# ==============================================================================
# AFR
# ==============================================================================
afr_results <- data.frame()

if (nrow(afr_matching) > 0) {
    for (i in 1:nrow(afr_matching)) {
        variant <- afr_matching[i, ]
        variant_genes <- clean_gene_names(variant$Gene.refGene)
        matching_genes <- variant_genes[variant_genes %in% ref_genes_expanded]

        for (gene in matching_genes) {
            ref_matches <- ref[ref$MAPPED_GENE == gene |
                                   str_detect(ref$MAPPED_GENE, paste0("\\b", gene, "\\b")), ]

            if (nrow(ref_matches) > 0) {
                for (j in 1:nrow(ref_matches)) {
                    result_row <- data.frame(
                        Chr = variant$Chr,
                        Start = variant$Start,
                        End = variant$End,
                        Ref = variant$Ref,
                        Alt = variant$Alt,
                        Gene.refGene = gene,
                        SNPS = ref_matches$SNPS[j],
                        CHR_ID = ref_matches$CHR_ID[j],
                        CHR_POS = ref_matches$CHR_POS[j],
                        MAPPED_GENE = ref_matches$MAPPED_GENE[j],
                        Dataset = "AFR",
                        stringsAsFactors = FALSE
                    )
                    afr_results <- rbind(afr_results, result_row)
                }
            }
        }
    }
}
# ==============================================================================
# STEP 5: CREATE SUMMARIES
# ==============================================================================
# EUR summary
if (nrow(eur_results) > 0) {
    eur_summary <- eur_results %>%
        group_by(MAPPED_GENE) %>%
        summarise(
            variant_count = n(),
            unique_variants = n_distinct(paste(Chr, Start, Ref, Alt)),
            dataset = "EUR",
            .groups = 'drop'
        ) %>%
        arrange(desc(variant_count))

    cat("EUR summary:\n")
    print(eur_summary)
} else {
    eur_summary <- data.frame()
}

# AFR summary
if (nrow(afr_results) > 0) {
    afr_summary <- afr_results %>%
        group_by(MAPPED_GENE) %>%
        summarise(
            variant_count = n(),
            unique_variants = n_distinct(paste(Chr, Start, Ref, Alt)),
            dataset = "AFR",
            .groups = 'drop'
        ) %>%
        arrange(desc(variant_count))

    cat("\nAFR summary:\n")
    print(afr_summary)
} else {
    afr_summary <- data.frame()
}

### SAVE ###
output_dir <- here("breast", "data", "2")

eur_dir <- file.path(output_dir, "EUR")
write_csv(eur_results, file.path(eur_dir, "EUR_gene_matching_results.csv"))
write_csv(eur_summary, file.path(eur_dir, "EUR_gene_matching_summary.csv"))

afr_dir <- file.path(output_dir, "AFR")
write_csv(afr_results, file.path(afr_dir, "AFR_gene_matching_results.csv"))
write_csv(afr_summary, file.path(afr_dir, "AFR_gene_matching_summary.csv"))

# # Combine and save
# if (nrow(eur_results) > 0 || nrow(afr_results) > 0) {
#     combined_results <- rbind(eur_results, afr_results)
#     combined_summary <- rbind(eur_summary, afr_summary)
#
#     write_csv(combined_results, file.path(output_dir, "combined_gene_matching_results.csv"))
#     write_csv(combined_summary, file.path(output_dir, "combined_gene_matching_summary.csv"))
#     cat("Combined results saved to:", output_dir, "\n")
# }
