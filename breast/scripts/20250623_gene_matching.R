# Gene Matching Script for GWAS Annotations
# Matches genes from annotated variants with reference GWAS hits

library(dplyr)
library(stringr)
library(readr)

# Function to clean and split gene names
clean_gene_names <- function(gene_string) {
    if (is.na(gene_string) || gene_string == "" || gene_string == ".") {
        return(character(0))
    }

    # Split by comma and clean whitespace
    genes <- str_split(gene_string, ",")[[1]]
    genes <- str_trim(genes)

    # Remove empty strings
    genes <- genes[genes != "" & genes != "."]

    return(genes)
}

# Read the annotated file
cat("Reading annotated file...\n")
annotated <- read_csv("your_annotated_file.csv",
                      col_types = cols(.default = "c")) # Read all as character to preserve formatting

# Read the reference file
cat("Reading reference file...\n")
reference <- read_csv("your_reference_file.csv",
                      col_types = cols(.default = "c"))

# Clean column names (remove any special characters)
colnames(annotated) <- make.names(colnames(annotated))
colnames(reference) <- make.names(colnames(reference))

# Display file information
cat("Annotated file dimensions:", nrow(annotated), "x", ncol(annotated), "\n")
cat("Reference file dimensions:", nrow(reference), "x", ncol(reference), "\n")

# Preview column names
cat("Annotated file columns:", paste(head(colnames(annotated), 10), collapse = ", "), "...\n")
cat("Reference file columns:", paste(colnames(reference), collapse = ", "), "\n")

# Extract unique genes from reference file
cat("Extracting reference genes...\n")
reference_genes <- unique(reference$MAPPED_GENE[!is.na(reference$MAPPED_GENE)])
reference_genes <- reference_genes[reference_genes != "" & reference_genes != "."]

# Handle genes with special characters (like LINC01488 - CCND1)
reference_genes_expanded <- c()
for (gene in reference_genes) {
    if (str_detect(gene, " - ")) {
        # Split genes connected by " - " and add both
        split_genes <- str_split(gene, " - ")[[1]]
        split_genes <- str_trim(split_genes)
        reference_genes_expanded <- c(reference_genes_expanded, split_genes, gene)
    } else {
        reference_genes_expanded <- c(reference_genes_expanded, gene)
    }
}

reference_genes_expanded <- unique(reference_genes_expanded)
cat("Total unique reference genes (including split):", length(reference_genes_expanded), "\n")

# Function to check if any genes in a comma-separated list match reference
check_gene_match <- function(gene_string, ref_genes) {
    variant_genes <- clean_gene_names(gene_string)
    if (length(variant_genes) == 0) return(FALSE)

    # Check if any variant gene matches any reference gene
    return(any(variant_genes %in% ref_genes))
}

# Find matching variants
cat("Finding gene matches...\n")
annotated$has_match <- sapply(annotated$Gene.refGene,
                              function(x) check_gene_match(x, reference_genes_expanded))

# Filter to only matching variants
matching_variants <- annotated[annotated$has_match, ]

cat("Found", nrow(matching_variants), "matching variants\n")

if (nrow(matching_variants) > 0) {
    # For each matching variant, find which specific genes match
    results <- data.frame()

    for (i in 1:nrow(matching_variants)) {
        variant <- matching_variants[i, ]
        variant_genes <- clean_gene_names(variant$Gene.refGene)

        # Find which genes match
        matching_genes <- variant_genes[variant_genes %in% reference_genes_expanded]

        for (gene in matching_genes) {
            # Find corresponding reference entries
            # Handle both exact matches and partial matches for complex gene names
            ref_matches <- reference[reference$MAPPED_GENE == gene |
                                         str_detect(reference$MAPPED_GENE, paste0("\\b", gene, "\\b")), ]

            if (nrow(ref_matches) > 0) {
                for (j in 1:nrow(ref_matches)) {
                    result_row <- data.frame(
                        # From annotated file
                        Chr = variant$Chr,
                        Start = variant$Start,
                        End = variant$End,
                        Ref = variant$Ref,
                        Alt = variant$Alt,
                        Gene.refGene = gene,  # The specific matching gene

                        # From reference file
                        SNPS = ref_matches$SNPS[j],
                        CHR_ID = ref_matches$CHR_ID[j],
                        CHR_POS = ref_matches$CHR_POS[j],
                        MAPPED_GENE = ref_matches$MAPPED_GENE[j],
                        stringsAsFactors = FALSE
                    )
                    results <- rbind(results, result_row)
                }
            }
        }
    }

    # Display results
    cat("\n=== MATCHING RESULTS ===\n")
    cat("Total matches found:", nrow(results), "\n\n")

    if (nrow(results) > 0) {
        # Show summary by gene
        gene_summary <- results %>%
            group_by(MAPPED_GENE) %>%
            summarise(
                variant_count = n(),
                unique_variants = n_distinct(paste(Chr, Start, Ref, Alt)),
                .groups = 'drop'
            ) %>%
            arrange(desc(variant_count))

        cat("Summary by gene:\n")
        print(gene_summary)

        cat("\nFirst 10 matches:\n")
        print(head(results, 10))

        # Write results to file
        output_file <- "gene_matching_results.csv"
        write_csv(results, output_file)
        cat("\nResults written to:", output_file, "\n")

        # Optional: Write summary file
        summary_file <- "gene_matching_summary.csv"
        write_csv(gene_summary, summary_file)
        cat("Summary written to:", summary_file, "\n")

    }
} else {
    cat("No matching variants found.\n")
    cat("Reference genes available:", paste(head(reference_genes_expanded, 20), collapse = ", "), "...\n")
    cat("Sample annotated genes:", paste(head(unique(annotated$Gene.refGene), 20), collapse = ", "), "...\n")
}

cat("\nScript completed!\n")
