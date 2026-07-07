library(biomaRt)
library(readxl)

library(biomaRt)
library(readxl)

# Read the Excel file
ref <- read_excel("~/Downloads/REF_GWAS.xlsx")

# Extract SNP IDs
rs_probes <- ref$`Best published SNP`

# Create SNP mart for hg38/GRCh38
snpmart <- useEnsembl(biomart = "snp", dataset = "hsapiens_snp")

# Remove any NA values and ensure rs IDs are clean
valid_snps <- rs_probes[!is.na(rs_probes) & rs_probes != "" & rs_probes != "NA"]

# Clean up rs IDs (remove any extra spaces or characters)
valid_snps <- gsub("^rs", "", valid_snps)  # Remove 'rs' prefix if present
valid_snps <- paste0("rs", valid_snps)      # Add it back consistently

print(paste("Querying", length(valid_snps), "SNPs for hg38 coordinates"))

# Query biomaRt to get hg38 positions from rs IDs
snp_results <- getBM(
    attributes = c('refsnp_id', 'chr_name', 'chrom_start', 'chrom_end', 'chrom_strand', 'allele'),
    filters = 'snp_filter',
    values = valid_snps,
    mart = snpmart
)

print(paste("Found", nrow(snp_results), "results"))

# Merge results back with original data
ref_with_hg38 <- merge(ref, snp_results,
                       by.x = "Best published SNP",
                       by.y = "refsnp_id",
                       all.x = TRUE)

write_xlsx(ref_with_hg38, "~/Documents/Nathanson/breast/data/REF_GWAS.xlsx")


ref1 <- ref[order(ref$Chr1), ]

head(ref)



# Create hg38 ID column
ref_with_hg38$hg38_ID <- paste(ref_with_hg38$chr_name, ref_with_hg38$chrom_start, sep = ":")

# Display comparison
comparison <- ref_with_hg38[, c("Best published SNP", "Chr1", "Position2",
                                "chr_name", "chrom_start", "ID", "hg38_ID")]
colnames(comparison) <- c("SNP_ID", "hg37_Chr", "hg37_Pos",
                          "hg38_Chr", "hg38_Pos", "hg37_ID", "hg38_ID")

print("Coordinate comparison (hg37 vs hg38):")
head(comparison, 10)









ref <- read_excel("~/Downloads/REF_GWAS.xlsx")

rs_probes <- ref$`Best published SNP`

snpmart <- useEnsembl(biomart = "snp", dataset = "hsapiens_snp")

ref$ID <- paste(ref$Chr1, ref$Position2, sep = ":")

valid_snps <- rs_probes[!is.na(rs_probes) & rs_probes != ""]

snp_results <- getBM(
    attributes = c('refsnp_id', 'allele', 'chr_name', 'chrom_start', 'chrom_strand'),
    filters = 'snp_filter',  # Use snp_filter for rs IDs
    values = valid_snps,
    mart = snpmart
)

head(snp_results)
ref$SNP <- snp_results
write.csv(snp_results)

# Method 2: Alternative - Query by chromosomal positions
# If you want to query by chromosome and position instead:
unique_positions <- unique(ref[!is.na(ref$Chr1) & !is.na(ref$Position2), c("Chr1", "Position2")])

if(nrow(unique_positions) > 0) {
    position_results <- getBM(
        attributes = c('refsnp_id', 'allele', 'chr_name', 'chrom_start', 'chrom_strand'),
        filters = c('chr_name', 'start', 'end'),
        values = list(
            chr_name = unique_positions$Chr1,
            start = unique_positions$Position2,
            end = unique_positions$Position2
        ),
        mart = snpmart
    )

    print("Position Query Results:")
    head(position_results)
}

# Check available filters and attributes if needed
print("Available filters:")
head(listFilters(snpmart), 10)

print("Available attributes:")
head(listAttributes(snpmart), 10)





