source(file.path("R", "config.R"))

# ========================
# GENE COORDINATES W/ MANE TRANSCRIPT
# ========================
mane <- read.delim(file.path("ch", "data", "MANE.GRCh38.v1.5.summary.txt"), header = TRUE)
genes <- read_excel(file.path("ch", "ss", "ch_genes.xlsx"), col_names = "Gene")
gene_list <- unique(genes$Gene)
# 60 genes
mane_coords$Gene[!(mane_coords$Gene %in% gene_list)]
# none so that's good

refseq_to_chr <- function(refseq_acc) {
    chr_num <- as.integer(sub("NC_0*([0-9]+)\\..*", "\\1", refseq_acc))
    ifelse(chr_num == 23, "X",
           ifelse(chr_num == 24, "Y", chr_num))
}

mane_coords <- mane %>%
    filter(symbol %in% gene_list) %>%
    mutate(Chromosome = sapply(GRCh38_chr, refseq_to_chr)) %>%  # Apply function to each element of GRCh38_chr
    group_by(symbol, Chromosome) %>%
    summarize(
        Min_Start = min(as.numeric(chr_start)),
        Max_End = max(as.numeric(chr_end)),
        N_Transcripts = n(),
        RefSeq_IDs = paste(RefSeq_nuc, collapse = ", "),
        .groups = "drop"
    ) %>%
    rename(Gene = symbol)
print(mane_coords)

pad_bp <- 50000
mane_pad <- mane_coords %>%
    mutate(
        Start_Padded = pmax(1L, as.integer(Min_Start) - pad_bp),
        End_Padded   = as.integer(Max_End) + pad_bp,
        Region       = paste0("chr", Chromosome, ":", Start_Padded, "-", End_Padded)
    ) %>%
    select(Gene, Min_Start, Max_End, Start_Padded, End_Padded, Region)

mane_final <- mane_pad %>%
    arrange(Gene) %>%
    mutate(line = paste0("    ", Gene, ': "', Region, '"')) %>%
    pull(line)

writeLines(mane_final, file.path("ch", "data", "ch_genes_dict.txt"))






