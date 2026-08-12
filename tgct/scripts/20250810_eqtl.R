# ============================================================
# CNA ANALYSIS - BATCH PROCESSING
# ============================================================
source("config.R")
setwd(here("tgct"))

# ========================
# LOAD DATA & CONVERT
# ========================
meta  <- read_xlsx(file.path("ss", "meta-hits-062526.xlsx"), sheet = "NovelHits")
# tgct <- read_xlsx(file.path("tgct_pooled.xlsx"))
# ss <- read_tsv(file.path("pennsamples_masterslist_062226.txt"))
# flags <- fread(file.path("rgcname_pmbbid_metadata_flags_freeze4.csv"), header = TRUE)

# ========================
# SELECT IDs
# ========================
new  <- read_xlsx(file.path("ss", "meta-hits-062526.xlsx"), sheet = "NovelHits")
old  <- read.table(file.path("data", "emilio-scores.txt"), header = TRUE)

old_ids <- paste0("chr", old$MarkerName)
new_ids <- paste0("chr", new$MarkerName)

writeLines(new_ids, file.path("data", "novel_hits_ids.txt"))
writeLines(old_ids, file.path("data", "old_hits_ids.txt"))


# ========================
# GTExr
# ========================
install.packages("gtexr")
get_eqtl_genes("Whole_Blood")
get_significant_single_tissue_eqtls(gencodeId = c(
    "ENSG00000132693.12",
    "ENSG00000203782.5"
))

get_significant_single_tissue_eqtls_by_location(
    tissueSiteDetailId,
    start,
    end,
    chromosome,
    datasetId = "gtex_v8",
    .return_raw = FALSE
)

