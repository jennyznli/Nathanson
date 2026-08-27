#!/usr/bin/env Rscript
# Snakemake-callable wrapper around qc_germline_gene()
# (R/germline_qc.R) for a single gene's VEP report CSV.
#
# Usage:
#   Rscript run_germline_qc_cli.R \
#     --input <gene>.vep.report.csv \
#     --gene <GENE> \
#     --output_csv <out>_qc.csv \
#     --log_dir <log dir> \
#     [--af_lower 0.18] [--af_upper 0.80] [--min_depth 14] [--min_altdepth 5]

suppressMessages(library(here))
suppressMessages(library(dplyr))
source(here("gray", "scripts", "germline_qc.R"))

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NA_character_) {
    idx <- match(flag, args)
    if (is.na(idx) || idx == length(args)) return(default)
    args[idx + 1]
}

input_csv <- get_arg("--input")
gene_name <- get_arg("--gene")
output_csv <- get_arg("--output_csv")
log_dir <- get_arg("--log_dir", file.path("gray", "log"))

af_lower <- as.numeric(get_arg("--af_lower", "0.18"))
af_upper <- as.numeric(get_arg("--af_upper", "0.80"))
min_depth <- as.numeric(get_arg("--min_depth", "14"))
min_altdepth <- as.numeric(get_arg("--min_altdepth", "5"))

if (any(is.na(c(input_csv, gene_name, output_csv)))) {
    stop("Usage: Rscript run_germline_qc_cli.R --input <csv> --gene <GENE> ",
         "--output_csv <out.csv> [--log_dir <dir>] ",
         "[--af_lower] [--af_upper] [--min_depth] [--min_altdepth]")
}

data <- read.csv(input_csv, stringsAsFactors = FALSE)

res <- qc_germline_gene(
    data, gene_name, dataset_label = gene_name,
    af_lower = af_lower, af_upper = af_upper,
    min_depth = min_depth, min_altdepth = min_altdepth,
    log_dir = log_dir
)

dir.create(dirname(output_csv), showWarnings = FALSE, recursive = TRUE)
write.csv(res$step3, output_csv, row.names = FALSE)

message("Wrote ", nrow(res$step3), " QC'd ", gene_name, " variants -> ", output_csv)
