# Plot mean depth as a function of distance from the nearest target boundary,
# to decide empirically how much slop/padding a targets BED actually needs.
#
# Inputs (per sample, from run_mosdepth_edge.sh):
#   <outdir>/<sample>.edge_window.bed  -- chrom, start, end, depth (mosdepth per-base, RLE)
#   targets.merged.bed                 -- chr-fixed, merged, UNPADDED targets (from process_targets_bed.sh)
#
# Usage: Rscript edge_effect_plot.R <outdir> <sample1> [<sample2> ...]

library(data.table)
library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)
outdir  <- args[1]
samples <- args[-1]

targets <- fread("/home/jennyzli/jenny/gray/bed/IDT-Target-HYB-PreCancer.hg38.targets.merged.bed", col.names = c("chrom", "start", "end", "name"))

edge_distance <- function(pos, chrom_edges) {
  edges <- sort(chrom_edges)
  idx <- findInterval(pos, edges)
  left  <- ifelse(idx >= 1, edges[pmax(idx, 1)], NA_real_)
  right <- ifelse(idx + 1 <= length(edges), edges[pmin(idx + 1, length(edges))], NA_real_)
  dist_left  <- pos - left
  dist_right <- right - pos
  dist <- pmin(dist_left, dist_right, na.rm = TRUE)
  inside <- (idx %% 2L) == 1L   # odd number of edges below pos => inside a target
  ifelse(inside, dist, -dist)
}

all_binned <- list()

for (sample in samples) {
  f <- file.path(outdir, paste0(sample, ".edge_window.bed"))
  cov <- fread(f, col.names = c("chrom", "start", "end", "depth"))

  # Expand RLE intervals to 1bp rows (fine at this scale -- targets +-100bp only).
  width <- cov$end - cov$start
  pos   <- unlist(mapply(function(s, e) s:(e - 1L), cov$start, cov$end, SIMPLIFY = FALSE))
  chrom <- rep(cov$chrom, width)
  depth <- rep(cov$depth, width)
  long  <- data.table(chrom, pos, depth)

  long[, dist := {
    d <- numeric(.N)
    for (ch in unique(chrom)) {
      ix <- chrom == ch
      ch_edges <- c(targets[chrom == ch, start], targets[chrom == ch, end])
      d[ix] <- edge_distance(pos[ix], ch_edges)
    }
    d
  }]

  long[, sample := sample]
  all_binned[[sample]] <- long
}

combined <- rbindlist(all_binned)

# Bin distance into 5bp buckets and summarize depth per bin per sample.
combined[, bin := floor(dist / 5) * 5]
summary_dt <- combined[, .(mean_depth = mean(depth), n = .N), by = .(sample, bin)]

p <- ggplot(summary_dt, aes(x = bin, y = mean_depth, color = sample)) +
  geom_line() +
  geom_point(size = 0.8) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  labs(
    x = "Distance from nearest target edge (bp, negative = outside target)",
    y = "Mean depth",
    title = "Capture edge-effect: depth vs. distance from target boundary"
  ) +
  theme_minimal()

sample_tag <- sprintf("n%d_samples", length(samples))
png_file <- file.path(outdir, paste0(sample_tag, ".edge_effect.png"))
csv_file <- file.path(outdir, paste0(sample_tag, ".edge_effect_summary.csv"))

ggsave(png_file, p, width = 8, height = 5, dpi = 150)
fwrite(summary_dt, csv_file)