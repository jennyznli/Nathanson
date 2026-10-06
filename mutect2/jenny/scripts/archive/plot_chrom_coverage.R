# Whole-chromosome coverage overview, plus a full-resolution zoom into the
# highest-depth window(s).
#
# Inputs:
#   <outdir>/<sample>.<chrom>.regions.bed.gz  -- windowed mean depth (from run_mosdepth_chrom.sh)
#   <outdir>/<sample>.per-base.bed.gz         -- full-resolution per-base depth (from run_mosdepth_edge.sh)
#
# Usage: Rscript plot_chrom_coverage.R <outdir> <sample> <chrom> [n_peaks] [flank_bp]

library(data.table)
library(ggplot2)

args     <- commandArgs(trailingOnly = TRUE)
outdir   <- args[1]
sample   <- args[2]
chrom    <- args[3]
n_peaks  <- if (length(args) >= 4) as.integer(args[4]) else 1L
flank_bp <- if (length(args) >= 5) as.integer(args[5]) else 2000L

regions <- fread(file.path(outdir, paste0(sample, ".", chrom, ".regions.bed.gz")),
                  col.names = c("chrom", "start", "end", "depth"))

# --- whole-chromosome overview ---
p_overview <- ggplot(regions, aes(x = (start + end) / 2, y = depth)) +
  geom_line(linewidth = 0.3) +
  labs(x = paste(chrom, "position (bp)"), y = "Mean depth (windowed)",
       title = paste("Coverage overview:", sample, chrom)) +
  theme_minimal()

ggsave(file.path(outdir, paste0(sample, ".", chrom, ".overview.png")),
       p_overview, width = 10, height = 4, dpi = 150)

# --- find the highest-depth window(s) to zoom into ---
peaks <- regions[order(-depth)][seq_len(n_peaks)]
cat("Top", n_peaks, "highest-depth window(s) in", chrom, ":\n")
print(peaks)

# --- zoomed, full-resolution view around the single highest peak ---
top        <- peaks[1]
zoom_start <- max(0, top$start - flank_bp)
zoom_end   <- top$end + flank_bp

perbase <- fread(file.path(outdir, paste0(sample, ".per-base.bed.gz")),
                  col.names = c("chrom", "start", "end", "depth"))
zoom <- perbase[chrom == ..chrom & end > zoom_start & start < zoom_end]

p_zoom <- ggplot(zoom, aes(x = start, xend = end, y = depth, yend = depth)) +
  geom_segment(linewidth = 0.6, color = "firebrick") +
  geom_vline(xintercept = c(top$start, top$end), linetype = "dashed", color = "grey40") +
  labs(x = paste(chrom, "position (bp)"), y = "Depth (per-base)",
       title = paste0("Zoom: ", sample, " ", chrom, ":", zoom_start, "-", zoom_end)) +
  theme_minimal()

ggsave(file.path(outdir, paste0(sample, ".", chrom, ".zoom_top1.png")),
       p_zoom, width = 10, height = 4, dpi = 150)

cat("Wrote", file.path(outdir, paste0(sample, ".", chrom, ".overview.png")), "\n")
cat("Wrote", file.path(outdir, paste0(sample, ".", chrom, ".zoom_top1.png")), "\n")
