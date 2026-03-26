library(here)
setwd("/Users/jennyzli/Documents/Nathanson")
source(here("R", "config.R"))
library(data.table, quietly=T)

# ========================
# read in data
# ========================
vars <- read.csv(file.path("ch", "data", "pilot_mutect2_ch_lof1.csv")) %>%
    select("vep_sample_id", "Chr", "Start", "REF", "ALT", "Gene", "Variant.Class", "Variant.Consequence", "Bio.type",
           "HGVSc", "HGVSp", "EXON", "INTRON", "Sample.Zyg", "Sample.Depth", "Sample.AltDepth", "Sample.AltFrac",
           "caller")

gList<-fread(file.path("ch", "data", "whitelist_filter_20230531", "NEJM_2017_genes_01262020.txt"))
whitelist.mis<-fread(file.path("ch", "data", "whitelist_filter_20230531", "CHIP_missense_vars_cv_04102022.txt"))
whitelist.splice<-fread(file.path("ch", "data", "whitelist_filter_20230531", "CHIP_splice_vars_agb_01262020.txt"))
whitelist.LoF<-fread(file.path("ch", "data", "whitelist_filter_20230531", "CHIP_nonsense_FS_vars_agb_01262020.txt"))

# ========================
# PRELIM QC - Whitelist Style
# ========================
process_ch <- function(data, dataset_label, AF_LOWER, AF_UPPER, DP, AD) {

    # Generate datetime stamp and log file
    datetime_stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
    log_file <- file.path("ch", "log", paste0(datetime_stamp, "_ch_qc_", dataset_label, ".txt"))
    dir.create(file.path("ch", "log"), showWarnings = FALSE, recursive = TRUE)

    # Initialize flags
    data <- data %>%
        mutate(
            pass_protein_coding = FALSE,
            pass_vaf = FALSE,
            pass_depth = FALSE,
            qc_pass = FALSE
        )

    # Step 1: Protein coding filter
    data$pass_protein_coding <- (data$Bio.type == "protein_coding") &
        (data$Sample.AltFrac != ".") &
        (data$Sample.Depth != ".") &
        (data$Sample.AltDepth != ".")

    # Convert to numeric
    data <- data %>%
        mutate(
            Sample.AltFrac = as.numeric(Sample.AltFrac),
            Sample.Depth = as.numeric(Sample.Depth),
            Sample.AltDepth = as.numeric(Sample.AltDepth)
        )

    # Remove -1 values from passing
    data$pass_protein_coding <- data$pass_protein_coding &
        !is.na(data$Sample.AltFrac) &
        !is.na(data$Sample.Depth) &
        !is.na(data$Sample.AltDepth) &
        (data$Sample.AltFrac != -1) &
        (data$Sample.Depth != -1) &
        (data$Sample.AltDepth != -1)

    # Step 2: VAF filter (only apply to those that passed step 1)
    data$pass_vaf <- data$pass_protein_coding &
        (data$Sample.AltFrac >= AF_LOWER) &
        (data$Sample.AltFrac <= AF_UPPER)

    # Step 3: Depth filter (only apply to those that passed step 2)
    data$pass_depth <- data$pass_vaf &
        (data$Sample.Depth >= DP) &
        (data$Sample.AltDepth >= AD)

    # Final QC pass
    data$qc_pass <- data$pass_depth

    # Log summary
    sink(log_file)
    cat("=== POST-VEP GENE QC LOG ===\n")
    cat("Dataset:", dataset_label, "\n\n")
    cat("=== FILTERING SUMMARY ===\n")
    cat("Initial variants:", nrow(data), "\n")
    cat("Pass protein coding:", sum(data$pass_protein_coding),
        sprintf("(%.1f%%)\n", 100 * sum(data$pass_protein_coding) / nrow(data)))
    cat("Pass VAF filter:", sum(data$pass_vaf),
        sprintf("(%.1f%%)\n", 100 * sum(data$pass_vaf) / nrow(data)))
    cat("Pass depth filter:", sum(data$pass_depth),
        sprintf("(%.1f%%)\n", 100 * sum(data$pass_depth) / nrow(data)))
    cat("\nFinal QC pass:", sum(data$qc_pass),
        sprintf("(%.1f%% of initial)\n\n", 100 * sum(data$qc_pass) / nrow(data)))
    sink()

    cat("Log saved to:", log_file, "\n")

    return(data)
}

# Run QC
vars <- process_ch(vars, "pilot_mutect2_ch", 0.02, 0.45, 15, 2)

# Examine variants that failed at each step
failed_protein_coding <- vars %>% filter(!pass_protein_coding)
failed_vaf <- vars %>% filter(pass_protein_coding & !pass_vaf)
failed_depth <- vars %>% filter(pass_vaf & !pass_depth)

# Work with passing variants
all_ch <- vars %>% filter(qc_pass)

cat("\nVariants removed at each step:\n")
cat("Failed protein coding:", nrow(failed_protein_coding), "\n")
cat("Failed VAF:", nrow(failed_vaf), "\n")
cat("Failed depth:", nrow(failed_depth), "\n")
cat("Passed all QC:", nrow(all_ch), "\n")

# ========================
# Extract protein change from HGVSp
# ========================
# VEP format: "p.Arg573Gln" or "p.R573Q"
# Extract everything after "p."
extractProteinChange <- function(hgvsp) {
    if (is.na(hgvsp) || hgvsp == "" || hgvsp == ".") {
        return(NA)
    }
    # Remove "p." prefix and decode URL encoding
    prot_change <- gsub("^p\\.", "", hgvsp)
    prot_change <- gsub("%3D", "=", prot_change)
    return(prot_change)
}

all_ch$ProteinChange <- sapply(all_ch$HGVSp, extractProteinChange)

# ========================
# Initialize whitelist flags
# ========================
all_ch <- all_ch %>%
    mutate(
        whitelist = FALSE,
        wl.mis = FALSE,
        wl.lof = FALSE,
        wl.splice = FALSE,
        wl.exception = FALSE,
        manualreview = FALSE
    )

# ========================
# 1) Handle MISSENSE variants
# ========================
vmis <- grepl("missense", all_ch$Variant.Consequence)

convert_1to3 <- function(aa_change) {
    if (is.na(aa_change) || aa_change == "" || aa_change == ".") {
        return(NA)
    }

    aa_map <- c(
        "A" = "Ala", "C" = "Cys", "D" = "Asp", "E" = "Glu", "F" = "Phe",
        "G" = "Gly", "H" = "His", "I" = "Ile", "K" = "Lys", "L" = "Leu",
        "M" = "Met", "N" = "Asn", "P" = "Pro", "Q" = "Gln", "R" = "Arg",
        "S" = "Ser", "T" = "Thr", "V" = "Val", "W" = "Trp", "Y" = "Tyr",
        "*" = "Ter", "=" = "="
    )

    # Pattern: one letter + digits + one letter/symbol
    pattern <- "^([A-Z*])(\\d+)([A-Z*=])$"

    if (grepl(pattern, aa_change)) {
        parts <- str_match(aa_change, pattern)
        ref_aa <- parts[2]
        position <- parts[3]
        alt_aa <- parts[4]

        # Convert
        ref_3 <- ifelse(ref_aa %in% names(aa_map), aa_map[ref_aa], ref_aa)
        alt_3 <- ifelse(alt_aa %in% names(aa_map), aa_map[alt_aa], alt_aa)

        return(paste0(ref_3, position, alt_3))
    } else {
        return(aa_change)
    }
}

whitelist.mis$AAChange_3letter <- sapply(whitelist.mis$AAChange, convert_1to3)

# Check if in missense whitelist
vmis_wl <- paste(all_ch$Gene, all_ch$ProteinChange, sep = "_") %in%
    paste(whitelist.mis$Gene, whitelist.mis$AAChange_3letter, sep = "_")

# Flag missense whitelist variants
all_ch$whitelist[vmis & vmis_wl] <- TRUE
all_ch$wl.mis[vmis & vmis_wl] <- TRUE

cat(sprintf("Missense variants in whitelist: %d\n", sum(vmis & vmis_wl)))
# 3

# ========================
# 2) Handle LoF and frameshift variants
# ========================
# Stop gain/loss/start lost OR frameshift
vlof <- grepl("Ter|\\*|X", all_ch$ProteinChange) |
    grepl("stop_gain|stop_lost|start_lost", all_ch$Variant.Consequence) |
    grepl("fs", all_ch$ProteinChange, fixed = TRUE) |
    grepl("frameshift", all_ch$Variant.Consequence)

# Check if gene is in whitelist
vLOFgene <- all_ch$Gene %in% whitelist.LoF$Gene

# Flag LoF variants
all_ch$whitelist[vlof & vLOFgene] <- TRUE
all_ch$wl.lof[vlof & vLOFgene] <- TRUE

cat(sprintf("LoF/frameshift variants in whitelist: %d\n",
            sum(vlof & vLOFgene, na.rm = TRUE)))
# 8

# ========================
# 3) Handle SPLICE variants
# ========================
vSplice <- grepl("splice", all_ch$Variant.Consequence)

# Check if gene is in splice whitelist
vSplicegene <- all_ch$Gene %in% whitelist.splice$Gene

# have to confirm correct transcript after we get that information -- should do so beforehand
# #confirm that the splicing refers to the correct transcript
# vSpliceCorrectTranscript<-apply(varsOI.func[,c('GeneDetail.refGene','Accession')],
#                                 1, function(x) {grepl(x[2],x[1],fixed=T)})
#
# # Flag splice variants
# all_ch$whitelist[vSplice & vSplicegene & vSpliceCorrectTranscript] <- TRUE
# all_ch$wl.splice[vSplice & vSplicegene & vSpliceCorrectTranscript] <- TRUE

# Flag splice variants
all_ch$whitelist[vSplice & vSplicegene] <- TRUE
all_ch$wl.splice[vSplice & vSplicegene] <- TRUE

# # Flag for manual review if splice gene but wrong transcript (shouldn't happen after canonical filter)
# all_ch$manualreview[(vSplice & vSplicegene) & (!vSpliceCorrectTranscript)] <- TRUE

cat(sprintf("Splice variants in whitelist: %d\n",
            sum(vSplice & vSplicegene)))
# 0

# ========================
# Summary of whitelist filtering
# ========================
cat("\n=== WHITELIST FILTERING SUMMARY ===\n")
cat(sprintf("Variants in ANY whitelist: %d (%.1f%%)\n",
            sum(all_ch$whitelist),
            100 * sum(all_ch$whitelist) / nrow(all_ch)))
cat(sprintf("  - Missense whitelist: %d\n", sum(all_ch$wl.mis)))
cat(sprintf("  - LoF whitelist: %d\n", sum(all_ch$wl.lof)))
cat(sprintf("  - Splice whitelist: %d\n", sum(all_ch$wl.splice)))
cat(sprintf("  - Flagged for manual review: %d\n", sum(all_ch$manualreview)))
# Variants in ANY whitelist: 11 (73.3%)
# > cat(sprintf("  - Missense whitelist: %d\n", sum(all_ch$wl.mis)))
# - Missense whitelist: 3
# > cat(sprintf("  - LoF whitelist: %d\n", sum(all_ch$wl.lof)))
# - LoF whitelist: 9
# > cat(sprintf("  - Splice whitelist: %d\n", sum(all_ch$wl.splice)))
# - Splice whitelist: 0
# > cat(sprintf("  - Flagged for manual review: %d\n", sum(all_ch$manualreview)))
# - Flagged for manual review: 0

# ========================
# GENE SPECIFIC RULES
# ========================
# Extract exon information from EXON column
# VEP format: "exon_number|total_exons" e.g., "16|22"
extract_exon_number <- function(exon_field) {
    if (is.na(exon_field) || exon_field == "" || exon_field == ".") {
        return(NA)
    }
    # Extract first number before the pipe
    exon_num <- as.numeric(gsub("\\|.*", "", exon_field))
    return(exon_num)
}

all_ch$ExonNumber <- sapply(all_ch$EXON, extract_exon_number)

# Extract amino acid position from ProteinChange for range-based filters
extract_aa_position <- function(prot_change) {
    if (is.na(prot_change) || prot_change == "" || prot_change == ".") {
        return(NA)
    }
    # Extract digits from protein change (e.g., "Arg573Gln" -> 573)
    position <- as.numeric(gsub("\\D", "", prot_change))
    return(position)
}

all_ch$AAPosition <- sapply(all_ch$ProteinChange, extract_aa_position)

# Define LoF and frameshift variants for exceptions
vlof <- grepl("Ter|\\*|X", all_ch$ProteinChange) |
    grepl("stop_gain|stop_lost|start_lost", all_ch$Variant.Consequence) |
    grepl("fs", all_ch$ProteinChange, fixed = TRUE) |
    grepl("frameshift", all_ch$Variant.Consequence)
vFS <- grepl("fs", all_ch$ProteinChange, fixed = TRUE) |
    grepl("frameshift_variant", all_ch$Variant.Consequence)
vSplice <- grepl("splice", all_ch$Variant.Consequence)

# need to check this elater!!
# vSpliceCorrectTranscript <- all_ch$RefSeqID %in% gList$Accession
vSpliceCorrectTranscript = TRUE # temporary


# ========================
# 4) Handle specific gene exceptions
# ========================
#### ASXL1 ###
# Frameshift/nonsense/splice-site in exon 11-12
vexon11 <- all_ch$ExonNumber == 11
vexon12 <- all_ch$ExonNumber == 12
asxl1Exception <- (all_ch$Gene == "ASXL1") & (vlof | vFS) & (vexon11 | vexon12)

all_ch$whitelist[asxl1Exception] <- TRUE
all_ch$wl.lof[asxl1Exception] <- TRUE
all_ch$wl.exception[asxl1Exception] <- TRUE

# #Splicing
asxl1ExceptionSplice <- (all_ch$Gene == "ASXL1") &
    vSplice & vSpliceCorrectTranscript &
    (vexon11 | vexon12)
all_ch$whitelist[asxl1ExceptionSplice] <- TRUE
all_ch$wl.splice[asxl1ExceptionSplice] <- TRUE
all_ch$wl.exception[asxl1ExceptionSplice] <- TRUE

cat(sprintf("ASXL1 exceptions: %d\n", sum(asxl1Exception | asxl1ExceptionSplice, na.rm = TRUE)))

### ASXL2: Frameshift/nonsense/splice-site in exon 11-12
# LOF
asxl2Exception <- (all_ch$Gene == "ASXL2") & (vlof | vFS) & (vexon11 | vexon12)
all_ch$whitelist[asxl2Exception] <- TRUE
all_ch$wl.lof[asxl2Exception] <- TRUE
all_ch$wl.exception[asxl2Exception] <- TRUE

# Splicing
asxl2ExceptionSplice <- (all_ch$Gene == "ASXL2") &
    vSplice & vSpliceCorrectTranscript &
    (vexon11 | vexon12)
all_ch$whitelist[asxl2ExceptionSplice] <- TRUE
all_ch$wl.splice[asxl2ExceptionSplice] <- TRUE
all_ch$wl.exception[asxl2ExceptionSplice] <- TRUE

cat(sprintf("ASXL2 exceptions: %d\n", sum(asxl2Exception | asxl2ExceptionSplice, na.rm = TRUE)))

#### PPM1D: Frameshift/nonsense in exon 5 or 6
vexon5 <- all_ch$ExonNumber == 5
vexon6 <- all_ch$ExonNumber == 6

ppm1dException <- (all_ch$Gene == "PPM1D") & (vlof | vFS) & (vexon5 | vexon6)
all_ch$whitelist[ppm1dException] <- TRUE
all_ch$wl.lof[ppm1dException] <- TRUE
all_ch$wl.exception[ppm1dException] <- TRUE

cat(sprintf("PPM1D exceptions: %d\n", sum(ppm1dException, na.rm = TRUE)))

### TET2: Missense mutations in catalytic domains (p.1104-1481 and 1843-2002)
vmis <- grepl("missense", all_ch$Variant.Consequence)
TETidx <- which(all_ch$Gene == "TET2" &
                    vmis &
                    !is.na(all_ch$AAPosition))

for(i in TETidx) {
    AApos <- all_ch$AAPosition[i]
    if(!is.na(AApos) && ((AApos >= 1104 & AApos <= 1481) | (AApos >= 1843 & AApos <= 2002))) {
        all_ch$whitelist[i] <- TRUE
        all_ch$wl.mis[i] <- TRUE
        all_ch$wl.exception[i] <- TRUE
    }
}

cat(sprintf("TET2 catalytic domain exceptions: %d\n",
            sum(all_ch$wl.exception & all_ch$Gene == "TET2", na.rm = TRUE)))

# CBL: RING finger missense p.381-421
CBLidx <- which(all_ch$Gene == "CBL" &
                    vmis &
                    !is.na(all_ch$AAPosition))

for(i in CBLidx) {
    AApos <- all_ch$AAPosition[i]
    if(!is.na(AApos) && AApos >= 381 & AApos <= 421) {
        all_ch$whitelist[i] <- TRUE
        all_ch$wl.mis[i] <- TRUE
        all_ch$wl.exception[i] <- TRUE
    }
}

cat(sprintf("CBL RING finger exceptions: %d\n",
            sum(all_ch$wl.exception & all_ch$Gene == "CBL", na.rm = TRUE)))

# CBLB: RING finger missense p.372-412
CBLBidx <- which(all_ch$Gene == "CBLB" &
                     vmis &
                     !is.na(all_ch$AAPosition))

for(i in CBLBidx) {
    AApos <- all_ch$AAPosition[i]
    if(!is.na(AApos) && AApos >= 372 & AApos <= 412) {
        all_ch$whitelist[i] <- TRUE
        all_ch$wl.mis[i] <- TRUE
        all_ch$wl.exception[i] <- TRUE
    }
}

cat(sprintf("CBLB RING finger exceptions: %d\n",
            sum(all_ch$wl.exception & all_ch$Gene == "CBLB", na.rm = TRUE)))

# ========================
# 5) Flag remaining exceptions for manual review
# could ptoentially code these in later if there's too many...
# ========================
# Non-frameshift indels (inframe insertions/deletions)
vNonFrameshiftDel <- grepl("inframe_deletion", all_ch$Variant.Consequence)
vNonFrameshiftIns <- grepl("inframe_insertion", all_ch$Variant.Consequence)

# GATA3: Frameshift/nonsense/splice-site ZNF domain
all_ch$manualreview[(vlof | vFS | vSplice) & (all_ch$Gene == "GATA3")] <- TRUE

# CREBBP: S1680del
all_ch$manualreview[vNonFrameshiftDel & (all_ch$Gene == "CREBBP")] <- TRUE

# CSF3R: truncating c.741-791
all_ch$manualreview[(vlof | vFS | vSplice) & (all_ch$Gene == "CSF3R")] <- TRUE

# DNMT3A: F732del, F752del
all_ch$manualreview[vNonFrameshiftDel & (all_ch$Gene == "DNMT3A")] <- TRUE

# EP300: VF1148_1149del
all_ch$manualreview[vNonFrameshiftDel & (all_ch$Gene == "EP300")] <- TRUE

# FLT3: FY590-591GD, del835
all_ch$manualreview[(vNonFrameshiftDel | vNonFrameshiftIns) & (all_ch$Gene == "FLT3")] <- TRUE

# JAK2: del/ins537-539L, del/ins538-539L, del/ins540-543MK, del/ins540-544MK, del542-543, del543-544, ins11546-547
all_ch$manualreview[(vNonFrameshiftDel | vNonFrameshiftIns) & (all_ch$Gene == "JAK2")] <- TRUE

# KDM6A: del419
all_ch$manualreview[vNonFrameshiftDel & (all_ch$Gene == "KDM6A")] <- TRUE

# KIT: ins503, del560, del579, del551-559
all_ch$manualreview[(vNonFrameshiftDel | vNonFrameshiftIns) & (all_ch$Gene == "KIT")] <- TRUE

# MPL: del513 W515-518KT
all_ch$manualreview[(vNonFrameshiftDel | vNonFrameshiftIns) & (all_ch$Gene == "MPL")] <- TRUE

# NPM1: Frameshift p.W288fs (insertion at c.859_860, 860_861, 862_863, 863_864)
vFrameshiftIndel <- grepl("frameshift_variant", all_ch$Variant.Consequence)
all_ch$manualreview[vFrameshiftIndel & (all_ch$Gene == "NPM1")] <- TRUE

cat("\n=== MANUAL REVIEW FLAGS ===\n")
cat(sprintf("Total flagged for manual review: %d\n", sum(all_ch$manualreview, na.rm = TRUE)))

# ========================
# Final summary and output
# ========================

cat("\n=== FINAL WHITELIST + EXCEPTION SUMMARY ===\n")
cat(sprintf("Total variants: %d\n", nrow(all_ch)))
cat(sprintf("Variants in whitelist (including exceptions): %d (%.1f%%)\n",
            sum(all_ch$whitelist),
            100 * sum(all_ch$whitelist) / nrow(all_ch)))
cat(sprintf("  - Missense whitelist: %d\n", sum(all_ch$wl.mis)))
cat(sprintf("  - LoF whitelist: %d\n", sum(all_ch$wl.lof)))
cat(sprintf("  - Splice whitelist: %d\n", sum(all_ch$wl.splice)))
cat(sprintf("  - Exception rules: %d\n", sum(all_ch$wl.exception)))
cat(sprintf("  - Flagged for manual review: %d\n", sum(all_ch$manualreview)))

# ========================
# Write output files
# ========================
datetime_stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
dir.create(file.path("ch", "output"), showWarnings = FALSE, recursive = TRUE)

check_vars <- data.frame(
    total_num_variants = nrow(all_ch),
    total_num_whitelist = sum(all_ch$whitelist),
    total_num_manualreview = sum(all_ch$manualreview)
)

# Write output files
write.csv(check_vars,
          file.path("ch", "data", paste0(datetime_stamp, "_pilot_mutect2_varcount.csv")),
          row.names = FALSE)

write.csv(all_ch,
          file.path("ch", "data", paste0(datetime_stamp, "_pilot_mutect2_allvariants.csv")),
          row.names = FALSE)

write.csv(all_ch[all_ch$whitelist, ],
          file.path("ch", "output", paste0(datetime_stamp, "_pilot_mutect2_wl.csv")),
          row.names = FALSE)

write.csv(all_ch[all_ch$manualreview, ],
          file.path("ch", "output", paste0(datetime_stamp, "_pilot_mutect2_manualreview.csv")),
          row.names = FALSE)


# # ========================
# # PRELIM QC
# # ========================
# process_ch <- function(data, dataset_label, AF_LOWER, AF_UPPER,
#                        DP, AD) {
#
#     # Generate datetime stamp
#     datetime_stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
#
#     # Set up gene-specific log file with datetime
#     log_file <- file.path("ch", "log",
#                           paste0(datetime_stamp, "_ch_qc_", dataset_label, ".txt"))
#     dir.create(file.path("ch", "log"), showWarnings = FALSE, recursive = TRUE)
#
#     # Function to write to both console and log
#     log_cat <- function(...) {
#         message <- paste0(...)
#         cat(message)
#         cat(message, file = log_file, append = TRUE)
#     }
#
#     # Function to calculate and log stats (excluding -1 values)
#     log_stats <- function(data, label) {
#         log_cat(label, ":\n")
#
#         # VAF stats (exclude -1)
#         vaf_clean <- data$Sample.AltFrac[data$Sample.AltFrac != -1]
#         if (length(vaf_clean) > 0) {
#             vaf_stats <- summary(vaf_clean)
#             log_cat("  VAF - Min: ", round(vaf_stats[1], 3),
#                     ", Median: ", round(vaf_stats[3], 3),
#                     ", Max: ", round(vaf_stats[6], 3), "\n")
#         } else {
#             log_cat("  VAF - No valid data\n")
#         }
#
#         # Total Depth stats (exclude -1)
#         depth_clean <- data$Sample.Depth[data$Sample.Depth != -1]
#         if (length(depth_clean) > 0) {
#             depth_stats <- summary(depth_clean)
#             log_cat("  Total Depth - Min: ", round(depth_stats[1], 1),
#                     ", Median: ", round(depth_stats[3], 1),
#                     ", Max: ", round(depth_stats[6], 1), "\n")
#         } else {
#             log_cat("  Total Depth - No valid data\n")
#         }
#
#         # Alt Depth stats (exclude -1)
#         altdepth_clean <- data$Sample.AltDepth[data$Sample.AltDepth != -1]
#         if (length(altdepth_clean) > 0) {
#             altdepth_stats <- summary(altdepth_clean)
#             log_cat("  Alt Depth - Min: ", round(altdepth_stats[1], 1),
#                     ", Median: ", round(altdepth_stats[3], 1),
#                     ", Max: ", round(altdepth_stats[6], 1), "\n")
#         } else {
#             log_cat("  Alt Depth - No valid data\n")
#         }
#     }
#
#     # Start fresh log
#     cat("", file = log_file, append = FALSE)
#     log_cat("=== POST-VEP GENE QC LOG ===\n")
#     log_cat("Dataset: ", dataset_label, "\n\n")
#
#     # Initial file info
#     log_cat("=== INITIAL DATA ===\n")
#     log_cat("Total rows: ", nrow(data), "\n")
#     log_cat("Total individuals: ", length(unique(data$Sample.ID)), "\n")
#     log_cat("Genomic range: ", min(data$Start, na.rm = TRUE), " - ", max(data$Start, na.rm = TRUE), "\n\n")
#
#     remove_pattern <- c("synonymous", "upstream", "3_prime_UTR_variant",
#                         "5_prime_UTR_variant", "inframe")
#
#     # Step 1: Gene and protein coding filter
#     log_cat("=== FILTERING STEPS ===\n\n")
#
#     step1 <- data %>%
#         filter(Bio.type == "protein_coding") %>%
#         filter(Sample.AltFrac != ".", Sample.Depth != ".", Sample.AltDepth != ".") %>%
#         mutate(
#             Sample.AltFrac = as.numeric(Sample.AltFrac),
#             Sample.Depth = as.numeric(Sample.Depth),
#             Sample.AltDepth = as.numeric(Sample.AltDepth),
#         ) %>%
#         filter(!is.na(Sample.AltFrac), !is.na(Sample.Depth), !is.na(Sample.AltDepth)) %>%
#         filter(Sample.AltFrac != -1, Sample.Depth != -1, Sample.AltDepth != -1)  # Remove -1 values
#
#     log_cat("Step 1 (Gene + Protein Coding): ",
#             nrow(step1), " variants, ",
#             length(unique(step1$Sample.ID)), " individuals\n")
#     log_stats(step1, "  Quality Metrics")
#     log_cat("\n")
#
#     # Step 2: VAF filter
#     step2 <- step1 %>%
#         filter(Sample.AltFrac >= AF_LOWER, Sample.AltFrac <= AF_UPPER)
#
#     log_cat("Step 2 (+ VAF Filter): ",
#             nrow(step2), " variants, ",
#             length(unique(step2$Sample.ID)), " individuals")
#     log_cat(" (removed ", nrow(step1) - nrow(step2), " variants)\n")
#     log_stats(step2, "  Quality Metrics")
#     log_cat("\n")
#
#     # Step 3: DP & AD filter
#     step3 <- step2 %>%
#         filter(Sample.Depth >= DP, Sample.AltDepth >= AD)
#
#     log_cat("Step 3 (+ Depth/AD Filter): ",
#             nrow(step3), " variants, ",
#             length(unique(step3$Sample.ID)), " individuals")
#     log_cat(" (removed ", nrow(step2) - nrow(step3), " variants)\n")  # FIXED: correct subtraction order
#     log_stats(step3, "  Quality Metrics")  # FIXED: log step3 stats instead of step2
#     log_cat("\n")
#
#     # Final summary
#     log_cat("=== FINAL RESULTS ===\n")
#     log_cat("Final genomic range: ", min(step3$Start, na.rm = TRUE), " - ", max(step3$Start, na.rm = TRUE), "\n\n")
#     log_cat("All consequences: ", paste(unique(step3$Variant.Consequence), collapse = ", "), "\n\n")
#
#     log_cat("\nLog saved to: ", log_file, "\n")
#
#     return(list(step1 = step1, step2 = step2, step3 = step3))
# }
#
# vars_qc1 <- process_ch(vars, "pilot_mutect2_ch", 0.02, 0.45, 18, 2)
# all_ch <- vars_qc1$step3


