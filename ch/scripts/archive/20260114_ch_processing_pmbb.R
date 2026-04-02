library(here)
setwd("/Users/jennyzli/Documents/Nathanson")
source(here("R", "config.R"))
source(here("ch", "config.R"))
library(data.table, quietly=T)

# ========================
# read in data
# ========================
gList<-fread(file.path("ch", "data", "whitelist_filter_20230531", "NEJM_2017_genes_01262020.txt"))
whitelist.mis<-fread(file.path("ch", "data", "whitelist_filter_20230531", "CHIP_missense_vars_cv_04102022.txt"))
whitelist.splice<-fread(file.path("ch", "data", "whitelist_filter_20230531", "CHIP_splice_vars_agb_01262020.txt"))
whitelist.LoF<-fread(file.path("ch", "data", "whitelist_filter_20230531", "CHIP_nonsense_FS_vars_agb_01262020.txt"))

all_ch <- read.csv(file.path("ch", "data", "ch_vep.all.qc1.csv")) %>% filter(MANE_SELECT != ".")
# should join by refseq...
# 4363119      62
unique(all_ch$Gene)
# [1] "ASXL1"   "BRAF"    "CBL"     "CREBBP"  "CSF3R"   "CTCF"    "CUX1"    "DNMT3A"  "EP300"   "GNAS"    "GNB1"
# [12] "MPL"     "NF1"     "PRPF40B" "PRPF8"   "SETD2"   "SETDB1"  "SRSF2"   "SUZ12"   "TET2"    "TP53"    "U2AF2"
length(unique(all_ch$Gene)) #22 genes

x <- all_ch %>% filter(MANE_SELECT == ".")
unique(x$Gene)
# "BRAF" "GNAS" -- these two don't have canonical?

table(all_ch$MANE_SELECT)
# .    NM_000546.6    NM_000760.4 NM_001031698.3 NM_001042492.3 NM_001127208.3 NM_001195427.2
# 308896         352330         206966           9498         687677         205387           7523
# NM_001366418.1    NM_001429.4    NM_002074.5    NM_004380.3    NM_005188.4    NM_005373.3    NM_006445.4
# 17468         313619          15895         118268         136381          55973         258486
# NM_006565.4    NM_007279.3    NM_014159.7    NM_015338.6    NM_015355.4    NM_022552.5    NM_181552.4
# 17420         310008         120436         175832         130739         238878         675439

# so two oft hem are missing...

x <- all_ch %>% filter(Gene == "BRAF")
x <- all_ch %>% filter(Gene == "GNAS")
# these don't have anything in canonical but have it in maneplusclinical?

# FILTER FUNCTIONAL VARIANTS - exonic and splicing...function Func.refGene

# ========================
# merge in canonical transcript
# ========================
# Extract RefSeq ID by removing version number (everything after the period)
all_ch$RefSeqID <- sapply(strsplit(all_ch$MANE_SELECT, "\\."), `[`, 1)

# Check extraction
table(all_ch$RefSeqID)

all_ch_canonical <- all_ch %>%
    filter(!is.na(RefSeqID)) %>%
    filter(RefSeqID %in% gList$Accession)
length(unique(all_ch_canonical$RefSeqID))

z <- all_ch %>%
    filter(!is.na(RefSeqID)) %>%
    filter(!(RefSeqID %in% gList$Accession))
unique(z$Gene)
# "BRAF"   "GNAS"   "NF1"    "SETDB1" "SRSF2"  "TP53"

x <- all_ch %>% filter(Gene == "NF1")
unique(x$RefSeqID)

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


all_ch_canonical$ProteinChange <- sapply(all_ch_canonical$HGVSp, extractProteinChange)


# ========================
# Initialize whitelist flags
# ========================
all_ch_canonical <- all_ch_canonical %>%
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
# VEP Consequence examples: "missense_variant", "missense_variant&splice_region_variant"
vmis <- grepl("missense_variant", all_ch_canonical$Variant.Consequence)

convert_1to3 <- function(aa_change) {
    if (is.na(aa_change) || aa_change == "" || aa_change == ".") {
        return(NA)
    }

    # Amino acid conversion
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
        # Handle special cases like duplications, deletions
        return(aa_change)
    }
}

# Apply
whitelist.mis$AAChange_3letter <- sapply(whitelist.mis$AAChange, convert_1to3)

# # Test
# head(data.frame(
#     original = whitelist.mis$AAChange,
#     converted = whitelist.mis$AAChange_3letter
# ), 10)

# Check if in missense whitelist
vmis_wl <- paste(all_ch_canonical$Gene, all_ch_canonical$ProteinChange, sep = "_") %in%
    paste(whitelist.mis$Gene, whitelist.mis$AAChange_3letter, sep = "_")

# Flag missense whitelist variants
all_ch_canonical$whitelist[vmis & vmis_wl] <- TRUE
all_ch_canonical$wl.mis[vmis & vmis_wl] <- TRUE

cat(sprintf("Missense variants in whitelist: %d\n", sum(vmis & vmis_wl)))
# 358...?

# ========================
# 2) Handle LoF and frameshift variants
# ========================
# Stop gain/loss: look for "X" or "*" in protein change
vlof <- grepl("X|\\*", all_ch_canonical$ProteinChange)

# OR use VEP consequence annotations
vlof_consequence <- grepl("stop_gain|stop_lost|start_lost", all_ch_canonical$Variant.Consequence)
sum(vlof_consequence)
# 284

# Frameshift: look for "fs" in protein change OR frameshift consequence
vFS <- grepl("fs", all_ch_canonical$ProteinChange, fixed = TRUE) |
    grepl("frameshift_variant", all_ch_canonical$Variant.Consequence)
# 548

# Check if gene is in LoF whitelist
vLOFgene <- all_ch_canonical$Gene %in% whitelist.LoF$Gene

# Flag LoF variants
all_ch_canonical$whitelist[(vlof | vlof_consequence | vFS) & vLOFgene] <- TRUE
all_ch_canonical$wl.lof[(vlof | vlof_consequence | vFS) & vLOFgene] <- TRUE

cat(sprintf("LoF/frameshift variants in whitelist: %d\n",
            sum((vlof | vlof_consequence | vFS) & vLOFgene, na.rm = TRUE)))

# 663 ?

# ========================
# 3) Handle SPLICE variants
# ========================
# Examples: "splice_acceptor_variant", "splice_donor_variant", "splice_region_variant"
vSplice <- grepl("splice_acceptor_variant|splice_donor_variant",
                 all_ch_canonical$Variant.Consequence)

# Check if gene is in splice whitelist
vSplicegene <- all_ch_canonical$Gene %in% whitelist.splice$Gene

# Since we already filtered for canonical transcripts (MANE_SELECT),
# all splice variants should be on the correct transcript
# But we can double-check using Feature.Accession or MANE_SELECT
vSpliceCorrectTranscript <- all_ch_canonical$RefSeqID %in% gList$Accession

# Flag splice variants
all_ch_canonical$whitelist[vSplice & vSplicegene & vSpliceCorrectTranscript] <- TRUE
all_ch_canonical$wl.splice[vSplice & vSplicegene & vSpliceCorrectTranscript] <- TRUE

# Flag for manual review if splice gene but wrong transcript (shouldn't happen after canonical filter)
all_ch_canonical$manualreview[(vSplice & vSplicegene) & (!vSpliceCorrectTranscript)] <- TRUE

cat(sprintf("Splice variants in whitelist: %d\n",
            sum(vSplice & vSplicegene & vSpliceCorrectTranscript)))
# 76

# ========================
# Summary of whitelist filtering
# ========================
cat("\n=== WHITELIST FILTERING SUMMARY ===\n")
cat(sprintf("Total variants after canonical filter: %d\n", nrow(all_ch_canonical)))
cat(sprintf("Variants in ANY whitelist: %d (%.1f%%)\n",
            sum(all_ch_canonical$whitelist),
            100 * sum(all_ch_canonical$whitelist) / nrow(all_ch_canonical)))
cat(sprintf("  - Missense whitelist: %d\n", sum(all_ch_canonical$wl.mis)))
cat(sprintf("  - LoF whitelist: %d\n", sum(all_ch_canonical$wl.lof)))
cat(sprintf("  - Splice whitelist: %d\n", sum(all_ch_canonical$wl.splice)))
cat(sprintf("  - Flagged for manual review: %d\n", sum(all_ch_canonical$manualreview)))

### GENE SPECIFIC RULES
#LoF OR SPLICE VARIANTS RESTRICTED TO EXONS 11 AND 12
#4) Handle the following exceptions
#ASXL1	Frameshift/nonsense/splice-site in exon 11-12
vlof<-grepl("X",varsOI.func$NonsynOI, fixed=T) #stop gain or stop loss
vFS<-grepl("fs",varsOI.func$NonsynOI, fixed=T) #frameshift
vexon11<-grepl("exon11",varsOI.func$transcriptOI, fixed=T)
vexon12<-grepl("exon12",varsOI.func$transcriptOI, fixed=T)
asxl1Exception<-(varsOI.func$Gene.refGene=="ASXL1")&(vlof|vFS)&(vexon11|vexon12)
varsOI.func[asxl1Exception,"whitelist"]=T
varsOI.func[asxl1Exception,"wl.lof"]=T
varsOI.func[asxl1Exception,"wl.exception"]=T

#Splicing
vSplice<-grepl("splicing",varsOI.func$Func.refGene, fixed=T)
vSpliceCorrectTranscript<-apply(varsOI.func[,c('GeneDetail.refGene','Accession')],
                                1, function(x) {grepl(x[2],x[1],fixed=T)})
vexon11splice<-grepl("exon11",varsOI.func$GeneDetail.refGene, fixed=T)
vexon12splice<-grepl("exon12",varsOI.func$GeneDetail.refGene, fixed=T)
asxl1ExceptionSplice<-(varsOI.func$Gene.refGene=="ASXL1")&
    vSplice&vSpliceCorrectTranscript&
    (vexon11splice|vexon12splice)
varsOI.func[asxl1ExceptionSplice,"whitelist"]=T
varsOI.func[asxl1ExceptionSplice,"wl.splice"]=T
varsOI.func[asxl1ExceptionSplice,"wl.exception"]=T
#ASXL2	Frameshift/nonsense/splice-site in exon 11-12
#Lof
asxl2Exception<-(varsOI.func$Gene.refGene=="ASXL2")&(vlof|vFS)&(vexon11|vexon12)
varsOI.func[asxl2Exception,"whitelist"]=T
varsOI.func[asxl2Exception,"wl.lof"]=T
varsOI.func[asxl2Exception,"wl.exception"]=T
#Splice
asxl2ExceptionSplice<-(varsOI.func$Gene.refGene=="ASXL2")&
    vSplice&vSpliceCorrectTranscript&
    (vexon11splice|vexon12splice)
varsOI.func[asxl2ExceptionSplice,"whitelist"]=T
varsOI.func[asxl2ExceptionSplice,"wl.splice"]=T
varsOI.func[asxl2ExceptionSplice,"wl.exception"]=T


#PPM1D	Frameshift/nonsense in exon 5 or 6
vlof<-grepl("X",varsOI.func$NonsynOI, fixed=T) #stop gain or stop loss
vFS<-grepl("fs",varsOI.func$NonsynOI, fixed=T) #frameshift
vexon5<-grepl("exon5",varsOI.func$transcriptOI, fixed=T)
vexon6<-grepl("exon6",varsOI.func$transcriptOI, fixed=T)
ppm1dException<-(varsOI.func$Gene.refGene=="PPM1D")&(vlof|vFS)&(vexon5|vexon6)
varsOI.func[ppm1dException,"whitelist"]=T
varsOI.func[ppm1dException,"wl.lof"]=T
varsOI.func[ppm1dException,"wl.exception"]=T

#TET2	missense mutations in catalytic domains (p.1104-1481 and 1843-2002)
TETidx<-which(varsOI.func$Gene.refGene=="TET2"&
                  varsOI.func$ExonicFunc.refGene=="nonsynonymous SNV"&
                  nchar(varsOI.func$NonsynOI)==6)
for(i in TETidx){
    AApos<-as.numeric(substr(varsOI.func$NonsynOI[i],2,5))
    if((AApos>=1104&AApos<=1481)|(AApos>=1843&AApos<=2002))
    {
        varsOI.func[i,"whitelist"]=T
        varsOI.func[i,"wl.mis"]=T
        varsOI.func[i,"wl.exception"]=T
    }
}

#CBL	RING finger missense p.381-421
CBLidx<-which(varsOI.func$Gene.refGene=="CBL"&
                  varsOI.func$ExonicFunc.refGene=="nonsynonymous SNV"&
                  nchar(varsOI.func$NonsynOI)==5)
for(i in CBLidx){
    AApos<-as.numeric(substr(varsOI.func$NonsynOI[i],2,4))
    if(AApos>=381&AApos<=421)
    {
        varsOI.func[i,"whitelist"]=T
        varsOI.func[i,"wl.mis"]=T
        varsOI.func[i,"wl.exception"]=T
    }
}

#CBLB	RING finger missense p.372-412
CBLBidx<-which(varsOI.func$Gene.refGene=="CBLB"&
                   varsOI.func$ExonicFunc.refGene=="nonsynonymous SNV"&
                   nchar(varsOI.func$NonsynOI)==5)
for(i in CBLBidx){
    AApos<-as.numeric(substr(varsOI.func$NonsynOI[i],2,4))
    if(AApos>=372&AApos<=412)
    {
        varsOI.func[i,"whitelist"]=T
        varsOI.func[i,"wl.mis"]=T
        varsOI.func[i,"wl.exception"]=T
    }
}

#5) flag remaining exceptions for manual review

vlof<-grepl("X|\\*",varsOI.func$NonsynOI, fixed=T) #stop gain or stop loss
vFS<-grepl("fs",varsOI.func$NonsynOI, fixed=T) #frameshift
vSplice<-grepl("splicing",varsOI.func$Func.refGene)

#GATA3	Frameshift/nonsense/splice-site ZNF domain
varsOI.func[(vlof|vFS|vSplice)&(varsOI.func$Gene.refGene=="GATA3"),"manualreview"]=T
#CREBBP	S1680del
varsOI.func[varsOI.func$ExonicFunc.refGene=="nonframeshift deletion" &
                (varsOI.func$Gene.refGene=="CREBBP"),"manualreview"]=T
#CSF3R	truncating c.741-791
varsOI.func[(vlof|vFS|vSplice)&(varsOI.func$Gene.refGene=="CSF3R"),"manualreview"]=T
#DNMT3A	F732del,	F752del
varsOI.func[varsOI.func$ExonicFunc.refGene=="nonframeshift deletion" &
                (varsOI.func$Gene.refGene=="DNMT3A"),"manualreview"]=T
#EP300	VF1148_1149del
varsOI.func[varsOI.func$ExonicFunc.refGene=="nonframeshift deletion" &
                (varsOI.func$Gene.refGene=="EP300"),"manualreview"]=T
#FLT3	FY590-591GD,	del835
varsOI.func[(varsOI.func$ExonicFunc.refGene=="nonframeshift deletion"|
                 varsOI.func$ExonicFunc.refGene=="nonframeshift insertion")
            &(varsOI.func$Gene.refGene=="FLT3"),"manualreview"]=T
#JAK2	del/ins537-539L, del/ins538-539L, del/ins540-543MK, del/ins540-544MK, del542-543, del543-544,	ins11546-547
varsOI.func[(varsOI.func$ExonicFunc.refGene=="nonframeshift deletion"|
                 varsOI.func$ExonicFunc.refGene=="nonframeshift insertion")
            &(varsOI.func$Gene.refGene=="JAK2"),"manualreview"]=T
#KDM6A	del419
varsOI.func[(varsOI.func$ExonicFunc.refGene=="nonframeshift deletion")
            &(varsOI.func$Gene.refGene=="KDM6A"),"manualreview"]=T
#KIT	ins503,	del560,	del579,	del551-559
varsOI.func[(varsOI.func$ExonicFunc.refGene=="nonframeshift deletion"|
                 varsOI.func$ExonicFunc.refGene=="nonframeshift insertion")
            &(varsOI.func$Gene.refGene=="KIT"),"manualreview"]=T
#MPL	del513 W515-518KT
varsOI.func[(varsOI.func$ExonicFunc.refGene=="nonframeshift deletion"|
                 varsOI.func$ExonicFunc.refGene=="nonframeshift insertion")
            &(varsOI.func$Gene.refGene=="MPL"),"manualreview"]=T
#NPM1	Frameshift p.W288fs (insertion at c.859_860, 860_861, 862_863, 863_864)
varsOI.func[(varsOI.func$ExonicFunc.refGene=="frameshift deletion"|
                 varsOI.func$ExonicFunc.refGene=="frameshift insertion")
            &(varsOI.func$Gene.refGene=="NPM1"),"manualreview"]=T

# Removed length of protein code

# Columns to be removed if they exist
if ("aalen" %in% colnames(varsOI.func)) {
    varsOI.func <- varsOI.func[,-c("aalen", "aapos", "aafirst10pctPeptide", "aalast10pctPeptide", "LOFfirst10pct", "LOFlast10pct")]
}

#"aalen", "aapos", "aafirst10pctPeptide", "aalast10pctPeptide", "LOFfirst10pct", "LOFlast10pct"



check_vars <- data.frame(Sample=sample_id,
                         total_num_variants=length(varsOI.func$Sample),
                         total_num_whitelist=length(varsOI.func[whitelist==T,]$Sample),
                         total_num_manualreview=length(varsOI.func[manualreview==T,]$Sample))

#write out files
write.csv(check_vars,paste(sample_id, ".annovar.varsOI.varcount.csv", sep=""), row.names=FALSE)
write.csv(varsOI.func,paste(sample_id, ".annovar.varsOI.allvariants.csv", sep=""), row.names=FALSE)
write.csv(varsOI.func[whitelist==T,],paste(sample_id, ".annovar.varsOI.wl.csv", sep=""), row.names=FALSE)
write.csv(varsOI.func[manualreview==T,],paste(sample_id, ".annovar.varsOI.manualreview.csv", sep=""), row.names=FALSE)
#
# # ========================
# # see multiallelics
# # ========================
# multi_df <- all_ch_canonical
# df <- identify_decomposed_multiallelics(multi_df)
# x <- (df %>% filter(is_multiallelic_site == TRUE))
# unique(x$Gene)
# # "CREBBP" "CUX1"   "DNMT3A" "EP300"  "PRPF8"  "SETD2"  "SUZ12"  "TET2"
#
# identify_decomposed_multiallelics <- function(df) {
#
#     # Create a position key
#     df <- df %>%
#         mutate(pos_key = paste0(Chr, ":", Start))
#
#     # Find positions with multiple alleles in the same sample
#     multiallelic_counts <- df %>%
#         group_by(Sample.ID, pos_key) %>%
#         summarise(n_alleles = n(), .groups = 'drop') %>%
#         filter(n_alleles > 1)
#
#     cat(sprintf("Found %d sample-position pairs with multiple alleles\n",
#                 nrow(multiallelic_counts)))
#
#     # Flag these in the dataframe
#     df <- df %>%
#         left_join(multiallelic_counts, by = c("Sample.ID", "pos_key")) %>%
#         mutate(
#             is_multiallelic_site = !is.na(n_alleles) & n_alleles > 1,
#             multiallelic_group = ifelse(is_multiallelic_site,
#                                         paste0(Sample.ID, "_", pos_key),
#                                         NA)
#         ) %>%
#         select(-n_alleles)
#
#     # Summary
#     cat(sprintf("\nTotal variants: %d\n", nrow(df)))
#     cat(sprintf("Variants at multiallelic sites: %d\n", sum(df$is_multiallelic_site)))
#     cat(sprintf("Unique multiallelic groups: %d\n",
#                 length(unique(df$multiallelic_group[!is.na(df$multiallelic_group)]))))
#
#     return(df)
# }
