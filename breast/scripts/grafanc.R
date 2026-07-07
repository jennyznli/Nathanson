# ========================
# Packages
# ========================
packages <- c("tidyr", "dplyr", "plotly", "readr", "readxl", "here",
              "stringr", "ggplot2",  "impute", "pals", "geneplotter", "plotly"
)
purrr::walk(packages, ~ require(.x, character.only = TRUE))
here()

DATE <- format(Sys.Date(), "%Y%m%d")

icd10_9 <- read.csv("icd10cmtoicd9gem.csv")

# ============================================================
# FUNCTIONS
# ============================================================

preprocessGrafAnc <- function(INFILE)
{
    dat <- read.table(INFILE, header=T)
    dat$AncGroupCont <- as.factor(substr(dat$AncGroupID, 1, 1))
    dat$AncGroupID <- as.factor(dat$AncGroupID)
    dat$AncGroupCont <- recode(dat$AncGroupCont,
                               "1" = "AFR",
                               "2" = "MEN",
                               "3" ="EUR",
                               "4" = "SAS",
                               "5" = "EAS",
                               "6" = "AMR",
                               "7" = "OCN",
                               "8" = "MIX")
    dat$AncGroupID <- recode(dat$AncGroupID,
                             "101" = "Nigeria",
                             "102" = "Western Africa",
                             "103" = "Central Africa",
                             "104" = "Kenya",
                             "105" = "Southern Africa",
                             "106" = "Northeastern Africa",
                             "107" = "African American",
                             "108" = "Other Africa",
                             "201" = "Northern Africa",
                             "202" = "Middle East 1",
                             "203" = "Middle East 2",
                             "301" = "Finland",
                             "302" = "Northern Europe",
                             "303" = "Western Europe",
                             "304" = "Southern Europe",
                             "305" = "Northeastern Europe",
                             "306" = "Southeastern Europe",
                             "307" = "Balkans",
                             "308" = "Other Europe",
                             "401" = "Asian India",
                             "402" = "Gujarati in India",
                             "403" = "Northern South Asia",
                             "404" = "Sri Lanka",
                             "405" = "Bangladesh",
                             "501" = "Japan Ryukyu",
                             "502" = "Japan Main Islands",
                             "503" = "Korea",
                             "504" = "Northern Asia",
                             "505" = "Northern China 1",
                             "506" = "Northern China 2",
                             "507" = "Southern China 1",
                             "508" = "Southern China 2",
                             "509" = "Southeast Asia",
                             "510" = "Thailand",
                             "511" = "Other East Asia",
                             "601" = "Latin American 1",
                             "602" = "Latin American 2",
                             "603" = "Native American",
                             "700" = "Oceania",
                             "800" = "Multiracial")
    return(dat)
}

# ============================================================
# GRAFANC 3D PLOT
# ============================================================

## GD1, GD2, GD3
all <- preprocessGrafAnc(here("all_pmbb_ancestry_pops.txt"))
v3 <- preprocessGrafAnc(here("pmbb_controls_v3.grafanc"))
v3_men <- preprocessGrafAnc(here("pmbb_ctrls_men60_v3.grafanc"))
v3_men <- as.table(here("pmbb_ctrls_men60_v3.grafanc"))
v3_men <- read.delim("pmbb_ctrls_men60_v3.grafanc", header = TRUE)

v2 <- preprocessGrafAnc(here("pmbb_all_v2.grafanc"))

# pmbb_all_v2.grafanc
# v2 <- preprocessGrafAnc(here("pmbb_all_v2.grafanc"))
# v2 <- as.table(here("pmbb_all_v2.grafanc"))
# df$hover_text <- paste("Sample:", df$Sample, "<br>AncGroupCont:", df$AncGroupCont)

v3_a <- v3[,c("Sample", "AncGroupCont")]
all_a <- all[,c("Sample", "AncGroupCont")]

all2_a <- all_a[all_a$Sample %in% v3_a$Sample, ]
comb_a <- all2_a[match(v3_a$Sample, all2_a$Sample),]
comb_a$AncGroupCont_John <- v3_a$AncGroupCont
mismatch <- comb_a[(comb_a$AncGroupCont != comb_a$AncGroupCont_John),]
# 113 mismatches

p <- plot_ly(
    data = df,
    x = ~GD1, y = ~GD2, z = ~GD3,
    color = ~AncGroupCont,
    # colors = colors,
    type = "scatter3d",
    mode = "markers",
    text = ~hover_text,
    hoverinfo = "text",
    marker = list(size = 4)
)

p

htmlwidgets::saveWidget(p, here("breast", paste0(DATE, "_grafanc_pmbb_3d_gd1_gd2_gd3.html")))

## EU1, IC2, IC3
p <- plot_ly(
    data = df,
    x = ~EU1, y = ~IC2, z = ~IC3,
    color = ~AncGroupCont,
    # colors = colors,
    type = "scatter3d",
    mode = "markers",
    text = ~hover_text,
    hoverinfo = "text",
    marker = list(size = 4)
)

p
htmlwidgets::saveWidget(p, here("breast", paste0(DATE, "_grafanc_pmbb_3d_eu1_ic2_ic3.html")))

# ============================================================
# SUBSET BY ANCESTRY
# ============================================================
# Create directories if they don't exist
dir.create("breast/ancestry_female", showWarnings = FALSE, recursive = TRUE)

# > table(df$AncGroupCont)
#
# AFR   MEN   EUR   SAS   EAS   AMR   MIX
# 7419   107 19178   301   402  1319   204

# Filter and export each ancestry group as a .txt file with two columns of the same ID
eur <- df %>% filter(AncGroupCont == "EUR")
write.table(cbind(0, eur$Sample), file = "breast/ancestry_female/eur.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

afr <- df %>% filter(AncGroupCont == "AFR")
write.table(cbind(0, afr$Sample), file = "breast/ancestry_female/afr.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

amr <- df %>% filter(AncGroupCont == "AMR")
write.table(cbind(0, amr$Sample), file = "breast/ancestry_female/amr.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

eas <- df %>% filter(AncGroupCont == "EAS")
write.table(cbind(0, eas$Sample), file = "breast/ancestry_female/eas.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

sas <- df %>% filter(AncGroupCont == "SAS")
write.table(cbind(0, sas$Sample), file = "breast/ancestry_female/sas.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

men <- df %>% filter(AncGroupCont == "MEN")
write.table(cbind(0, men$Sample), file = "breast/ancestry_female/men.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
#
# ocn <- df %>% filter(AncGroupCont == "OCN")
# write.table(cbind(0, ocn$Sample), file = "breast/ancestry_female/ocn.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

mix <- df %>% filter(AncGroupCont == "MIX")
write.table(cbind(0, mix$Sample), file = "breast/ancestry_female/mix.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

# then copy over...






# Create directories if they don't exist
dir.create("breast", showWarnings = FALSE)
dir.create("breast/ancestry", showWarnings = FALSE, recursive = TRUE)

# Function to write ancestry files with FID=0 and IID=Sample
write_ancestry_file <- function(df_subset, filename) {
    write.table(
        data.frame(FID = 0, IID = df_subset$Sample),
        file = filename,
        row.names = FALSE,
        col.names = FALSE,
        quote = FALSE
    )
}

# Write files for each ancestry group
write_ancestry_file(df %>% filter(AncGroupCont == "EUR"), "breast/ancestry/eur.txt")
write_ancestry_file(df %>% filter(AncGroupCont == "AFR"), "breast/ancestry/afr.txt")
write_ancestry_file(df %>% filter(AncGroupCont == "AMR"), "breast/ancestry/amr.txt")
write_ancestry_file(df %>% filter(AncGroupCont == "EAS"), "breast/ancestry/eas.txt")
write_ancestry_file(df %>% filter(AncGroupCont == "SAS"), "breast/ancestry/sas.txt")
write_ancestry_file(df %>% filter(AncGroupCont == "MEN"), "breast/ancestry/men.txt")
write_ancestry_file(df %>% filter(AncGroupCont == "MIX"), "breast/ancestry/mix.txt")




