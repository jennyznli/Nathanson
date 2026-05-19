# ========================
# Packages
# ========================
packages <- c("tidyr", "dplyr", "plotly", "readr", "readxl", "here",
              "stringr", "ggplot2",  "impute", "pals", "geneplotter", "plotly"
)
purrr::walk(packages, ~ require(.x, character.only = TRUE))
here()

DATE <- format(Sys.Date(), "%Y%m%d")

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

# output <- read.table(here("grafanc", "output"), header = TRUE)
df <- preprocessGrafAnc(here("tecac", "output"))

colors <-  c("Finland" = "#E31A1C",
             "Northern Europe" = "#33A02C",
             "Western Europe" = "#1F78B4",
             "Southern Europe" = "#FF7F00",
             "Northeastern Europe" = "#6A3D9A",
             "Southeastern Europe" = "#B15928",
             "Balkans" = "#A6CEE3",
             "Other Europe" = "#B2DF8A")

df$hover_text <- paste("Sample:", df$Sample, "<br>AncGroupCont:", df$AncGroupCont)

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

htmlwidgets::saveWidget(p, "grafanc.html")


# ============================================================
# LATIN/HISPANIC
# ============================================================

df_latin <- df %>% filter(AncGroupID == "Latin American 1"| AncGroupID == "Latin American 2")
# native american
# replication data set

# p_latin <- plot_ly(
#     data = df_latin,
#     x = ~GD1, y = ~IC1, z = ~GD3,
#     color = ~AncGroupCont,
#     # colors = colors,
#     type = "scatter3d",
#     mode = "markers",
#     text = ~hover_text,
#     hoverinfo = "text",
#     marker = list(size = 4)
# )
#
# p_latin
#
# htmlwidgets::saveWidget(p_latin, "latin.html")

# match and get phenotypes in replication dataset
# 0 control, 1  case, 2 mothers, 3 fathers
ss <- as.data.frame(read.csv(here("tecac", "nci-db-071323.csv"), header = TRUE))
ss <- ss[(ss$ID_GWAS %in% df_latin$Sample), ]
ss <- ss[match(df_latin$Sample, ss$ID_GWAS),]
df_latin$Phenotype <- as.numeric(ss$Phenotype)

write.table(cbind(df_latin_noparents$Sample, df_latin_noparents$Sample, df_latin_noparents$Phenotype),
            file = here("tecac", "20250610_latin_phenotypes.txt"),
            row.names = FALSE,
            col.names = FALSE,
            sep = "\t")

write.csv(df_latin_noparents,
            file = here("tecac", "20250610_latin_grafanc.csv"),
            row.names = FALSE)


