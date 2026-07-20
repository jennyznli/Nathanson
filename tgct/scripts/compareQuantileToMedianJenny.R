library(dplyr)
library(ggplot2)
library(readxl)
library(here)

# ------------------- convertBetaToOR ------------------------------ #
convertBetaToOR <- function( B, B.se )
  # input: B (numeric), the beta weight from regression; the logOR of cases v controls
  #        B.se (numeric), standard error of beta
  # output: OR (numeric), the odds ratio
  #         CI.hi/lo (numeric), 95% confidence interval of the odds ratio
{
  OR <- exp(B)
  CI.hi <- exp( B + 1.96 * B.se )
  CI.lo <- exp( B - 1.96 * B.se )
  return( list(OR, CI.lo, CI.hi) )
}
# ------------------------------------------------------------------ #

# ------------------- compareQuantileToMedian --------------------- #
# combined cohort: no ancestry PCs, just case/control status
compareQuantileToMedian <- function( dat, median, quantile )
  # input: dat (data.frame), the df containing the PRS scores, phenotype, quantile, and site
  #        median (integer vector), the range of quantiles to use as the median
  #        quantile (integer vector), the range of quantiles to compare against the median
  # output: a list with the Beta coefficient, which is the logOR of quantile vs median
  #     also returns standard error of beta
{
  dat <- dat[ dat$quantile %in% c(median, quantile),]

  if( dim(dat)[1] == 0 )
  {
    stop("dat is empty after subsetting, check your median/quantile definitions")
  }

  dat$PHENO.q <- 0
  dat$PHENO.q[ dat$quantile %in% quantile ] <- 1

  fit <- glm( data = dat, PHENO ~ PHENO.q, family = "binomial")
  B <- summary(fit)$coefficients[2,1]
  B.se <- summary(fit)$coefficients[2,2]

  return( list(B, B.se) )
}

# Replication only: ancestry PCs (EV1-10)
compareQuantileToMedianCov <- function( dat, median, quantile )
{
  dat <- dat[ dat$quantile %in% c(median, quantile),]

  if( dim(dat)[1] == 0 )
  {
    stop("dat is empty after subsetting, check your median/quantile definitions")
  }

  dat$PHENO.q <- 0
  dat$PHENO.q[ dat$quantile %in% quantile ] <- 1

  fit <- glm( data = dat, PHENO ~ PHENO.q + site + EV1+EV2+EV3+EV4+EV5+EV6+EV7+EV8+EV9+EV10, family = "binomial")
  B <- summary(fit)$coefficients[2,1]
  B.se <- summary(fit)$coefficients[2,2]

  return( list(B, B.se) )
}
# ------------------------------------------------------------------ #

r.df <- data.frame(range.lo = seq(from = 1, to = 96, by = 5),
                   range.hi = seq(from = 5, to = 100, by = 5) )
median.range <- c(45:55)

# -- helper: run the 20-bin quantile-vs-median loop with a given compare fn
runQuantileLoop <- function(dat, compareFn)
{
  out <- c()
  for( i in 1:20 )
  {
    out <- rbind(out, unlist(compareFn( dat, median.range, seq(from = r.df[i,1], to = r.df[i,2], by = 1))))
  }
  OR.out <- convertBetaToOR( out[,1], out[,2] )
  data.frame( OR = OR.out[[1]], CI.lo = OR.out[[2]], CI.hi = OR.out[[3]])
}

plotOR <- function(OR.df)
{
  ggplot( data = OR.df ) +
    geom_segment(aes( x = seq(1:20), xend = seq(1:20), y = CI.lo, yend = CI.hi), color = "blue") +
    geom_point(aes(x = seq(1:20), y = OR), color = "blue") +
    theme_minimal() + xlab("Quantiles") + ylab("Odds Ratio") +
    scale_x_continuous(breaks = seq(1:20), labels = r.df[,2]) +
    geom_hline(yintercept=1, color = "black", linetype = "dashed")
}
# ------------------------------------------------------------------ #

# ------------------- load data ------------------------------------ #
ss <- read_excel(here("tgct", "ss", "20260710_tgct_master.xlsx"))
ss <- ss %>% mutate(Rep_ID = paste0(Final_ID, "_", Final_ID))
ss <- ss %>%
    mutate(Final_ID2 = if_else(GenoSource2 == "Replication", Rep_ID, Final_ID))

# combined case/control list
casectl <- read.table(here("tgct", "data", "tgct_case_control_combined.txt"), header=F,
                      col.names = c("FID", "IID", "STATUS"))

# ------------------- prs scores  ------------------------------------ #
adj <- read.table(here("tgct", "data", "tgct_pooled_adj.profile"), header=T)
org <- read.table(here("tgct", "data", "tgct_pooled.profile"), header=T)
dim(org)

# ------------------- replication  ------------------------------------ #
# replication ancestry PCs
rep_anc <- read.table(here("tgct", "data", "tecac.sample"), header=T)
rep_anc <- rep_anc %>%
    filter(ID_1 != 0) %>%
    mutate(across(starts_with("EV"), as.numeric))

# ------------------- PMBB  ------------------------------------ #
cov <- fread(file.path(pmbb4, "PMBB-Release-2026-4.0_phenotype_covariates.txt"), header = TRUE)
person <- fread(file.path(pmbb4, "PMBB-Release-2026-4.0_phenotype_person.txt"), header = TRUE)

# ------------------- combine  ------------------------------------ #
prs <- org
prs$SCORESUM_ADJ <- adj$SCORESUM
prs <- prs %>% left_join(ss, by = c("IID" = "Final_ID2"))

# ------------------- replication  ------------------------------------ #
# rep_prs <- prs %>% filter(GenoSource2 == "Replication")

rep_anc <- rep_anc %>% left_join(prs, by = c("ID_2" = "IID"))
rep_anc$PHENO <- as.numeric(rep_anc$phenotype)
rep_anc$SCORESUM <- as.numeric(rep_anc$SCORESUM)
rep_anc$SCORESUM_ADJ <- as.numeric(rep_anc$SCORESUM_ADJ)
rep_anc$site <- substr(rep_anc$ID_2,1,2)
rep_anc$site <- as.factor(rep_anc$site)

table(rep_anc$phenotype)
dim(rep_anc)
# 0    1
# 4491 5508

### PLOT UNADJ
rep_anc$quantile <- as.numeric(ntile(rep_anc$SCORESUM, 100))
OR.df.rep <- runQuantileLoop(rep_anc, compareQuantileToMedianCov)
p2 <- plotOR(OR.df.rep)
print(p2)
ggsave(here("tgct", "figures", "tgct_replication_dose_response.png"), plot = p2, width = 7, height = 5)

### PLOT ADJ
rep_anc$quantile <- as.numeric(ntile(rep_anc$SCORESUM_ADJ, 100))
OR.df.rep <- runQuantileLoop(rep_anc, compareQuantileToMedianCov)
p2 <- plotOR(OR.df.rep)
print(p2)
ggsave(here("tgct", "figures", "tgct_replication_dose_response_adj.png"), plot = p2, width = 7, height = 5)

# ------------------- PMBB  ------------------------------------ #
pmbb_prs <- prs %>% filter(GenoSource2 == "PMBB")

# ------------------- pooled  ------------------------------------ #
pooled <- prs %>%
    left_join(casectl %>% select(IID, STATUS), by = "IID")
pooled$PHENO <- pooled$STATUS
pooled$site <- pooled$GenoSource2

# don't have ancestry yet....

### PLOT
pooled$quantile <- as.numeric(ntile(pooled$SCORESUM, 100))
OR.df <- runQuantileLoop(pooled, compareQuantileToMedian)
p1 <- plotOR(OR.df)
print(p1)
ggsave(here("tgct", "figures", "tgct_pooled_dose_response.png"), plot = p1, width = 7, height = 5)


pooled$quantile <- as.numeric(ntile(pooled$SCORESUM_ADJ, 100))
OR.df <- runQuantileLoop(pooled, compareQuantileToMedian)
p1 <- plotOR(OR.df)
print(p1)
ggsave(here("tgct", "figures", "tgct_pooled_dose_response_adj.png"), plot = p1, width = 7, height = 5)


# ------------------------------------------------------------------ #


