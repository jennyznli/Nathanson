library(dplyr)

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
compareQuantileToMedian <- function( dat, median, quantile )
  # input: dat (data.frame), the df containing the PRS scores, covariates, pheno type, and
  #     quantiles
  #        median (integer vector), the range of quantiles to use as the median
  #        quantile (integer vector), the range of quantiles to compare against the median
  # output: a list with the Beta coefficient, which is the logOR of quantile vs median
  #     also returns standard error of beta
{
  # truncate data to just the subjects in the median and quantile groups
  dat <- dat[ dat$quantile %in% c(median, quantile),]
  
  if( dim(dat)[1] == 0 )
  {
    stop("dat is empty after subsetting, check your median/quantile definitions")
  }
  
  # subjects in the median quantiles are coded 0, subjects in the comparison quantile are 1
  dat$PHENO.q <- 0
  dat$PHENO.q[ dat$quantile %in% quantile ] <- 1
  
  # same model as in tecac association
  fit <- glm( data = dat, PHENO ~ PHENO.q  + EV1+EV2+EV3+EV4+EV5+EV6+EV7+EV8+EV9+EV10+site, family = "binomial")
  B <- summary(fit)$coefficients[2,1]
  B.se <- summary(fit)$coefficients[2,2]
  
  return( list(B, B.se) )
}
# ------------------------------------------------------------------ #

dat <- read.table("gwas2021_rep_scores_b38.txt", header=T, sep = ",")


# prs 
#prs <- read.table("~/Documents/NathansonLab/GWAS2026_Manuscript/tecac.prs.all.b38.sscore", header=F)
prs <- read.table("/project/knathans_tecac/TGCT_PRS_jenny/pooled/tecac_prs_pooled_score.profile", header=T)
colnames(prs)[6] <- "SCORE"
cov <- read.table("~/Documents/NathansonLab/GWAS2026_Manuscript/casectl_qc2.eigenvec", header=F)
colnames(cov) <- c("FID", "IID", "EV1", "EV2", "EV3", "EV4", "EV5" ,"EV6", "EV7", "EV8", "EV9", "EV10")
#colnames(prs)[2] <- "IID"
#colnames(prs)[5] <- "SCORE"
prs$IID <- substr( prs$IID, 1 , (nchar(prs$IID) - 1)/2)

cov <- cov[ match(prs$IID, cov$IID),]

prs <- merge(prs, cov, by = "IID")
nci <- read.table("~/Documents/NathansonLab/TGCT_GWAS/nci-db-071323.csv", header=T, sep =",")
nci <- nci[ nci$ID_GWAS %in% prs$IID,]
nci <- nci[ match(prs$IID, nci$ID_GWAS),]
prs$PHENO <- nci$Phenotype
prs$site <- substr(prs$IID,1,2)

prs$quantile <- ntile(prs$SCORE,100)
median.range <- c(45:55)

r.df <- data.frame(range.lo = seq(from = 1, to = 96, by = 5),
                   range.hi = seq(from = 5, to = 100, by = 5) )

# -- compare each quantile to the median and convert to OR
out <- c()
for( i in 1:20 )
{
  out <- rbind(out, unlist(compareQuantileToMedian( prs, median.range, seq(from = r.df[i,1], to = r.df[i,2], by = 1))))
}

OR.out <- convertBetaToOR( out[,1], out[,2] )
OR.df <- data.frame( OR = OR.out[[1]], CI.lo = OR.out[[2]], CI.hi = OR.out[[3]])
# --


# plot
p1 <- ggplot( data = OR.df ) + 
  geom_segment(aes( x = seq(1:20), xend = seq(1:20), y = CI.lo, yend = CI.hi), color = "blue") +
  geom_point(aes(x = seq(1:20), y = OR), color = "blue") + 
  theme_minimal() + xlab("Quantiles") + ylab("Odds Ratio") +
  scale_x_continuous(breaks = seq(1:20), labels = r.df[,2]) +
  geom_hline(yintercept=1, color = "black", linetype = "dashed")
print(p1)


