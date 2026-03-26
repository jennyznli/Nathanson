# ========================
# MODELING
# ========================
library(here)
library(splines)
library(sandwich)
library(lmtest)

setwd("/Users/jennyzli/Documents/Nathanson")
source(here("R", "config.R"))

# ========================
# READ DATA
# ========================
b1_case <- read.csv(file.path("ch", "data", "pmbb_brca1_case.csv"))$x
b2_case <- read.csv(file.path("ch", "data", "pmbb_brca2_case.csv"))$x

m.out4 <- readRDS(file.path("ch", "data", "ch_psm_matched4.rds"))
weights <- m.out4$weights
m.data <- match_data(m.out4)

# ========================
# REGRESSION MODELS
# ========================


fit_glm <- glm(
    CHIP_Binary ~ BRCA12_Case + distance + Batch + Smoke_History + Sequenced_gender +
        PC1 + PC2 + PC3 + PC4 + PC5 + PC6
    data = m.data,
    weights = weights,
    family = quasibinomial()
)

# ========================
# SPLINES
# ========================
model <- glm(
    outcome ~ ns(age, df = 4) + sex + bmi,
    family = binomial,
    data = mydata
)
# test if splines are needed
# fit two models
# if sig-> non linear relationships
anova(model_linear, model_spline, test="LRT")


# ========================
# FIRTH
# ========================
# Firth logistic regression if separation occurs or events are sparse
library(logistf)
fit_firth <- logistf(
    chip_any ~ carrier_any + ns(age, df=4) + smoking + log(mean_depth) + callability_chip +
        PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 +
        cancer_pretherapy + freeze,
    data = matched_df
)
summary(fit_firth)

# ========================
# G-COMPUTATION
# ========================
# Firth logistic regression if separation occurs or events are sparse
library(logistf)
fit_firth <- logistf(
    chip_any ~ carrier_any + ns(age, df=4) + smoking + log(mean_depth) + callability_chip +
        PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 +
        cancer_pretherapy + freeze,
    data = matched_df
)
summary(fit_firth)

marginaleffects::avg_comparisons()

# ========================
# G-COMPUTATION
# ========================
library("marginaleffects")
# g computation? to estimate the ATT and estimate SE and CI


# get marginal log RR - risk ratio/relative risk
avg_comparisons(fit_glm,
                variables = "BRCA12_Case",
                vcov = ~subclass, # cluster robust SE
                newdata = subset(treat == 1),  # version of matched datsets with only treated units
                comparison = "lnratioavg" # compute marginal log RR
)

# get marginal RR - risk ratio/relative risk
avg_comparisons(fit_glm,
                variables = "BRCA12_Case",
                vcov = ~subclass, # cluster robust SE
                newdata = subset(treat == 1),  # version of matched datsets with only treated units
                comparison = "lnratioavg", # compute marginal log RR
                transform = "exp" # exponentiates marginal log RR and CI
)

# get marginal OR - odds ratio
avg_comparisons(fit_glm,
                variables = "BRCA12_Case",
                vcov = ~subclass, # cluster robust SE
                newdata = subset(treat == 1),  # version of matched datsets with only treated units
                comparison = "lnoravg", # compute marginal log RR
                transform = "exp" # exponentiate marginal log RR and CI
)

# get marginal RD - risk difference
avg_comparisons(fit_glm,
                variables = "BRCA12_Case",
                vcov = ~subclass, # cluster robust SE
                newdata = subset(treat == 1)  # version of matched datsets with only treated units
)

# PS subclassification - omit weights, use subclasses directly
# only appropriate when small subclasses relativce to sample size n
avg_comparisons(fit_glm,
                variables = "BRCA12_Case",
                vcov = "HC3", # cluster robust SE
                newdata = subset(treat == 1)) # version of matched datsets with only treated units

fit2 <- glm(
    CHIP_Binary ~ BRCA12_Case * (distance + Batch + Smoke_History + Sequenced_gender +
                                     PC1 + PC2 + PC3 + PC4 + PC5 + PC6)
    data = m.data,
    weights = weights
)

# average estimated potential outcomes
# interperetation of estimates as expectedd potential only valid if
# all covarites are interacted with tratment
avg_predictions(fit_glm,
                variables = "BRCA12_Case",
                vcov = ~subclass, # cluster robust SE
                newdata = subset(treat == 1)) # version of matched datsets with only treated units


library("sandwich") # - cluster robust SE

