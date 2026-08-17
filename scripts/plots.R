## Plots for NICE Technical Forum

## Load packages
library(tidyverse)
library(MASS)
library(ggplot2)
library(glm2)
library(glmnet)
library(randomForest)
library(SuperLearner)
library(survey)
library(boot)
library(viridis)
library(viridisLite)

## Source auxiliary scripts
source("./scripts/utils.R")
source("./scripts/dgp.R")
source("./scripts/estimators.R")

## Simulation 1a 
# Plot OLS_lin vs OLS_int vs IPW_glin vs TMLE_Qlin_glin
true_value <- 1

# Formulas
data <- sim1a()
covariates <- colnames(data)[grepl("^W", colnames(data))]
Q_lin <- as.formula(paste("Y ~ A +", paste(covariates, collapse = " + ")))
Q_int <- as.formula(paste("Y ~", "A*(", paste(covariates, collapse = " + "), ")"))
g_lin <- as.formula(paste("A ~", paste(covariates, collapse = " + ")))
g_int <- as.formula(paste("A ~", paste0("(", paste(covariates, collapse = " + "), ")^2")))

# List of estimators
est_sim1a <- list(
  OLS = function(data){ estimate_OLS(data, Qform = Q_lin)},
  OLS_int = function(data){ estimate_OLS(data, Qform = Q_int)},
  IPW = function(data){ estimate_IPW(data, gform = g_lin)},
  TMLE = function(data){ estimate_TMLE(data, Qform = Q_lin, gform = g_lin)}
)

# Estimation & plot
res <- MC_est(k = 1, n_sample = 2000, simulation = sim1a, est_sim1a)
res[,2:8] <- signif(res[,2:8], digits = 2)
res

# Save results
p_title <- "Miss-specified Q, well-specified g \nSample size n = 2000"
p <- plot_est_CI(res, custom_title = p_title)
p
ggsave("./figs/testplot.pdf", p, width=6, height=4, units="in")

