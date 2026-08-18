# ============================================================
# dependencies.R
# ============================================================
# This file is NOT meant to be run as a script.
#
# renv's default ("implicit") snapshot mode only records packages it can
# detect being used, by scanning the repo's R files for literal
# library()/require()/pkg::fun() patterns. Since this repo's actual
# analysis scripts (which use glmnet, glm2, etc.) live in a separate
# project, renv had no evidence these packages were needed and silently
# dropped some of them from renv.lock during the last snapshot.
#
# This file exists purely to give renv's scanner that evidence, so every
# package in the required list gets captured correctly on the next
# renv::snapshot() -- regardless of what analysis code is or isn't
# checked into this repo yet.
# ============================================================

# Data manipulation 
library(tidyverse)
library(MASS)
library(survey)
library(boot)
library(ggrepel)
library(gridExtra)
library(scales)
library(viridis)

# Plotting and visualisation
library(ggplot2)
library(skimr)

# Machine learning / ensemble methods
library(glm2)
library(gam)
library(glmnet)
library(randomForest)
library(SuperLearner)
library(lightgbm)
library(ranger)
library(earth)
library(caret)

# Targeted learning packages
library(tmle)
library(ctmle)

# Miscallaneous
library(future.apply)
library(parallelly)
library(writexl)