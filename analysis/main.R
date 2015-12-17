rm(list=ls())

## chage this to your working directory
setwd('~/Dropbox/org_vs_conv/analysis/writeup/analysis')
library(rjags)
library(R2jags)
library(runjags)
library(ggmcmc)
source('src/prep.R')
load('data/metadat.RData')
## ************************************************************
## model with a study random effect and a study-specific within year
## effect
source('src/make_data_three_level.R')
source('src/models/study_obs_study.R')

JAGS.data <- makeData(meta.dat)
out <- runAnalysis(JAGS.data, scale=1e1)
