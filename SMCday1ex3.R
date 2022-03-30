# PACKAGES AND DIRECORIES------------------------------------------------------
rm(list=ls())
thiswd <- "/Users/alice/Library/CloudStorage/OneDrive-UniversityofWarwick/COURSESandCONF/SMCmasterclass"
setwd(thiswd)
library(MASS)
library(mvtnorm)

# 0 - read the data in ----------------------------------------------------

data <- read.csv("parus.csv")
