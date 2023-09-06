# Chapter 13 
## Pedigree models

# Clean the working environment
rm(list = ls())

# Load packages
library("tidyverse")
library("nadiv")
# install.packages("tidyverse")
# install.packages("nadiv")

# Example 13.1 ------------------------------------------------------------

# Prepare pedigree
a = seq(1, 12)
s = c(NA, NA, NA, NA, 1, 3, 6, NA, 3, 3, 6, 6)
d = c(NA, NA, NA, NA, 2, 4, 5, 5, 8, 8, 8, 8)
ped = data.frame(a, s, d)

# Create A-inverse
Ainv = nadiv::makeAinv(ped)$Ainv

# Create D-inverse
Dinv = nadiv::makeD(ped)$Dinv

# Prepare data
pig = seq(5, 12)
ww = c(17.0, 20.0, 18.0, 13.5, 20.0, 15.0, 25.0, 19.5)
pen = c(1, 1, 1, 1, 2, 2, 2, 2)

data = data.frame(pig, pen, ww)

data$pen = as.factor(pen)
data$pig= factor(x = data$pig, levels = ped$a)

# Variances
varA = 90 
varD = 80
varE = 120

# Variance ratios 
alpha1 = varE / varA
alpha2 = varE / varD

# Setting up the incidence matrices for the MME
# Model: y = Xb + Za + Wd + e

X = model.matrix(ww ~ -1 + pen, data = data)

Z = model.matrix(ww ~ -1 + pig, data = data)

W = Z

XpX = crossprod(X)
ZpZ = crossprod(Z)
WpW = crossprod(W)

Xpy = crossprod(X, data$ww)
Zpy = crossprod(Z, data$ww)
Wpy = crossprod(W, data$ww)

XpZ = crossprod(X, Z)
XpW = crossprod(X, W)

ZpX = crossprod(Z, X)
ZpW = crossprod(Z, W)

WpX = crossprod(W, X)
WpZ = crossprod(W, Z)

LHS = rbind(cbind(XpX, XpZ, XpW), 
            cbind(ZpX, ZpZ + (Ainv * alpha1), ZpW),
            cbind(WpX, WpZ, WpW + (Dinv * alpha2)))

RHS = rbind(Xpy,
            Zpy,
            Wpy)

solutions = solve(LHS, RHS)

round(solutions, 3)

