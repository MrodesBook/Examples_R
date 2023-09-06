# Chapter 5 
## Repeatability Model

# Clean the working environment
rm(list = ls())

# Load packages
library("pedigreemm")
# install.packages("pedigreemm")

# Example 5.1 -------------------------------------------------------------

# Prepare pedigree
a = seq(1, 8)
s = c(NA, NA, NA, 1, 3, 1, 3, 1)
d = c(NA, NA, NA, 2, 2, 5, 4, 7)
ped = data.frame(a, s, d)

pedX = pedigree(label = a,
                sire  = s,
                dam   = d)

Ainv = getAInv(ped = pedX)

# Prepare data
cow = rep(seq(4, 8), each = 2)
fy = c(201, 280, 150, 200, 160, 190, 180, 250, 285, 300)
hys = c(1, 3, 1, 4, 2, 3, 1, 3, 2, 4)
pairity = rep(c(1,2), 5)

data = data.frame(cow, pairity, hys, fy)

data$pairity = as.factor(pairity)
data$hys = as.factor(hys)
data$hys = relevel(factor(hys), ref = "1")
data$cow = factor(x = data$cow, levels = pedX@label)

# Variances
varA = 20
varP = 12
varE = 28
varY = 60

# Variance ratios 
alpha1 = varE / varA
alpha2 = varE / varP
rep = (varA + varP) / varY 

# Setting up the incidence matrices for the MME
X = model.matrix(fy ~ -1 + pairity + hys, data = data)
X = X[,-4]

Z = model.matrix(fy ~ cow - 1, data = data)

W = model.matrix(fy ~ droplevels(cow) - 1, data = data)

XpX = crossprod(X)
ZpZ = crossprod(Z)
WpW = crossprod(W)

Xpy = crossprod(X, data$fy)
Zpy = crossprod(Z, data$fy)
Wpy = crossprod(W, data$fy)

XpZ = crossprod(X, Z)
XpW = crossprod(X, W)

ZpX = crossprod(Z, X)
ZpW = crossprod(Z, W)

WpX = crossprod(W, X)
WpZ = crossprod(W, Z)

LHS = rbind(cbind(XpX, XpZ, XpW), 
            cbind(ZpX, ZpZ + (Ainv * alpha1), ZpW),
            cbind(WpX, WpZ, WpW + (diag(nrow = ncol(W)) * alpha2)))

RHS = rbind(Xpy,
            Zpy,
            Wpy)

solutions = solve(LHS, RHS)

round(solutions, 3)

