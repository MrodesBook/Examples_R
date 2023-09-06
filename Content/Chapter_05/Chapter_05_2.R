# Chapter 5
## Model with Common Environmental Effects

# Clean the working environment
rm(list = ls())

# Load packages
library("pedigreemm")
# install.packages("pedigreemm")

# Example 5.2 -------------------------------------------------------------

# Prepare pedigree
a = seq(1, 15)
s = c(NA, NA, NA, NA, NA, 1, 1, 1, 3, 3, 3, 3, 1, 1, 1)
d = c(NA, NA, NA, NA, NA, 2, 2, 2, 4, 4, 4, 4, 5, 5, 5)
ped = data.frame(a, s, d)

pedX = pedigree(label = a,
                sire  = s,
                dam   = d)

Ainv = getAInv(ped = pedX)

# Prepare data
piglet = seq(6, 15)
ww = c(90, 70, 65, 98, 106, 60, 80, 100, 85, 68)
sex = c("Male", "Female", "Female", "Female", "Male", "Female", "Female", "Male", "Female", "Male")
fsfamily = c(1, 1, 1, 2, 2, 2, 2, 3, 3, 3) 

data = data.frame(piglet, fsfamily, sex, ww)

data$sex = as.factor(sex)
data$sex = relevel(factor(sex), ref = "Male")
data$fsfamily = as.factor(fsfamily)
data$piglet = factor(x = data$piglet, levels = pedX@label)

# Variances 
varP = 100
varA = 20
varC = 15
varE = 65
varY = 100

# Variance ratios 

alpha1 = varE / varA
alpha2 = varE / varC

rep = (varA + varP) / varY 

# Setting up the incidence matrices for the MME
X = model.matrix(ww ~ -1 + sex, data = data)

Z = model.matrix(ww ~ piglet - 1, data = data)

W = model.matrix(ww ~ fsfamily - 1, data = data)

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
            cbind(WpX, WpZ, WpW + (diag(nrow = ncol(W)) * alpha2)))

RHS = rbind(Xpy,
            Zpy,
            Wpy)

solutions = solve(LHS, RHS)

round(solutions, 3)

