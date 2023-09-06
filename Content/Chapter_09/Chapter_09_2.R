# Chapter 9
## Model with no associative effects

# Clean the working environment
rm(list = ls())

# Load packages
library("pedigreemm")
# install.packages("pedigreemm")

# Example 9.1 Model with no associative effects ---------------------------

# Prepare pedigree
a = seq(1, 15)
s = c(NA, NA, NA, NA, NA, NA, 1, 1, 2, 1, 2, 3, 2, 3, 3)
d = c(NA, NA, NA, NA, NA, NA, 4, 4, 5, 4, 5, 6, 5, 6, 6)

ped = data.frame(a, s, d)

pedX = pedigree(label = a,
                sire  = s,
                dam   = d)

Ainv = getAInv(ped = pedX)

# Prepare data
animal = seq(7, 15)
gr = c(5.50, 9.80, 4.90, 8.23, 7.50, 10.0, 4.50, 8.40, 6.40)
sex = c("Male", "Female", "Female", "Male", "Female", "Female", "Male", "Female", "Male")
pen = c(1, 1, 1, 2, 2, 2, 3, 3, 3)
litter = c(1, 1, 2, 1, 2, 3, 2, 3, 3)

data = data.frame(animal, sex, pen, litter, gr)

data$sex = as.factor(sex)
data$sex = relevel(factor(sex), ref = "Male")
data$pen = as.factor(pen)
data$litter = as.factor(litter)
data$animal = factor(x = data$animal, levels = pedX@label)

# (Co)variances
varGD = 25.70 
varC = 12.5 
varED = 40.6
varES = 10.0

varE = varED + (3-1)*varES

# Variance ratios
alpha1 = varE / varGD
alpha2 = varE / varC


# Setting up the incidence matrices for the MME

X = model.matrix(gr ~ -1 + pen + sex, data = data)

Zd = model.matrix(gr ~ -1 + animal, data = data)

W = model.matrix(gr ~ -1 + litter, data = data)

XpX = crossprod(X)
ZdpZd = crossprod(Zd)
WpW = crossprod(W)

Xpy = crossprod(X, data$gr)
Zdpy = crossprod(Zd, data$gr)
Wpy = crossprod(W, data$gr)

XpZd = crossprod(X, Zd)
XpW = crossprod(X, W)

ZdpX = crossprod(Zd, X)
ZdpW = crossprod(Zd, W)

WpX = crossprod(W, X)
WpZd = crossprod(W, Zd)

LHS = rbind(cbind(XpX, XpZd, XpW), 
            cbind(ZdpX, ZdpZd + (Ainv * alpha1), ZdpW),
            cbind(WpX, WpZd, WpW + diag(1, nrow = nrow(WpW), ncol = ncol(WpW)) * alpha2))

RHS = rbind(Xpy,
            Zdpy,
            Wpy)

solutions = solve(LHS, RHS)

round(solutions, 3)

