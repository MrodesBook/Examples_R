# Chapter 9
## Animal Model with Social Interaction Effects

# Clean the working environment
rm(list = ls())

# Load packages
library("pedigreemm")
# install.packages("pedigreemm")

# Example 9.1 -------------------------------------------------------------

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
varGS = 3.60
covGDS = 2.25
varC = 12.5 
varED = 40.6
varES = 10.0
rho = 0.2

varE = varED + (3-1)*varES

varEstar = varE - (rho * varE) 

G = matrix(c(varGD, covGDS,
             covGDS, varGS),
           nrow = 2, byrow = TRUE)

Ginv = solve(G)

# Variance ratios 
alpha = Ginv * varEstar
round(alpha, 4)

alpha4 = varEstar / (rho * varE) 

alpha5 = varEstar / varC

# Setting up the incidence matrices for the MME
# Model: y = Xb + Z_Du_D +Z_Su_S + Vg + Wc +e

X = model.matrix(gr ~ -1 + sex, data = data)

Zd = model.matrix(gr ~ -1 + animal, data = data)

Zs = model.matrix(animal ~ -1 + pen, data = data)
Zs = Zs %*% t(Zs)
diag(Zs) = 0
Zsp = matrix(0, nrow = nrow(Zd), ncol = 6)
Zs = cbind(Zsp, Zs)
colnames(Zs) = seq(1, 15)

V = model.matrix(gr ~ -1 + pen, data = data)

W = model.matrix(gr ~ -1 + litter, data = data)

XpX = crossprod(X)
ZdpZd = crossprod(Zd)
ZspZs = crossprod(Zs)
VpV = crossprod(V)
WpW = crossprod(W)

Xpy = crossprod(X, data$gr)
Zdpy = crossprod(Zd, data$gr)
Zspy = crossprod(Zs, data$gr)
Vpy = crossprod(V, data$gr)
Wpy = crossprod(W, data$gr)

XpZd = crossprod(X, Zd)
XpZs = crossprod(X, Zs)
XpV = crossprod(X, V)
XpW = crossprod(X, W)

ZdpX = crossprod(Zd, X)
ZdpZs = crossprod(Zd, Zs)
ZdpV = crossprod(Zd, V)
ZdpW = crossprod(Zd, W)

ZspX = crossprod(Zs, X)
ZspZd = crossprod(Zs, Zd)
ZspV = crossprod(Zs, V)
ZspW = crossprod(Zs, W)

VpX = crossprod(V, X)
VpZd = crossprod(V, Zd)
VpZs = crossprod(V, Zs)
VpW = crossprod(V, W)

WpX = crossprod(W, X)
WpZd = crossprod(W, Zd)
WpZs = crossprod(W, Zs)
WpV = crossprod(W, V)

LHS = rbind(cbind(XpX, XpZd, XpZs, XpV, XpW), 
            cbind(ZdpX, ZdpZd + (Ainv * alpha[1,1]), ZdpZs + (Ainv * alpha[1,2]), ZdpV, ZdpW),
            cbind(ZspX, ZspZd + (Ainv * alpha[1,2]), ZspZs + (Ainv * alpha[2,2]), ZspV, ZspW),
            cbind(VpX, VpZd, VpZs, VpV + diag(1, nrow = nrow(VpV), ncol = ncol(VpV)) * alpha4, VpW),
            cbind(WpX, WpZd, WpZs, WpV, WpW + diag(1, nrow = nrow(WpW), ncol = ncol(WpW)) * alpha5))

# isSymmetric(LHS, check.attributes = FALSE)

RHS = rbind(Xpy,
            Zdpy,
            Zspy,
            Vpy,
            Wpy)

solutions = solve(LHS, RHS)

round(solutions, 3)

TBV = solutions[3:17] + (3-1)*solutions[18:32]

round(TBV, 3)

