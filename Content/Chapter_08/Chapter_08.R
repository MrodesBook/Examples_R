# Chapter 8
## Animal Model for a Maternal Trait

# Clean the working environment
rm(list = ls())

# Load packages
library("pedigreemm")
# install.packages("pedigreemm")


# Example 8.1 -------------------------------------------------------------

# Prepare pedigree
a = seq(1, 14)
s = c(NA, NA, NA, NA,  1, 3, 4, 3, 1, 3, 3, 8, 9, 3)
d = c(NA, NA, NA, NA,  2, 2, 6, 5, 6, 2, 7, 7, 2, 6)
ped = data.frame(a, s, d)

pedX = pedigree(label = a,
                sire  = s,
                dam   = d)

Ainv = getAInv(ped = pedX)

# Prepare data
calf = seq(5, 14)
bw = c(35.0, 20.0, 25.0, 40.0, 42.0, 22.0, 35.0, 34.0, 20.0, 40.0)
herds = c(1, 1, 1, 1, 2, 2, 2, 3, 3, 3)
pen = c(1, 2, 2, 1, 1, 2, 2, 2, 1, 2)
dam = c(2, 2, 6, 5, 6, 2, 7, 7, 2, 6)

data = data.frame(calf, herds, pen, bw, dam)

data$herds = as.factor(herds)
data$herds = relevel(factor(herds), ref = "1")
data$pen = as.factor(pen)
data$calf = factor(x = data$calf, levels = pedX@label)
data$dam = factor(x = data$dam, levels = pedX@label)

# (Co)variances
varGA = 150 
varMA = 90
covMA = -40
varP = 40
varE = 350

G = matrix(c(varGA, covMA,
             covMA, varMA),
           nrow = 2, byrow = TRUE)

Ginv = solve(G)

# Variance ratios 
alpha = Ginv * varE

alpha4 = varE / varP

# Setting up the incidence matrices for the MME
X = model.matrix(bw ~ -1 + pen + herds, data = data)

Z = model.matrix(bw ~ -1 + calf, data = data)

W = model.matrix(bw ~ -1 + dam, data = data)

S = model.matrix(bw ~ -1 + droplevels(dam), data = data)

XpX = crossprod(X)
ZpZ = crossprod(Z)
WpW = crossprod(W)
SpS = crossprod(S)

Xpy = crossprod(X, data$bw)
Zpy = crossprod(Z, data$bw)
Wpy = crossprod(W, data$bw)
Spy = crossprod(S, data$bw)

XpZ = crossprod(X, Z)
XpW = crossprod(X, W)
XpS = crossprod(X, S)

ZpX = crossprod(Z, X)
ZpW = crossprod(Z, W)
ZpS = crossprod(Z, S)

WpX = crossprod(W, X)
WpZ = crossprod(W, Z)
WpS = crossprod(W, S)

SpX = crossprod(S, X)
SpZ = crossprod(S, Z)
SpW = crossprod(S, W)

LHS = rbind(cbind(XpX, XpZ, XpW, XpS), 
            cbind(ZpX, ZpZ + (Ainv * alpha[1,1]), ZpW + (Ainv * alpha[1,2]), ZpS),
            cbind(WpX, WpZ + (Ainv * alpha[1,2]), WpW + (Ainv * alpha[2,2]), WpS),
            cbind(SpX, SpZ, SpW, SpS + diag(1, nrow = nrow(SpS), ncol = ncol(SpS)) * alpha4))

RHS = rbind(Xpy,
            Zpy,
            Wpy,
            Spy)

solutions = solve(LHS, RHS)

round(solutions, 3)

