# Chapter 4
## A Model for an Animal Evaluation (Animal Model)
## Accuracy of evaluations
## A Sire Model

# Clean the working environment
rm(list = ls())

# Load packages
library("pedigreemm")
library("tidyverse")
# install.packages("pedigreemm")
# install.packages("tidyverse")

# Example 4.1 -------------------------------------------------------------

# Prepare pedigree
a = seq(1, 8)
s = c(NA, NA, NA, 1, 3, 1, 4, 3)
d = c(NA, NA, NA, NA, 2, 2, 5, 6)
ped = data.frame(a, s, d)

pedX = pedigree(label = a,
                sire  = s,
                dam   = d)

Ainv = getAInv(ped = pedX)

# Prepare data
wwg = c(NA, NA, NA, 4.5, 2.9, 3.9, 3.5, 5.0)
sex = c("Male", "Female", "Male", "Male", "Female", "Female", "Male", "Male")

data = data.frame(a, s, d, sex, wwg)

data$sex = as.factor(sex)
data$sex = relevel(factor(sex), ref = "Male")
data$a = factor(x = data$a, levels = pedX@label)

# Variance components 
varA = 20
varE = 40

# Variance ratios
alpha = varE/varA

# Setting up the incidence matrices for the MME
X = model.matrix(wwg ~ -1 + sex, data = data)

Z = model.matrix(wwg ~ a - 1, data = data)

y = na.omit(data$wwg) 

XpX = crossprod(X)
XpZ = crossprod(X, Z)
ZpX = crossprod(Z, X)
ZpZ = crossprod(Z)
Xpy = crossprod(X, y)
Zpy = crossprod(Z, y)

LHS = rbind(cbind(XpX, XpZ), 
            cbind(ZpX, ZpZ + Ainv*alpha))

RHS = rbind(Xpy, 
            Zpy)

solutions = solve(LHS, RHS)

round(solutions, 3)

# Accuracy of evaluations

CM = solve(LHS)

diagCM = diag(CM[3:10, 3:10])
r_square = 1 - diagCM*alpha
r = sqrt(r_square)
SEP = sqrt(diagCM*varE)

accuracy = round(data.frame(a, diagCM, r_square, r, SEP), 3)

# Example 4.2 -------------------------------------------------------------

# Prepare pedigree
a = seq(1, 8)
s = c(NA, NA, NA, 1, 3, 1, 4, 3)
d = c(NA, NA, NA, NA, 2, 2, 5, 6)
ped = data.frame(a, s, d)

pedsires = ped %>%
  filter(ped$a %in% ped$s)

pedS = pedigree(label = pedsires$a,
                sire  = pedsires$s,
                dam   = pedsires$d)

Asinv = getAInv(ped = pedS)

# Variance components 
varS = 0.25 * 20
varE = 60 - 5

# Variance ratios
alpha = varE/varS

# Setting up the incidence matrices for the MME
data$s = factor(x = data$s, levels = unique(pedX@sire))
Z = model.matrix( wwg ~ s - 1, data = data)

XpX = crossprod(X)
XpZ = crossprod(X, Z)
ZpX = crossprod(Z, X)
ZpZ = crossprod(Z)
Xpy = crossprod(X, y)
Zpy = crossprod(Z, y)

LHS = rbind(cbind(XpX, XpZ), 
            cbind(ZpX, ZpZ + Asinv*alpha))

RHS = rbind(Xpy, 
            Zpy)

solutions_sire = solve(LHS, RHS)

round(solutions_sire, 3)

