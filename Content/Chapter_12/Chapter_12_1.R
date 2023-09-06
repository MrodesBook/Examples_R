# Chapter 12
## SS-GBLUP

# Clean the working environment
rm(list = ls())

# Load packages
library("tidyverse")
library("pedigreemm")
# install.packages("tidyverse")
# install.packages("pedigreemm")

# Example 12.1 ------------------------------------------------------------

# Prepare pedigree and data
a = seq(1, 26)
s = c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 13, 15, 15, 14, 14, 14, 1, 14, 14, 14, 14, 14)
d = c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 4, 2, 5, 6, 9, 9, 3, 8, 11, 10, 7, 12)
dyd = c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 9.0, 13.4, 12.7, 15.4, 5.9, 7.7, 10.2, 4.8, NA, NA, NA, NA, NA, NA)  

data = data.frame(a, s, d, dyd)

pedX = pedigree(label = data$a,
                sire  = data$s,
                dam   = data$d)

A = getA(ped = pedX)

Ainv = getAInv(ped = pedX)

# Genotypes
g13 = c(2, 0, 1, 1, 0, 0, 0, 2, 1, 2)
g14 = c(1, 0, 0, 0, 0, 2, 0, 2, 1, 0)
g15 = c(1, 1, 2, 1, 1, 0, 0, 2, 1, 2)
g16 = c(0, 0, 2, 1, 0, 1, 0, 2, 2, 1)
g17 = c(0, 1, 1, 2, 0, 0, 0, 2, 1, 2)
g18 = c(1, 1, 0, 1, 0, 2, 0, 2, 2, 1)
g19 = c(0, 0, 1, 1, 0, 2, 0, 2, 2, 0)
g20 = c(0, 1, 1, 0, 0, 1, 0, 2, 2, 0)
g21 = c(2, 0, 0, 0, 0, 1, 2, 2, 1, 2)
g22 = c(0, 0, 0, 1, 1, 2, 0, 2, 0, 0)
g23 = c(0, 1, 1, 0, 0, 1, 0, 2, 2, 1)
g24 = c(1, 0, 0, 0, 1, 1, 0, 2, 0, 0)
g25 = c(0, 0, 0, 1, 1, 2, 0, 2, 1, 0)
g26 = c(1, 0, 1, 1, 0, 2, 0, 1, 0, 0)

geno = rbind(g13, g14, g15, g16, g17, g18, g19, g20, g21, g22, g23, g24, g25, g26)

rownames(geno) = seq(13, 26)

# Prepare genomic relationship matrix G
M = as.matrix(geno)
p = colMeans(M)/2
q = 1-p

Za = sweep(M,2,2*p,"-")

Gm = tcrossprod(Za) / (2 * (sum(p * q)))
round(Gm, 3)

# A22 and A22 inverse
A22 = A[13:26,13:26]
A22inv = solve(A22)

# Blended G
Gblended = 0.95*Gm + 0.05*A22

# Tuned G
a = ((mean(A22)*mean(diag(Gblended))) - (mean(Gblended)*mean(diag(A22)))) / ((1*mean(diag(Gblended)))-(mean(Gblended)*1))
b = ((1*mean(diag(A22))) - (mean(A22)*1)) / ((1*mean(diag(Gblended)))-(mean(Gblended)*1))

Gtuned1 = a + Gblended*b 

Ginv = solve(Gtuned1)

GiAi = Ginv - A22inv

# H-matrix
Hinv = Ainv + rbind(cbind((rep(0, 12) %*% t(rep(0, 26)))), 
                    cbind((rep(0, 14) %*% t(rep(0, 12))), GiAi))

# Setting up the incidence matrices for the MME
X = model.matrix(dyd ~ 1, data = data)

data$a = factor(x = data$a, levels = data$a)
Z = model.matrix(dyd ~ a - 1, data = data)

alpha = 245/35.241

y = na.omit(data$dyd)

XpX = crossprod(X)
XpZ = crossprod(X, Z)
ZpX = crossprod(Z, X)
ZpZ = crossprod(Z)
Xpy = crossprod(X, y)
Zpy = crossprod(Z, y)

LHS = rbind(cbind(XpX, XpZ), 
            cbind(ZpX, ZpZ + Hinv*alpha))
RHS = rbind(Xpy, 
            Zpy)

solutions = solve(LHS, RHS)

round(solutions, 3)

