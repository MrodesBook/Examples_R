# Chapter 6
## Equal Design Matrices and No Missing Records

# Clean the working environment
rm(list = ls())

# Load packages
library("pedigreemm")
# install.packages("pedigreemm")

# Example 6.1 -------------------------------------------------------------

# Prepare data
wwg = c(4.5, 2.9, 3.9, 3.5, 5.0)
pwg = c(6.8, 5.0, 6.8, 6.0, 7.5)
sex = c("Male", "Female", "Female", "Male", "Male" )
a = seq(4, 8)
s = c(1, 3, 1, 4, 3)
d = c(NA, 2, 2, 5, 6)

data = data.frame(a, s, d, sex, wwg, pwg)
data$sex = factor(data$sex)
data$sex = relevel(factor(sex), ref = "Male")

y2 = c(wwg, pwg)
trait = c("WWG", "WWG", "WWG", "WWG", "WWG", "PWG", "PWG", "PWG", "PWG", "PWG")
sex2 = c(sex, sex)
data2 = data.frame(y2, trait, sex2)
data2$sex2 = factor(data2$sex2)
data2$sex2 = relevel(factor(sex2), ref = "Male")

# Prepare pedigree
a = seq(1, 8)
s = c(NA, NA, NA, 1, 3, 1, 4, 3)
d = c(NA, NA, NA, NA, 2, 2, 5, 6)

ped = data.frame(a, s, d)

pedX = pedigree(label = a,
                sire  = s,
                dam   = d)

Ainv = getAInv(ped = pedX)

# (Co)variances
G_0 = matrix(c(20, 18, 18, 40), nrow = 2)
R_0 = matrix(c(40, 11, 11, 30), nrow = 2)

G_0inv = solve(G_0)
R_0inv = solve(R_0)

# Setting up the incidence matrices for the MME
X = model.matrix(wwg ~ -1 + sex, data = data)

I = diag(1, 2)
IX = kronecker(I, X)

data$a = factor(x = data$a, levels = pedX@label)

Z = model.matrix(wwg ~ a - 1, data = data)

IZ = kronecker(I, Z)

XRX = crossprod(IX, kronecker(R_0inv, diag(1, 5))) %*% IX
XRZ = crossprod(IX, kronecker(R_0inv, diag(1, 5))) %*% IZ

ZRZ = crossprod(IZ, kronecker(R_0inv, diag(1, 5))) %*% IZ
ZRX = crossprod(IZ, kronecker(R_0inv, diag(1, 5))) %*% IX

XRy = crossprod(IX, kronecker(R_0inv, diag(1, 5))) %*% y2
ZRy = crossprod(IZ, kronecker(R_0inv, diag(1, 5))) %*% y2

AiGi = kronecker(G_0inv, Ainv) 

LHS = rbind(cbind(XRX, XRZ), 
            cbind(ZRX, ZRZ + AiGi))

RHS = rbind(XRy, 
            ZRy)

solutions = solve(LHS, RHS)

round(solutions, 3)

