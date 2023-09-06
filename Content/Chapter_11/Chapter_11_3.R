# Chapter 11
## Multivariate genomic models

# Clean the working environment
rm(list = ls())

# Example 11.13 -----------------------------------------------------------

# Genotypes
g4 = c(0, 2, 1, 1, 1, 2, 0, 1, 1, 1)
g5 = c(2, 1, 2, 1, 1, 0, 1, 0, 2, 0)
g6 = c(1, 2, 1, 0, 0, 1, 1, 1, 2, 0)
g7 = c(1, 2, 2, 1, 0, 1, 1, 1, 2, 1)
g8 = c(1, 1, 1,	1, 1, 0, 1, 1, 2, 0)

geno = rbind(g4, g5, g6, g7, g8)

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

# (Co)variances
G_0 = matrix(c(20, 18, 18, 40), nrow = 2)
R_0 = matrix(c(40, 11, 11, 30), nrow = 2)

G_0inv = solve(G_0)
R_0inv = solve(R_0)

# Prepare genomic relationship matrix G
M = as.matrix(geno)
p = colMeans(M)/2
q = 1-p

Za = sweep(M,2,2*p,"-")

G = tcrossprod(Za) / (2 * (sum(p * q)))
round(G, 3)

# Make matrix invertible
fG = G + (diag(0.01, nrow(G)))
Ginv = solve(fG)

# Setting up the incidence matrices for the MME
X = model.matrix(wwg ~ -1 + sex, data = data)

I = diag(1, 2)
IX = kronecker(I, X)

data$a = factor(x = data$a, levels = data$a)

Z = model.matrix(wwg ~ -1 + a, data = data)

IZ = kronecker(I, Z)

XRX = crossprod(IX, kronecker(R_0inv, diag(1, 5))) %*% IX
XRZ = crossprod(IX, kronecker(R_0inv, diag(1, 5))) %*% IZ

ZRZ = crossprod(IZ, kronecker(R_0inv, diag(1, 5))) %*% IZ
ZRX = crossprod(IZ, kronecker(R_0inv, diag(1, 5))) %*% IX

XRy = crossprod(IX, kronecker(R_0inv, diag(1, 5))) %*% y2
ZRy = crossprod(IZ, kronecker(R_0inv, diag(1, 5))) %*% y2

AiGi = kronecker(G_0inv, Ginv) 

LHS = rbind(cbind(XRX, XRZ), 
            cbind(ZRX, ZRZ + AiGi))

RHS = rbind(XRy, 
            ZRy)

solutions = solve(LHS, RHS)

round(solutions, 3)

