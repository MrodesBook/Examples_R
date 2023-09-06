# Chapter 11
## SNP-BLUP
## GBLUP
## Computing SNP solutions from GBLUP
## Computing base population allele frequencies
 
# Clean the working environment
rm(list = ls())

# Load packages
library("pedigreemm")
# install.packages("pedigreemm")

# Example 11.2 ------------------------------------------------------------

# Prepare data
a = seq(13, 26)
s = c(NA, NA, 13, 15, 15, 14, 14, 14, 1, 14, 14, 14, 14, 14)
d = c(NA, NA, 4, 2, 5, 6, 9, 9, 3, 8, 11, 10, 7, 12)
dyd = c(9.0, 13.4, 12.7, 15.4, 5.9, 7.7, 10.2, 4.8, NA, NA, NA, NA, NA, NA)  
# Animals 13 to 20 are assumed as the reference population 
# Animals 21 to 26 are assumed as the selection candidates
data = data.frame(a, s, d, dyd)

data$a = factor(x = data$a, levels = data$a)

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
  
geno_ref = geno[c(1:8),]

geno_cand = geno[c(9:14),]

pm = colMeans(geno) / 2
round(pm, 3)

# Variance ratios 
varA = 35.241 
varE = 245
k = 2 * (sum(pm * (1 - pm)))
alpha = k*(varE / varA)

Z1 = sweep(geno_ref,2,2*pm,"-")

# Setting up the incidence matrices for the MME
y = na.omit(data$dyd)
X = matrix(rep(1, length(y)))

XpX = crossprod(X)
ZpZ = crossprod(Z1)

Xpy = crossprod(X, y)
Zpy = crossprod(Z1, y)

XpZ = crossprod(X, Z1)
ZpX = crossprod(Z1, X)

LHS = rbind(cbind(XpX, XpZ), 
            cbind(ZpX, ZpZ + (diag(nrow = ncol(Z1)) * alpha)))

RHS = rbind(Xpy,
            Zpy)

sol_SNPBLUP = solve(LHS, RHS)

# SNP effects
round(sol_SNPBLUP, 3)

# Compute DGV for the reference population 
round(Z1%*%sol_SNPBLUP[-1], 3)

# Compute DGV for the selection candidates
Z2 = sweep(geno_cand,2,2*pm,"-")

round(Z2%*%sol_SNPBLUP[-1], 3)

# Re-analyse using EDCs as weights
edc = c(558, 722, 300, 73, 52, 87, 64, 103)

dii = (1/edc)
D = diag(x=dii)

Dinv = diag(x=edc)
  
XpX = t(X)%*%(Dinv)%*%X
ZpZ = t(Z1)%*%(Dinv)%*%Z1

Xpy = t(X)%*%(Dinv)%*%y
Zpy = t(Z1)%*%(Dinv)%*%y

XpZ = t(X)%*%(Dinv)%*%Z1
ZpX = t(Z1)%*%(Dinv)%*%X

LHS = rbind(cbind(XpX, XpZ), 
            cbind(ZpX, ZpZ + (diag(nrow = ncol(Z1)) * alpha)))

RHS = rbind(Xpy,
            Zpy)

sol_SNPBLUPwt = solve(LHS, RHS)

round(sol_SNPBLUPwt, 3)

# Compute DGV for the reference population 
round(Z1%*%sol_SNPBLUPwt[-1], 3)

# Compute DGV for the selection candidates
Z2 = sweep(geno_cand,2,2*pm,"-")

round(Z2%*%sol_SNPBLUPwt[-1], 3)

# Example 11.3 ------------------------------------------------------------

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
X = model.matrix(dyd ~ 1, data = data)

Z = model.matrix(dyd ~ -1 + a, data = data)

# Variance ratios 
alpha = 245/35.241

y = na.omit(data$dyd)

XpX = crossprod(X)
XpZ = crossprod(X, Z)

ZpX = crossprod(Z, X)
ZpZ = crossprod(Z)

Xpy = crossprod(X, y)
Zpy = crossprod(Z, y)

LHS = rbind(cbind(XpX, XpZ), 
            cbind(ZpX, ZpZ + Ginv*alpha))

RHS = rbind(Xpy, 
            Zpy)

solutions_GBLUP = solve(LHS, RHS)

round(solutions_GBLUP, 3)

# Example 11.4 ------------------------------------------------------------

a_hat = solutions_GBLUP[-1]
# To get exactly the same g as in the book, use a_hat from the book (rounded values)
# a_hat = matrix(c(0.069, 0.116, 0.049, 0.260, -0.500, -0.359, 0.146, 
#                 -0.231, 0.028, 0.115, -0.240, 0.143, 0.054, 0.353))

k = 2 * (sum(p * q))

g = 1/k * t(Za) %*% Ginv %*% a_hat

round(g, 3)

# Example 11.6 ------------------------------------------------------------

# Prepare pedigree
a_full = seq(1, 26)
s_full = c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 13, 15, 15, 14, 14, 14, 1, 14, 14, 14, 14, 14)
d_full = c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 4, 2, 5, 6, 9, 9, 3, 8, 11, 10, 7, 12)

pedX = pedigree(label = a_full,
                sire  = s_full,
                dam   = d_full)

Ainv = getAInv(ped = pedX)

# Setting up the incidence matrices for the MME

# Vector of gene content for SNP #1
y = geno[,1]

X = model.matrix(y ~ 1)

aid = factor(x = seq(13, 26), levels = Ainv@Dimnames[[1]])

M = model.matrix(y ~ -1 + aid)

eta = 0.01

XpX = crossprod(X)
XpM = crossprod(X, M)
MpX = crossprod(M, X)
MpM = crossprod(M)
Xpy = crossprod(X, y)
Mpy = crossprod(M, y)

LHS = rbind(cbind(XpX, XpM), 
            cbind(MpX, MpM + Ainv*eta))

RHS = rbind(Xpy, 
            Mpy)

solutions = solve(LHS, RHS)

round(solutions, 3)

round(solutions[2:27] + solutions[1], 3)

round(solutions[2:27] + solutions[1])

