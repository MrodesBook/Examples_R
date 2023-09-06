# Chapter 13 
## Genomic models

# Clean the working environment
rm(list = ls())

# Load packages
library("tidyverse")
# install.packages("tidyverse")

# Example 13.3 ------------------------------------------------------------

# Read pedigree and SNP genotypes
ped = read_table2(file = "ped_snp.txt", col_names = FALSE)

colnames(ped) = c("Animal", "Sire", "Dam", "Sex",
                   paste("Snp_", 1:(length(ped)-4), sep = ""))

geno = ped[,-c(1:4)]

# Prepare data
pig = seq(5, 15)
ww = c(17.0, 20.0, 18.0, 13.5, 20.0, 15.0, 25.0, 19.5, 22.5, 16.0, 24.5)
pen = c(1, 1, 1, 1, 2, 2, 2, 2, 1, 1, 2)

data = data.frame(pig, pen, ww)

data$pen = as.factor(pen)
data$pig = factor(x = data$pig, levels = ped$Animal)

# Prepare genomic relationship matrices G and D
M = as.matrix(geno)
p = colMeans(M)/2
q = 1-p

Ma = sweep(M,2,2*p,"-")

Md2 = matrix(0,nrow=nrow(M),ncol=ncol(M))
tmp = numeric(nrow(M))
for(i in 1:ncol(Md2)){
  tmp[M[,i]==0] = -2*p[i]^2 
  tmp[M[,i]==1] = 2*p[i]*q[i] 
  tmp[M[,i]==2] = -2*q[i]^2 
  Md2[,i] = tmp
}

GG = tcrossprod(Ma) / (2 * (sum(p * q)))
round(GG, 3)

DD = tcrossprod(Md2) / (sum((2 * p * q)^2))  
round(DD, 3)

# Make matrices invertible
fG = GG + (diag(0.01, nrow(GG)))
Ginv = solve(fG)

fD = DD + (diag(0.01, nrow(DD)))
Dinv = solve(fD)

# Variances
varA = 90 
varD = 80
varE = 120

# Variance ratios 
alpha1 = varE / varA
alpha2 = varE / varD

# Setting up the incidence matrices for the MME
# Model: y = Xb + Za + Wd + e

X = model.matrix(ww ~ -1 + pen, data = data)

Z = model.matrix(ww ~ -1 + pig , data = data)

W = Z

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
            cbind(ZpX, ZpZ + (Ginv * alpha1), ZpW),
            cbind(WpX, WpZ, WpW + (Dinv * alpha2)))

RHS = rbind(Xpy,
            Zpy,
            Wpy)

solutions = solve(LHS, RHS)

round(solutions, 3)

# Example 13.5 ------------------------------------------------------------

# Assume genetic variance for epistatic effect 
varAA = 6

# Variance ratios
alpha3 = varE / varAA

# Prepare epistatic additive by additive genomic relationship matrix
HadamardGG = GG*GG 
k = sum(diag(HadamardGG))/15
GAA = HadamardGG / k
round(GAA, 3)

# Make matrix invertible
fGAA = GAA + (diag(0.01, nrow(GAA)))
GAAinv = solve(fGAA)

# Setting up the incidence matrices for the MME
# Model: y = Xb + Za + Wd + Saa + e

S = Z

XpS = crossprod(X, S)
ZpS = crossprod(Z, S)
WpS = crossprod(W, S)

SpX = crossprod(S, X)
SpZ = crossprod(S, Z)
SpW = crossprod(S, W)
SpS = crossprod(S)

Spy = crossprod(S, data$ww)

rm(LHS, RHS)

LHS = rbind(cbind(XpX, XpZ, XpW, XpS), 
            cbind(ZpX, ZpZ + (Ginv * alpha1), ZpW, ZpS),
            cbind(WpX, WpZ, WpW + (Dinv * alpha2), WpS),
            cbind(SpX, SpZ, SpW, SpS + (GAAinv * alpha3)))

RHS = rbind(Xpy,
            Zpy,
            Wpy,
            Spy)

solutions_epi = solve(LHS, RHS)

round(solutions_epi, 3)

