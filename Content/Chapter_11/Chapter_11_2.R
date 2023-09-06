# Chapter 11
## Haplotype models

# Clean the working environment
rm(list = ls())

# Example 11.8 ------------------------------------------------------------

# Prepare data
wwg = c(4.5, 2.9, 3.9, 3.5, 5.0)
sex = c("Male", "Female", "Female", "Male", "Male" )
a = seq(4, 8)
s = c(1, 3, 1, 4, 3)
d = c(NA, 2, 2, 5, 6)

data = data.frame(a, s, d, sex, wwg)

data$a = factor(x = data$a, levels = data$a)

y = data$wwg

# Genotypes
g1 = c(2, 0, 1, 0, 1, 0, 2, 2, 0)
g2 = c(0, 2, 1, 1, 1, 1, 0, 2, 1)
g3 = c(2, 1, 1, 0, 1, 2, 0, 0, 1)
g4 = c(1, 1, 1, 1, 1, 1, 1, 2, 0)
g5 = c(1, 0, 0, 0, 0, 1, 0, 1, 2) 

geno = rbind(g1, g2, g3, g4, g5)

rownames(geno) = seq(1, 5)

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
X = model.matrix(wwg ~ 1, data = data)

Z = model.matrix(wwg~ -1 + a, data = data)

# Variance ratios 
alpha = 40/20

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

# Pseudo-SNPs
h1 = c(1, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 2, 0, 0, 0, 0)
h2 = c(0, 0, 1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 1, 1, 0, 0)
h3 = c(0, 1, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 1)
h4 = c(0, 1, 1, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0)
h5 = c(0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0) 
      
hap = rbind(h1, h2, h3, h4, h5)

rownames(hap) = seq(1, 5)

# Prepare haplotype-based genomic relationship matrix G
M = as.matrix(hap)
p = colMeans(M)/2
q = 1-p

Za = sweep(M,2,2*p,"-")

Ghap = tcrossprod(Za) / (2 * (sum(p * q)))
round(Ghap, 3)

# Make matrix invertible
fGhap = Ghap + (diag(0.01, nrow(Ghap)))
Ginvhap = solve(fGhap)

# Setting up the incidence matrices for the MME
rm(LHS, RHS)

LHS = rbind(cbind(XpX, XpZ), 
            cbind(ZpX, ZpZ + Ginvhap*alpha))

RHS = rbind(Xpy, 
            Zpy)

solutions_hap = solve(LHS, RHS)

round(solutions_hap, 3)

round(cor(solutions_hap, solutions_GBLUP), 3)

