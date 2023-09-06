# Chapter 12
## APY

# Clean the working environment
rm(list = ls())

# Load packages
library("tidyverse")
# install.packages("tidyverse")

# Source function
source("APY_fja.R")

# Example 12.2 ------------------------------------------------------------

# Prepare data
a = seq(13, 26)
s = c(NA, NA, 13, 15, 15, 14, 14, 14, 1, 14, 14, 14, 14, 14)
d = c(NA, NA, 4, 2, 5, 6, 9, 9, 3, 8, 11, 10, 7, 12)
dyd = c(9.0, 13.4, 12.7, 15.4, 5.9, 7.7, 10.2, 4.8, NA, NA, NA, NA, NA, NA)  

data = data.frame(a, s, d, dyd)

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

G = tcrossprod(Za) / (2 * (sum(p * q)))
round(G, 3)

# Make matrix invertible
fG = G + (diag(0.01, nrow(G)))
Ginv = solve(fG)

geno_ref = tibble(geno_id = rownames(geno), geno_order = seq(1:length(geno_id)) )

# Set number of core animals 
ncore = 4

# Set seed to get the same random core animals as in the book
set.seed(12345)

index_core = sample_n(geno_ref[2], size = ncore)

geno_ref2 = geno_ref %>%
  dplyr::mutate(core_status = if_else(geno_order %in% index_core$geno_order, true = 1, false = 0))

# Run APY inverse function
apyginvlist = APY_inverse(G, geno_ref2)

round(apyginvlist$APY_Ginv, 3)

APYinv = apyginvlist$APY_Ginv

# Setting up the incidence matrices for the MME
apyginvlist$index_file$geno_id = as.integer(apyginvlist$index_file$geno_id)

pheno_apy = inner_join(data, apyginvlist$index_file[, c(1,4)], by = c("a" = "geno_id")) %>%
  arrange(apy_order)

X = model.matrix(dyd ~ 1, data = pheno_apy)

pheno_apy$a = factor(x = pheno_apy$a, levels = pheno_apy$a)
Z = model.matrix(dyd ~ -1 + a, data = pheno_apy)

# Variance ratios 
alpha = 245/35.241

y = na.omit(pheno_apy$dyd)

XpX = crossprod(X)
XpZ = crossprod(X, Z)

ZpX = crossprod(Z, X)
ZpZ = crossprod(Z)

Xpy = crossprod(X, y)
Zpy = crossprod(Z, y)

LHS = rbind(cbind(XpX, XpZ), 
            cbind(ZpX, ZpZ + APYinv*alpha))

RHS = rbind(Xpy, 
            Zpy)

solutions_APY = solve(LHS, RHS)

# Reorder solutions 
solutions_APY_sorted = tibble(solutions_APY) %>%
  mutate(aid = rownames(solutions_APY)) %>%
  arrange(aid)

round(solutions_APY_sorted[-1,1], 3)

# Solve GBLUP

# Setting up the incidence matrices for the MME
X = model.matrix(dyd ~ 1, data = data)

data$a = factor(x = data$a, levels = data$a)
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

# Correlations between DGV form direct G inverse and APY
round(cor(solutions_APY_sorted[-1,1], solutions_GBLUP[-1,]), 3)

