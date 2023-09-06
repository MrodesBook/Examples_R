# Chapter 18
## Gibbs sampling

# Clean the working environment
rm(list = ls())

# Load packages
library("pedigreemm")
# install.packages("pedigreemm")

# Example 18.1 ------------------------------------------------------------

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
calf = seq(4,8)
sex = c("Male", "Female", "Female", "Male", "Male")
wwg = c(4.5, 2.9, 3.9, 3.5, 5.0)

data = data.frame(calf, sex, wwg)

data$calf = factor(x = data$calf, levels = pedX@label)
data$sex = as.factor(sex)
data$sex = relevel(factor(sex), ref = "Male")

# Variance components 
varA = 20
varE = 40

# Variance ratios 
alpha = varE/varA

# Degree of belief 
ve = -2
vu = -2

# Prior
se = 0
su = 0

# Setting up the incidence matrices for the MME
X = model.matrix(wwg ~ -1 + sex, data = data)
Z = model.matrix(wwg ~ -1 + calf, data = data)

XpX = crossprod(X)
XpZ = crossprod(X, Z)
ZpX = crossprod(Z, X)
ZpZ = crossprod(Z)
Xpy = crossprod(X, data$wwg)
Zpy = crossprod(Z, data$wwg)

LHS = rbind(cbind(XpX, XpZ), 
            cbind(ZpX, ZpZ + Ainv * alpha))
RHS = rbind(Xpy, 
            Zpy)

solutions = solve(LHS, RHS)

round(solutions, 3)

# Set number of parameters
obs = length(data$wwg)
lev = ncol(X) + ncol(Z)

# Set seed to get the same random numbers as in the book
set.seed(12345)
# rnorm(lev)

theta = rep(0, lev)
mumu = rep(0, lev)

for(i in 1:lev){
  mu = (RHS[i] - LHS[i,]%*%theta - LHS[i,i]*theta[i])/LHS[i,i]
  mumu[i] = mu
  var = varE/LHS[i,i]
  theta[i] = rnorm(1, mean = mu, sd = sqrt(var)) 
  }

b = theta[1:ncol(X)]
u = theta[(ncol(X) + 1):(lev)]
e = data$wwg - X%*%b - Z%*%u

round(e, 3)

epe = t(e) %*% e
round(epe, 3)

# Draw sigma2e
df = obs + ve
S = epe + se
sigma2e = S / rchisq(1, df=df)

# Draw sigma2u
df = ncol(Z) + vu
S = t(u) %*% Ainv %*% u + su
sigma2u = S / rchisq(1, df=df)

round(sigma2e, 3)
round(sigma2u, 3)

# The next round of iteration is then commenced using the updated values computed for the parameters.

