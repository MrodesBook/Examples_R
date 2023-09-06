# Chapter 17
## REML

# Clean the working environment
rm(list = ls())

# Load packages
library("pedigreemm")
# install.packages("pedigreemm")

# Example 17.7 ------------------------------------------------------------

# Prepare pedigree
a = seq(1, 8)
s = c(NA, NA, NA, 1, 3, 1, 4, 3)
d = c(NA, NA, NA, NA, 2, 2, 5, 6)
ped = data.frame(a, s, d)

pedX = pedigree(label = a,
                sire  = s,
                dam   = d)

A = getA(ped = pedX)

Ainv = getAInv(ped = pedX)

# Prepare data
calf = seq(4,8)
sex = c("Male", "Female", "Female", "Male", "Male")
wwg = c(2.6, 0.1, 1.0, 3.0, 1.0)

data = data.frame(calf, sex, wwg)

data$sex = as.factor(sex)
data$sex = relevel(factor(sex), ref = "Male")
data$calf = factor(x = data$calf, levels = pedX@label)

# Variance components 
varA = 0.2
varE = 0.4

# Variance ratios
alpha = varE/varA

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
            cbind(ZpX, ZpZ + Ainv*alpha))
RHS = rbind(Xpy, 
            Zpy)

solutions = solve(LHS, RHS)

round(solutions, 3)

res = data$wwg - X%*%solutions[1:2] - Z%*%solutions[3:10]

round(res,4)

CM = solve(LHS)
C22 = CM[3:10, 3:10]
traceAC = sum(diag(Ainv %*% C22))

apAa = t(solutions[3:10]) %*% Ainv %*% solutions[3:10]

V = Z %*% A %*% t(Z) * varA  + diag(5) * varE

Vinv = solve(V)

P = Vinv - Vinv %*% X %*% solve(t(X) %*% Vinv %*% X) %*% t(X) %*% Vinv

ypPy = t(data$wwg) %*% P %*% data$wwg

round(ypPy, 4)

L = 0.5*(-ypPy - log(det(V)) - log(det(t(X) %*% Vinv %*% X)))
round(L, 4)

D1E = 0.5*(t(res)%*%res/varE^2 - (5-2-8)/varE - traceAC/varA)
D1A = 0.5*(apAa/varA^2 - 8/varA + traceAC*varE/varA^2)

A11 = -0.5*((t(res)%*%P%*%res)/varE^2)
A12 = -0.5*(t(solutions[3:10])%*%t(Z)%*%P%*%res) / (varE*varA)
A22 = -0.5*(t(solutions[3:10])%*%t(Z)%*%P%*%Z%*%solutions[3:10]) / varA^2

Ainf = rbind(cbind(-A11, -A12), 
            cbind(-A12, -A22))

Ainf_inv = solve(Ainf)

round(Ainf, 4)
round(Ainf_inv, 4)

varAE_01 = rbind(varE, varA) + Ainf_inv %*% rbind(D1E, D1A)
round(varAE_01, 4)

# Expectation maximization (EM) REML
varA_01 = (apAa + traceAC*varE) / 8
varE_01 = t(res)%*%data$wwg / (5 - 2)

round(varA_01, 4)
round(varE_01, 4)


