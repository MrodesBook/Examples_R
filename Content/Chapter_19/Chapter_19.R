# Chapter 19 
## Gauss-Seidel

# Clean the working environment
rm(list = ls())

# Load packages
library("pedigreemm")
# install.packages("pedigreemm")

# Example 19.2 ------------------------------------------------------------

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
alpha = varE / varA

# Setting up the incidence matrices for the MME
X = model.matrix(wwg ~ -1 + sex, data = data)
Z = model.matrix(wwg ~ calf - 1, data = data)

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

# Initial values of sex and animal effects are set as in the book
solutions = c(4.333, 3.400, 0.000, 0.000, 0.000, 0.167, -0.500, 0.500, -0.833, 0.667)
convergence = 1 

# Set the number of iterations 
# Here we use 20 as in the book
iterations = 20

solutions_all = data.frame(solutions)
solutions_all = rbind(solutions_all, convergence)   

for (i in 1:iterations){
  solutions_old = solutions
  for (j in 1:ncol(LHS)){
    solutions[j]  = (RHS[j] - LHS[j,-j]%*%solutions[-j])/LHS[j,j]
  }
  convergence = (sum((solutions - solutions_old)^2)) / sum(solutions^2)
  solutions_new = data.frame(solutions)
  solutions_new = rbind(solutions_new, convergence)  
  solutions_all = cbind(solutions_all, solutions_new)
}

round(solutions_all, 3)  
  
