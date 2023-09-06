# Chapter 4
## Reduced Animal Model

# Clean the working environment
rm(list = ls())

# Load packages
library("pedigreemm")
library("tidyverse")
# install.packages("pedigreemm")
# install.packages("tidyverse")


# Example 4.3 -------------------------------------------------------------

# Prepare data
wwg = c(4.5, 2.9, 3.9, 3.5, 5.0)
sex = c("Male", "Female", "Female", "Male", "Male")
an = seq(4, 8)
sn = c(1, 3, 1, 4, 3)
dn = c(NA, 2, 2, 5, 6)

data = data.frame(an, sn, dn, sex, wwg)

data$sex = as.factor(sex)
data$sex = relevel(factor(sex), ref = "Male")
data$an = factor(x = data$a, levels = seq(1, 8))

# Prepare pedigree
a = seq(1, 6)
s = c(NA, NA, NA, 1, 3, 1)
d = c(NA, NA, NA, NA, 2, 2)
ped = data.frame(a, s, d)

pedX = pedigree(label = a,
                sire  = s,
                dam   = d)

Ainv = getAInv(ped = pedX)

# Variance components 
varA = 20
varE = 40

# Variance ratios
alpha = varE/varA

# Setting up the incidence matrices for the MME
X = model.matrix(wwg ~ -1 + sex, data = data)

W = matrix(0, nrow = nrow(data), ncol = nrow(Ainv))

rownames(W) = as.character(data$an)
colnames(W) = as.character(pedX@label)

for (i in 1:nrow(data)) {
  a_i = as.character(data[i, "an"])
  if (a_i %in% colnames(W) ) {
    W[i, a_i] = 1
  }
  else {
    s_i = as.character(data[i, "sn"])
    d_i = as.character(data[i, "dn"])
    W[i, s_i] = 0.5
    W[i, d_i] = 0.5
  }
}

R = matrix(0, nrow = nrow(data), ncol = nrow(data))

rownames(R) = as.character(data$an)
colnames(R) = as.character(data$an)

for (i in 1:nrow(data)) {
  a_i = as.character(data[i, "an"])
  if (a_i %in% unique(c(data$sn, data$dn)) ) {
    R[i, a_i] = varE
  }
  else {
    R[i, a_i] = varE + varA * 0.5
  }
}

Ri = solve(R)

XpRiX = crossprod(X, Ri) %*% X
XpRiW = crossprod(X, Ri) %*% W
WpRiX = crossprod(W, Ri) %*% X
WpRiW = crossprod(W, Ri) %*% W

y = na.omit(data$wwg) 
XpRiy = crossprod(X, Ri) %*% y
WpRiy = crossprod(W, Ri) %*% y

LHS = rbind(cbind(XpRiX, XpRiW), 
            cbind(WpRiX, WpRiW + Ainv*(1/varA) ))

RHS = rbind(XpRiy, 
            WpRiy)

solutions = solve(LHS, RHS)

round(solutions, 3)

#Solutions for non-parents

i = diag(2)
d = diag(2, x = 0.5)
k = solve(i+(solve(d))*alpha)

W_n = W[4:5,]
X_n = X[4:5,]
y_n = y[4:5]

b = solutions[1:2,]
a_p = solutions[-(1:2),]

a_n = W_n %*% a_p + k %*% (y_n - X_n %*% b - W_n %*% a_p)

round(a_n, 3)

