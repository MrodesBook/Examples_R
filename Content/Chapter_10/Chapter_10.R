# Chapter 10
## Fixed Regression Model
## Random Regression Model

# Clean the working environment
rm(list = ls())

# Load packages
library("pedigreemm")
library("tidyverse")
library("orthopolynom")
# install.packages("pedigreemm")
# install.packages("tidyverse")
# install.packages("orthopolynom")


# Example 10.1 -------------------------------------------------------------

# Prepare the matrix of Legendre polynomials evaluated at different DIM
# Appendix G

dim = c(4, 38, 72, 106, 140, 174, 208, 242, 276, 310)

dmin = min(dim)
dmax = max(dim)
adim = -1 + 2*(dim - dmin) / (dmax - dmin)

k = 5

M = matrix(0, nrow = length(dim), ncol = k)
for (i in 1:k) {
  M[,i] = adim^(i-1)
}

legendre = legendre.polynomials((k-1), normalized = TRUE)

Lambda = matrix(0, nrow = k, ncol = k)
Lambda[1,1] = unlist(legendre[1])
Lambda[1:2,2] = unlist(legendre[2])
Lambda[1:3,3] = unlist(legendre[3])
Lambda[1:4,4] = unlist(legendre[4])
Lambda[1:5,5] = unlist(legendre[5])

Phi = M %*% Lambda

# Prepare data
htd = seq(1:10)

tdy_4 = c(17.0, 18.6, 24.0, 20.0, 20.0, 15.6, 16.0, 13.0, 8.2, 8.0)
tdy_5 = c(23.0, 21.0, 18.0, 17.0, 16.2, 14.0, 14.2, 13.4, 11.8, 11.4)
tdy_6 = c(NA, NA, NA, NA, NA, 10.4, 12.3, 13.2, 11.6, 8.4)
tdy_7 = c(NA, NA, NA, 22.8, 22.4, 21.4, 18.8, 18.3, 16.2, 15.0)
tdy_8 = c(22.2, 20.0, 21.0, 23.0, 16.8, 11.0, 13.0, 17.0, 13.0, 12.6)

data = data.frame(c(tdy_4, tdy_5, tdy_6, tdy_7, tdy_8), rep(1:10, 5), rep(4:8, each = 10))

colnames(data) = c("tdy", "htd", "animals")

data$htd = as.factor(htd)
data$htd = relevel(factor(htd), ref = "10")
data$animals = factor(x = data$animals, levels = pedX@label)

data = drop_na(data)

colnames(Phi) = c("Leg0", "Leg1", "Leg2", "Leg3", "Leg4")

data = cbind(data, rbind(Phi, Phi, Phi[1:5,], Phi[1:7,], Phi))

# Variances
varA = 5.521
varP = 8.470
varE = 3.710

# Variance ratios 
alpha1 = varE / varA
alpha2 = varE / varP

# Prepare pedigree
a = seq(1, 8)
s = c(NA, NA, NA, 1, 3, 1, 3, 1)
d = c(NA, NA, NA, 2, 2, 5, 4, 7)
ped = data.frame(a, s, d)

pedX = pedigree(label = a,
                sire  = s,
                dam   = d)

Ainv = getAInv(ped = pedX)

# Setting up the incidence matrices for the MME
X1 = model.matrix(tdy ~ -1 + htd, data = data)
X1pX1 = crossprod(X1)

X2 = rbind(Phi, Phi, Phi[1:5,], Phi[1:7,], Phi)
X2pX2 = crossprod(X2)

Q = model.matrix(tdy ~ -1 + animals, data = data)

Z = model.matrix(tdy ~ -1 + droplevels(animals), data = data)

QpQ = crossprod(Q)
ZpZ = crossprod(Z)

y = drop_na(data)

X1py = crossprod(X1, data$tdy)
X2py = crossprod(X2, data$tdy)
Qpy = crossprod(Q, data$tdy)
Zpy = crossprod(Z, data$tdy)

X1pX2 = crossprod(X1, X2)
X1pQ = crossprod(X1, Q)
X1pZ = crossprod(X1, Z)

X2pX1 = crossprod(X2, X1)
X2pQ = crossprod(X2, Q)
X2pZ = crossprod(X2, Z)

QpX1 = crossprod(Q, X1)
QpX2 = crossprod(Q, X2)
QpZ = crossprod(Q, Z)

ZpX1 = crossprod(Z, X1)
ZpX2 = crossprod(Z, X2)
ZpQ = crossprod(Z, Q)

LHS = rbind(cbind(X1pX1, X1pX2, X1pQ, X1pZ), 
            cbind(X2pX1, X2pX2, X2pQ, X2pZ),
            cbind(QpX1, QpX2, QpQ + (Ainv * alpha1), QpZ),
            cbind(ZpX1, ZpX2, ZpQ, ZpZ + diag(1, nrow = nrow(ZpZ), ncol = ncol(ZpZ)) * alpha2))

RHS = rbind(X1py,
            X2py,
            Qpy,
            Zpy)

# Constraints
LHS = LHS[,-c(1)]
LHS = LHS[-c(1),]
RHS = RHS[-c(1),]

solutions = solve(LHS, RHS)

round(solutions, 4)

# Example 10.2 -------------------------------------------------------------

# Covariance matrices for the random regression coefficients for animal effect and pe effects
G = matrix(c(3.297, 0.594, -1.381,
             0.594, 0.921, -0.289,
             -1.381, -0.289, 1.005),
           nrow = 3, byrow = TRUE)

P = matrix(c(6.872, -0.254, -1.101,
             -0.254, 3.171, 0.167,
             -1.101, 0.167, 2.457),
           nrow = 3, byrow = TRUE)

Ginv = solve(G)
Pinv = solve(P)

t(Phi[4,1:3]) %*% G %*% Phi[4,1:3]

t(Phi[4,1:3]) %*% G %*% Phi[5,1:3]

# Setting up the incidence matrices for the MME
X1 = model.matrix(tdy ~ htd, data = data)

X1RX1 = t(X1) %*% diag(1/varE, nrow(data)) %*% X1

X2 = rbind(Phi, Phi, Phi[1:5,], Phi[1:7,], Phi)

X2RX2 = t(X2) %*% diag(1/varE, nrow(data)) %*% X2

Qp = matrix(0, nrow = 42, ncol = 15)
Qp[1:10, 1:3] = Phi[ ,1:3]
Qp[11:20, 4:6] = Phi[ ,1:3]
Qp[21:25, 7:9] = Phi[1:5,1:3]
Qp[26:32, 10:12] = Phi[1:7,1:3]
Qp[33:42, 13:15] = Phi[ ,1:3]

# For animal 6, Q6' is
t(Phi[1:5,1:3])

t(Qp[21:25, 7:9])

QpRQp = t(Qp) %*% diag(1/varE, nrow(data)) %*% Qp

Z = Qp
ZRZ = QpRQp

# Add nr rows and columns for 3 animals without records and 3 regressions 3*3

QRQ = rbind(matrix(0, nrow = 3*3, ncol = 8*3), 
            cbind(matrix(0, nrow = nrow(QpRQp), ncol = 3*3), QpRQp) )

Q = cbind(matrix(0, ncol = 3*3, nrow = nrow(data)), Qp)
     
X1Ry = t(X1) %*% diag(1/varE, nrow(data)) %*% data$tdy
X2Ry = t(X2) %*% diag(1/varE, nrow(data)) %*% data$tdy
QRy = t(Q) %*% diag(1/varE, nrow(data)) %*% data$tdy
ZRy = t(Z) %*% diag(1/varE, nrow(data)) %*% data$tdy

X1RX2 = t(X1) %*% diag(1/varE, nrow(data)) %*% X2
X1RQ = t(X1) %*% diag(1/varE, nrow(data)) %*% Q
X1RZ = t(X1) %*% diag(1/varE, nrow(data)) %*% Z

X2RX1 = t(X2) %*% diag(1/varE, nrow(data)) %*% X1
X2RQ = t(X2) %*% diag(1/varE, nrow(data)) %*% Q
X2RZ = t(X2) %*% diag(1/varE, nrow(data)) %*% Z

QRX1 = t(Q) %*% diag(1/varE, nrow(data)) %*% X1
QRX2 = t(Q) %*% diag(1/varE, nrow(data)) %*% X2
QRZ = t(Q) %*% diag(1/varE, nrow(data)) %*% Z

ZRX1 = t(Z) %*% diag(1/varE, nrow(data)) %*% X1
ZRX2 = t(Z) %*% diag(1/varE, nrow(data)) %*% X2
ZRQ = t(Z) %*% diag(1/varE, nrow(data)) %*% Q

LHS = rbind(cbind(X1RX1, X1RX2, X1RQ, X1RZ), 
            cbind(X2RX1, X2RX2, X2RQ, X2RZ),
            cbind(QRX1, QRX2, QRQ + kronecker(Ainv, Ginv), QRZ),
            cbind(ZRX1, ZRX2, ZRQ, ZRZ + kronecker(diag(1, 5),Pinv)) )

RHS = rbind(X1Ry,
            X2Ry,
            QRy,
            ZRy)

# Constraints
LHS = LHS[,-c(1)]
LHS = LHS[-c(1),]
RHS = RHS[-c(1),]

solutions = solve(LHS, RHS)

round(solutions, 4)

