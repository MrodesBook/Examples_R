# Chapter 14 
## Genomic models

# Clean the working environment
rm(list = ls())

# Load packages
library("pedigreemm")
library("tidyverse")
library("MASS")
# install.packages("pedigreemm")
# install.packages("tidyverse")
# install.packages("MASS")

# Example 14.4 ------------------------------------------------------------

# Prepare data
animal = seq(1, 12)
sire = c(NA, NA, 1, 1, NA, NA, 5, 5, 1, 1, 5, 5)
dam = c(NA, NA, 2, 3, NA, NA, 6, 7, 6, 7, 2, 3)
records = c(11.0, 12.00, 13.00, 14.0, 15.00, 16.00, 17.0, 18.00, 19.00, 20.0, 21.00, 22.0)
herd = c(2, 1, 2, 1, 1, 2, 1, 2, 1, 2, 1, 2)
breed = c(1, 1, 1, 1, 2, 2, 2, 2, 12, 12, 12, 12)

data = data.frame(animal, sire, dam, breed, herd, records)

# Genotypes
g1 = c(2, 0, 1, 1, 0, 0, 2, 2, 1, 2)
g2 = c(1, 1, 0, 0, 1, 2, 0, 1, 1, 0)
g3 = c(2, 1, 1, 0, 0, 1, 1, 1, 1, 1)
g4 = c(2, 0, 2, 1, 0, 0, 2, 1, 2, 1)
g5 = c(1, 1, 2, 1, 1, 0, 2, 2, 1, 2)
g6 = c(0, 0, 1, 1, 0, 1, 0, 1, 2, 1)
g7 = c(1, 1, 2, 0, 1, 0, 1, 1, 2, 2)
g8 = c(2, 0, 2, 1, 0, 0, 2, 1, 1, 2)
g9 = c(1, 0, 1, 1, 0, 0, 1, 2, 2, 1)
g10 = c(2, 1, 2, 1, 0, 0, 2, 2, 2, 2)
g11 = c(2, 2, 1, 0, 1, 1, 1, 1, 2, 1)
g12 = c(2, 1, 2, 1, 0, 1, 1, 2, 1, 2)

geno = rbind(g1, g2, g3, g4, g5, g6, g7, g8, g9, g10, g11, g12)

rownames(geno) = seq(1, 12)
  
# Genotype matrices for purebred 1 and 2
M1 = geno[c(1:4),]
M2 = geno[c(5:8),]

# Matrix of marker alleles for the crossbreds derived from the purebred 1 and 2
Q1 = matrix(c(1, 0, 1, 0, 0, 0, 1, 1, 1, 1,
              1, 0, 1, 1, 0, 0, 1, 1, 1, 1,
              1, 1, 1, 0, 1, 0, 1, 1, 1, 1,   
              1, 1, 1, 1, 0, 0, 1, 1, 1, 1), byrow = TRUE, nrow = 4)

Q2 = matrix(c(0, 0, 0, 1, 0, 0, 0, 1, 1, 0,
              1, 1, 1, 0, 0, 0, 1, 1, 1, 1,
              1, 1, 0, 0, 0, 1, 0, 0, 1, 0,   
              1, 0, 1, 0, 0, 1, 0, 1, 0, 1), byrow = TRUE, nrow = 4)

# SNP frequencies
p1 = colSums(rbind(M1, Q1)) / (2*nrow(M1) + nrow(Q1)) 

p2 = colSums(rbind(M2, Q2)) / (2*nrow(M2) + nrow(Q2)) 

round(p1, 3)
round(p2, 3)

# Prepare genomic matrices  
Z1 = sweep(M1,2,2*p1,"-")
round(Z1, 3)

W1 = sweep(Q1,2,p1,"-")
round(W1, 3)

Z2 = sweep(M2,2,2*p2,"-")

W2 = sweep(Q2,2,p2,"-")

G11 = tcrossprod(Z1) / (2 * (sum(p1 * (1 - p1))))
G13 = tcrossprod(Z1, W1) / (2 * (sum(p1 * (1 - p1))))
G331 = tcrossprod(W1) / (2 * (sum(p1 * (1 - p1))))

G22 = tcrossprod(Z2) / (2 * (sum(p2 * (1 - p2))))  
G23 = tcrossprod(Z2, W2) / (2 * (sum(p2 * (1 - p2))))
G332 = tcrossprod(W2) / (2 * (sum(p2 * (1 - p2))))

G1 = rbind(cbind(G11, G13), 
            cbind(t(G13), G331))

round(G1, 3)

G2 = rbind(cbind(G22, G23), 
           cbind(t(G23), G332))

round(G2, 3)

# Make matrices invertible
fG1 = G1 + (diag(0.01, nrow(G1)))
G1inv = solve(fG1)

fG2 = G2 + (diag(0.01, nrow(G2)))
G2inv = solve(fG2)

# Add zeros for the opposite purebred
G12 = cbind(G1inv[,1:4], matrix(0, 8, 4), G1inv[,5:8])
G10inv = rbind(G12[1:4,], matrix(0, 4, 12), G12[5:8,])

G22 = cbind(matrix(0, 8, 4), G2inv)
G20inv = rbind(matrix(0, 4, 12), G22)

# (Co)variances
r11 = 0.25
r22 = 0.25
r33 = 1/4.75

S1 = matrix(c(1, 0.92, 0.92, 1.5), nrow = 2)
S2 = matrix(c(2, 1.385, 1.385, 1.5), nrow = 2)

S1inv = solve(S1)
S2inv = solve(S2)

# Setting up the incidence matrices for the MME
data$herd = as.factor(herd)
data$animal = factor(animal, levels = c(1:12))
data$sire = factor(sire, levels = unique(data$sire))
data$dam = factor(dam, levels = unique(data$dam))

X = model.matrix(records ~ -1 + herd, data = data)

X1 = rbind(X[1:4, 1:2], matrix(0, 8, 2))
X2 = rbind(matrix(0, 4, 2), X[5:8, 1:2], matrix(0, 4, 2))
X3 = rbind(matrix(0, 8, 2), X[9:12, 1:2]) 

y1 = c(data$records[1:4], rep(0,8))  
y2 = c(rep(0,4), data$records[5:8], rep(0,4))
y3 = c(rep(0,8), data$records[9:12])

Zd0 = matrix(0, nrow(data), nrow(data))

Zd1 = Zd0
diag(Zd1) = c(rep(1,4), rep(0,8))

Zd2 = Zd0
diag(Zd2) = c(rep(0,4), rep(1,4), rep(0,4))

Zd13 = Zd0
Zd13[9,1] = 0.5
Zd13[10,1] = 0.5
Zd13[11,2] = 0.5
Zd13[12,3] = 0.5

Zd23 = Zd0
Zd23[9,6] = 0.5
Zd23[10,7] = 0.5
Zd23[11,5] = 0.5
Zd23[12,5] = 0.5

X1pX1 = (t(X1) * r11) %*% X1
X2pX2 = (t(X2) * r22) %*% X2
X3pX3 = (t(X3) * r33) %*% X3

Zd1pX1 = (t(Zd1) * r11) %*% X1
Zd2pX2 = (t(Zd2) * r22) %*% X2
Zd13pX3 = (t(Zd13) * r33) %*% X3
Zd23pX3 = (t(Zd23) * r33) %*% X3

Zd1pZ1 = (t(Zd1) * r11) %*% Zd1
Zd2pZ2 = (t(Zd2) * r22) %*% Zd2

Zd13pZ13 = (t(Zd13) * r33) %*% Zd13
Zd23pZ23 = (t(Zd23) * r33) %*% Zd23
Zd23pZ13 = (t(Zd23) * r33) %*% Zd13

ZM2 = matrix(0, 2, 2)
ZM12 = matrix(0, 12, 2)
ZM1212 = matrix(0, 12, 12)

LHS = rbind(cbind(X1pX1, ZM2, ZM2, t(Zd1pX1), t(ZM12), t(ZM12), t(ZM12)),
            cbind(ZM2, X2pX2, ZM2, t(ZM12), t(ZM12), t(Zd2pX2), t(ZM12)),      
            cbind(ZM2, ZM2, X3pX3, t(ZM12), t(Zd13pX3), t(ZM12), t(Zd23pX3)),
            cbind(Zd1pX1, ZM12, ZM12, Zd1pZ1 + S1inv[1,1] * G10inv, t(S1inv[1,2] * G10inv), t(ZM1212), t(ZM1212)),
            cbind(ZM12, ZM12, Zd13pX3, S1inv[1,2] * G10inv, Zd13pZ13 + S1inv[2,2] * G10inv, t(ZM1212), t(Zd23pZ13)),
            cbind(ZM12, Zd2pX2, ZM12, ZM1212, ZM1212, Zd2pZ2 + S2inv[1,1] * G20inv, t(S2inv[1,2] * G20inv)),
            cbind(ZM12, ZM12, Zd23pX3, ZM1212, Zd23pZ13, S2inv[1,2] * G20inv, Zd23pZ23 + S2inv[2,2] * G20inv ))                                

# isSymmetric(LHS, check.attributes = FALSE)

RHS = rbind((t(X1) * r11) %*% y1,
            (t(X2) * r22) %*% y2,
            (t(X3) * r33) %*% y3,
            (t(Zd1) * r11) %*% y1,
            (t(Zd13) * r33) %*% y3,
            (t(Zd2) * r22) %*% y2,
            (t(Zd23) * r33) %*% y3)

solutions = MASS::ginv(LHS) %*% RHS

round(solutions, 3)

