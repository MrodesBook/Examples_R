# Chapter 15
## Joint Analysis of Quantitative and Binary Traits
## Analysis using a linear model

# Clean the working environment
rm(list = ls())

# Example 15.2.2 Linear model ---------------------------------------------

# Prepare data
bw = c(41, 37.5, 41.5, 40, 43, 42, 35, 46, 40.5, 39, 41.4, 43, 34, 47, 42, 44.5, 49, 41.6, 36, 42.7, 32.5, 44.4, 46, 47, 51, 39, 44.5, 40.5, 43.5, 42.5, 48.8, 38.5, 52, 48, 41, 50.5, 43.7, 51, 51.6, 45.3, 36.5, 50.5, 46, 45, 36, 43.5, 36.5)
cd = c(rep(0,11),1,0,1,rep(0,9),1,1, rep(0,5),1, 0,0,0,0,rep(1,5),0,0,1,0,0,0,0)
sex = c(1, 1, rep(2,8), 1, 1, 2, rep(1,6), 2,2,2,1,1,2,2,1,1,2,1,1,1,1,2,2,1,1,1,2,1,2,1,1,1,2,2,2)
origin = c(rep(1,7), rep(2,3), rep(1,5), rep(2,2), 1, rep(2,5), 1, 1, 1, 2, rep(1,7), 2, 2, 2, 2, rep(1,7), 2, 2)
season = c(1,1,1,2,2,2,2,1,1,2,1,1,2,2,2,2,2,1,1,1,2,2,2,2,2,2,1,1,1,2,2,2,2,2,1,1,2,2,1,1,1,2,2,2,2,1,1)
sire = c(rep(1,10), rep(2,7), rep(3,6), rep(4,4), rep(5,11), rep(6,9))

data = data.frame(origin, sire, season, sex, bw, cd)
data$sex = factor(data$sex)
data$sex = relevel(factor(sex), ref = "2")
data$season = factor(data$season)
data$season = relevel(factor(season), ref = "2")
data$origin = factor(data$origin)
data$sire = factor(data$sire)

y2 = c(bw, cd)

# Prepare pedigree
Ainv = matrix(c(1.424,  0.182, -0.667, -0.364,  0.000,  0.000,
                0.182,  1.818,  0.364, -0.727, -0.364, -0.727,
                -0.667,  0.364,  1.788,  0.000, -0.727, -0.364,
                -0.364, -0.727,  0.000,  1.455,  0.000,  0.000,
                0.000, -0.364, -0.727,  0.000,  1.455,  0.000,
                0.000, -0.727, -0.364,  0.000,  0.000,  1.455), 
              byrow = TRUE, 
              ncol = 6)

# (Co)variances
G_0 = matrix(c(0.7178, 0.1131, 0.1131, 0.0466), nrow = 2)
R_0 = matrix(c(20, 2.089, 2.089, 1.036), nrow = 2)

G_0inv = solve(G_0)
R_0inv = solve(R_0)

# Setting up the incidence matrices for the MME
X = model.matrix(bw ~ -1 + origin + season + sex, data = data)

I = diag(1, 2)
IX = kronecker(I, X)

Z = model.matrix(bw ~ -1 + sire, data = data)

IZ = kronecker(I, Z)

XRX = crossprod(IX, kronecker(R_0inv, diag(1, nrow(data)))) %*% IX
XRZ = crossprod(IX, kronecker(R_0inv, diag(1, nrow(data)))) %*% IZ

ZRZ = crossprod(IZ, kronecker(R_0inv, diag(1, nrow(data)))) %*% IZ
ZRX = crossprod(IZ, kronecker(R_0inv, diag(1, nrow(data)))) %*% IX

XRy = crossprod(IX, kronecker(R_0inv, diag(1, nrow(data)))) %*% y2
ZRy = crossprod(IZ, kronecker(R_0inv, diag(1, nrow(data)))) %*% y2

AiGi = kronecker(G_0inv, Ainv) 

LHS = rbind(cbind(XRX, XRZ), 
            cbind(ZRX, ZRZ + AiGi))

RHS = rbind(XRy, 
            ZRy)

solutions = solve(LHS, RHS)

round(solutions, 4)

# Please note ordering of solutions is different compared to the book;
# Fixed effects for both traits first, then the sire effects for both traits.

