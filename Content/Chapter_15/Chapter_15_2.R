# Chapter 15
## Joint Analysis of Quantitative and Binary Traits

# Clean the working environment
rm(list = ls())

# Example 15.2 ------------------------------------------------------------

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
R_0 = matrix(c(20, 2.089, 2.089, 1.036), nrow = 2)

r_12 = R_0[1,2]/sqrt(R_0[1,1]*R_0[2,2])

b = (r_12/sqrt(R_0[1,1])) * (1/(sqrt(1-r_12^2))) 

G = matrix(c(0.7178, 0.1131, 0.1131, 0.0466), nrow = 2)

Gc_0 = matrix(c(1, 0, -b, 1), byrow = TRUE, ncol = 2)

Gc = Gc_0 %*% G %*% t(Gc_0)

Gc_inv = solve(Gc)

# Setting up the incidence matrices for the MME
X1 = model.matrix(bw ~ -1 + origin + season + sex, data = data)
Z1 = model.matrix(bw ~ -1 + sire, data = data)

X2 = X1
Z2 = Z1

W = diag(nrow(data))

IR1 = diag(nrow(data)) * (1/R_0[1,1])

XpX = t(X1) %*% IR1 %*% X1
ZpZ = t(Z1) %*% IR1 %*% Z1

XpZ = t(X1) %*% IR1 %*% Z1
ZpX = t(Z1) %*% IR1 %*% X1

X2pX2 = t(X2) %*% W %*% X2
Z2pZ2 = t(Z2) %*% W %*% Z2

X2pZ2 = t(X2) %*% W %*% Z2
Z2pX2 = t(Z2) %*% W %*% X2

q = data$cd
v = rep(0,6)

Xpy = t(X1) %*% IR1 %*% data$bw
Zpy = t(Z1) %*% IR1 %*% data$bw
X2py2 = t(X2) %*% q 
Z2py2 = t(Z2) %*% q 

LHS = rbind(cbind(XpX, XpZ, matrix(0,4,4), matrix(0,4,6)), 
            cbind(ZpX, ZpZ + (Ainv * Gc_inv[1,1]), matrix(0,6,4), (Ainv * Gc_inv[1,2])),
            cbind(matrix(0,4,4), matrix(0,4,6), X2pX2, X2pZ2),
            cbind(matrix(0,6,4), (Ainv * Gc_inv[2,1]), Z2pX2, Z2pZ2 + (Ainv * Gc_inv[2,2])))

RHS = rbind(Xpy,
            Zpy - (Ainv * Gc_inv[1,2]) %*% v,
            X2py2,
            Z2py2 - (Ainv * Gc_inv[2,2]) %*% v)

solutions = solve(LHS, RHS)

round(solutions, 4)

# Presented solutions correspond to iteration 0. 

