# Chapter 4
# Animal Model with Groups

# Clean the working environment
rm(list = ls())

# Load packages
library("tidyverse")
library("nadiv")
# install.packages("tidyverse")
# install.packages("nadiv")

# Example 4.4 -------------------------------------------------------------

# Original pedigree
a = seq(1, 8)
s = c(NA, NA, NA, 1, 3, 1, 4, 3)
d = c(NA, NA, NA, NA, 2, 2, 5, 6)
ped = data.frame(a, s, d)

# Extend the pedigree with unknown parent groups
g1 = max(a) + 1
g2 = g1 + 1

pedEx = ped %>%
  dplyr::mutate(s = ifelse(is.na(s), g1, s)) %>%
  dplyr::mutate(d = ifelse(is.na(d), g2, d)) %>%
  add_row(a = g1, d = NA, s = NA, .before = 1) %>%
  add_row(a = g2, d = NA, s = NA, .before = 2) 

# Create Ainv
Ainv = nadiv::makeAinv(pedEx, ggroups = 2)$Ainv

Ann = Ainv[1:8, 1:8]
Anp = Ainv[1:8, 9:10]
Apn = Ainv[9:10, 1:8]
App = Ainv[9:10, 9:10]

# Prepare data
calves = seq(4, 8)
wwg = c(4.5, 2.9, 3.9, 3.5, 5.0)
sex = c("Male", "Female", "Female", "Male", "Male" )

data = data.frame(calves, sex, wwg)

data$sex = as.factor(sex)
data$sex = relevel(factor(sex), ref = "Male")

data$calves = factor(x = data$calves, levels = ped$a)

# Variance components 
varA = 20
varE = 40

# Variance ratios
alpha = varE/varA

# Model: y = Xb + ZQg + Za + e
# Q = TQ*

X = model.matrix(wwg ~ -1 + sex, data = data)
Z = model.matrix(wwg ~ calves - 1, data = data)

XpX = crossprod(X)
ZpZ = crossprod(Z)
Xpy = crossprod(X, data$wwg)
Zpy = crossprod(Z, data$wwg)
XpZ = crossprod(X, Z)
ZpX = crossprod(Z, X)

LHS = rbind(cbind(XpX, XpZ, 0, 0), 
            cbind(ZpX, ZpZ + (Ann * alpha),  Anp * alpha),
            cbind(0, 0, Apn * alpha, App * alpha))

RHS = rbind(Xpy,
            Zpy,
            0,
            0)

# Constraints
LHS = LHS[,-c(11)]
LHS = LHS[-c(11),]
RHS = RHS[-c(11),]

solutions_QP = solve(LHS, RHS)

rownames(solutions_QP) = c("Males", "Females", 
                           "1", "2", "3", "4", "5", "6", "7", "8", "G10")

round(solutions_QP, 3)

# Without Q-P transformation
Ainv = nadiv::makeAinv(ped)$Ainv

Q = nadiv::ggcontrib(pedEx)

X = model.matrix(wwg ~ -1 + sex, data = data)

data$calves = factor(x = data$calves, levels = rownames(Ainv))
Z = model.matrix(wwg ~ calves - 1, data = data)

XpX = crossprod(X)
XpZ = crossprod(X, Z)
XpZQ = XpZ %*% Q

ZpX = crossprod(Z, X)
ZpZ = crossprod(Z)
ZpZQ = ZpZ %*% Q

QpZpX = crossprod(Q, ZpX)
QpZpZ = crossprod(Q, ZpZ)
QpZpZQ = QpZpZ %*% Q

Xpy = crossprod(X, data$wwg)
Zpy = crossprod(Z, data$wwg)
QpZpy = crossprod(Q, Zpy)

LHS = rbind(cbind(XpX, XpZ, XpZQ), 
            cbind(ZpX, ZpZ + (Ainv * alpha), ZpZQ),
            cbind(QpZpX, QpZpZ, QpZpZQ))

RHS = rbind(Xpy,
            Zpy,
            QpZpy)

# Constraints
LHS = LHS[,-c(11)]
LHS = LHS[-c(11),]
RHS = RHS[-c(11),]

solutions = solve(LHS, RHS)
rownames(solutions) = c("Males", "Females", 
                           "1", "2", "3", "4", "5", "6", "7", "8", "G10")

round(solutions, 3)

a_star = solutions[c(3:10)] + Q %*% c(0, solutions_QP[11])

round(a_star, 3)

# Compare from QP absorbed 
cbind(round(solutions_QP[3:10], 3), round(a_star, 3))

# To calculate from the solutions obtained by animal model without groups
XpX = crossprod(X)
ZpZ = crossprod(Z)
Xpy = crossprod(X, data$wwg)
Zpy = crossprod(Z, data$wwg)
XpZ = crossprod(X, Z)
ZpX = crossprod(Z, X)

LHS = rbind(cbind(XpX, XpZ), 
            cbind(ZpX, ZpZ + (Ainv * alpha)))

RHS = rbind(Xpy,
            Zpy)

solutions_am = solve(LHS, RHS)

round(solutions_am, 3)

a_star_am = solutions_am[c(3:10)] + Q %*% c(0, solutions_QP[11])

round(a_star_am, 3)

# Compare
cbind(round(solutions_QP[3:10], 3), round(a_star, 3), round(a_star_am, 3))

