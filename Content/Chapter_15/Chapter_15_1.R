# Chapter 15 
## Threshold model

# Clean the working environment
rm(list = ls())

# Load packages
library("tidyverse")
library("pedigreemm")
# install.packages("tidyverse")
# install.packages("pedigreemm")

# Source function
source("threshold_fja.R")

# Example 15.1 ------------------------------------------------------------

# Prepare pedigree
a = seq(1, 4)
s = c(NA, NA, 1, 3)
d = c(NA, NA, NA, NA)
ped = data.frame(a, s, d)

pedX = pedigree(label = a,
                sire  = s,
                dam   = d)

Ainv = getAInv(ped = pedX)

# Prepare data
herd = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2)
sex = c("Male", "Female", "Male", "Female", "Male", "Female",  "Male", 
        "Female", "Male", "Female", "Male", "Male", "Female", "Male", 
        "Female", "Male", "Male", "Female", "Female", "Male")
sire = c(1, 1, 1, 2, 2, 2, 3, 3, 3, 1, 1, 1, 2, 2, 3, 3, 4, 4, 4, 4)
c1 = c(1, 1, 1, 0, 1, 3, 1, 0, 1, 2, 1, 0, 1, 1, 0, 0, 0, 1, 2, 2)
c2 = c(0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0)
c3 = c(0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 0)

data = data.frame(herd, sex, sire, c1, c2, c3)
data$total = data$c1 + data$c2 + data$c3

# Variance components 
varS = 1/19
varE = 40

# Proportion before t1
sum(data$c1)/sum(data$total)
# Proportion before t2
(sum(data$c1) + sum(data$c2))/sum(data$total)

t1_init = 0.468
t2_init = 1.080

inits = matrix(c(t1_init, t2_init, rep(0, 6)))

# Setting up the incidence matrices for the MME
data$sex = as.factor(sex)
data$sex = relevel(factor(sex), ref = "Male")
data$herd = as.factor(herd)

X = model.matrix(total ~ -1 + herd + sex, data = data)
X = X[,-1]

# Constraints: Herd "1" and Sex "Male" are set to zero

data$sire = factor(x = data$sire, levels = pedX@label)
Z = model.matrix(total ~ -1 + sire, data = data)

# Solve using function "thr"
solutions = thr(Ainv, X, Z, ncats = 3, 
                cat = as.matrix(cbind(data$c1, data$c2, data$c3)), 
                inits, 
                ratio = 1/varS, 
                disp = TRUE)

