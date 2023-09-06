## The non-linear (threshold) model for categorical traits.

## Arguments
## Ainv: inverse of the additive relationship matrix
## X: incidence matrix relating data to fixed effects
## Z: incidence matrix relating data to individuals
## ncats: number of categories in the data
## cat: matrix of category of response
## inits: vecotor of initial values for the threshold effects, fixed effects and random effects
## ratio: ratio of the variance components
## disp: a logical value. If true, show the solutions at each iteration


## Literature1 : Mrode, R.A. 2005. Linear Models for the Prediction of Animal Breeding Values. CAB International, Oxon, UK.
## Literature 2: Gianola, D. and Foulley, J.L. 1983. Sire Evaluation for Ordered Categorical Data with a Threshold Model. Genetic Selection Evolution.

## Author: Gota Morota <morota at wisc dot edu>
## Create: 16-Apr-2010
## Last-Modified: 18-Apr-2010
## License: GPLv3 or later

`thr` <-
function(Ainv, X, Z, ncats, cat, inits, ratio, disp){

	# starting values
	B <- inits
	# total
	total <- apply(cat, 1, sum)
	# number of records
	n <- length(total)
	# number of category
	m <- ncats
	# number of threshold
	t <- m-1

	diff <- 1
	it = 0
	while (diff > 10E-9){
		it = it + 1
		# nannja kore
		tmp <- cbind(X, Z)%*%B[(t+1):dim(B)[1]]
		d <- matrix(0, nrow=n, ncol=t)
		for (i in 1:n){
			for (j in 1:t){
				d[i, j] <- B[j] - tmp[i]
			}
		}

		phi <- dnorm(d)
		Phi <- pnorm(d)

		P <- matrix(0, nrow=n, ncol=m)
		for (i in 1:n){
			for (j in 2:m){
				if (j == m){
					P[i,j] <- 1.0 - Phi[i, j-1]
				}
				else {
					P[i,j] <- Phi[i,j] - Phi[i,j-1]
				}
			}
		}

		P[,1] <- Phi[,1]


		#  W matrix
		w <- as.numeric()
		for (i in 1:n){
			# first element
			tmp1 <- ((0-phi[i,1])^2/(Phi[i,1])  )
			# last element
			tmp2 <- ((phi[i,t] - 0)^2/(P[i,m])  )
			tmp3 <- 0
			for (j in 2:t){
				for (k in 2:(m-1)){
					#print(k)
					tmp3 <- tmp3 + (   (phi[i,j-1]-phi[i,j])^2/(P[i,k])       )
				}
			}
	
			w[i] <- total[i]*(tmp1 + tmp2 + tmp3)
		}

		W <- diag(w)




		# v vector
		v <- as.numeric()
		for (i in 1:n){
			# first element
			tmp1 <-   cat[i,1]*((0-phi[i,1])/(Phi[i,1])  )
			# last element
			tmp2 <- cat[i,m]*((phi[i,t] - 0)/(P[i,m])  )
			tmp3 <- 0
			for (j in 2:t){
				for (k in 2:(m-1)){
					tmp3 <- tmp3 + cat[i,k]*(   (phi[i,j-1]-phi[i,j])/(P[i,k])       )
				}
			}
	
			v[i] <-(tmp1 + tmp2 + tmp3)
		}



		# L matrix
		L <- matrix(0, ncol=m-1, nrow=n)
		for (i in 1:n){
	
			for (j in 1:t){
				if (j == 1){
					L[i,j] <- -total[i]*phi[i,j]* ((phi[i,j] - 0   )/P[i,j] - (phi[i,j+1] - phi[i,j]   )/P[i,j+1]  )  

				}
				else if (j == t){
					L[i,j] <- -total[i]*phi[i,j]* ((phi[i,j] - phi[i,j-1]   )/P[i,j] - (0 - phi[i,j]   )/P[i,j+1]  )  

				}
				else {
					L[i,j] <- -total[i]*phi[i,j]* ((phi[i,j] - phi[i,j-1]   )/P[i,j] - (phi[i,j+1] - phi[i,j]   )/P[i,j+1]  ) 
			
				}
		
			}
		}


		# Q matrix
		Q <- matrix(0, ncol=t, nrow=t)
		for (k in 1:t){
			tmp1 <- 0
			for (j in 1:n){
				tmp1 <- tmp1 + total[j]*phi[j,k]^2*(  (P[j,k]+P[j,k+1])/(P[j,k]*P[j,k+1])  )	
			}
			Q[k,k] <- tmp1 
			for (l in (k+1):t){
				tmp2 <- 0
				if (l > t) break
				for (j in 1:n){
					tmp2 <- tmp2 + -total[j]*(  (phi[j,k]*phi[j,k+1])/(P[j,k+1])  )	
				}

			Q[k,l] <- tmp2
			Q[l,k] <- Q[k,l]	
			}
		
		}


		# p matrix
		p <- as.numeric()


		for (k in 1:t) {
			tmp <- 0
			for (j in 1:n ){
				tmp <- tmp + phi[j,k]*(  (cat[j,k])/(P[j,k]) - (cat[j,k+1])/(P[j,k+1])      )
		
			}	
			p[k] <- tmp
		
		}


		mme1 <- cbind(Q, t(L)%*%X, t(L)%*%Z)
		mme2 <- cbind(t(X)%*%L, t(X)%*%W%*%X, t(X)%*%W%*%Z)
		mme3 <- cbind(t(Z)%*%L, t(Z)%*%W%*%X, t(Z)%*%W%*%Z + Ainv*ratio)
		LHS <- rbind(mme1, mme2, mme3)

		mme4 <- matrix(p)
		mme5 <- t(X)%*%v
		mme6 <- t(Z)%*%v - ratio*Ainv%*%matrix(B[(dim(B)[1]-length(s)+1):dim(B)[1],])
		RHS <- rbind(mme4, mme5, mme6)

		newB <- solve(LHS)%*%RHS
		diff <-  (  sum((B + newB - B)^2) ) /sum((B + newB)^2)
		B <- B + newB 
		if (disp == TRUE){
			cat('\n')
			cat("iteration ", it, '\n')
			print(B)
		}


	}


}