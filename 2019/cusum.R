# CUSUM
BB_stat <- function(Y, base = exp(1)) {
  Y2 <- Y*Y
  T <- length(Y) # number of observations
  Y2mean <- sum(Y2)/T

  KL <- 0
  for (k in 1:T){
    KL <- c(KL, (sum(Y2[1:k]) - k*Y2mean)/sqrt(T))
  }
  KL <- abs(KL[-1])
  tau <- min(which(KL == max(KL)))

  r <- floor(log(T, base)) # here we can take log or log10, as it doesn't change convergence or r/T to 0
  w <- c((1:r)/(r+1), ((r+1):1)/(r+1)) # Bartlett weights
  
  C <- 0
  for (j in 0:r){
    Cj <- 0
    for (i in 1:(T-j)){
      Cj <- Cj + (Y2[i]-Y2mean)*(Y2[i+j]-Y2mean) #sample covariances
    }
    C <- c(C, Cj/T)
  }
  C <- C[-1]
  C <- c(rev(C[-1]), C) # symmetry of covariance
  v <- sqrt(sum(w*C)) # triangular kernel
  
  return(c(tau, KL[tau]/v))
}

Y <- c(1, 2, 3, 5, 8, 1, 2, 3, 5, 8)

BB_test <- function(Y, base = exp(1), p = 0.05) {
  if (p == 0.05)
    BB_krit <- 1.358
  if (p == 0.01)
    BB_krit <- 1.628
  BB <- BB_stat(Y, base)
  return(BB)
}

print(crit(0.05))

# ICSS
print(BB_stat(Y))
