BB_stat_l <- function(Y) {
  X <- Y*Y
  T <- length(Y) # number of observations in our sample
  Xmean <- sum(X)/T # mean value of Y^2 series
  
  KL <- 0
  for (k in 1:T){
    KL <- c(KL, (sum(X[1:k]) - k*Xmean)/sqrt(T)) # calculating KL statistic
  }
  KL <- abs(KL[-1])
  tau <- min(which(KL == max(KL))) # suspicious time moment
  
  q <- log10(T) # using square root function (logarithm finds more breaks, but also gives more mistakes)
  w <- c((1:q)/(q+1), ((q+1):1)/(q+1)) # Bartlett weights (triangular kernel with window [-q; q])
  
  C <- c()
  for (j in 0:q){
    Cj <- 0
    for (i in 1:(T-j)){
      Cj <- Cj + (X[i]-Xmean)*(X[i+j]-Xmean) # sample covariances in [-q; q] window
    }
    C <- c(C, Cj/T)
  }
  C <- c(rev(C[-1]), C) # here we use symmetry of covariance to simplify the code
  s <- sqrt(sum(w*C)) # triangular kernel
  return(c(tau, KL[tau]/s))
}  
BB_test_l <- function(Y, p = 0.99) {
  BB <- BB_stat_l(Y)
  if (BB[2] > BB_crit(p)){ # check whether statistic is larger than asymptotic quantile or not
    return(c('Break', BB))}
  else
    return(c('No break', BB))
}

w <- 0.1
b <- 0.7
g <- 0.2
logg <- c()
sqrtt <- c()

for (n in 1:100){
  ksi <- rnorm(1999, mean = 0, sd = 1)
  sigma <- c(0)
  Y <- c(rnorm(1, mean = 0, sd = 1))
  
  for (i in 2:1000){
    sigma[i] <- sqrt(w + b*sigma[i-1]^2 + g*Y[i-1]^2)
    Y[i] <- ksi[i-1]*sigma[i]
  }
  
  w <- 0.3 # external shock w_2 = w_1 + 0.2
  
  for (i in 1001:2000){
    sigma[i] <- sqrt(w + b*sigma[i-1]^2 + g*Y[i-1]^2)
    Y[i] <- ksi[i-1]*sigma[i]
  }
  if (BB_test(Y)[1] == 'Break')
    sqrtt <- c(sqrtt, BB_test(Y, 0.90)[2])
  if (BB_test_l(Y)[1] == 'Break')
    logg <- c(logg, BB_test_l(Y, 0.90)[2])
}
print(sqrtt[!(900<sqrtt) | !(sqrtt<1100)])
print(logg[!(900<logg) | !(logg<1100)])

plot(Y, xlab = "Observations")
plot(KL(Y), xlab = "Observations")

print(paste0('ƒл€ данных ',length(Y), ' наблюдений с помощью KL-статистики вы€влен структурный сдвиг в момент времени ', 't_break - ', BB_stat(Y)[1]))
