# libraries
{
library("forecast")
library("lubridate")
library("xts")
library("readxl")
library("dplyr")
library("readr")
library("rugarch")
library("boot") # not needed?
library("fGarch")
}

# functions
{
  
}

Y <- read_excel("C:/Users/novle/Documents/Programming/R/coursework/N2.xlsx", col_names = FALSE, range = "Ћист1!A2:A1001")
dats <- seq.Date(from = as.Date("2012-05-11"),
                                  by = "day",
                                  length.out = 2000) #length(p1)

#Step 1: preprocess data, estimate coefficients for GARCH(1, 1),
#build innovation series and calculate T_n
Y <- read_xlsx("C:/Users/novle/Documents/Programming/R/coursework/trial.xlsx", col_names = FALSE, range = "CUSUM (1 break)!A2:A1001")
dats <- read_xlsx("C:/Users/novle/Documents/Programming/R/coursework/shares 2015-2020.xlsx", col_names = FALSE, range = "Sheet1!B3:B1002")

Preprocess_data <- function(Y){
  Y <- as.matrix(Y)
  storage.mode(Y) <- "numeric" 
  Y <- as.vector(Y)
  return(Y)
}
Preprocess_dates <- function(dats){
  dats <- as.matrix(dats)
  storage.mode(dats) <- "character"
  dats <- as.vector(dats)
  return(dats)
}

plot(Preprocess_data(Y), xlab = "Observations")
plot(T_k(Preprocess_data(Y)), xlab = "Observations")

#Get_GARCH_params <- function(Y, dats, show_plot=FALSE){
  df <- data.frame(Preprocess_data(Y))
  colnames(df) <- c("Log.yield")
  rownames(df) <- Preprocess_dates(dats)
  
  values <- df$Log.yield
  log_yields <- xts(values, order.by=ymd(rownames(df)))
  if (show_plot == TRUE)
    tsdisplay(log_yields)
   
  model <- ugarchspec(
    variance.model = list(garchOrder = c(1, 1)), 
    mean.model = list(armaOrder = c(0, 0)))
  model_est <- ugarchfit(spec = model, data = log_yields)
  
  #get coefficients [1] - mu, [2] - omega, [3] - alpha1, [4] - beta1
  coefs <- as.vector(coef(model_est))
  return(coefs)
}
get_params <- function(Y){
  est_model <- garchFit(formula = ~ garch(1, 1), data = Y,
                        include.mean = FALSE, trace = FALSE)
  params <- est_model@fit$par[1:3] # We are getting omega, alpha1 and beta1
  return(params)
}
get_sigmas <- function(Y){
  est_model <- garchFit(formula = ~ garch(1, 1), data = Y,
                        include.mean = FALSE, trace = FALSE)
  sigmas <- est_model@sigma.t # We are getting sigmas series
  return(sigmas)
}

#Get_GARCH_params(Preprocess_data(Y), dats)
get_params(Y)

# We are getting mu, omega, alpha1, beta1

innov <- Y/get_sigmas(Y)
res <- innov[,1]

Obs_boot <- function(res, coefs) {
  #Step 2: build bootstrap innovation series
  innov_boot <- sample(res, replace = TRUE) # bootstrapping residuals
  sigmas_boot <- c(var(res)*(length(res)-1)/length(res)) # Bootstraped sigmas - start with sample variance
  #Step 3: build bootstrapped observations recursively
  obs_boot <- c(sqrt(var(res)*(length(res)-1)/length(res))*innov_boot[1]) # Bootstraped sigmas - start with sample variance
  for (i in 2:length(res)) {
    new_sigma_boot <- coefs[1] + coefs[2]*(obs_boot[i-1])^2 + 
      coefs[3]*sigmas_boot[i-1]
    sigmas_boot <- c(sigmas_boot, new_sigma_boot)
    obs_boot <- c(obs_boot, sqrt(new_sigma_boot)*innov_boot[i])
  }
  return(obs_boot)
}

x_new <- Obs_boot(res, get_params(Y)) # save this

#Step 4: calculate bootstrapped T_n statistic for x_new

params <- get_params(x_new)
innov = Y/get_sigmas(x_new)
innov[,1]

T_k <- function(Y, uncorrelated=TRUE) {
  X <- Y*Y
  N <- length(Y) # number of observations in our sample
  Xmean <- mean(X) # mean value of Y^2 series

  KL <- numeric(0) # for Kokoszka-Leipus statistics
  for (k in 1:N){
    KL <- c(KL, abs(sum(X[1:k]) - k*Xmean)/sqrt(N)) # calculating KL statistic
  }
  
  # estimation of variance
  # for innovation series, which are uncorrelated obs.
  if (uncorrelated == TRUE){ 
    Z <- X*X
    s <- sqrt(mean(Z)-mean(X)*mean(X))
  }
  # old version: for series, which are correlated - Bartlett weights cause distortion
  if (uncorrelated == FALSE){
    q <- floor(sqrt(N)) # using square root function (logarithm finds more breaks, but also gives more mistakes)
    w <- c((1:q)/(q+1), ((q+1):1)/(q+1)) # Bartlett weights (triangular kernel with window [-q; q])
    C <- c()
    for (j in 0:q){
      Cj <- 0
      for (i in 1:(N-j)){
        Cj <- Cj + (X[i]-Xmean)*(X[i+j]-Xmean) # sample covariances in [-q; q] window
      }
      C <- c(C, Cj/N)
    }
    C <- c(rev(C[-1]), C) # here we use symmetry of covariance to simplify the code
    s <- sqrt(sum(w*C)) # triangular kernel
  }
  
  return(KL/s)
} # calculates KL statistic for each observation from the sample

BB_crit <- function(p = 0.99) {
  if (p == 0.99) { # significance level 0.99
    return(1.627)
  } else if (p == 0.95) { # significance level 0.95
    return(1.358)
  } else if (p == 0.90) { # significance level 0.90
    return(1.224)
  } else if (p == 0.75) { # significance level 0.75
    return(1.019)
  } else if (p == 0.50) { # significance level 0.50
    return(0.828)
  } else if (p == 0.25) { # significance level 0.25
    return(0.677)
  } else if (p == 0.10) { # significance level 0.10
    return(0.571)
  } else if (p == 0.05) { # significance level 0.05
    return(0.520)
  } #p + 2*(-1)^i*exp(-2*i*i*x*x)
} # asymptotic quantiles for supremum of Brownian bridge's absolute value (KL)

BB_stat <- function(Y, uncorrelated=TRUE) {
  X <- Y*Y
  N <- length(Y) # number of observations in our sample
  Xmean <- mean(X) # mean value of Y^2 series
  
  KL <- numeric(0)
  for (k in 1:N){
    KL <- c(KL, abs(sum(X[1:k]) - k*Xmean)/sqrt(N)) # calculating KL statistic
  }
  tau <- min(which(KL == max(KL))) # supremum for statistics
  
  # estimation of variance
  # for innovation series, which are uncorrelated obs.
  if (uncorrelated == TRUE){ 
    Z <- X*X
    s <- sqrt(mean(Z)-Xmean*Xmean)
  }
  # for series, which are correlated
  if (uncorrelated == FALSE){
    q <- floor(sqrt(N)) # using square root function (logarithm finds more breaks, but also gives more mistakes)
    w <- c((1:q)/(q+1), ((q+1):1)/(q+1)) # Bartlett weights (triangular kernel with window [-q; q])
    C <- c()
    for (j in 0:q){
      Cj <- 0
      for (i in 1:(N-j)){
        Cj <- Cj + (X[i]-Xmean)*(X[i+j]-Xmean) # sample covariances in [-q; q] window
      }
      C <- c(C, Cj/N)
    }
    C <- c(rev(C[-1]), C) # here we use symmetry of covariance to simplify the code
    s <- sqrt(sum(w*C)) # triangular kernel
  }
  
  return(c(tau, KL[tau]/s))
}  # returns suspicious time moment and value of KL statistic

BB_test <- function(Y, uncorrelated=TRUE, p = 0.99) {
  BB <- BB_stat(Y, uncorrelated)
  if (BB[2] > BB_crit(p)){ # check whether statistic is larger than asymptotic quantile or not
    return(c('Break', BB))}
  else
    return(c('No break', BB))
} # returns most suspicious time moment and value of statistic

plot(innov[,1])
T_boot <- c(BB_test(innov[,1])[3])
T_boot
plot(T_k(innov[,1]))
plot(T_k(Preprocess_data(Y)))


options(max.print=100000)

plot(b)
plot(boot(b, stat, R=10, l=90, sim="fixed")$t[1,])

typeof(Y)

b

size2 = floor(1.1*length(b))

stat <- function(x) {
  return(cumsum(x))
}
a <- c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20)
dmean <- 1
sum((a-dmean)*(a-dmean))/(length(a)-1)
var(a)
start_time <- Sys.time()

end_time <- Sys.time()
end_time - start_time
set.seed(626)
dats <- as.xts(dataf)
plot(b)



GARCH_1break <- function(w = 0.1, b = 0.7, g = 0.2, dw = 0.2, db = 0, dg = 0) {
  
  ksi <- rnorm(1999, mean = 0, sd = 1)
  sigma <- c(0)
  Y <- c(rnorm(1, mean = 0, sd = 1))
  
  for (i in 2:1000){
    sigma[i] <- sqrt(w + b*sigma[i-1]^2 + g*Y[i-1]^2)
    Y[i] <- ksi[i-1]*sigma[i]
  }
  
  w <- w + dw # external shock w_2 = w_1 + dw (0.2 by default)
  b <- b + db # external shock b_2 = b_1 + db (0 by default)
  g <- g + dg # external shock g_2 = g_1 + dg (0 by default)
  
  for (i in 1001:2000){
    sigma[i] <- sqrt(w + b*sigma[i-1]^2 + g*Y[i-1]^2)
    Y[i] <- ksi[i-1]*sigma[i]
  }
  return(Y)
}
GARCH_2breaks <- function(w = 0.1, b = 0.7, g = 0.2, dw1 = 0.2, db1 = 0, dg1 = 0, dw2 = 0.1, db2 = -0.1, dg2 = 0) {
  
  ksi <- rnorm(2999, mean = 0, sd = 1)
  sigma <- c(0)
  Y <- c(rnorm(1, mean = 0, sd = 1))
  
  for (i in 2:1000){
    sigma[i] <- sqrt(w + b*sigma[i-1]^2 + g*Y[i-1]^2)
    Y[i] <- ksi[i-1]*sigma[i]
  }
  
  w <- w + dw1 # external shock w_2 = w_1 + dw1 (0.2 by default)
  b <- b + db1 # external shock b_2 = b_1 + db1 (0 by default)
  g <- g + dg1 # external shock g_2 = g_1 + dg1 (0 by default)
  
  for (i in 1001:2000){
    sigma[i] <- sqrt(w + b*sigma[i-1]^2 + g*Y[i-1]^2)
    Y[i] <- ksi[i-1]*sigma[i]
  }
  
  w <- w + dw2 # external shock w_3 = w_2 + dw2 (0.1 by default)
  b <- b + db2 # external shock b_3 = b_2 + db2 (-0.1 by default)
  g <- g + dg2 # external shock g_3 = g_2 + dg2 (0 by default)
  
  for (i in 2001:3000){
    sigma[i] <- sqrt(w + b*sigma[i-1]^2 + g*Y[i-1]^2)
    Y[i] <- ksi[i-1]*sigma[i]
  }
  return(Y)
}

Y <- GARCH_2breaks(dw2 = -0.2)
plot(Y, xlab = "Observations")
plot(KL(Y), xlab = "Observations")

print(paste0('ƒл€ данных ',length(Y), ' наблюдений с помощью KL-статистики вы€влен структурный сдвиг в момент времени ', 't_break - ', BB_stat(Y)[1]))
print(paste0('ƒл€ данных ', length(Y), ' наблюдений с помощью KL-статистики вы€влены следующие моменты структурных сдвигов: ', 't_break - ', ICSS_refinement(Y, p = 0.99)))