# libraries
{
library("readxl")
library("readr")
library("fGarch")
}

# functions
{
Preprocess_data <- function(Y){
  Y <- as.matrix(Y)
  storage.mode(Y) <- "numeric" 
  Y <- as.vector(Y)
  return(Y)
} # preprocess observations
Preprocess_dates <- function(dats){
  dats <- as.matrix(dats)
  storage.mode(dats) <- "character"
  dats <- as.vector(dats)
  return(dats)
} # preprocess dates

# estimate parameters or sigmas of GARCH(1,1)
get_params <- function(Y){
  suppressWarnings(est_model <- garchFit(formula = ~ garch(1, 1), data = Y,
                        include.mean = FALSE, trace = FALSE))
  suppressWarnings(params <- est_model@fit$par[1:3]) # We are getting omega, alpha1 and beta1
  return(params) # comment: warnings are suppressed, because garchFit sometimes have problems with std errors
}
get_sigmas <- function(Y){
  suppressWarnings(est_model <- garchFit(formula = ~ garch(1, 1), data = Y,
                        include.mean = FALSE, trace = FALSE))
  suppressWarnings(sigmas <- est_model@sigma.t) # We are getting sigmas series
  return(sigmas) # comment: warnings are suppressed, because garchFit sometimes have problems with std errors

}

# builds bootstrapped observations using bootrapped innovation series
Obs_boot <- function(res, innov_boot, coefs) {
  sigmas_boot <- c(var(res)*(length(res)-1)/length(res)) # Bootstraped sigmas - start with sample variance
  obs_boot <- c(sqrt(var(res)*(length(res)-1)/length(res))*innov_boot[1]) # Bootstraped sigmas - start with sample variance
  
  for (i in 2:length(res)) {
    new_sigma_boot <- coefs[1] + coefs[2]*(obs_boot[i-1])^2 + 
      coefs[3]*sigmas_boot[i-1]
    sigmas_boot <- c(sigmas_boot, new_sigma_boot)
    obs_boot <- c(obs_boot, sqrt(new_sigma_boot)*innov_boot[i])
  }
  return(obs_boot)
}

# calculates KL statistic for each observation from the sample
T_k <- function(Y, uncorrelated=TRUE) {
  X <- Y*Y
  N <- length(Y) # number of observations in our sample
  Xmean <- mean(X) # mean value of Y^2 series
  
  KL <- numeric(0) # for Kokoszka-Leipus statistics
  for (k in 1:N){
    KL <- c(KL, abs(sum(X[1:k]) - k*Xmean)/sqrt(N)) # calculating KL statistic
  }
  
  # estimation of variance (two options)
  # new version: for innovation series, which are uncorrelated observations
  if (uncorrelated == TRUE){ 
    Z <- X*X
    s <- sqrt(mean(Z)-mean(X)*mean(X))
  }
  # old version: for series, which are correlated - but Bartlett weights cause size distortion
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
}

# calculates k_hat - esimated moment of structural break
k_hat <- function(Y, uncorrelated=TRUE) {
  X <- Y*Y
  N <- length(Y) # number of observations in our sample
  
  k_hats <- numeric(0)
  for (k in 1:(N-1)){
    k_hats <- c(k_hats, abs((N-k)*sum(X[1:k]) - k*sum(X[(k+1):N])))
  }
  k <- min(which(k_hats == max(k_hats))) # maximum
  return(k)
} 

# returns asymptotic quantiles for supremum of Brownian bridge's absolute value
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
  } #p = 1 + 2*sum((-1)^k*exp(-2*k*k*x*x))
}
# returns suspicious time moment and value of KL statistic
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
    s <- sqrt(mean(Z)-mean(X)*mean(X))
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
}
# returns suspicious time moment, value of statistic, and decision
BB_test <- function(Y, uncorrelated=TRUE, p = 0.99) {
  BB <- BB_stat(Y, uncorrelated=TRUE)
  if (BB[2] > BB_crit(p)){ # check whether statistic is larger than asymptotic quantile or not
    return(c('Break', BB))}
  else
    return(c('No break', BB))
}

# testing for 0 structural breaks with bootstrap
test_for_no_breaks <- function(Y, B=100){
  #Step 1: take data and estimate coefficients for GARCH(1,1)
  params0 <- get_params(Y)
  # We are getting omega, alpha1, beta1
  
  #Step 2: build innovation series
  innov <- Y/get_sigmas(Y)
  
  #Step 3: calculate residual-based CUSUM test statistic
  T_start <- BB_stat(innov)[2]
  
  #Step 4: use bootstrap
  T_boot <- numeric(0)
  for (i in 1:B){
    #Step 4.1: build bootstrapped observations recursively
    innov_boot <- sample(innov, replace = TRUE) # bootstrapping innovation series
    x_new <- Obs_boot(innov, innov_boot, params0)
    
    #Step 4.2: calculate residual-based bootstrapped T_n statistic
    innov_boot <- x_new/get_sigmas(x_new)
    T_boot <- c(T_boot, BB_stat(innov_boot)[2])
  }
  
  #Step 5: calculate p-value for the test
  p_val <- mean(T_boot>T_start)
  return(c(p_val, T_start, T_boot)) # returns p-value and bootstrapped statistic
}

# testing for presence of 1 structural break with bootstrap
test_for_one_str_break <- function(Y, B=100){
  #Step 1: get the estimated change point and estimate coefficients for two GARCH(1,1) models
  K <- k_hat(Y)
  # get estimated change point
  Y1 <- Y[1:K]
  Y2 <- Y[(K+1):length(Y)]
  params1 <- get_params(Y1)
  params2 <- get_params(Y2)
  # get omega, alpha1, beta1 for both models
  
  #Step 2: build innovation series for both models
  innov1 <- Y1/get_sigmas(Y1)
  innov2 <- Y2/get_sigmas(Y2)
  
  #Step 3: calculate residual-based CUSUM test statistic
  T_st_1 <- BB_stat(innov1)[2]
  T_st_2 <- BB_stat(innov2)[2]
  M_start <- max(T_st_1, T_st_2)
  
  #Step 4: use bootstrap
  M_boot <- numeric(0)
  for (i in 1:B){
    #Step 4.1: build bootstrapped observations recursively
    innov_boot_1 <- sample(innov1, replace = TRUE)
    innov_boot_2 <- sample(innov2, replace = TRUE)
    # bootstrapping innovation series
    x_new_1 <- Obs_boot(innov1, innov_boot_1, params1)
    x_new_2 <- Obs_boot(innov2, innov_boot_2, params2)
    
    #Step 4.2: calculate bootstrapped M_n statistic for x_new
    innov_boot_1 <- x_new_1/get_sigmas(x_new_1)
    innov_boot_2 <- x_new_2/get_sigmas(x_new_2)
    T_boot_1 <- BB_stat(innov_boot_1)[2]
    T_boot_2 <- BB_stat(innov_boot_2)[2]
    M_boot <- c(M_boot, max(T_boot_1, T_boot_2))
  }
  
  #Step 5: calculate p-value
  p_val <- mean(M_boot>M_start)
  return(c(p_val, M_start, K, M_boot))
}

# returns all structural breaks on a given significance level
ICSS <- function(Y, B=100, p = 0.05) {
  T <- length(Y) # number of observations in our sample
  t1 <- 1
  t2 <- T # starting endpoints
  tfirst <- 1
  tlast <- T
  tau <- numeric(0) # for all possible points of structural break
  
  while (tfirst < tlast - 1) {
    if (test_for_no_breaks(Y[t1:t2], B)[1] > p){
      tau <- sort(unique(tau))
      return(tau) # stop the procedure, if there are no breaks at all
    }
    if (test_for_one_str_break(Y[t1:t2], B)[1] > p){
      tcus <- k_hat(Y[t1:t2]) + t1-1 # taking into account that k_hat returns relative position
      tau <- c(tau, tcus)
      tau <- sort(unique(tau))
      return(tau) # stop the procedure, if there is only 1 break on this interval
    }
    #if (BB_test(Y[t1:t2], p)[1] == 'No break'){
    #print('t2:')
    #print(t2)
    #  tau <- sort(unique(tau))
    #  return(tau) # stop the procedure, if there are no more breaks
    #} 
    else {
      tcus <- k_hat(Y[t1:t2]) + t1-1 # taking into account that k_hat returns relative position
      tau <- c(tau, tcus)
      t2 <- tcus-1 # looking at left interval
      
      repeat {
        if (test_for_no_breaks(Y[t1:t2], B)[1] > p){
          tau <- c(tau, tcus)
          break # go further, if there are no breaks at all on this interval
        }
        if (test_for_one_str_break(Y[t1:t2], B)[1] > p){
          tcus <- k_hat(Y[t1:t2]) + t1-1 # taking into account that k_hat returns relative position
          tau <- c(tau, tcus)
          break # go further, if there are no more breaks to the left
        } 
        else {
          tcus <- k_hat(Y[t1:t2]) + t1-1 # remember for the case if there are no breaks to the left
          t2 <- k_hat(Y[t1:t2])-1 + t1-1
        }
      }
      
      tfirst <- tcus
      t1 <- tfirst
      t2 <- T # new endpoints for next part of an algorithm
      
      repeat {
        if (test_for_no_breaks(Y[t1:t2], B)[1] > p){
          tau <- c(tau, tcus)
          break # go further, if there are no breaks at all on this interval
        }
        if (test_for_one_str_break(Y[t1:t2], B)[1] > p){
          tcus <- k_hat(Y[t1:t2]) + t1-1 # taking into account that k_hat returns relative position
          tau <- c(tau, tcus)
          break # go further, if there are no more breaks to the right
        } 
        else {
          tcus <- k_hat(Y[t1:t2]) + t1-1 # remember for the case if there are no breaks to the right
          t1 <- k_hat(Y[t1:t2]) + t1-1
        }
      }
      
      tlast <- t1
      t1 <- tfirst
      t2 <- tlast-1 # new endpoints for next iteration of an algorithm
    }
  }
  tau <- sort(unique(tau))
  return(tau) # stop the procedure, if the remaining interval is too small
}

# refining the moments of structural breaks
ICSS_refinement <- function(Y, B = 100, p = 0.05) {
  tau <- ICSS(Y, B, p) # getting series of possible structural breaks
  if ((length(tau) == 0)) { # in case if there are no structural breaks at all
    #tau <- BB_stat(Y)[1] 
    #D <- BB_stat(Y)[2] # if info about p-value of structural break is needed
    return('отсутствуют')
  }
  
  tau_ref <- tau # series of possible structural breaks without 0 and T
  iteration <- 1
  while (iteration < 5) { # setting limit to number of iterations; though, it coincides quickly
    
    tau <- c(1, tau_ref, length(Y)+1)
    tau_ref <- c()
    K <- length(tau) # remember total number of structural breaks
    
    for (n in 2:(K-1)){
      # check if the moment is really a structural break
      # for interval between two adjacent potential breaks 
      tprev <- tau[n-1]
      tnext <- tau[n+1]-1
      
      if (test_for_no_breaks(Y[tprev:tnext], B)[1] < p){ 
        tau_ref <- c(tau_ref, tau[n-1]-1 + k_hat(Y[tprev:tnext]))
      }
    }
    
    iteration <- iteration + 1
    
    if ((length(tau_ref) == length(tau[2:(length(tau)-1)])) & ((is.null(tau_ref)) | (max(abs(tau_ref-tau[2:(length(tau)-1)]))<5)))
      break # end refinement process if there is no change in tau or if change is very small
  }
  return(tau_ref)
}
}

#take simulated data
#Y <- read_xlsx("GARCH_with_breaks.xlsx", col_names = FALSE, range = cell_cols(1))
#Y <- Preprocess_data(Y)

#Y <- read_xlsx("GARCH_with_2breaks.xlsx", col_names = FALSE, range = cell_cols(1))
#Y <- Preprocess_data(Y)

#set p-value:
p <- 0.05

jnj <- read_xlsx("JNJ.xlsx", col_names = FALSE, range = cell_cols(5))
jnj <- Preprocess_data(jnj) # convert to numeric vector
jnj <- jnj[2:length(jnj)] #delete first row with <ticket>
plot(jnj) # checking for splits
Y <- diff(log(jnj))

#plot(Y, xlab = "Observations")
plot(T_k(Y/get_sigmas(Y)), xlab = "Observations", ylab = "Statistic for innovations")
plot(T_k(Y, uncorrelated = FALSE), xlab = "Observations", ylab = "Statistic for original observations")

ICSS(Y, p=0.05)

start_time <- Sys.time()
ICSS_refinement(Y, p=0.05)
end_time <- Sys.time()
end_time - start_time

print(paste0("ћоменты структурного сдвига - наблюдени€ ", K))
ind <- c(1,1029,1979,1990,3000)

for (i in 1:4){
  print(if_any_str_break(Y[ind[i]:ind[i+1]])[1])
}
for (i in 1:3){
  print(if_any_str_break(Y[ind[i]:ind[i+2]])[1])
}
for (i in 1:3){
  print(str_break(Y[ind[i]:ind[i+2]]))[1]
}

a<-if_any_str_break(Y[1990:3000])
b<-if_any_str_break(Y[600:1400],1000)
c<-if_any_str_break(Y[K:2000])
if_any_str_break(Y)[1]

d<-str_break(Y[600:1400],1000)

A0<-if_any_str_break(Y[1:1000],5000)
A00<-if_any_str_break(Y,5000)
A000<-if_any_str_break(Y[1001:2000],5000)

A1<-str_break(Y,5000)
A11<-str_break(Y[1:1000],5000)
A111<-str_break(Y[1001:2000],5000)

start_time <- Sys.time()
B <- test_for_one_str_break(Y, 1000)
end_time <- Sys.time()
end_time - start_time
C <- B[4:length(B)]


plot(Y)
plot(T_k(Y/get_sigmas(Y)), xlab = "Observations", ylab = "Statistic for innovations")
plot(get_sigmas(Y), xlab = "Observations", ylab = "Conditional variances before detecting breaks")
renewed <- c(get_sigmas(Y[1:(k_hat(Y)-1)]), get_sigmas(Y[k_hat(Y):5000]))
plot(renewed, xlab = "Observations", ylab = "Conditional variances after detecting breaks")

goog <- read_xlsx("C:/Users/novle/Downloads/Goog.xlsx", col_names = FALSE, range = cell_cols(5))
goog <- Preprocess_data(goog) # convert to numeric vector
goog <- goog[2:length(goog)] #delete first row with <ticket>
Y <- diff(log(goog))

plot(Y, xlab = "Observations")
plot(T_k(Y/get_sigmas(Y)), xlab = "Observations")
plot(T_k(Y))

Y <- GARCH_2breaks()
plot(Y, xlab = "Observations")
plot(T_k(Y/get_sigmas(Y)), xlab = "Observations")

start_time <- Sys.time()
ICSS_refinement(Y, p=0.05)
end_time <- Sys.time()
end_time - start_time
test_for_no_breaks(Y[1:1117])
test_for_no_breaks(Y[1117:3000])

k_hat(Y)

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
GARCH_2breaks <- function(w = 0.1, b = 0.7, g = 0.2, dw1 = 0.2, db1 = 0, dg1 = 0, dw2 = 0.1, db2 = -0.1, dg2 = 0.05) {
  
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