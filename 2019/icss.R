library(readxl)
library(readr)


# CUSUM test realisation

KL <- function(Y){
  X <- Y*Y
  T <- length(Y) # number of observations in our sample
  Xmean <- sum(X)/T # mean value of Y^2 series
  
  KL <- 0
  for (k in 1:T){
    KL <- c(KL, (sum(X[1:k]) - k*Xmean)/sqrt(T)) # calculating KL statistic
  }
  KL <- abs(KL[-1])
  q <- floor(sqrt(T)) # using square root function (logarithm finds more breaks, but also gives more mistakes)
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
  return(KL/s)
} # calculates KL statistic for each observation from the sample

BB_crit <- function(p = 0.01) {
  if (p == 0.95) # significance level 0.05
    BB_cr <- 0.520
  if (p == 0.90) # significance level 0.10
    BB_cr <- 0.571
  if (p == 0.75) # significance level 0.25
    BB_cr <- 0.677
  if (p == 0.50) # significance level 0.50
    BB_cr <- 0.828
  if (p == 0.25) # significance level 0.75
    BB_cr <- 1.019
  if (p == 0.10) # significance level 0.90
    BB_cr <- 1.224
  if (p == 0.05) # significance level 0.95
    BB_cr <- 1.358
  if (p == 0.01) # significance level 0.99
    BB_cr <- 1.627
  return(BB_cr)
} # asymptotic quantiles for supremum of Brownian bridge's absolute value (KL)

BB_st <- function(Y) {
  X <- Y*Y
  T <- length(Y) # number of observations in our sample
  Xmean <- sum(X)/T # mean value of Y^2 series

  KL <- 0
  for (k in 1:T){
    KL <- c(KL, (sum(X[1:k]) - k*Xmean)/sqrt(T)) # calculating KL statistic
  }
  KL <- abs(KL[-1])
  tau <- min(which(KL == max(KL))) # suspicious time moment
  
  q <- floor(sqrt(T)) # using square root function (logarithm finds more breaks, but also gives more mistakes)
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
}  # returns suspicious time moment and value of KL statistic

BB_test <- function(Y, p = 0.01) {
  BB <- BB_stat(Y)
  if (BB[2] > BB_crit(p)){ # check whether statistic is larger than asymptotic quantile or not
    return(c('Break', BB))}
  else
    return(c('No break', BB))
} # returns most suspicious time moment and value of statistic

# ICSS procedure

ICSS_iter <- function(Y, p = 0.01) {
  
  T <- length(Y) # number of observations in our sample
  t1 <- 1
  t2 <- T # starting endpoints
  tfirst <- 0
  tlast <- T
  tau <- c(tfirst, tlast) # all possible points of structural break (including the borders)
  tcus <- 0
  
  while (tfirst < tlast - 1) {
    if (BB_test(Y[t1:t2], p)[1] == 'No break'){
      #print('t2:')
      #print(t2)
      tau <- sort(unique(tau))
      return(tau) # stop the procedure, if there are no more breaks
    } 
    else {
      tcus <- as.numeric(BB_test(Y[t1:t2], p)[2]) + t1-1 # taking into account that BB_test returns relative position
      tau <- c(tau, tcus) # remember first break-like point
      t2 <- tcus # looking at left interval
      
      repeat {
        if (BB_test(Y[t1:t2], p)[1] == 'No break'){
          break # go further, if there are no more breaks to the left
        } 
        else {
          t2 <- as.numeric(BB_test(Y[t1:t2], p)[2])+t1-1
        }
      }

      tfirst <- t2
      tau <- c(tau, tfirst) # remember break-like point
      t1 <- tcus + 1
      t2 <- T # new endpoints for next part of an algorithm
      
      repeat {
        if (BB_test(Y[t1:t2], p)[1] == 'No break'){
          break # go further, if there are no more breaks to the right
        } 
        else {
          t1 <- as.numeric(BB_test(Y[t1:t2], p)[2])+1 + t1-1
        }
      }
      
      tlast <- t1 - 1
      tau <- c(tau, tlast) # remember break-like point
      t1 <- tfirst + 1 
      t2 <- tlast # new endpoints for next iteration of an algorithm
    }
  }
  tau <- sort(unique(tau))
  return(tau) # stop the procedure, if the remaining interval is too small
} # returns series of suspicious moments

ICSS_refinement <- function(Y, p = 0.01) {
  tau <- ICSS_iter(Y, p) # getting series of possible structural breaks
  if ((length(tau) <= 2)) { # in case there are no structural breaks at all
    #tau <- BB_stat(Y)[1] 
    #D <- BB_stat(Y)[2] # if info about p-value of structural break is needed
    return('отсутствуют')
  }
    
  tau_ref <- tau[2:(length(tau)-1)] # series of possible structural breaks without 0 and T
  iteration <- 1
  while (iteration < 20) { # setting limit to number of iterations; though, it coincides quickly
    
    tau <- c(0, tau_ref, length(Y))
    tau_ref <- c()
    K <- length(tau)
    
    for (n in 2:(K-1)){
      tprev <- tau[n-1]+1
      tnext <- tau[n+1]
      if (BB_test(Y[tprev:tnext], p)[1] == 'Break') # check if the moment is really a structural break for interval from two adjacent moments 
        tau_ref <- c(tau_ref, tau[n-1] + as.numeric(BB_test(Y[tprev:tnext], p)[2]))
    }
    
    iteration <- iteration + 1
    
    if ((length(tau_ref) == length(tau[2:(length(tau)-1)])) & ((is.null(tau_ref)) | (max(abs(tau_ref-tau[2:(length(tau)-1)]))<3)))
      break # end refinement process if there is no change in tau or if change is very small
  }
  return(tau_ref)
} # refining the moments of structural breaks

GARCH_QUAK <- function(N = 1000, d = 0.25){
  pi <- numeric(0)
  pi <- c(pi, 1)
  d <- 0.25
  for (j in 1:100){
    pi <- c(pi, pi[length(pi)]*(j-d-1)/j)
  }
  ksi <- rnorm(N-1, mean = 0, sd = 1)
  sigma <- numeric(0)
  Y <- c(rnorm(1, mean = 0, sd = 1))
  pi
  for (i in 2:N){
    sigma_new <- sqrt(w + b*sigma[i-1]^2 + g*Y[i-1]^2)
    sigma <- c(sigma, sigma_new)
    Y[i] <- ksi[i-1]*sigma[i]
  } 
  return(Y) # generates GARCH(1,1) with given parameters and length N
}

# Algorithm based on sum of KS statistics (one break case)

KS_sum_stat <- function(Y, h = 4) {
  T <- length(Y)
  D <- c()
  if (T >= 8) {
    for (k in h:(T-h)){
      l <- floor(k/2)
      r <- floor((T+k)/2)
      Y_L1 <- Y[1:l]
      Y_L2 <- Y[(l+1):(k-1)]
      Y_R1 <- Y[k:r]
      Y_R2 <- Y[(r+1):T] # generates 4 samples of size not less than 2 each
      
      D_left <- as.numeric(ks.test(Y_L1, Y_L2)[1])
      D_right <- as.numeric(ks.test(Y_R1, Y_R2)[1]) # KS statistic for left and right pair
      D <- c(D, D_left+D_right)
    }
    return(D) # returns vector of statistics for each observation
  }
  else
    return('More observations are needed!')
} # calculating sum of KS statistics for each observation

tau_KS_validation <- function(Y, h = 4, p = 0.99) { 
  T <- length(Y)
  delta <- 200
  KSstats <- KS_sum_stat(Y)
  tau_susp <- which(KSstats == min(KSstats)) # search minimum in KS statistic vector
  tau <- c()
  
  Y_1 <- Y[1:max((tau_susp-delta), h)]
  Y_2 <- Y[min((tau_susp+delta), (T-h)):T]
  if (as.numeric(ks.test(Y_1, Y_2, exact = FALSE)[2]) < 1-p) { # testing on significance level p
    tau <- c(tau, tau_susp)
    return(tau)
  }
  else
    return('отсутствуют')
} # checking the possible t with KS test

# Data preprocessing

#Y <- read_excel("C:/Users/novle/Documents/Programming/R/coursework/GARCH_with_breaks.xlsx", col_names = FALSE, range = "CUSUM (1 break)!A1:A2000")
#Y <- read_excel("C:/Users/novle/Documents/Programming/R/coursework/GARCH_with_breaks.xlsx", col_names = FALSE, range = "ICSS (2 breaks)!A1:A3000")

# Microsoft
#Y <- read_excel("C:/Users/novle/Documents/Programming/R/coursework/US - 2007-2018.xlsx", col_names = FALSE, range = "US.1999-2018!D3:D5066")
#dates <- read_excel("C:/Users/novle/Documents/Programming/R/coursework/US - 2007-2018.xlsx", col_names = FALSE, range = "US.1999-2018!B3:B5066")

# Amazon
#Y <- read_excel("C:/Users/novle/Documents/Programming/R/coursework/US - 2007-2018.xlsx", col_names = FALSE, range = "US.1999-2018!I3:I4743")
#dates <- read_excel("C:/Users/novle/Documents/Programming/R/coursework/US - 2007-2018.xlsx", col_names = FALSE, range = "US.1999-2018!G3:G4743")

# Apple
#Y <- read_excel("C:/Users/novle/Documents/Programming/R/coursework/US - 2007-2018.xlsx", col_names = FALSE, range = "US.1999-2018!N3:N5066")
#dates <- read_excel("C:/Users/novle/Documents/Programming/R/coursework/US - 2007-2018.xlsx", col_names = FALSE, range = "US.1999-2018!L3:L5066")

# Verizon
#Y <- read_excel("C:/Users/novle/Documents/Programming/R/coursework/US - 2007-2018.xlsx", col_names = FALSE, range = "US.2007-2018!D3:D2946")
#dates <- read_excel("C:/Users/novle/Documents/Programming/R/coursework/US - 2007-2018.xlsx", col_names = FALSE, range = "US.2007-2018!B3:B2946")

# Boeing
#Y <- read_excel("C:/Users/novle/Documents/Programming/R/coursework/US - 2007-2018.xlsx", col_names = FALSE, range = "US.2007-2018!I3:I2953")
#dates <- read_excel("C:/Users/novle/Documents/Programming/R/coursework/US - 2007-2018.xlsx", col_names = FALSE, range = "US.2007-2018!G3:G2953")

# Exxon Mobil
#Y <- read_excel("C:/Users/novle/Documents/Programming/R/coursework/US - 2007-2018.xlsx", col_names = FALSE, range = "US.2007-2018!N3:N2956")
#dates <- read_excel("C:/Users/novle/Documents/Programming/R/coursework/US - 2007-2018.xlsx", col_names = FALSE, range = "US.2007-2018!L3:L2956")

# Walmart
#Y <- read_excel("C:/Users/novle/Documents/Programming/R/coursework/US - 2007-2018.xlsx", col_names = FALSE, range = "US.2007-2018!S3:S2957")
#dates <- read_excel("C:/Users/novle/Documents/Programming/R/coursework/US - 2007-2018.xlsx", col_names = FALSE, range = "US.2007-2018!Q3:Q2957")

# Walt Disney
#Y <- read_excel("C:/Users/novle/Documents/Programming/R/coursework/US - 2007-2018.xlsx", col_names = FALSE, range = "US.2007-2018!X3:X2986")
#dates <- read_excel("C:/Users/novle/Documents/Programming/R/coursework/US - 2007-2018.xlsx", col_names = FALSE, range = "US.2007-2018!V3:V2986")

# Johnson & Johnson
#Y <- read_excel("C:/Users/novle/Documents/Programming/R/coursework/US - 2007-2018.xlsx", col_names = FALSE, range = "US.2007-2018!AC3:AC3018")
#dates <- read_excel("C:/Users/novle/Documents/Programming/R/coursework/US - 2007-2018.xlsx", col_names = FALSE, range = "US.2007-2018!AA3:AA3018")

# Coca Cola
#Y <- read_excel("C:/Users/novle/Documents/Programming/R/coursework/US - 2007-2018.xlsx", col_names = FALSE, range = "US.2007-2018!AH3:AH3018")
#dates <- read_excel("C:/Users/novle/Documents/Programming/R/coursework/US - 2007-2018.xlsx", col_names = FALSE, range = "US.2007-2018!AF3:AF3018")

# Google (Alphabet)
#Y <- read_excel("C:/Users/novle/Documents/Programming/R/coursework/US - 2007-2018.xlsx", col_names = FALSE, range = "US.2007-2018!AM3:AM2948")
#dates <- read_excel("C:/Users/novle/Documents/Programming/R/coursework/US - 2007-2018.xlsx", col_names = FALSE, range = "US.2007-2018!AK3:AK2948")

###

# SBER (Sberbank)
Y <- read_excel("C:/Users/novle/Documents/Programming/R/coursework/shares 2015-2020.xlsx", col_names = FALSE, range = "Sheet1!D3:D1316")
dates <- read_excel("C:/Users/novle/Documents/Programming/R/coursework/shares 2015-2020.xlsx", col_names = FALSE, range = "Sheet1!B3:B1316")

Y <- read_excel("C:/Users/novle/Documents/Programming/R/coursework/shares 2015-2020.xlsx", col_names = FALSE, range = "Sheet1!H3:H11834")
dates <- read_excel("C:/Users/novle/Documents/Programming/R/coursework/shares 2015-2020.xlsx", col_names = FALSE, range = "Sheet1!F3:F11834")

Y <- read_excel("C:/Users/novle/Documents/Programming/R/coursework/shares 2015-2020.xlsx", col_names = FALSE, range = "Sheet1!L3:L13591")
dates <- read_excel("C:/Users/novle/Documents/Programming/R/coursework/shares 2015-2020.xlsx", col_names = FALSE, range = "Sheet1!J3:J13591")
Z <- read_excel("C:/Users/novle/Documents/Programming/R/coursework/shares 2015-2020.xlsx", col_names = FALSE, range = "Sheet1!K3:K13591")

Y <- read_excel("C:/Users/novle/Documents/Programming/R/coursework/shares 2015-2020.xlsx", col_names = FALSE, range = "Sheet1!Q3:Q1314")
dates <- read_excel("C:/Users/novle/Documents/Programming/R/coursework/shares 2015-2020.xlsx", col_names = FALSE, range = "Sheet1!O3:O1314")

Y <- read_excel("C:/Users/novle/Documents/Programming/R/coursework/shares 2015-2020.xlsx", col_names = FALSE, range = "Sheet1!U3:U23635")
dates <- read_excel("C:/Users/novle/Documents/Programming/R/coursework/shares 2015-2020.xlsx", col_names = FALSE, range = "Sheet1!S3:S23635")

Y <- read_excel("C:/Users/novle/Documents/Programming/R/coursework/shares 2015-2020.xlsx", col_names = FALSE, range = "Sheet1!Y3:Y16136")
dates <- read_excel("C:/Users/novle/Documents/Programming/R/coursework/shares 2015-2020.xlsx", col_names = FALSE, range = "Sheet1!W3:W16136")

Y <- as.matrix(Y)
storage.mode(Y) <- "numeric" 
Y <- as.vector(Y)
#Z <- as.matrix(Z) 
#storage.mode(Z) <- "numeric" 
#Z <- as.vector(Z) 

dates <- as.matrix(dates)
storage.mode(dates) <- "character"
dates <- as.vector(dates)

#plot(Z)

plot(Y, xlab = "Observations")
plot(KL(Y), xlab = "Observations")
#plot(KS_sum_stat(Y))

print(paste0('ƒл€ данных ', length(Y), ' наблюдений с помощью KL-статистики вы€влены следующие моменты структурных сдвигов: ', 't_break - ', ICSS_refinement(Y, p = 0.99)))
print(paste0('ƒл€ данных ', length(Y), ' наблюдений с помощью KL-статистики вы€влены следующие моменты структурных сдвигов: ', 't_break - ', ICSS_refinement(Y, p = 0.95)))
print(dates[19756])
print(KL(Y)[19755])
print(KL(Y)[19756])
print(KL(Y)[19757])
# Microsoft
print(dates[1118])
#print(paste0('ƒл€ данных ', length(Y), ' наблюдений с помощью суммы KS-статистик вы€влен следующий момент структурного сдвига: ', 't_break - ', tau_KS_validation(Y, p = 0.99)))
#print(dates[986])

# Amazon
print(dates[893])
print(dates[2715])

# Apple
#print(dates[901])
#print(dates[2610])

# Verizon (only for p = 0.95) !! approximately one year difference
#print(BB_stat(Y))
#print(paste0('ƒл€ данных ', length(Y), ' наблюдений с помощью KL-статистики вы€влены следующие моменты структурных сдвигов: ', 't_break - ', ICSS_refinement(Y, p = 0.95)))
#print(dates[333])
#print(dates[498])

# Boeing
#print(dates[1181])
#print(paste0('ƒл€ данных ', length(Y), ' наблюдений с помощью суммы KS-статистик вы€влен следующий момент структурного сдвига: ', 't_break - ', tau_KS_validation(Y, p = 0.99)))
#print(dates[1100])

# Exxon Mobil (only for p = 0.95)
#print(BB_stat(Y))
#print(dates[492])
#print(paste0('ƒл€ данных ', length(Y), ' наблюдений с помощью суммы KS-статистик вы€влен следующий момент структурного сдвига: ', 't_break - ', tau_KS_validation(Y, p = 0.99)))
#print(dates[756])

# Walmart (only for p = 0.95)
#print(paste0('ƒл€ данных ', length(Y), ' наблюдений с помощью KL-статистики вы€влены следующие моменты структурных сдвигов: ', 't_break - ', ICSS_refinement(Y, p = 0.95)))
#print(dates[533])
#print(dates[2106])

# Walt Disney (only for p = 0.95)
#print(BB_stat(Y))
#print(dates[1227])
#print(paste0('ƒл€ данных ', length(Y), ' наблюдений с помощью суммы KS-статистик вы€влен следующий момент структурного сдвига: ', 't_break - ', tau_KS_validation(Y, p = 0.99)))
#print(dates[2008])

# Johnson & Johnson (No structural breaks even for p = 0.90, but for KS-method one break found exactly at point of maximum of KS-statistic)
#print(plot(KL(Y)))
#print(paste0('ƒл€ данных ', length(Y), ' наблюдений с помощью суммы KS-статистик вы€влен следующий момент структурного сдвига: ', 't_break - ', tau_KS_validation(Y, p = 0.95)))
#print(dates[1407])

# Coca Cola (only for p = 0.95) !!! one year difference
#print(BB_stat(Y))
#print(dates[624])
#print(paste0('ƒл€ данных ', length(Y), ' наблюдений с помощью суммы KS-статистик вы€влен следующий момент структурного сдвига: ', 't_break - ', tau_KS_validation(Y, p = 0.99)))
#print(dates[879])

# Google (Alphabet) (only for p = 0.95)
print(dates[520])
#print(paste0('ƒл€ данных ', length(Y), ' наблюдений с помощью суммы KS-статистик вы€влен следующий момент структурного сдвига: ', 't_break - ', tau_KS_validation(Y, p = 0.99)))
#print(dates[602])