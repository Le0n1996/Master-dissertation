library(readxl)

# CUSUM

BB_crit <- function(p = 0.99) {
  if (p == 0.05) # significance level 0.05
    BB_cr <- 0.520
  if (p == 0.10) # significance level 0.10
    BB_cr <- 0.571
  if (p == 0.25) # significance level 0.25
    BB_cr <- 0.677
  if (p == 0.50) # significance level 0.50
    BB_cr <- 0.828
  if (p == 0.75) # significance level 0.75
    BB_cr <- 1.019
  if (p == 0.90) # significance level 0.90
    BB_cr <- 1.224
  if (p == 0.95) # significance level 0.95
    BB_cr <- 1.358
  if (p == 0.99) # significance level 0.99
    BB_cr <- 1.628
  return(BB_cr)
}

BB_stat <- function(Y, base = exp(1)) {
  Y2 <- Y*Y
  T <- length(Y) # number of observations
  Y2mean <- sum(Y2)/T
  
  KL <- 0
  for (k in 1:T){
    KL <- c(KL, (sum(Y2[1:k]) - k*Y2mean)/sqrt(T)) # calcu
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

BB_test <- function(Y, base = exp(1), p = 0.99) {
  BB <- BB_stat(Y, base)
  if (BB[2] > BB_crit(p))
    return(c('Break', BB))
  else
    return(c('No break', BB))
}

# ICSS

#Y <- read_excel("C:/Users/Лев/Documents/R/coursework/ICSS_easy.xls", col_names = FALSE)
#Y <- read_excel("C:/Users/Лев/Documents/R/coursework/ICSS_check.xls", col_names = FALSE, range = "CUSUM_KL!A1:A2000")
#Y <- read_excel("C:/Users/Лев/Documents/R/coursework/ICSS_check.xls", col_names = FALSE, range = "ICSS_KL!A1:A3000")
#Y <- read_excel("C:/Users/Лев/Documents/R/coursework/ICSS_check.xls", col_names = FALSE, range = "CUSUM_KL!B1:B2000")
Y <- read_excel("C:/Users/Лев/Documents/R/coursework/amzn.xlsx", col_names = TRUE) 

Y=as.matrix(Y)
storage.mode(Y) <- "numeric" 
Y=as.vector(Y)

ICSS_iter <- function(Y, base = exp(1)) {
  
  T <- length(Y) # number of observations
  t1 <- 1
  t2 <- T # starting endpoints
  tfirst <- 0
  tlast <- T
  tau <- c(tfirst, tlast) # all possible points of structural break (including the borders)
  tcus <- 0
  
  while (tfirst < tlast - 1) {
    if (BB_test(Y[t1:t2], base)[1] == 'No break'){
      tau <- sort(unique(tau))
      return(tau)
    } 
    else {
      tcus <- as.numeric(BB_test(Y[t1:t2], base)[2])+t1-1 # take into account that BB_test()[2] returns relative position
      tau <- c(tau, tcus) # remember first break-like point
      t2 <- tcus
      
      repeat {
        if (BB_test(Y[t1:t2], base)[1] == 'No break'){
          break
        } 
        else {
          t2 <- as.numeric(BB_test(Y[t1:t2], base)[2])+t1-1
        }
      }
      
      tfirst <- t2
      tau <- c(tau, tfirst) # remembering break-like point
      t1 <- tcus + 1
      t2 <- T # new endpoints for next iteration of an algorithm
      
      repeat {
        if (BB_test(Y[t1:t2], base)[1] == 'No break'){
          break
        } 
        else {
          t1 <- as.numeric(BB_test(Y[t1:t2], base)[2])+1 + t1-1
        }
      }
      
      tlast <- t1 - 1
      tau <- c(tau, tlast) # remembering break-like point
      t1 <- tfirst + 1 
      t2 <- tlast # new endpoints for next iteration of an algorithm
    }
  }
  tau <- sort(unique(tau))
  return(tau)
}

print(ICSS_iter(Y))

BB_test(Y[0:1989])
BB_test(Y[1099:3000])