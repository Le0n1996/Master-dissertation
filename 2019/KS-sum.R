library(readxl)


# Algorithm based on sum of KS statistics

KS_sum_stat <- function(Y, h = 4) {  # calculating sum of KS statistics for each observation
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
} 
tau_KS_validation <- function(Y, h = 4, p = 0.99) {   # checking the possible t with KS test
  T <- length(Y)
  delta <- 200
  KSstats <- KS_sum_stat(Y)
  tau_susp <- which(KSstats (Y) == min(KSstats(Y))) # search minimum in KS statistic vector
  tau <- c()
  
  Y_1 <- Y[1:max((tau_susp-delta), h)]
  Y_2 <- Y[min((tau_susp+delta), (T-h)):T]
  if (as.numeric(ks.test(Y_1, Y_2)[2]) < 1-p) { # testing on significance level p
    tau <- c(tau, tau_susp)
    return(tau)
  }
  else
    return('отсутствуют')
}

# Data preprocessing

#Y <- read_excel("C:/Users/Ћев/Documents/R/coursework/ICSS_easy.xls", col_names = FALSE)
Y <- read_excel("C:/Users/Ћев/Documents/R/coursework/ICSS_check.xls", col_names = FALSE, range = "CUSUM_KL!A1:A2000")
#Y <- read_excel("C:/Users/Ћев/Documents/R/coursework/ICSS_check.xls", col_names = FALSE, range = "ICSS_KL!A1:A3000")
#Y <- read_excel("C:/Users/Ћев/Documents/R/coursework/amzn.xlsx", col_names = TRUE) 

Y <- read_excel("C:/Users/Ћев/Documents/R/coursework/GARCH_with_breaks.xlsx", col_names = FALSE, range = "CUSUM (1 break)!A1:A2000")
#Y <- read_excel("C:/Users/Ћев/Documents/R/coursework/GARCH_with_breaks.xlsx", col_names = FALSE, range = "ICSS (2 breaks)!A1:A3000")


Y <- as.matrix(Y) 
storage.mode(Y) <- "numeric" 
Y <- as.vector(Y) 

plot(KS_sum_stat(Y))
print(paste0('ƒл€ данных ', length(Y), ' наблюдений вы€влен следующий момент структурного сдвига: ', 't_break - ', tau_KS_validation(Y, p = 0.99)))