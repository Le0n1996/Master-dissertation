library(readxl)

tau_KS <- function(Y, h = 4) {
  T <- length(Y)
  D <- c()
  for (k in h:(T-h)){
    Y_L1 <- Y[1:floor(k/2)]
    Y_L2 <- Y[(floor(k/2)+1):(k-1)]
    Y_R1 <- Y[k:(floor((T+k)/2))]
    Y_R2 <- Y[(floor((T+k)/2)+1):T] # generates 4 samples
    
    D_left <- as.numeric(ks.test(Y_L1, Y_L2)[1])
    D_right <- as.numeric(ks.test(Y_R1, Y_R2)[1]) # KS statistic for left and right pair
    D <- c(D, D_left+D_right)
  }
  return(which(D == min(D))) # search for all minimums
} # calculate tau_KS

# Data preprocessing

#Y <- read_excel("C:/Users/Лев/Documents/R/coursework/ICSS_easy.xls", col_names = FALSE)
#Y <- read_excel("C:/Users/Лев/Documents/R/coursework/ICSS_check.xls", col_names = FALSE, range = "CUSUM_KL!A1:A2000")
Y <- read_excel("C:/Users/Лев/Documents/R/coursework/ICSS_check.xls", col_names = FALSE, range = "ICSS_KL!A1:A3000")
#Y <- read_excel("C:/Users/Лев/Documents/R/coursework/amzn.xlsx", col_names = TRUE) 

Y <- as.matrix(Y) 
storage.mode(Y) <- "numeric" 
Y <- as.vector(Y) 

#Y <- c(1,2,3,1,2,3,1,2,3,1,2,3,1,2,3,1,2,3,1,2,3)

#for (i in 1:(k-1))
#  prob_diff_left <- c(prob_diff_left, abs(ecdf(Y_L1)(Y[i])-ecdf(Y_L2)(Y[i])))
#for (j in k:T)
#  prob_diff_right <- c(prob_diff_right, abs(ecdf(Y_R1)(Y[j])-ecdf(Y_R2)(Y[j])))
#gg <- c(1, 2, 0.1, 3, 2, 0.01, 4, 2.5, 10, 20, 100, 4, 4, 32, 5, 12, 23, 1, 2, 9)
#print(sort(gg))
#print(ecdf(gg)(30))

#plot(D)

tau_susp <- tau_KS(Y)

delta <- 400
h <- 4
T <- length(Y)
tau <- c()

Y_1 <- Y[1:max((tau_susp-delta), h)]
Y_2 <- Y[min((tau_susp+delta), (T-h)):T]
if (as.numeric(ks.test(Y_1, Y_2)[2]) < 0.05)
  tau <- c(tau, tau_susp)
Y <- Y[(tau+1):T]

print(tau)
tau_susp <- tau_KS(Y)

delta <- 400
h <- 4
T <- length(Y)

Y_1 <- Y[1:max((tau_susp-delta), h)]
Y_2 <- Y[min((tau_susp+delta), (T-h)):T]
if (as.numeric(ks.test(Y_1, Y_2)[2]) < 0.05){
  tau <- c(tau, tau_susp)
}
print(ks.test(Y_1, Y_2))

print(tau)
print(tau_susp)

Y_1 <- Y[805:max((tau_susp-delta), h)]
Y_2 <- Y[min((tau_susp+delta), (T-h)):T]
if (as.numeric(ks.test(Y_1, Y_2)[2]) < 0.05){
  tau <- c(tau, tau_susp)
}
print(ks.test(Y_1, Y_2))