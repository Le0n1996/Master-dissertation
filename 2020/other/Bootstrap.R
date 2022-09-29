# libraries
{
  library("readxl")
  library("fGarch")
  library("ggplot2")
  library("hrbrthemes")
  library("latex2exp")
}

# functions
{
  # Converting and plotting
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
  # import data and show plot for splits detection
  take_data <- function(path, show_plot = TRUE){
    data <- read_xlsx(path, col_names = FALSE, range = cell_cols(5)) # data downloaded from finam.ru
    data <- Preprocess_data(data) # convert to numeric vector
    data <- data[2:length(data)] # delete first row with <ticket>
    if (show_plot == TRUE)
      plot(data) # checking for splits
    return(data)
  }
  dates <- function(path){
    dats <- read_xlsx(path, col_names = FALSE, range = cell_cols(3))
    dats <- Preprocess_dates(dats) # convert to character vector
    dats <- dats[2:length(dats)] # delete first row with <ticket>
    for (i in 1:length(dats)){
      dats[i] <- paste(substr(dats[i], 1, 4), substr(dats[i], 5, 6), substr(dats[i], 7, 8), sep="-")
    }
    return(dats)
  }
  # detect split (if necessary)
  detect_split <- function(data){
    mod1 <- data[2:length(data)]
    mod2 <- data[1:(length(data)-1)]
    ratio <- mod1/mod2
    k1 <- which(ratio == min(ratio)) # candidate for usual split
    k2 <- which(ratio == max(ratio)) # candidate for reversed split
    
    print(paste0('Возможный момент для сплита акций: ', data[k2], " до и ", data[k2+1], " после."))
    print(paste0('Возможный момент для обратного сплита акций: ', data[k1], " до и ", data[k1+1], " после."))
    return(c(k1, k2)) # two candidates for split moment
  }
  
  # Estimate parameters or sigmas of GARCH(1,1)
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
  Obs_boot <- function(inn, innov_boot, coefs) {
    sigmas_boot <- c(var(inn)*(length(inn)-1)/length(inn)) # Bootstraped sigmas - start with sample variance
    obs_boot <- c(sqrt(var(inn)*(length(inn)-1)/length(inn))*innov_boot[1])
    
    for (i in 2:length(inn)) {
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
  
  # Tests:
  # For 0 structural breaks with bootstrap
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
  # For presence of 1 structural break with bootstrap
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
  
  # Iterated algorithm
  ICSS <- function(Y, B=100, p = 0.05) {
    T <- length(Y) # number of observations in our sample
    t1 <- 1
    t2 <- T # starting endpoints
    tfirst <- 1
    tlast <- T
    tau <- c(1, T+1) # for all possible points of structural break
    
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
  # Algorithm, refining the moments of structural breaks
  ICSS_refinement <- function(Y, B = 100, p = 0.05) {
    tau <- ICSS(Y, B, p) # getting series of possible structural breaks
    if ((length(tau) <= 2)) { # in case if there are no structural breaks at all
      #tau <- BB_stat(Y)[1] 
      #D <- BB_stat(Y)[2] # if info about p-value of structural break is needed
      return(NULL)
    }
    #print(tau) # uncomment if we want to see all steps
    tau_ref <- tau[2:(length(tau)-1)] # series of possible structural breaks without 0 and T
    iteration <- 1
    while (iteration <= 5) { # setting limit to number of iterations; though, it coincides quickly
      
      tau <- c(1, tau_ref, length(Y)+1)
      tau_ref <- c()
      K <- length(tau) # remember total number of structural breaks
      
      for (n in 2:(K-1)){
        # check if the moment is really a structural break
        # for interval between two adjacent potential breaks 
        tprev <- tau[n-1]
        tnext <- tau[n+1]-1
        #print(c(tprev, tnext)) # uncomment if we want to see all steps
        #print(iteration) # uncomment if we want to see all steps
        if (test_for_no_breaks(Y[tprev:tnext], B)[1] < p){ 
          tau_ref <- c(tau_ref, tau[n-1]-1 + k_hat(Y[tprev:tnext]))
        }
      }
      #print(tau_ref) # uncomment if we want to see all steps
      iteration <- iteration + 1
      #print(length(tau_ref)) # uncomment if we want to see all steps
      #print(length(tau[2:(length(tau)-1)])) # uncomment if we want to see all steps
      #print(tau) # uncomment if we want to see all steps
      if (is.null(tau_ref)) # end refinement process if there are no breaks after all
        break
      if (length(tau_ref) == length(tau[2:(K-1)])){
        if (max(abs(tau_ref-tau[2:(K-1)]))<5)
          break # end refinement process if there is no change in tau or if change is very small
      }
    }
    return(sort(unique(tau_ref)))
  }
  
  # Functions for GARCH simulations
  GARCH_0breaks <- function(N = 1000, w = 0.1, b = 0.7, g = 0.2){
    ksi <- rnorm(N-1, mean = 0, sd = 1)
    sigma <- c(0)
    Y <- c(rnorm(1, mean = 0, sd = 1))
    
    for (i in 2:N){
      sigma[i] <- sqrt(w + b*sigma[i-1]^2 + g*Y[i-1]^2)
      Y[i] <- ksi[i-1]*sigma[i]
    } 
    return(Y) # generates GARCH(1,1) with given parameters and length N
  }
  GARCH_1break <- function(N = 2000, w = 0.1, b = 0.7, g = 0.2, dw = 0.2, db = 0, dg = 0) {
    
    ksi <- rnorm(N-1, mean = 0, sd = 1)
    sigma <- c(0)
    Y <- c(rnorm(1, mean = 0, sd = 1))
    
    for (i in 2:(N/2)){
      sigma[i] <- sqrt(w + b*sigma[i-1]^2 + g*Y[i-1]^2)
      Y[i] <- ksi[i-1]*sigma[i]
    }
    
    w <- w + dw # external shock w_2 = w_1 + dw (0.2 by default)
    b <- b + db # external shock b_2 = b_1 + db (0 by default)
    g <- g + dg # external shock g_2 = g_1 + dg (0 by default)
    
    for (i in (N/2+1):N){
      sigma[i] <- sqrt(w + b*sigma[i-1]^2 + g*Y[i-1]^2)
      Y[i] <- ksi[i-1]*sigma[i]
    }
    return(Y) # generates two GARCH(1,1) with given parameters and length N in total
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
}

Test_T200 <- GARCH_0breaks(200)
Test_T500 <- GARCH_0breaks(500)
Test_T1000 <- GARCH_0breaks(1000)
Test_T2000 <- GARCH_0breaks(2000)
Test_T5000 <- GARCH_0breaks(5000)

Test_M200 <- GARCH_1break(200)
Test_M500 <- GARCH_1break(500)
Test_M1000 <- GARCH_1break(1000)
Test_M2000 <- GARCH_1break(2000)
Test_M5000 <- GARCH_1break(5000)

plot(Test_T200)
plot(Test_T500)
plot(Test_T1000)
plot(Test_T2000)
plot(Test_T5000)

plot(Test_M200)
plot(Test_M500)
plot(Test_M1000)
plot(Test_M2000)
plot(Test_M5000)

start_time <- Sys.time()
t200 <- test_for_no_breaks(Test_T200, B=10000)
end_time <- Sys.time()
end_time - start_time

start_time <- Sys.time()
t500 <- test_for_no_breaks(Test_T500, B=10000)
end_time <- Sys.time()
end_time - start_time

start_time <- Sys.time()
t1000 <- test_for_no_breaks(Test_T1000, B=10000)
end_time <- Sys.time()
end_time - start_time

start_time <- Sys.time()
t2000 <- test_for_no_breaks(Test_T2000, B=10000)
end_time <- Sys.time()
end_time - start_time

start_time <- Sys.time()
t5000 <- test_for_no_breaks(Test_T5000, B=10000)
end_time <- Sys.time()
end_time - start_time

start_time <- Sys.time()
m200 <- test_for_one_str_break(Test_M200, B=10000)
end_time <- Sys.time()
end_time - start_time

start_time <- Sys.time()
m500 <- test_for_one_str_break(Test_M500, B=10000)
end_time <- Sys.time()
end_time - start_time

start_time <- Sys.time()
m1000 <- test_for_one_str_break(Test_M1000, B=10000)
end_time <- Sys.time()
end_time - start_time

start_time <- Sys.time()
m2000 <- test_for_one_str_break(Test_M2000, B=10000)
end_time <- Sys.time()
end_time - start_time

start_time <- Sys.time()
m5000 <- test_for_one_str_break(Test_M5000, B=10000)
end_time <- Sys.time()
end_time - start_time

t200[3:10002]
t500[3:10002]
t1000[3:10002]
t2000[3:10002]
t5000[3:10002]
m200[4:10003]
m500[4:10003]
m1000[4:10003]
m2000[4:10003]
m5000[4:10003]

sort(b[3:202])[200*0.99]
sort(c[3:202])[200*0.99]

dats <- dates("DIS.xlsx")[1:1000]

vec <- sort(t500[3:10002])
Tstat <- vec
q_09 <- Tstat[10000*0.9]
q_095 <- Tstat[10000*0.95]
q_099 <- Tstat[10000*0.99]


q_09 <- Mstat[10000*0.9]
q_095 <- Mstat[10000*0.95]
q_099 <- Mstat[10000*0.99]
# now you need to ggplot it
plot(get_sigmas(cvx_log_yields), xlab = "Observations", ylab = "Chevron: Conditional variances before detecting breaks")
renewed <- c(get_sigmas(cvx_log_yields[1:2538]), get_sigmas(cvx_log_yields[2539:length(cvx_log_yields)]))
plot(renewed, xlab = "Observations", ylab = "Chevron: Conditional variances after detecting breaks")

df = data.frame(value = 1:length(cvx_log_yields), p = renewed)
plot(df)

ggplot(data = df, mapping = aes(x = value, y = p)) + 
  geom_point() + labs(x = "Наблюдения", y = "Условная дисперсия") +
  ggtitle("Пересчитанная условная дисперсия") + theme_ipsum() + theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5),
axis.title.x = element_text(color = "black", face = "bold.italic", size = 12, hjust = 0.5, vjust = 1.5),
axis.title.y = element_text(color = "black", size = 12, hjust = 0.5, vjust = 4.5)) + scale_x_continuous(breaks=seq(0, 2586, by = 500)) + scale_y_continuous(breaks=seq(0, 0.18, by = 0.03)) +
  geom_segment(aes(x=2533, y=renewed[2533]-0.002, xend=2533, yend=renewed[2533]+0.002), size = 2, color = "red")
  
save.image(file = "my_work_space_26_05.RData")
getwd()


graph <- ggplot(data = df, mapping = aes(x = value, y = p)) + 
  geom_point() + labs(x = TeX("Значение статистики $T_n$"), y = "Эмпирическая функция распределения") +
  xlim(vec[1],vec[1000]) + scale_x_continuous(breaks=seq(0.5, 2.0, by = 0.5)) +
  scale_y_continuous(breaks=c(0, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99)) +
  ggtitle(TeX("Распределение статистики $T_n$ (n = 500)")) + theme_ipsum()
graph <- graph + geom_segment(aes(x=vec[1], y=0.9, xend=q_09, yend=0.9), linetype="dashed", color = "blue") +
  geom_segment(aes(x=vec[1], y=0.95, xend=q_095, yend=0.95), linetype="dashed", color = "blue") +
  geom_segment(aes(x=vec[1], y=0.99, xend=q_099, yend=0.99), linetype="dashed", color = "blue") +
  geom_segment(aes(x=q_09, y=0, xend=q_09, yend=0.9), linetype="dashed", color = "blue") +
  geom_segment(aes(x=q_095, y=0, xend=q_095, yend=0.95), linetype="dashed", color = "blue") +
  geom_segment(aes(x=q_099, y=0, xend=q_099, yend=0.99), linetype="dashed", color = "blue")
graph <- graph + theme(
  axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5),
  axis.title.x = element_text(color = "black", face = "bold.italic", size = 12, hjust = 0.5, vjust = 1.5),
  axis.title.y = element_text(color = "black", size = 12, hjust = 0.5, vjust = 4.5))
graph + annotate(geom="text", x=q_09-0.005, y=-0.04, label=round(q_09, digits = 3), color = 'blue') +
  annotate(geom="text", x=q_095+0.005, y=-0.04, label=round(q_095, digits = 3), color = 'blue') +
  annotate(geom="text", x=q_099, y=-0.04, label=round(q_099, digits = 3), color = 'blue')

ggplot(stat_function(fun = asymp, color = "red") )
asymp <- function(x){
  k <- 1:20
  k2 <- k*k
  val <- 2*(-1)^k*exp(-2*k2*x^2)
  1+sum(val)
}
asymp(0.35809863)

t500[3:10002]
t200[1:2]
t500[1:2]
t1000[1:2]
t2000[1:2]
t5000[1:2]

plot(Test_T200)

plot(sort(t500[3:10002]))


plot(new_one)

start_time <- Sys.time()
kl_stats <- numeric(0)
boot_p_vals <- numeric(0)
for (i in 1:200){
  new_one <- GARCH_1break(N=200, w=0.1, b=0.8, g=0.1)
  #new_one[250] <- 5 #different test
  kl_stats <- c(kl_stats, BB_stat(new_one, uncorrelated = FALSE)[2])
  boot_p_vals <- c(boot_p_vals, test_for_no_breaks(new_one, B=500)[1])
}
end_time <- Sys.time()
end_time - start_time

results_kl_200 <- sort(kl_stats)
results_boot_200 <- sort(boot_p_vals)
plot(results_kl_200)
plot(results_boot_200)
length(results_kl_200[results_kl_200>=1.358])/200
length(results_kl_200[results_kl_200>=1.628])/200
length(results_boot_200[results_boot_200<=0.05])/200
length(results_boot_200[results_boot_200<=0.01])/200
results
