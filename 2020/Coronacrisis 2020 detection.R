# Libraries:
{
  library("readxl")
  library("fGarch")
}

# Functions:
{
  # Converting and plotting
  Preprocess_data <- function(X){
    X <- as.matrix(X)
    storage.mode(X) <- "numeric" 
    X <- as.vector(X)
    return(X)
  } # preprocess observations
  Preprocess_dates <- function(dats){
    dats <- as.matrix(dats)
    storage.mode(dats) <- "character"
    dats <- as.vector(dats)
    return(dats)
  } # preprocess dates
  # import data and show plot for splits detection
  take_data <- function(path, show_plot = TRUE){
    data <- read_xlsx(path, col_names = FALSE, range = cell_cols(5)) # <CLOSE> for data downloaded from finam.ru
    data <- Preprocess_data(data) # convert to numeric vector
    data <- data[2:length(data)] # delete first row with <ticker>
    if (show_plot == TRUE)
      plot(data) # checking for splits
    return(data)
  }
  dates <- function(path){
    dats <- read_xlsx(path, col_names = FALSE, range = cell_cols(3)) # <DATE> for data downloaded from finam.ru
    dats <- Preprocess_dates(dats) # convert to character vector
    dats <- dats[2:length(dats)] # delete first row with <ticket>
    for (i in 1:length(dats)){ # converting date to type yyyy-mm-dd
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
    
    print(paste0('Возможный момент для сплита акций: Стоимость акции ', data[k2], " до и ", data[k2+1], " после."))
    print(paste0('Возможный момент для обратного сплита акций: Стоимость акции ', data[k1], " до и ", data[k1+1], " после."))
    return(c(k1, k2)) # two candidates for split moment
  }
  
  # Estimate parameters of GARCH(1,1) and get sigmas
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
  
  # Builds bootstrapped observations using bootrapped innovation series
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
  
  # Calculates sum of squares for statistic for each observation from the sample
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
  # Calculates k_hat - estimated moment of structural break
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
  # Returns suspicious time moment and value of KL statistic
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
  # For 0 structural breaks with bootstrap (returns p-value)
  # Null hypothesis: 0 breaks, alternative: 1 or more breaks
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
  # For presence of 1 structural break with bootstrap (returns p-value)
  # Null hypothesis: exactly 1 structural break, alternative: 2 or more breaks
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
  
  # Iterated algorithm, which uses two previous tests
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
  # Algorithm, refining the moments of structural breaks given by ICSS procedure
  ICSS_refinement <- function(Y, B = 100, p = 0.05) {
    tau <- ICSS(Y, B, p) # getting series of possible structural breaks
    if ((length(tau) <= 2)) { # in case if there are no structural breaks at all
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
      #print(tau) # uncomment if we want to see all steps
      #print(tau_ref) # uncomment if we want to see all steps
      iteration <- iteration + 1
      
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

# Significance level
p <- 0.01 # the default one for this work

# Importing data:
{
  # Data from January 2010 to April 2020 (~2500 obs.)
  # technology
  {
  msft <- take_data("MSFT.xlsx", TRUE)
  msft_dates <- dates("MSFT.xlsx")
  msft_log_yields <- diff(log(msft))
  
  aapl <- take_data("AAPL.xlsx", TRUE)
  aapl_dates <- dates("AAPL.xlsx")
  moments <- detect_split(aapl)
  aapl <- c(aapl[1:moments[1]], 7*aapl[(moments[1]+1):length(aapl)])
  plot(aapl) # checking for result
  aapl_log_yields <- diff(log(aapl))
  
  ibm <- take_data("IBM.xlsx", TRUE)
  ibm_dates <- dates("IBM.xlsx")
  ibm_log_yields <- diff(log(ibm))
  
  intc <- take_data("INTC.xlsx", TRUE)
  intc_dates <- dates("INTC.xlsx")
  intc_log_yields <- diff(log(intc))
  }
  
  # entertainment
  {
  dis <- take_data("DIS.xlsx", TRUE)
  dis_dates <- dates("DIS.xlsx")
  dis_log_yields <- diff(log(dis))
  }
  
  # financial
  {
  jpm <- take_data("JPM.xlsx", TRUE)
  jpm_dates <- dates("JPM.xlsx")
  jpm_log_yields <- diff(log(jpm))
  }
  
  # healthcare
  {
  jnj <- take_data("JNJ.xlsx", TRUE)
  jnj_dates <- dates("JNJ.xlsx")
  jnj_log_yields <- diff(log(jnj)) # get logarithmic yields 
  }

  # retail
  {
  wmt <- take_data("WMT.xlsx", TRUE)
  wmt_dates <- dates("WMT.xlsx")
  wmt_log_yields <- diff(log(wmt))
  }
  
  # energy
  {
  xom <- take_data("XOM.xlsx", TRUE)
  xom_dates <- dates("XOM.xlsx")
  xom_log_yields <- diff(log(xom))
  
  cvx <- take_data("CVX.xlsx", TRUE)
  cvx_dates <- dates("CVX.xlsx")
  cvx_log_yields <- diff(log(cvx))
  }
  
  # aerospace & defence
  {
  ba <- take_data("BA.xlsx", TRUE)
  ba_dates <- dates("BA.xlsx")
  ba_log_yields <- diff(log(ba))
  }
  
  # consumers - food and beverages
  {
  ko <- take_data("KO.xlsx", TRUE)
  ko_dates <- dates("KO.xlsx")
  moments <- detect_split(ko)
  ko <- c(ko[1:moments[1]], 2*ko[(moments[1]+1):length(ko)])
  plot(ko) # checking for result
  ko_log_yields <- diff(log(ko))
  
  mcd <- take_data("MCD.xlsx", TRUE)
  mcd_dates <- dates("MCD.xlsx")
  mcd_log_yields <- diff(log(mcd))
  }
  
##########
  
  # Data from April 2019 to April 2020 (~ 250 obs.)
  # communication services and entertainment
  {
  googl <- take_data("GOOGL.xlsx", TRUE)
  googl_dates <- dates("GOOGL.xlsx")
  googl_log_yields <- diff(log(googl))
    
  fb <- take_data("FB.xlsx", TRUE)
  fb_dates <- dates("FB.xlsx")
  fb_log_yields <- diff(log(fb))
  
  nflx <- take_data("NFLX.xlsx", TRUE)
  nflx_dates <- dates("NFLX.xlsx")
  nflx_log_yields <- diff(log(nflx))
  }
  
  #retail
  {
  amzn <- take_data("AMZN.xlsx", TRUE)
  amzn_dates <- dates("AMZN.xlsx")
  amzn_log_yields <- diff(log(amzn))
  }
  
  # consumers - food and beverages
  {
  pep <- take_data("PEP.xlsx", TRUE)
  pep_dates <- dates("PEP.xlsx")
  pep_log_yields <- diff(log(pep))
  
  sbux <- take_data("SBUX.xlsx", TRUE)
  sbux_dates <- dates("SBUX.xlsx")
  sbux_log_yields <- diff(log(sbux))
  }

##########
  
  # Extra:
  # Zoom (startup)
  #zm <- take_data("ZM.xlsx", TRUE)
  #zm_dates <- dates("ZM.xlsx")
  #zm_log_yields <- diff(log(zm))
  
##########
}

# Results:

#estimate common GARCH(1,1) parameters
mean(get_params(aapl_log_yields)[2], get_params(ba_log_yields)[2],
get_params(cvx_log_yields)[2], get_params(dis_log_yields)[2],
get_params(ibm_log_yields)[2], get_params(intc_log_yields)[2],
get_params(jnj_log_yields)[2], get_params(jpm_log_yields)[2],
get_params(ko_log_yields)[2], get_params(mcd_log_yields)[2],
get_params(msft_log_yields)[2], get_params(wmt_log_yields)[2],
get_params(xom_log_yields)[2])

mean(get_params(aapl_log_yields)[3], get_params(ba_log_yields)[3],
     get_params(cvx_log_yields)[3], get_params(dis_log_yields)[3],
     get_params(ibm_log_yields)[3], get_params(intc_log_yields)[3],
     get_params(jnj_log_yields)[3], get_params(jpm_log_yields)[3],
     get_params(ko_log_yields)[3], get_params(mcd_log_yields)[3],
     get_params(msft_log_yields)[3], get_params(wmt_log_yields)[3],
     get_params(xom_log_yields)[3])


# healthcare
print("Johnson&Johnson: ")
plot(KL(jnj_log_yields))
print(ICSS_refinement(jnj_log_yields, p)) # one break
print(jnj_dates[2010]) # signals that volatility shift is 2018/01/19
plot(get_sigmas(jnj_log_yields), xlab = "Observations", ylab = "IBM: Conditional variances before detecting breaks")
renewed <- c(get_sigmas(jnj_log_yields[1:2009]), get_sigmas(jnj_log_yields[2010:length(jnj_log_yields)]))
plot(renewed, xlab = "Observations", ylab = "IBM: Conditional variances after detecting breaks")
# least affected by coronacrisis, as healthcare organisation

print("Disney: ")
print(ICSS_refinement(dis_log_yields, p)) # no breaks
print("Walmart: ")
print(ICSS_refinement(wmt_log_yields, p)) # one break
print(wmt_dates[1413]) # signals that volatility shift was long ago, 2015/08/14

# technology sector
print("Microsoft: ")
print(ICSS_refinement(msft_log_yields, p)) # no breaks
print("Apple: ")
print(ICSS_refinement(aapl_log_yields, p)) # no breaks
print("IBM: ")
print(ICSS_refinement(ibm_log_yields, p)) # one break
print(ibm_dates[2191]) # signals that volatility shift is 2018/10/08
print("Intel: ")
print(ICSS_refinement(intc_log_yields, p)) # one break
print(intc_dates[2012]) # signals that volatility shift is 2018/01/24
plot(get_sigmas(ibm_log_yields), xlab = "Observations", ylab = "IBM: Conditional variances before detecting breaks")
renewed <- c(get_sigmas(ibm_log_yields[1:2190]), get_sigmas(ibm_log_yields[2191:length(ibm_log_yields)]))
plot(renewed, xlab = "Observations", ylab = "IBM: Conditional variances after detecting breaks")
plot(get_sigmas(intc_log_yields), xlab = "Observations", ylab = "Intel: Conditional variances before detecting breaks")
renewed <- c(get_sigmas(intc_log_yields[1:2011]), get_sigmas(intc_log_yields[2012:length(intc_log_yields)]))
plot(renewed, xlab = "Observations", ylab = "Intel: Conditional variances after detecting breaks")

# financial sector
print("JPMorgan: ")
print(ICSS_refinement(jpm_log_yields, p)) # one break
print(jpm_dates[2533]) # signals that volatility shift is 2020/02/20
plot(get_sigmas(jpm_log_yields), xlab = "Observations", ylab = "JPMorgan: Conditional variances before detecting breaks")
renewed <- c(get_sigmas(jpm_log_yields[1:2532]), get_sigmas(jpm_log_yields[2533:length(jpm_log_yields)]))
plot(renewed, xlab = "Observations", ylab = "JPMorgan: Conditional variances after detecting breaks")

# energy sector
print("ExxonMobil: ")
print(ICSS_refinement(xom_log_yields, p)) # no breaks
print("Chevron: ")
print(ICSS_refinement(cvx_log_yields, p)) # one break
print(cvx_dates[2539]) # signals that volatility shift is 2020/02/20
plot(get_sigmas(cvx_log_yields), xlab = "Observations", ylab = "Chevron: Conditional variances before detecting breaks")
renewed <- c(get_sigmas(cvx_log_yields[1:2538]), get_sigmas(cvx_log_yields[2539:length(cvx_log_yields)]))
plot(renewed, xlab = "Observations", ylab = "Chevron: Conditional variances after detecting breaks")

# aerospace & defence
print("Boeing: ")
print(ICSS_refinement(ba_log_yields, p)) # one or three breaks
# 2533th observation is definitely a structural break
print(ba_dates[2533]) # signals that volatility shift is 2020/02/20
plot(get_sigmas(ba_log_yields), xlab = "Observations", ylab = "Boeing: Conditional variances before detecting breaks")
renewed <- c(get_sigmas(ba_log_yields[1:2532]), get_sigmas(ba_log_yields[2533:length(ba_log_yields)]))
plot(renewed, xlab = "Observations", ylab = "Boeing: Conditional variances after detecting breaks")

# consumers: food and beverages
print("Coca-cola: ")
print(ICSS_refinement(ko_log_yields, p)) # one break
print(ko_dates[2219]) # signals that volatility shift is 2018/11/16
plot(get_sigmas(ko_log_yields), xlab = "Observations", ylab = "Coca-cola: Conditional variances before detecting breaks")
renewed <- c(get_sigmas(ko_log_yields[1:2218]), get_sigmas(ko_log_yields[2219:length(ko_log_yields)]))
plot(renewed, xlab = "Observations", ylab = "Coca-cola: Conditional variances after detecting breaks")
print("McDonalds: ")
print(ICSS_refinement(mcd_log_yields, p)) # one break
print(mcd_dates[2533]) # signals that volatility shift is 2020/02/21
plot(get_sigmas(mcd_log_yields), xlab = "Observations", ylab = "McDonalds: Conditional variances before detecting breaks")
renewed <- c(get_sigmas(mcd_log_yields[1:2532]), get_sigmas(mcd_log_yields[2533:length(mcd_log_yields)]))
plot(renewed, xlab = "Observations", ylab = "McDonalds: Conditional variances after detecting breaks")

#For only last year
print("Google: ")
print(ICSS_refinement(googl_log_yields, p)) # no breaks
print("Netflix: ")
print(ICSS_refinement(nflx_log_yields, p)) # no breaks
print("Starbucks: ")
print(ICSS_refinement(sbux_log_yields, p)) # no breaks
print("Zoom: ")
print(ICSS_refinement(zm_log_yields, p)) # no breaks

print("Amazon: ")
print(ICSS_refinement(amzn_log_yields, p)) # one break
print(amzn_dates[193]) # signals that volatility shift is 2020/01/29
plot(get_sigmas(amzn_log_yields), xlab = "Observations", ylab = "Amazon: Conditional variances before detecting breaks")
renewed <- c(get_sigmas(amzn_log_yields[1:192]), get_sigmas(amzn_log_yields[193:length(amzn_log_yields)]))
plot(renewed, xlab = "Observations", ylab = "Amazon: Conditional variances after detecting breaks")

print("Facebook: ")
print(ICSS_refinement(fb_log_yields, p)) # two breaks
print(fb_dates[37]) # signals about volatility shift at 2019/06/14
print(fb_dates[208]) # signals that volatility shift is 2020/02/20
plot(get_sigmas(fb_log_yields), xlab = "Observations", ylab = "Facebook: Conditional variances before detecting breaks")
renewed <- c(get_sigmas(fb_log_yields[1:36]), get_sigmas(fb_log_yields[37:207]), get_sigmas(fb_log_yields[208:length(fb_log_yields)]))
plot(renewed, xlab = "Observations", ylab = "Facebook: Conditional variances after detecting breaks")
print(get_params(fb_log_yields[1:36])) # looks strange, but data is fine

print("Pepsi: ")
print(ICSS_refinement(pep_log_yields, p)) # three breaks most likely
print(pep_dates[59]) # signals about volatility shift at 2019/07/17
print(pep_dates[115]) # signals about volatility shift at 2019/10/04
print(pep_dates[211]) # signals that volatility shift is 2020/02/25
plot(get_sigmas(pep_log_yields), xlab = "Observations", ylab = "Pepsi: Conditional variances before detecting breaks")
renewed <- c(get_sigmas(pep_log_yields[1:58]), get_sigmas(pep_log_yields[59:114]), get_sigmas(pep_log_yields[115:210]), get_sigmas(pep_log_yields[211:length(pep_log_yields)]))
plot(renewed, xlab = "Observations", ylab = "Pepsi: Conditional variances after detecting breaks")
print(get_params(pep_log_yields[1:58]))
print(get_params(pep_log_yields[115:210]))

# Conclusion:

# plenty of observations show presence of structural break at February 20th
# (financial (JPM), energy (CVX), aerospace (BA), consumer-food (McD))
# on 0.01 significance level.

# Moreover, this remains even if sample size is small (~ 260 obs. for this one year):
# Amazon and Facebook shares show the presence of this structural break