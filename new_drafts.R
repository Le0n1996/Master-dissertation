#Step 1: take data and estimate coefficients for GARCH(1,1)
{
  SBERtr <- read_excel("C:/Users/novle/Documents/Programming/R/coursework/shares 2015-2020.xlsx", col_names = FALSE, range = "Sheet1!D3:D1316")
  SBERdates <- read_excel("C:/Users/novle/Documents/Programming/R/coursework/shares 2015-2020.xlsx", col_names = FALSE, range = "Sheet1!B3:B1316")
  SBERtr <- as.matrix(SBERtr)
  storage.mode(SBERtr) <- "numeric" 
  SBERtr <- as.vector(SBERtr)
  SBERdates <- as.matrix(SBERdates)
  storage.mode(SBERdates) <- "character"
  SBERdates <- as.vector(SBERdates)  
  
  Y <- read.delim("C:/Users/novle/Documents/Programming/R/coursework/GARCH_with_breaks.csv",
                  sep = ";",  dec = ",", colClasses=c("numeric", "numeric"))
  
  p1 <- Y$X1break[!is.na(Y$X1break)]
  #p2 <- na.omit(Y$X2breaks)
  dats <- seq.Date(from = as.Date("2012-05-11"),
                   by = "day",
                   length.out = 1000) #length(p1)
  #Y <- read_excel("C:/Users/novle/Documents/Programming/R/coursework/N2.xlsx", col_names = FALSE, range = "Лист1!A2:A1001")
  
  p1 <- SBERtr
  dats <- SBERdates
  
  data1 <- data.frame(p1)
  colnames(data1) <- c("Log.yield")
  data1 <- slice(data1, 1:1000)
  rownames(data1) <- dats
  colnames(data1) <- c("Log.yield") # data is ready
  
  y <- data1$Log.yield
  t <- ymd(dats)
  logyields <- xts(y, order.by=t)
  tsdisplay(logyields)
  
  model <- ugarchspec(
    variance.model = list(garchOrder = c(1, 1)), 
    mean.model = list(armaOrder = c(0, 0)))
  model_est <- ugarchfit(spec = model,
                         data=logyields)
  
  #get coefficients coef(model_est)[1]
  coefs <- as.vector(coef(model_est))
}
coefs

save.image(file = "my_work_space_06.RData")
load("my_work_space_23_05.RData")

Get_GARCH_params_another_draft <- function(TS, tm=seq.Date(from = as.Date("2000-01-01"), by = "day", length.out = length(TS))) {
  #Y <- as.matrix(TS)
  #storage.mode(Y) <- "numeric" 
  #Y <- as.vector(Y)
  #dats <- as.matrix(tm)
  #storage.mode(dats) <- "character"
  #dats <- as.vector(dats)
  
  #Y <- read.delim("C:/Users/novle/Documents/Programming/R/coursework/GARCH_with_breaks.csv",
  #                sep = ";",  dec = ",", colClasses=c("numeric", "numeric"))
  
  #p1 <- Y$X1break[!is.na(Y$X1break)]
  #p2 <- na.omit(Y$X2breaks)
  #dats <- seq.Date(from = as.Date("2012-05-11"),
  #                 by = "day",
  #                 length.out = 1000) #length(p1)
  #data1 <- slice(df, 1:1000)
  
  df <- data.frame(Y)
  rownames(df) <- dats
  colnames(df) <- c("Log.yield") # data is ready
  
  values <- df$Log.yield
  t <- ymd(dats)
  log_yields <- xts(values, order.by=t)
  tsdisplay(log_yields)
  
  model <- ugarchspec(
    variance.model = list(garchOrder = c(1, 1)), 
    mean.model = list(armaOrder = c(0, 0)))
  model_est <- ugarchfit(spec = model,
                         data=log_yields)
  
  #get coefficients [1] - mu, [2] - omega, [3] - alpha1, [4] - beta1
  coefs <- as.vector(coef(model_est))
}

Innov_series <- function(Y, coefs) {
  observations <- as.vector(Y)
  sq <- observations[,1]*observations[,1]
  sigma_start <- (sum(sq)/length(sq) - coefs[1]^2) #initialization with sample variance
  sigmas <- c(sigma_start)
  sigma_start
  
  for (i in 1:nrow(Y)){#
    new_sigma <- coefs[2] + coefs[3]*sq[i] + coefs[4]*sigmas[length(sigmas)]
    sigmas <- c(sigmas, new_sigma)
  }
  innov <- observations/sigmas[2:(nrow(Y)+1)]
  return(innov)
}

tsboot(a, stat, R=1, l=5, sim="fixed", orig.t=FALSE)$t[1,]

options(max.print=100000)

start_time <- Sys.time()

end_time <- Sys.time()
end_time - start_time

Innov_series <- function(Y, coefs) {
  Y <- Preprocess_data(Y)
  obs <- as.vector(Y)
  sq <- obs*obs
  sigma_start <- (sum(sq)/length(sq) - coefs[1]^2) #initialization with sample variance
  sigmas <- c(sigma_start)
  
  for (i in 1:length(obs)) {
    new_sigma <- coefs[2] + coefs[3]*sq[i] + 
      coefs[4]*sigmas[length(sigmas)] #i
    sigmas <- c(sigmas, new_sigma)
  }
  innov <- obs/sigmas[2:(length(obs)+1)]
  return(innov)
}

typeof(Y)

est_model <- garchFit(formula = ~ garch(1, 1), data = Y1,
                      include.mean = FALSE, trace = FALSE,
                      stepwise=FALSE, approx=FALSE)
params <- est_model@fit$par[1:3] # We are getting omega, alpha1 and beta1

tryCatch(
  {
    quak = b/c
  },
  error=function(cond) {
    message(paste("ZERO!", url))
    # Choose a return value in case of error
    return(NA)
  },
  warning=function(cond) {
    message(paste("URL caused a warning:", url))
    message("Here's the original warning message:")
    message(cond)
    # Choose a return value in case of warning
    return(NULL)
  },
  finally={
    message("Some other message at the end")
  }
)  

get_params <- function(Y){
  tryCatch(
    expr = {
      suppressWarnings(est_model <- garchFit(formula = ~ garch(1, 1), data = Y,
                                             include.mean = FALSE, trace = FALSE))
      suppressWarnings(params <- est_model@fit$par[1:3]) # We are getting omega, alpha1 and beta1
      return(params) # comment: warnings are suppressed, because garchFit sometimes have problems with std errors
    },
    error = function(e){
      suppressWarnings(est_model <- garchFit(formula = ~ garch(1, 1), data = Y,
                                             include.mean = FALSE, trace = FALSE, algorithm = "lbfgsb+nm")) # another method for hessian crush case
      suppressWarnings(params <- est_model@fit$par[1:3]) # We are getting omega, alpha1 and beta1
      return(params) # comment: warnings are suppressed, because garchFit sometimes have problems with std errors
    },
    warning = function(w){
      message('Another one!')
    },
    finally = {
    }
  )
}
get_sigmas <- function(Y){
  tryCatch(
    expr = {
      suppressWarnings(est_model <- garchFit(formula = ~ garch(1, 1), data = Y,
                                             include.mean = FALSE, trace = FALSE))
      suppressWarnings(sigmas <- est_model@sigma.t) # We are getting sigmas series
      return(sigmas) # comment: warnings are suppressed, because garchFit sometimes have problems with std errors
    },
    error = function(Y){
      message('However, we are trying!')
      suppressWarnings(est_model <- garchFit(formula = ~ garch(1, 1), data = Y,
                                             include.mean = FALSE, trace = FALSE, algorithm = "lbfgsb+nm")) # another method for hessian crush case
      suppressWarnings(sigmas <- est_model@sigma.t) # We are getting sigmas series
      return(sigmas) # comment: warnings are suppressed, because garchFit sometimes have problems with std errors
    },
    warning = function(w){
      message('Another one!')
    },
    finally = {
    }
  )
}

# returns asymptotic quantiles for supremum of Brownian bridge's absolute value
BB_crit <- function(p = 0.05) {
  if (p == 0.01) { # significance level 0.01
    return(1.627)
  } else if (p == 0.05) { # significance level 0.05
    return(1.358)
  } else if (p == 0.10) { # significance level 0.1
    return(1.224)
  } else if (p == 0.25) { # significance level 0.25
    return(1.019)
  } else if (p == 0.50) { # significance level 0.5
    return(0.828)
  } else if (p == 0.75) { # significance level 0.75
    return(0.677)
  } else if (p == 0.90) { # significance level 0.9
    return(0.571)
  } else if (p == 0.95) { # significance level 0.95
    return(0.520)
  } #p = 1 + 2*sum((-1)^k*exp(-2*k*k*x*x))
}
# returns suspicious time moment, value of statistic, and decision
BB_test <- function(Y, uncorrelated=TRUE, p = 0.05) {
  BB <- BB_stat(Y, uncorrelated=TRUE)
  if (BB[2] > BB_crit(p)){ # check whether statistic is larger than p-value or not
    return(c('Break', BB))}
  else
    return(c('No break', BB))
}