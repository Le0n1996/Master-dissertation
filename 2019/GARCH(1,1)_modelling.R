# GARCH(1,1) modelling

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