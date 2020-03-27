gendata <- function (size, beta, method, censor.control, sample) {
  ## size: sample size, should be the same for both case-cohort and cox model
  cosize <- size + 100  #cohort size
  z <- cbind (rbinom (cosize, 1, 0.5), rnorm (cosize, 0, 1))
  exp.zbeta <- exp (as.vector (z %*% beta))
  
  if (method == "Exponential") {
    Ftime <- rexp (cosize, exp.zbeta)
  }
  if (method == "Weibull") {
    wshape = 2
    wscale <- exp.zbeta ^ (- 1 / wshape)
    Ftime <- rweibull (cosize, wshape, wscale)
  }
  
  Ctime <- runif (cosize, 0, censor.control)       # censoring time
  delta <- ifelse (Ftime < Ctime, 1, 0)   # censoring indicator
  Otime <- apply (cbind (Ftime, Ctime), 1, min)  # observed time
  
  if (sample == "Cox"){
    cox.index <- sample(1:cosize, size)
    data <- data.frame(time = Otime, delta = delta, z = z)[cox.index, ]
  }
  if (sample == "cc"){
    noncase.index <- which(delta == 0)
    case.index <- which(delta == 1)
    case.size <- length(case.index)
    cssize <- (size - case.size)* mean(delta) / (1 - mean(delta)) ## same censoring rate for subcohort
    subcohort.index <- c(sample(noncase.index, size - case.size), sample(case.index, round(cssize)))
    xi <- rep(0, cosize)
    xi[subcohort.index] <- 1
    cc.index <- xi == 1 | delta == 1 
    data <- data.frame(time = Otime, delta = delta, xi = xi, z = z)[cc.index, ]
  }
  
  return(data)
}

CP.indicator <- function (real, est, SE.hat) {
  real.left <- real - SE.hat * 1.96
  real.right <- real + SE.hat * 1.96
  (est >= real.left) & (est <= real.right)
}

both.simulation <- function (nsim, all_size, coef, init, method, censor.control) {
  np <- length (coef)
  res.cc <- data.frame (matrix (0, nsim, 1 + np + 1 + 1 + np + np))
  res.cox <- data.frame (matrix (0, nsim, 1 + np + 1 + 1 + np + np))
  col_censor.rate <- 1
  col_est         <- 2 : (1 + np)
  col_singular    <- 2 + np
  col_success     <- 3 + np
  col_se          <- (4 + np) : (3 + np * 2)
  col_cp          <- (4 + np * 2) : (3 + np * 3)
  names (res.cc) <- c ("censor.rate", rep ("est", np), "singular", "success", rep ("se", np), rep ("cp", np))
  names (res.cox) <- c ("censor.rate", rep ("est", np), "singular", "success", rep ("se", np), rep ("cp", np))
  
  for (i in 1 : nsim) {
    cc.data <- gendata(all_size, coef, method, censor.control, "cc")    # data.frame (time, delta, xi, z)
    res.cc[i, col_censor.rate] <- 1 - mean(cc.data[cc.data$xi == 1, ]$delta)
    est <- cc.estimation (cc.data, init)
    res.cc [i, c(col_est,col_singular,col_success,col_se)] <- c (est $ est, est $ singular, est $ success, est $ se)
    res.cc [i, col_cp] <- CP.indicator (coef, est $ est, est $ se)
    
    cox.data <- gendata (all_size, coef, method, censor.control, "Cox")    # data.frame (time, delta, z)
    res.cox [i, col_censor.rate] <- 1 - mean (cox.data $ delta)
    est <- cox.estimation (cox.data, init)
    res.cox [i, c(col_est,col_singular,col_success,col_se)] <- c (est $ est, est $ singular, est $ success, est $ se)
    res.cox [i, col_cp] <- CP.indicator (coef, est $ est, est $ se)
  }
  
  mean.censor.rate <- list (cc = mean (res.cc $ censor.rate), cox = mean (res.cox $ censor.rate))
  nsingular <- list (cc = sum (res.cc $ singular), cox = sum (res.cox $ singular))
  
  valid_index.cc <- (res.cc $ singular == FALSE) & (res.cc $ success == TRUE)
  valid_index.cox <- (res.cox $ singular == FALSE) & (res.cox $ success == TRUE)
  res2.cc <- res.cc [valid_index.cc, ]
  res2.cox <- res.cox [valid_index.cox, ]
  
  summary.cc = data.frame (coef = coef,
                           bias = apply (res2.cc [ , col_est], 2, mean) - coef,
                           sd = apply (res2.cc [ , col_est], 2, sd),
                           se = apply (res2.cc [ , col_se], 2, mean),
                           cp = apply (res2.cc [ , col_cp], 2, mean))
  summary.cox = data.frame (coef = coef,
                            bias = apply (res2.cox [ , col_est], 2, mean) - coef,
                            sd = apply (res2.cox [ , col_est], 2, sd),
                            se = apply (res2.cox [ , col_se], 2, mean),
                            cp = apply (res2.cox [ , col_cp], 2, mean))
  
  #  list (results = list (res.cc, res.cox),
  list (summary = list (cc = summary.cc, cox = summary.cox),
        mean.censor.rate = mean.censor.rate,
        nsingular = nsingular)
}

both.simulation.wrap <- function (nsim, all_size, coef, init, method, censor.control) {
  t1 <- proc.time ()
  
  variable <- paste (all_size, paste (round(coef,1), collapse="_"), method, censor.control, sep = "_")
  
  res <- both.simulation (nsim, all_size, coef, init, method, censor.control)
  save (res, file = paste ("results/", variable, ".Rdata", sep = ""))
  
  t2 <- proc.time()
  t2 - t1
}
