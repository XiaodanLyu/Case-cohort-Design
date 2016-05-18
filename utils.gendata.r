size <- c (200, 300)
cssize <- list(c(47, 19), c(77, 29))
coef <- list (c (log(2), -0.5), c (log(2), 0), c (0, -0.5))
method <- c ("Exponential", "Weibull")
cr <- list (rbind (c (0.4924, 0.1298),
                   c (0.5196, 0.1155),
                   c (0.7151, 0.1948)),
            rbind (c (0.9180, 0.4501),
                   c (0.9310, 0.4705),
                   c (1.1002, 0.5470)))

gendata <- function (size, cssize, beta, method, censor.control, sample) {
  ## size: sample size, should be the same for both case-cohort and cox model
  allsize <- size + 100  #cohort size
  z <- cbind (rbinom (allsize, 1, 0.5), rnorm (allsize, 0, 1))
  exp.zbeta <- exp (as.vector (z %*% beta))
  
  if (method == "Exponential") {
    Ftime <- rexp (allsize, exp.zbeta)
  }
  if (method == "Weibull") {
    wshape = 2
    wscale <- exp.zbeta ^ (- 1 / wshape)
    Ftime <- rweibull (allsize, wshape, wscale)
  }
  
  Ctime <- runif (allsize, 0, censor.control)       # censoring time
  delta <- ifelse (Ftime < Ctime, 1, 0)   # censoring indicator
  Otime <- apply (cbind (Ftime, Ctime), 1, min)  # observed time
  
  if (sample == "Cox"){
    cox.index <- sample(1:allsize, size)
    data <- data.frame(time = Otime, delta = delta, z = z)[cox.index, ]
  }
  if (sample == "cc"){
    noncase.index <- which(delta == 0)
    case.index <- which(delta == 1)
    case.size <- length(case.index)
    subcohort.index <- c(sample(noncase.index, size - case.size), sample(case.index, cssize))
    xi <- rep(0, allsize)
    xi[subcohort.index] <- 1
    cc.index <- xi == 1 | delta == 1 
    data <- data.frame(time = Otime, delta = delta, xi = xi, z = z)[cc.index, ]
  }
  
  return(data)
}

censoring.test <- function(size, beta, method, censor.control) {
  nsim <- 500
  cc <- NULL
  for (i in 1:nsim)
    cc <- c(cc, 1-mean(gendata(size, beta, method, censor.control)$delta))
  mean(cc)
}

for (nsize in 1 : length (size))
  for (ncoef in 1 : length (coef))
    for (nmethod in 1:length(method))
      for (ncr in 1 : ncol (cr[[1]]))
        print(censoring.test (size [nsize], coef [[ncoef]], method [nmethod], cr [[nmethod]] [ncoef, ncr]))

