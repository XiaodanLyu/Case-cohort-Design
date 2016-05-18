## utilities
aat <- function (a) {
  a %*% t (a)
}
tailsum <- function (a) {
  sum (a) + a - cumsum (a)
}
risksum <- function (A) {
  apply (A, 2, tailsum)
}
var.sample <- function (sample.matrix) {
  sample.dim <- ncol (sample.matrix)
  m <- sample.matrix
  m1 <- t (t (m) - apply (m, 2, mean))
  m2 <- t (apply (m1, 1, aat))
  m3 <- apply (m2, 2, mean)
  v <- matrix (m3, sample.dim, sample.dim)
  v
}


## cox proportional hazard regression coef estimation
#' data = data.frame (time, delta, covariate)
coxsum <- function (data, coef) {
  data <- as.matrix (data [sort.list (data $ time), ])
  d <- as.vector (data [ , 2])
  z <- as.matrix (data [ , 3 : ncol (data)])
  b <- coef
  
  ezb <- exp (as.vector (z %*% b))
  zezb <- ezb * z
  z2 <- t (apply (z, 1, aat))
  z2ezb <- ezb * z2
  
  S0 <- tailsum (ezb)
  S1 <- risksum (zezb)
  S2 <- risksum (z2ezb)
  
  S1S0 <- S1 / S0
  S2S0 <- S2 / S0
  S1S02 <- t (apply (S1S0, 1, aat))
  
  Gradient <- apply ((z - S1S0) * d, 2, sum)
  Hessian <- matrix (apply ((S2S0 - S1S02) * d, 2, sum),
                     length (b))
  
  list (S1 = S1, S2 = S2, S1S0 = S1S0, S2S0 = S2S0, S1S02 = S1S02, Gradient = Gradient, Hessian = Hessian)
}
coxest.nr <- function (data, init, tolerance = 0.0001, maxiter = 200) {
  goodenough <- function (v1, v2) { sum (abs (v1 - v2)) <= tolerance }
  iterate <- function (guess, singular, success, iter) {
    while (singular == FALSE & success == FALSE & iter > 0) {
      SS <- coxsum (data, guess)
      f <- SS $ Gradient
      d <- SS $ Hessian
      
      if (is.na(rcond(d)) | rcond (d) <= .Machine $ double.eps) {
        singular <- TRUE
      } else {
        new <- guess + solve (d, f)
        if (goodenough (guess, new)) {
          success <- TRUE
        }
        guess <- new
      }
      iter <- iter - 1
    }
    list (guess = guess, singular = singular, success = success)
  }
  iterate (init, singular = FALSE, success = FALSE, maxiter)
}
coxseest <- function (data, est) {
  size <- nrow (data)
  d <- data [ , 2]
  z <- data [ , 3 : ncol (data)]
  
  SS <- coxsum (data, est)
  S1S0 <- SS $ S1S0
  
  H <- SS $ Hessian
  H.inv <- solve (H / size)
  
  var.matrix <- var.sample ((z - S1S0) * d)
  var <- diag (H.inv %*% var.matrix %*% H.inv) / size
  se <- sqrt (var)
  se
}
cox.bootstrap <- function(sample.data, init, B=200)
{
  p <- length(init)
  nn <- nrow(sample.data)
  betaboot <- matrix(0, B, p)
  ind <- logical(B)
  for (k in 1:B)
  {
    index <- sample(nn, replace = T)
    boot.data <- sample.data[index, ]
    best <- coxest.nr(boot.data, init)
    betaboot[k, ] <- best$guess
    ind[k] <- !best$singular & best$success
  }
  apply(betaboot[ind, ], 2, sd)
}
cox.estimation <- function (data, init) {
  est <- coxest.nr (data, init)    # list (guess, singular, success)
  if (est $ singular) {
    se <- rep (NA, length (init))
  } else {
    se <- rbind(coxseest (data, est$guess), cox.bootstrap(data, init))
  }
  list (est = est $ guess, singular = est $ singular, success = est $ success, se = se)
}
