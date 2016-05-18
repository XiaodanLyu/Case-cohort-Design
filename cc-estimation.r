## utilities
aat <- function (a) {
    a %*% t (a)
}
tailsum <- function (a) {
    sum (a) + a - cumsum (a)
}
tailsum.cc <- function (a, d, x) {
    tailsum (a * x) + a * d * (1 - x)
}
risksum <- function (A) {
    apply (A, 2, tailsum)
}
risksum.cc <- function (A, d, x) {
    apply (A * x, 2, tailsum) + A * d * (1 - x)
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


## case-cohort coef estimation
#' data = data.frame (time, delta, xi, covariate)
ccsum <- function (data, coef) {
    data <- as.matrix (data [sort.list (data $ time), ])
    d <- as.vector (data [ , 2])
    x <- as.vector (data [ , 3])
    z <- as.matrix (data [ , 4 : ncol (data)])
    b <- coef

    ezb <- exp (as.vector (z %*% b))
    zezb <- ezb * z
    z2 <- t (apply (z, 1, aat))
    z2ezb <- ezb * z2

    S0 <- tailsum.cc (ezb, d, x)
    S1 <- risksum.cc (zezb, d, x)
    S2 <- risksum.cc (z2ezb, d, x)

    S1S0 <- S1 / S0
    S2S0 <- S2 / S0
    S1S02 <- t (apply (S1S0, 1, aat))

    Gradient <- apply ((z - S1S0) * d, 2, sum)
    Hessian <- matrix (apply ((S2S0 - S1S02) * d, 2, sum),
                       length (b))

    list (S1 = S1, S2 = S2, S1S0 = S1S0, S2S0 = S2S0, S1S02 = S1S02, Gradient = Gradient, Hessian = Hessian)
}
ccest.nr <- function (data, init, tolerance = 0.0001, maxiter = 200) {
    goodenough <- function (v1, v2) { sum (abs (v1 - v2)) <= tolerance }
    iterate <- function (guess, singular, success, iter) {
        while (singular == FALSE & success == FALSE & iter > 0) {
            SS <- ccsum (data, guess)
            f <- SS $ Gradient
            d <- SS $ Hessian

            if (is.na(rcond(d)) | rcond (d) <= .Machine $ double.eps) {
                singular <- TRUE
            } else {
                new <- guess + solve (d, f)
                success <- goodenough (guess, new)
                guess <- new
            }
            iter <- iter - 1
        }
        list (guess = guess, singular = singular, success = success)
    }
    iterate (init, singular = FALSE, success = FALSE, maxiter)
}
ccseest <- function (data, est) {
    size <- nrow (data)
    d <- data [ , 2]
    z <- data [ , 4 : ncol (data)]

    SS <- ccsum (data, est)
    S1S0 <- SS $ S1S0
    H <- SS $ Hessian
    H.inv <- solve (H / size)

    var.matrix <- var.sample ((z - S1S0) * d)
    var <- diag (H.inv %*% var.matrix %*% H.inv) / size
    se <- sqrt (var)
    se
}
cc.bootstrap <- function(sample.data, init, B=200)
{
  p <- length(init)
  nn <- nrow(sample.data)
  betaboot <- matrix(0, B, p)
  ind <- logical(B)
  for (k in 1:B)
  {
    index <- sample(nn, replace = T)
    boot.data <- sample.data[index, ]
    best <- ccest.nr(boot.data, init)
    betaboot[k, ] <- best$guess
    ind[k] <- !best$singular & best$success
  }
  apply(betaboot[ind, ], 2, sd)
}
cc.estimation <- function (data, init) {
    est <- ccest.nr (data, init)    # list (guess, singular, success)
    if (est $ singular) {
        se <- NA
    } else {
        se <- rbind(ccseest (data, est$guess), cc.bootstrap(data, init))
    }
    list (est = est $ guess, singular = est $ singular, success = est $ success, se = se)
}
