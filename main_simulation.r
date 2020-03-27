## developed by @betamint
source ("functions/cox-estimation.r")
source ("functions/cc-estimation.r")
source ("functions/utils.simulation.r")
source ("functions/utils.gendata.r")

options (width = 150)

nsim <- 1000
all_size <- c(200, 300)
coef <- list (c (log(2), -0.5), c (log(2), 0), c (0, -0.5))
init <- c(0, 1)
method <- c ("Exponential", "Weibull")
cr <- list (rbind (c (1.1125, 0.4924, 0.1298),
                   c (1.1161, 0.5196, 0.1155),
                   c (1.5665, 0.7151, 0.1948)),
            rbind (c (1.4923, 0.9180, 0.4501),
                   c (1.4750, 0.9310, 0.4705),
                   c (1.7626, 1.1002, 0.5470)))

set.seed (2015)											  
for (nsize in 1 : length (all_size))
for (ncoef in 1 : length (coef))
for (nmethod in 1)
for (ncr in 1 : ncol (cr[[1]]))
  both.simulation.wrap(nsim = nsim, all_size = all_size[nsize],
                       coef = coef [[ncoef]], init = init, method = method [nmethod],
                       censor.control = cr [[nmethod]] [ncoef, ncr])						 						  