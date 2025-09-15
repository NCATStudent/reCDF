# Required Packages
require(tidyverse)
require(tidyr)
require(ggpmisc)
require(svMisc)
require(CBPS)
require(geiger)
require(lsa)
require(DescTools)
require(MVTests)
require(caret)
require(dbarts)
require(parallel)
require(gbm)
require(pbmcapply)
require(extraDistr)
require(plyr)
require(data.table)
RNGkind("L'Ecuyer-CMRG")
require(schoenberg)
require(lsr)
require(car)
require(psych)
require(quickmatch)
require(devtools)
require(Matching)
require(rgenoud)
require(optmatch)
require(Hmisc)
require(rio)
require(sjstats)
require(schoolmath)
require(MatchIt)
require(magrittr)

require(writexl)
require(ncvreg)
require(np)
require(dplyr)
require(tidyverse)
require(survey)
require(sampling)
require(KernSmooth)
require(caTools)
require(gridExtra)
require(lars)
require(glmnet)
require(Metrics)
Sample <- Vectorize(sample)
Dat <- function(x, p, size, seed) {
  set.seed(seed)
  return(x[x$part == p, ][sample(nrow(x[x$part == p, ]), size), ])
}
Dat.v <- Vectorize(Dat, vectorize.args = c("size", "p"))
require(MASS)
require(mgcv)
require(parallel)
require(ranger)
require(xgboost)
library(dplyr)
library(mgcv)
require(haven)
`%notin%` <- Negate(`%in%`)

require(survey)
require(svrep)
require(svMisc)
require(doParallel)
require(splines)
require(future.apply)
require(latex2exp)
require(ggplot2)

N = 1e05
nB = c(.01, .20)*N
mod = 'f3'
if (mod == "f1") {
  set.seed(seed)
  x1 <- rnorm(N, mean = 2, sd = 1)
  x2 <- rnorm(N, mean = 2, sd = 1)
  x3 <- rnorm(N, mean = 4, sd = 1)
  x4 <- rnorm(N, mean = 4, sd = 1)
  X <- cbind(x1, x2, x3, x4)
  p <- ncol(X)
  y <- rnorm(N, 4 * x1 + 4 * x2 + 2 * x3 + 2 * x4, 3)
} else if (mod == "f2") {
  set.seed(seed)
  x1 <- runif(N, min = 0, max = 4)
  x2 <- runif(N, min = 0, max = 4)
  x3 <- runif(N, min = 4, max = 8)
  x4 <- runif(N, min = 4, max = 8)
  X <- cbind(x1, x2, x3, x4)
  p <- ncol(X)
  y <- rnorm(N, 4 * x1^2 + 4 * x2^2 + 2 * x3^2 + 2 * x4^2 + (x1 + x2)^2 + (x3 + x4)^2, sd = 50)
} else if (mod == "f3") {
  set.seed(seed)
  x1 <- runif(N, min = -1, max = 1)
  x2 <- runif(N, min = -1, max = 1)
  x3 <- runif(N, min = -1, max = 1)
  x4 <- runif(N, min = -1, max = 1)
  X <- data.frame(x1, x2, x3, x4)
  p <- ncol(X)
  err <- rnorm(N, mean = 0, sd = sqrt(.5))
  y <- -sin(x1) + x2^2 + x3 - exp(-x4^2) + err
} else if (mod == "f4") {
  set.seed(seed)
  x1 <- rnorm(N)
  x2 <- rnorm(N)
  x3 <- rnorm(N)
  x4 <- rnorm(N)
  x5 <- rnorm(N)
  x6 <- rnorm(N)
  X <- data.frame(x1, x2, x3, x4, x5, x6)
  p <- ncol(X)
  err <- rnorm(N)
  y <- x1 + .707 * x2^2 + 2 * ifelse(x3 > 0, 1, 0) + .873 * log(abs(x1)) * abs(x3) +
    .894 * x2 * x4 + 2 * ifelse(x5 > 0, 1, 0) + .46 * exp(x6) + err
} else {
  stop("You must select one of the following models: (f1, f2, f3, f4)")
}


pop <- as.data.frame(cbind(X, y))


lm_pop = lm(y~., data = pop %>% dplyr::select(
  paste0('x', 1:p), 'y'
))

b_N1 = lm_pop$coefficients
b_N2 = lm_pop$coefficients
b_N3 =  
# population characteristics
mu <- mean(pop$y)
F_N <- c(.01, .10, .25, .50, .75, .90, .99)
q <- quantile(pop$y, probs = F_N, type = 1)

