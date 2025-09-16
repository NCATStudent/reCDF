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
require(RcppAlgos)

N <- 1e05
nB <- c(.01, .20) * N
nA <- 1000
seed <- 101
mod <- "f1"
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

# population characteristics
mu <- mean(pop$y)
F_N <- c(.01, .10, .25, .50, .75, .90, .99)
q <- quantile(pop$y, probs = F_N, type = 1)

# run for model f1
nBs <- nB[1]
ts <- q[1]
mX_bs <- 0 + 4 * x1 + 4 * x2 + 2 * x3 + 2 * x4
vX <- sqrt(3)
V2_const <- (1 / nBs) * (nA - 1) / (N * (N - 1)) * 1 / nA

# get variances 

## my beloved

V1 <- (1 - nA / N) / nA * var(pnorm((ts - mX_bs) / vX))

V2_2 <- -V2_const * sum(pnorm((ts - mX_bs) / vX))^2

## pain in the ass

raw_index <- expandGrid(i_it = 1:nA, h_it = 1:nA, nThreads = detectCores() - 1)

Rh_df <- A[raw_index["h_it"] %>% unlist(), ] %>%
  dplyr::select(-c(y, weight, Prob)) %>%
  mutate(Rh = qN_vals - predict(lm_B, newdata = .)) %>%
  mutate(h_it = raw_index["h_it"]) %>%
  dplyr::select(h_it, Rh)

Ri_df <- A[raw_index["i_it"] %>% unlist(), ] %>%
  dplyr::select(-c(y, weight, Prob)) %>%
  mutate(Ri = qN_vals - predict(lm_B, newdata = .)) %>%
  mutate(i_it = raw_index["i_it"]) %>%
  dplyr::select(i_it, Ri)

da_Ghats <- Rh_df %>%
  cbind(Ri_df) %>%
  mutate(da_min = pmin(Rh, Ri)) %>%
  dplyr::select(
    h_it, i_it, everything()
  ) %>%
  mutate(
    Ghat_min = e_ecdf(da_min)
  ) %>%
  dplyr::select(Ghat_min) %>%
  unlist()

v2_1 <- sum(da_Ghats) * (1 / (nB * nA^2))

# now the other one
v2_2 <- (1 / (nA^2 * nB)) * sum(A %>%
                                  dplyr::select(-c(weight, y, Prob)) %>%
                                  mutate(Rh = qN_vals - predict(lm_B, newdata = .)) %>%
                                  mutate(Ghat_h = e_ecdf(Rh)) %>%
                                  dplyr::select(Ghat_h) %>%
                                  unlist())^2
vars_to <- v1 + (v2_1 - v2_2)


pnorm(1.68)



lm_pop <- lm(y ~ ., data = pop %>% dplyr::select(
  paste0("x", 1:p), "y"
))

b_N1 <- lm_pop$coefficients # 1e+05
b_N2 <- lm_pop$coefficients # 1e+07
b_N3 <- lm_pop$coefficients # 1e+08


mlm_est_func <- function(qN_vals, nA) {
  m_lm <- (1 / sum(A$weight)) * sum(A$weight * e_ecdf(qN_vals - lm_pred_A))
  
  v1 <- (1 - nA / N) / nA * var(e_ecdf(qN_vals - lm_pred_A))
  
  # this part is gonna be computationally painful, fair warning :/
  raw_index <- expand.grid(i_it = 1:nA, h_it = 1:nA)
  
  Rh_df <- A[raw_index["h_it"] %>% unlist(), ] %>%
    dplyr::select(-c(y, weight, Prob)) %>%
    mutate(Rh = qN_vals - predict(lm_B, newdata = .)) %>%
    mutate(h_it = raw_index["h_it"]) %>%
    dplyr::select(h_it, Rh)
  
  Ri_df <- A[raw_index["i_it"] %>% unlist(), ] %>%
    dplyr::select(-c(y, weight, Prob)) %>%
    mutate(Ri = qN_vals - predict(lm_B, newdata = .)) %>%
    mutate(i_it = raw_index["i_it"]) %>%
    dplyr::select(i_it, Ri)
  
  da_Ghats <- Rh_df %>%
    cbind(Ri_df) %>%
    mutate(da_min = pmin(Rh, Ri)) %>%
    dplyr::select(
      h_it, i_it, everything()
    ) %>%
    mutate(
      Ghat_min = e_ecdf(da_min)
    ) %>%
    dplyr::select(Ghat_min) %>%
    unlist()
  
  v2_1 <- sum(da_Ghats) * (1 / (nB * nA^2))
  
  # now the other one
  v2_2 <- (1 / (nA^2 * nB)) * sum(A %>%
                                    dplyr::select(-c(weight, y, Prob)) %>%
                                    mutate(Rh = qN_vals - predict(lm_B, newdata = .)) %>%
                                    mutate(Ghat_h = e_ecdf(Rh)) %>%
                                    dplyr::select(Ghat_h) %>%
                                    unlist())^2
  vars_to <- v1 + (v2_1 - v2_2)
  
  return(data.frame(
    vars_to,
    m_lm
  ))
}