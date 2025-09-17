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

N <- 10000
nB <- c(.05, .20) * N
nA <- .05 * N
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

# get variances

## my beloved
pop_var_500 <- pop_var_2000 <- c()

for (i in 1:length(q)) {
  ts <- q[i]
  nBs <- 500

  mX_bs <- 0 + 4 * x1 + 4 * x2 + 2 * x3 + 2 * x4
  vX <- sqrt(3)
  V2_const <- (1 / nBs) * (nA - 1) / (N * (N - 1)) * 1 / nA


  V1 <- (1 - nA / N) / nA * var(pnorm((ts - mX_bs) / vX))

  V2_2 <- -V2_const * sum(pnorm((ts - mX_bs) / vX))^2

  ## pain in the ass

  raw_index <- expandGrid(
    ru_it = (ts - mX_bs) / vX,
    rv_it = (ts - mX_bs) / vX, nThreads = detectCores() - 1, return_df = TRUE
  ) %>%
    mutate(r_min = pmin(ru_it, rv_it)) %>%
    mutate(G_min = pnorm(r_min)) %>%
    dplyr::select(G_min) %>%
    sum()

  V2_1 <- V2_const * raw_index

  pop_var_500[i] <- V1 + V2_1 + V2_2

  print(i)
}

for (i in 1:length(q)) {
  ts <- q[i]
  nBs <- 2000

  mX_bs <- 0 + 4 * x1 + 4 * x2 + 2 * x3 + 2 * x4
  vX <- sqrt(3)
  V2_const <- (1 / nBs) * (nA - 1) / (N * (N - 1)) * 1 / nA


  V1 <- (1 - nA / N) / nA * var(pnorm((ts - mX_bs) / vX))

  V2_2 <- -V2_const * sum(pnorm((ts - mX_bs) / vX))^2

  ## pain in the ass

  raw_index <- expandGrid(
    ru_it = (ts - mX_bs) / vX,
    rv_it = (ts - mX_bs) / vX, nThreads = detectCores() - 1, return_df = TRUE
  ) %>%
    mutate(r_min = pmin(ru_it, rv_it)) %>%
    mutate(G_min = pnorm(r_min)) %>%
    dplyr::select(G_min) %>%
    sum()

  V2_1 <- V2_const * raw_index

  pop_var_2000[i] <- V1 + V2_1 + V2_2

  print(i)
}

pop_vars <-
  data.frame(
    perc = paste0(c(F_N, F_N) * 100, "%"),
    F_N = c(F_N, F_N),
    q = c(q, q),
    nB = rep(c(500, 2000),
      each = length(pop_var_500)
    ), pop_var = c(pop_var_500, pop_var_2000)
  ) # %>%
rownames_to_column("perc")

# grab MC variance

setwd("/Users/jeremyflood/Library/CloudStorage/OneDrive-Personal/Documents/Grad School/2024-2025/Fall 2025/reCDF/reCDF/Variance Estimation/Additional Requests/Confirming Variance/")
perf_df <- openxlsx::read.xlsx("modf1_results.xlsx", sheet = "summary")

cdf_pop_var = perf_df %>%
  dplyr::select(perc, nB, miss, est_type, perc, var_type, MC_var) %>%
  filter(
    est_type == "cdf",
    var_type == "asymp",
    miss == "MAR"
  ) %>%
  left_join(pop_vars, by = c("perc", "nB")) %>%
  dplyr::select(-c(var_type, miss, q)) %>%
  mutate(MC_var = signif(MC_var, digits = 5),
         pop_var = signif(pop_var, digits = 5),
         rb = round(100*(MC_var - pop_var)/pop_var, 2)) %>%
  arrange(nB) %>%
  dplyr::select(est_type, nB, perc, MC_var, pop_var, diff)

## now for q

ecdf_pop = ecdf(pnorm((y - mX_bs) / vX))

qN_dat = pop %>% dplyr::select(y) %>%
  mutate(mX = mX_bs) %>%
  mutate(ecdf_vals =purrr::map(y, function(x){mean(pnorm((x - mX) / vX))}) %>% as.numeric()) %>%  # why tf woud rowwise not work...?
  arrange(y) %>%
  ungroup()

search_func <- function(x, ...) {
  return(
qN_dat %>%
  filter(ecdf_vals >= x) %>%
  arrange(y) %>%
  dplyr::slice(1) %>%
  dplyr::select('y') %>%
  unlist() %>%
  unname() 
  )
}
  
pop_qvar_df = pop_vars %>%
  mutate(
    UL_crit = F_N + qnorm(1 - .10 / 2) * sqrt(pop_var),
    LL_crit = F_N - qnorm(1 - .10 / 2) * sqrt(pop_var)
  ) %>%
  mutate(LL = purrr::map(LL_crit, search_func) %>% as.numeric(),
         UL = purrr::map(UL_crit, search_func) %>% as.numeric()) %>%
  dplyr::select(-c(UL_crit, LL_crit)) %>%
  mutate(pop_q_var = ((UL-LL)/(2*qnorm(1-.10/2)))^2) %>%
  dplyr::select(
    perc,
    F_N, 
    q, 
    nB, pop_q_var
  )


perf_df %>%
  dplyr::select(perc, nB, miss, est_type, perc, var_type, MC_var) %>%
  filter(
    est_type == "t",
    var_type == "asymp",
    miss == "MAR"
  ) %>%
  left_join(pop_qvar_df, by = c("perc", "nB")) %>%
  dplyr::select(-c( var_type, miss, q)) %>%
  mutate(diff = round(MC_var - pop_q_var, 3)) %>%
  arrange(nB) %>%
  dplyr::select(est_type, nB, perc, MC_var, pop_q_var, diff)

Rh_df <- pop[raw_index["h_it"] %>% unlist(), ] %>%
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