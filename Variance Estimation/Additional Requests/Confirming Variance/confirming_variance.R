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

# get data

setwd('/Users/jeremyflood/Library/CloudStorage/OneDrive-Personal/Documents/Grad School/2024-2025/Fall 2025/reCDF/reCDF/Variance Estimation/Additional Requests/Confirming Variance/smaller N/Data/Cleaned Results')
total_dat <- c()

for (i in 1:length(list.files())) {
  total_dat[[i]] <-
      openxlsx::read.xlsx(list.files()[i], sheet = "summary")
    
  print(paste0("completed ", i, " out of ", length(list.files())))
}

MC_var_results <- total_dat %>%
  bind_rows() %>%
  mutate(N = nB / .20) %>%
  filter(var_type == "asymp", miss == "MAR") %>%
  mutate(group = ifelse(est_type == 't', 'q', 'cdf')) %>%
  arrange(N) %>%
  dplyr::select(
    perc, group, MC_var, nB, N
  )

N_vect <- MC_var_results$N %>%
  table() %>%
  names()
search_func <- function(x, ...) {
  return(
    qN_dat %>%
      filter(ecdf_vals >= x) %>%
      arrange(y) %>%
      dplyr::slice(1) %>%
      dplyr::select("y") %>%
      unlist() %>%
      unname()
  )
}


# now we need to get results from variance test
final_df <- c()
for (j in 1:length(N_vect)) {
  N <- as.numeric(N_vect[j])
  nB <- .20 * N
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
  nBs <- nB
  # get variances

  ## my beloved
  pop_var <- c()

  for (i in 1:length(q)) {
    ts <- q[i]

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

    pop_var[i] <- V1 + V2_1 + V2_2

    print(i)
  }

  pop_vars_cdf <-
    data.frame(
      perc = paste0(F_N*100, '%'),
      F_N = F_N,
      q = q,
      nB = nB,
      N = N,
      pop_var = c(pop_var),
      group = "cdf"
    )

  # now for q
  ecdf_pop <- ecdf(pnorm((y - mX_bs) / vX))

  qN_dat <- pop %>%
    dplyr::select(y) %>%
    mutate(mX = mX_bs) %>%
    mutate(ecdf_vals = purrr::map(y, function(x) {
      mean(pnorm((x - mX) / vX))
    }) %>% as.numeric()) %>% # why tf woud rowwise not work...?
    arrange(y) %>%
    ungroup()


  pop_vars_q <- pop_vars_cdf %>%
    mutate(
      UL_crit = F_N + qnorm(1 - .10 / 2) * sqrt(pop_var),
      LL_crit = F_N - qnorm(1 - .10 / 2) * sqrt(pop_var)
    ) %>%
    mutate(
      LL = purrr::map(LL_crit, search_func) %>% as.numeric(),
      UL = purrr::map(UL_crit, search_func) %>% as.numeric()
    ) %>%
    dplyr::select(-c(UL_crit, LL_crit)) %>%
    mutate(pop_var = ((UL - LL) / (2 * qnorm(1 - .10 / 2)))^2, group = "q", N = N) %>%
    dplyr::select(
      perc,
      F_N,
      q,
      nB, pop_var, N, group
    )

  final_df[[j]] <- rbind(pop_vars_q, pop_vars_cdf)
  print(paste0("completed ", j, " out of ", length(N_vect)))
}

actual_var_df = final_df %>% 
  bind_rows() %>%
  rownames_to_column() %>%
  dplyr::select(-rowname) %>%
  dplyr::select(nB, N, group, perc, everything()) 

actual_var_df %>% 
  left_join(MC_var_results, by = c('nB', 'N', 'group', 'perc')) %>%
  mutate(abs_diff = abs(pop_var - MC_var), nB= factor(nB)) %>%
  mutate(abs_diff = ifelse(group == 'q' & perc == '99%', NA, abs_diff)) %>%
  filter(nB %in% c(
   # 200,
    400,
    600,
    1000,
    1600,
    2000
  )) %>%
  ggplot(aes(x = perc, y = abs_diff, color = nB, group = nB)) +
  geom_line(linewidth = .5) +
  geom_point(size = 1) +
  facet_grid(
    scales = "free", rows = dplyr::vars(group),
    labeller = label_parsed
  ) #+
