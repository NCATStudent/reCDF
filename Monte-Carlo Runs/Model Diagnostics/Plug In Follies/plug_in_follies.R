# investigating why median and 75th for plug in

# requirements
require(tidyverse)
require(sampling)
require(openxlsx)
require(magrittr)
require(mcreplicate)
require(parallel)

# reading data
setwd("/Users/jeremyflood/Library/CloudStorage/OneDrive-Personal/Documents/Grad School/2024-2025/Fall 2025/reCDF/reCDF/Monte-Carlo Runs/r_15")

mods <- paste0("modf", 1:4, "_results.xlsx")
all_data <- c()
for (i in 1:length(mods)) {
  all_data[[i]] <- read.xlsx(mods[i], sheet = "perf")
}

all_data %<>% bind_rows()

# why lm_plug does so well at 50% under MAR
## confirm issue
all_data %>%
  filter(
    group == "cdf",
    est %in% c("lm_plug", "m_lm"),
    miss == "MNAR",
    perc == "75%"
  ) %>%
  arrange(nB) %>%
  dplyr::select(-c(rvar, r, group, miss)) %>%
  pivot_wider(names_from = "est", values_from = "rrmse")
## replicate
### q1. how good are the predicted mean / medians
samp_func <- function(
    mod,
    seed = 101,
    N = 1e05,
    nB = 20000,
    nA = 1000,
    miss,
    r, option = 1) {
  # get pop
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
  F_N <-c(.50, seq(0.001, .99, by =.001)) %>% sort()
  q <- quantile(pop$y, probs = F_N, type = 1)
  mu <- mean(pop$y)
  med <- median(pop$y)

  # get A

  A <- pop[sample(1:nA), ] %>%
    mutate(prob = nA / N, weight = N / nA)

  # get B

  B_generator <- function(miss, nB, r) {
    sizes_rounded <- round(c(r * nB, (1 - r) * nB),
      digits = 0
    )

    sizes <- ifelse(sizes_rounded < 1, 1, sizes_rounded)
    if (miss == "MAR") {
      strat <- paste0("x", c(abs(cor(X, y)) %>% which.max()))
      pop_restrat <- pop %>%
        mutate(restrat = ifelse(get(strat) <= median(get(strat)), 1,
          2
        )) %>%
        arrange(restrat)
      s.B <- sampling::strata(pop_restrat,
        stratanames = "restrat",
        size = sizes,
        method = "srswor"
      )
      B <- getdata(pop_restrat, s.B)
    } else if (miss == "MNAR") {
      pop_restrat <- pop %>%
        mutate(mnar_strat = ifelse(y < median(y), "0", "1") %>% as.factor()) %>%
        arrange(mnar_strat)
      s.B <- sampling::strata(pop_restrat,
        stratanames = "mnar_strat",
        size = sizes,
        method = "srswor"
      )
      B <- getdata(pop_restrat, s.B)
    }
    return(B)
  }
  B <- B_generator(miss, nB, r)

  # build regression model
  lm_mod <- lm(y ~ ., data = B %>% dplyr::select(
    "y", paste0("x", 1:p)
  ))

  lm_pred_A <- predict(lm_mod, newdata = A)

  # compare
  final_df <- data.frame(
    "mu" = mu,
    "hat_mu" = mean(lm_pred_A),
    "med" = med,
    "hat_med" = median(lm_pred_A),
    "mad" = mad(pop$y),
    "mad_hat" = mad(lm_pred_A, constant = 1),
    "mad_hat_pop_med" = mad(lm_pred_A, center = median(pop$y), constant = 1),
    "sd" = sd(pop$y),
    "sd_hat" = sd(lm_pred_A),
    "miss" = miss,
    r = r,
    mod = mod
  )
  if (option == 2) {
    output <- data.frame(F_N = F_N, y_N = q, lm_pred_A = quantile(lm_pred_A, F_N), "miss" = miss, "r" = r, "mod" = mod) %>%
      rownames_to_column("perc")
  } else {
    output <- final_df
  }
  return(list(output))
}

med_test <- list(
  mc_replicate(1000, expr = samp_func(mod = "f1", seed = NULL, miss = "MAR", r = 1-.15, option = 1)),
  mc_replicate(1000, expr = samp_func(mod = "f1", seed = NULL, miss = "MNAR", r = 1-.15, option = 1)),
  mc_replicate(1000, expr = samp_func(mod = "f2", seed = NULL, miss = "MAR", r = 1-.15, option = 1)),
  mc_replicate(1000, expr = samp_func(mod = "f2", seed = NULL, miss = "MNAR", r = 1-.15, option = 1)),
  mc_replicate(1000, expr = samp_func(mod = "f3", seed = NULL, miss = "MAR", r = 1-.15, option = 1)),
  mc_replicate(1000, expr = samp_func(mod = "f3", seed = NULL, miss = "MNAR", r = 1-.15, option = 1)),
  mc_replicate(1000, expr = samp_func(mod = "f4", seed = NULL, miss = "MAR", r = 1-.15, option = 1)),
  mc_replicate(1000, expr = samp_func(mod = "f4", seed = NULL, miss = "MNAR", r = 1-.15, option = 1)),
  
  mc_replicate(1000, expr = samp_func(mod = "f1", seed = NULL, miss = "MAR", r = .15, option = 1)),
  mc_replicate(1000, expr = samp_func(mod = "f1", seed = NULL, miss = "MNAR", r = .15, option = 1)),
  mc_replicate(1000, expr = samp_func(mod = "f2", seed = NULL, miss = "MAR", r = .15, option = 1)),
  mc_replicate(1000, expr = samp_func(mod = "f2", seed = NULL, miss = "MNAR", r = .15, option = 1)),
  mc_replicate(1000, expr = samp_func(mod = "f3", seed = NULL, miss = "MAR", r = .15, option = 1)),
  mc_replicate(1000, expr = samp_func(mod = "f3", seed = NULL, miss = "MNAR", r = .15, option = 1)),
  mc_replicate(1000, expr = samp_func(mod = "f4", seed = NULL, miss = "MAR", r = .15, option = 1)),
  mc_replicate(1000, expr = samp_func(mod = "f4", seed = NULL, miss = "MNAR", r = .15, option = 1)),
  mc_replicate(1000, expr = samp_func(mod = "f1", seed = NULL, miss = "MAR", r = 1-.30, option = 1)),
  mc_replicate(1000, expr = samp_func(mod = "f1", seed = NULL, miss = "MNAR", r = 1-.30, option = 1)),
  mc_replicate(1000, expr = samp_func(mod = "f2", seed = NULL, miss = "MAR", r = 1-.30, option = 1)),
  mc_replicate(1000, expr = samp_func(mod = "f2", seed = NULL, miss = "MNAR", r = 1-.30, option = 1)),
  mc_replicate(1000, expr = samp_func(mod = "f3", seed = NULL, miss = "MAR", r = 1-.30, option = 1)),
  mc_replicate(1000, expr = samp_func(mod = "f3", seed = NULL, miss = "MNAR", r = 1-.30, option = 1)),
  mc_replicate(1000, expr = samp_func(mod = "f4", seed = NULL, miss = "MAR", r = 1-.30, option = 1)),
  mc_replicate(1000, expr = samp_func(mod = "f4", seed = NULL, miss = "MNAR", r = 1-.30, option = 1)),
  
  mc_replicate(1000, expr = samp_func(mod = "f1", seed = NULL, miss = "MAR", r = .30, option = 1)),
  mc_replicate(1000, expr = samp_func(mod = "f1", seed = NULL, miss = "MNAR", r = .30, option = 1)),
  mc_replicate(1000, expr = samp_func(mod = "f2", seed = NULL, miss = "MAR", r = .30, option = 1)),
  mc_replicate(1000, expr = samp_func(mod = "f2", seed = NULL, miss = "MNAR", r = .30, option = 1)),
  mc_replicate(1000, expr = samp_func(mod = "f3", seed = NULL, miss = "MAR", r = .30, option = 1)),
  mc_replicate(1000, expr = samp_func(mod = "f3", seed = NULL, miss = "MNAR", r = .30, option = 1)),
  mc_replicate(1000, expr = samp_func(mod = "f4", seed = NULL, miss = "MAR", r = .30, option = 1)),
  mc_replicate(1000, expr = samp_func(mod = "f4", seed = NULL, miss = "MNAR", r = .30, option = 1))
) %>%
  bind_rows() %>%
  group_by(miss, mod, r) %>%
  dplyr::summarize_if(is.numeric, mean) %>%
  dplyr::select(miss, mod, r, med, hat_med)

dist_test <- list(
  mc_replicate(1000, expr = samp_func(mod = "f1", seed = NULL, miss = "MAR", r = 1-.15, option = 2)),
  mc_replicate(1000, expr = samp_func(mod = "f1", seed = NULL, miss = "MNAR", r = 1-.15, option = 2)),
  mc_replicate(1000, expr = samp_func(mod = "f2", seed = NULL, miss = "MAR", r = 1-.15, option = 2)),
  mc_replicate(1000, expr = samp_func(mod = "f2", seed = NULL, miss = "MNAR", r = 1-.15, option = 2)),
  mc_replicate(1000, expr = samp_func(mod = "f3", seed = NULL, miss = "MAR", r = 1-.15, option = 2)),
  mc_replicate(1000, expr = samp_func(mod = "f3", seed = NULL, miss = "MNAR", r = 1-.15, option = 2)),
  mc_replicate(1000, expr = samp_func(mod = "f4", seed = NULL, miss = "MAR", r = 1-.15, option = 2)),
  mc_replicate(1000, expr = samp_func(mod = "f4", seed = NULL, miss = "MNAR", r = 1-.15, option = 2)),
  
  mc_replicate(1000, expr = samp_func(mod = "f1", seed = NULL, miss = "MAR", r = .15, option = 2)),
  mc_replicate(1000, expr = samp_func(mod = "f1", seed = NULL, miss = "MNAR", r = .15, option = 2)),
  mc_replicate(1000, expr = samp_func(mod = "f2", seed = NULL, miss = "MAR", r = .15, option = 2)),
  mc_replicate(1000, expr = samp_func(mod = "f2", seed = NULL, miss = "MNAR", r = .15, option = 2)),
  mc_replicate(1000, expr = samp_func(mod = "f3", seed = NULL, miss = "MAR", r = .15, option = 2)),
  mc_replicate(1000, expr = samp_func(mod = "f3", seed = NULL, miss = "MNAR", r = .15, option = 2)),
  mc_replicate(1000, expr = samp_func(mod = "f4", seed = NULL, miss = "MAR", r = .15, option = 2)),
  mc_replicate(1000, expr = samp_func(mod = "f4", seed = NULL, miss = "MNAR", r = .15, option = 2)),
  mc_replicate(1000, expr = samp_func(mod = "f1", seed = NULL, miss = "MAR", r = 1-.30, option = 2)),
  mc_replicate(1000, expr = samp_func(mod = "f1", seed = NULL, miss = "MNAR", r = 1-.30, option = 2)),
  mc_replicate(1000, expr = samp_func(mod = "f2", seed = NULL, miss = "MAR", r = 1-.30, option = 2)),
  mc_replicate(1000, expr = samp_func(mod = "f2", seed = NULL, miss = "MNAR", r = 1-.30, option = 2)),
  mc_replicate(1000, expr = samp_func(mod = "f3", seed = NULL, miss = "MAR", r = 1-.30, option = 2)),
  mc_replicate(1000, expr = samp_func(mod = "f3", seed = NULL, miss = "MNAR", r = 1-.30, option = 2)),
  mc_replicate(1000, expr = samp_func(mod = "f4", seed = NULL, miss = "MAR", r = 1-.30, option = 2)),
  mc_replicate(1000, expr = samp_func(mod = "f4", seed = NULL, miss = "MNAR", r = 1-.30, option = 2)),
  
  mc_replicate(1000, expr = samp_func(mod = "f1", seed = NULL, miss = "MAR", r = .30, option = 2)),
  mc_replicate(1000, expr = samp_func(mod = "f1", seed = NULL, miss = "MNAR", r = .30, option = 2)),
  mc_replicate(1000, expr = samp_func(mod = "f2", seed = NULL, miss = "MAR", r = .30, option = 2)),
  mc_replicate(1000, expr = samp_func(mod = "f2", seed = NULL, miss = "MNAR", r = .30, option = 2)),
  mc_replicate(1000, expr = samp_func(mod = "f3", seed = NULL, miss = "MAR", r = .30, option = 2)),
  mc_replicate(1000, expr = samp_func(mod = "f3", seed = NULL, miss = "MNAR", r = .30, option = 2)),
  mc_replicate(1000, expr = samp_func(mod = "f4", seed = NULL, miss = "MAR", r = .30, option = 2)),
  mc_replicate(1000, expr = samp_func(mod = "f4", seed = NULL, miss = "MNAR", r = .30, option = 2)),
  mc_replicate(1000, expr = samp_func(mod = "f1", seed = NULL, miss = "MAR", r = 1-.001, option = 2)),
  mc_replicate(1000, expr = samp_func(mod = "f1", seed = NULL, miss = "MNAR", r = 1-.001, option = 2)),
  mc_replicate(1000, expr = samp_func(mod = "f2", seed = NULL, miss = "MAR", r = 1-.001, option = 2)),
  mc_replicate(1000, expr = samp_func(mod = "f2", seed = NULL, miss = "MNAR", r = 1-.001, option = 2)),
  mc_replicate(1000, expr = samp_func(mod = "f3", seed = NULL, miss = "MAR", r = 1-.001, option = 2)),
  mc_replicate(1000, expr = samp_func(mod = "f3", seed = NULL, miss = "MNAR", r = 1-.001, option = 2)),
  mc_replicate(1000, expr = samp_func(mod = "f4", seed = NULL, miss = "MAR", r = 1-.001, option = 2)),
  mc_replicate(1000, expr = samp_func(mod = "f4", seed = NULL, miss = "MNAR", r = 1-.001, option = 2)),
  
  mc_replicate(1000, expr = samp_func(mod = "f1", seed = NULL, miss = "MAR", r = .001, option = 2)),
  mc_replicate(1000, expr = samp_func(mod = "f1", seed = NULL, miss = "MNAR", r = .001, option = 2)),
  mc_replicate(1000, expr = samp_func(mod = "f2", seed = NULL, miss = "MAR", r = .001, option = 2)),
  mc_replicate(1000, expr = samp_func(mod = "f2", seed = NULL, miss = "MNAR", r = .001, option = 2)),
  mc_replicate(1000, expr = samp_func(mod = "f3", seed = NULL, miss = "MAR", r = .001, option = 2)),
  mc_replicate(1000, expr = samp_func(mod = "f3", seed = NULL, miss = "MNAR", r = .001, option = 2)),
  mc_replicate(1000, expr = samp_func(mod = "f4", seed = NULL, miss = "MAR", r = .001, option = 2)),
  mc_replicate(1000, expr = samp_func(mod = "f4", seed = NULL, miss = "MNAR", r = .001, option = 2))
) 

MAR_result = 
 dist_test %>%
  bind_rows() %>%
  group_by(perc, miss, mod, r) %>%
  dplyr::summarize_if(is.numeric, mean) %>%
  arrange(mod, perc) %>%
#  filter(perc == '50%') %>%
   
  mutate(err = y_N - lm_pred_A) %>%
   left_join( 
     med_test %>%
       dplyr::select(miss, mod, r, everything()) %>%
       pivot_longer(cols = 10:ncol(.)) %>%
       dplyr::select(-c(name, value)) %>%
       unique(), by = c('miss', 'mod', 'r'))


dist_test %>%
  bind_rows() %>%
  group_by(perc, miss, mod, r) %>%
  dplyr::summarize_if(is.numeric, mean) %>%
  filter(miss == 'MNAR') %>%
  mutate(err = y_N - lm_pred_A) %>%
  group_by(miss, mod, r) %>%
  filter(abs(err) == max(abs(err))) %>%
  filter(r %in% c(.15, 1-.15)) %>%
  dplyr::select(miss, mod, r, F_N, y_N, lm_pred_A, err) %>%
  rename(c('tN'= 'y_N', 'hat_tN' = 'lm_pred_A', 'alpha' = 'F_N')) %>%
  arrange(r, mod) %>%
  filter(r %in% c(
    .15,
    1-.15
  )) %>%
  filter(r == 1-.15) #%>%
  filter( %in% seq(.7, .8, .001))
  arrange(mod,r,  F_N) 
  filter(F_N %in% seq(.7, .8, .001), r == .15)
  print(n = 100)

  
  dist_test %>%
    bind_rows() %>%
    group_by(perc, miss, mod, r) %>%
    dplyr::summarize_if(is.numeric, mean) %>%
    filter(miss == 'MNAR') %>%
    mutate(err = y_N - lm_pred_A) %>%
    rename(c('tN'= 'y_N', 'hat_tN' = 'lm_pred_A', 'alpha' = 'F_N')) %>%
    ungroup() %>%
    #filter(alpha %in% c('0.1', '0.15', '0.2', '0.25', '0.3', '0.35'))%>%
    filter(r == .001) %>%
    arrange(mod, alpha) %>%
    print(n = 100)


#%>%
 mutate(F_N = gsub(x = perc, pattern = '%', replacement = '') %>% as.numeric() * 1/100) %>%
  arrange(mod, r) %>%
  group_by(mod, r) %>%
  filter(abs(err) == min(abs(err))) %>%
  mutate(rc = 1 - r) %>%
  mutate(diff_rc_F_N = abs(rc - F_N))

# %>%
bind_rows() %>%
  group_by(miss, mod) %>%
  dplyr::summarize_if(is.numeric, mean)

med_test %>%
  group_by(miss, mod) %>%
  dplyr::summarize_if(is.numeric, mean) %>%
  filter(miss != "MNAR")

# get B
all_info <- samp_func("f1", miss = "MAR", seed = NULL, r = .001)
B <- all_info[[1]]
A <- all_info[[2]]
mu <- all_info[[3]]
med <- all_info[[7]]
# build regression model
lm_mod <- lm(y ~ ., data = B %>% dplyr::select(
  "y", paste0("x", 1:p)
))

lm_pred_A <- predict(lm_mod, newdata = A)

# compare
data.frame(
  "mu" = mu,
  "hat_mu" = mean(lm_pred_A),
  "med" = med,
  "hat_med" = median(lm_pred_A)
)


lm_pred_B <- predict(lm_mod, newdata = B)

# calculate estimators
## B eCDF
B_plug <- ecdf(B$y)(q)
## Plug in
lm_plug <- ecdf(lm_pred_A)(q)
## reCDF
v_lm <- B$y - lm_pred_B
e_ecdf <- ecdf(v_lm)

m_lm <- c()

for (i in 1:length(q)) {
  m_lm[i] <- e_ecdf(q[i] - lm_pred_A) %>% mean()
}

# compile
names(B_plug) <- names(m_lm) <- names(lm_plug) <- names(q)

final_cdf <-
  data.frame(
    B_plug,
    lm_plug,
    m_lm
  ) %>%
  rownames_to_column("perc")

final_cdf
