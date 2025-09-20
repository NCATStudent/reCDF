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
require(openxlsx)
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
require(hablar)
require(pbapply)
require(pbmcapply)
require(mcreplicate)
require(foreach)
require(doParallel)
require(styler)
require(RcppAlgos)
require(gitcreds)

seed <- 101
nsim <- 5000 # if you wish to change L, make sure to also change it in the function below.
L <- 1
N <- 10000
nA <- .05*N
alp <- .10
mod <- "f3"
r <- .15

# Superpopulation Model
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
q_N <- quantile(pop$y, probs = F_N, type = 1)


As <- lapply(1:nsim, function(seed) {
  set.seed(seed)
  pop[sample(1:N, nA, replace = FALSE), ]
})

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
      mutate(restrat = ifelse(y < median(y), "0", "1") %>% as.factor()) %>%
      arrange(restrat)
    s.B <- sampling::strata(pop_restrat,
                            stratanames = "restrat",
                            size = sizes,
                            method = "srswor"
    )
    B <- getdata(pop_restrat, s.B)
  }
  return(B)
}


## MAR

Bs_MAR_1nA <- mclapply(1:nsim, function(seed) {
  B_generator(miss = "MAR", nB = .05*N, r = r)
}, mc.cores = 10)

Bs_MAR_20nA <- mclapply(1:nsim, function(seed) {
  B_generator(miss = "MAR", nB = .20*N, r = r)
}, mc.cores = 10)



## MNAR

Bs_MNAR_1nA <- mclapply(1:nsim, function(seed) {
  B_generator(miss = "MNAR", nB = .05*N, r = r)
}, mc.cores = 10)

Bs_MNAR_20nA <- mclapply(1:nsim, function(seed) {
  B_generator(miss = "MNAR", nB = .20*N, r = r)
}, mc.cores = 10)

# for testing
# B_perm_MAR_1nA= Bs_MAR_1nA[[1]]
# B_perm_MAR_20nA = Bs_MAR_20nA[[1]]
# B_perm_MNAR_1nA = Bs_MNAR_1nA[[1]]
# B_perm_MNAR_20nA = Bs_MNAR_20nA[[1]]
# A = As[[1]]
# B = B_perm_MNAR_20nA



samp.f <- function(samp, A, B_MAR_nA, B_MAR_20nA, B_MNAR_1nA, B_MNAR_20nA, r, pop, p, miss, L = 1, mc_cores = 15, alp = .10, seed = 101) {
  options(dplyr.summarise.inform = FALSE)
  
  var_cdf <- function(B, A, p, alp, da_miss, isboot, mc_cores = 15) {
    nA <- nrow(A)
    nB <- nrow(B)
    # defining A
    ecdf_A <- ecdf(A$y)
    ecdf_B <- ecdf(B$y)
    
    # lmod
    lm_B <- lm(y ~ .,
               data = B[, c(
                 paste0("x", 1:p),
                 "y"
               )]
    )
    lm_pred_B <- predict(lm_B, newdata = B)
    lm_pred_A <- predict(lm_B, newdata = A)
    
    pA_pred <- svydesign(
      id = ~1,
      weights = ~weight,
      fpc = ~fpc,
      data = A %>%
        mutate(
          fpc = nrow(A) / sum(A$weight),
          y = lm_pred_A
        )
    )
    ecdf_lm <- svycdf(~y, pA_pred)
    
    
    #
    # # plug-in
    
    A_plug <- ecdf_A(q_N)
    B_plug <- ecdf_B(q_N)
    lm_plug <- ecdf_lm[[1]](q_N)
    names(lm_plug) <- names(A_plug) <- names(lm_plug) <- names(q_N)
    
    
    # reCDF
    v_lm <- B$y - lm_pred_B
    e_ecdf <- ecdf(v_lm)
    
    mlm_est_func <- function(qN_vals, nA) {
      m_lm <- (1 / sum(A$weight)) * sum(A$weight * e_ecdf(qN_vals - lm_pred_A))
      
      v1 <- (1 - nA / N) / nA * var(e_ecdf(qN_vals - lm_pred_A))
      
      # this part is gonna be computationally painful, fair warning :/
      raw_index <- expandGrid(nThreads = detectCores() -1, return_df = TRUE,i_it = 1:nA, h_it = 1:nA)
      
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
    
    # idea here -- don't run asymp variance for boot reps
    if (isboot == TRUE) {
      m_lm <- vars_to <- c()
      for (i in 1:length(q_N)) {
        m_lm[i] <- (1 / sum(A$weight)) * sum(A$weight * e_ecdf(q_N[i] - lm_pred_A))
        vars_to[i] <- NA
      }
    } else {
      vars_to_results <- mcmapply(
        FUN = mlm_est_func,
        qN_vals = q_N,
        MoreArgs = (
          list(
            "nA" = nA
          )
        ), mc.cores = mc_cores
      )
      
      vars_to <-
        vars_to_results[1, ] %>%
        bind_rows() %>%
        unlist()
      
      
      m_lm <-
        vars_to_results[2, ] %>%
        bind_rows() %>%
        unlist()
    }
    names(vars_to) <- names(m_lm) <- names(q_N)
    
    final_cdf <- rbind(
      F_N,
      A_plug,
      B_plug,
      lm_plug,
      m_lm
    ) %>%
      as.data.frame() %>%
      rownames_to_column("name") %>%
      mutate(group = "cdf") %>%
      mutate(val = "mc")
    
    
    
    # organizing variances
    
    
    quant_est <- function(F_N, nA) {
      # sorry this is so messy, but ie_result is needed
      ### EXCUSE THE MESS
      
      # R
      q_mlm <- function(a, b, e_ecdf = e_ecdf) {
        return(e_ecdf(a - b) %>% mean())
      }
      
      vq <- Vectorize(q_mlm,
                      vectorize.args = "a"
      )
      
      mlm_df <- data.frame(
        probs = vq(lm_pred_A, lm_pred_A, e_ecdf = e_ecdf),
        t = lm_pred_A
      ) %>%
        arrange(probs)
      
      ie_fun <- function(dats, a) {
        return(dats %>%
                 lapply(dplyr::filter, probs >= a) %>%
                 lapply(function(x) {
                   if (nrow(x) == 0) {
                     data.frame(
                       probs = NA,
                       t = NA
                     )
                   } else {
                     x
                   }
                 }) %>%
                 lapply(setNames, c("probs", "t")) %>%
                 lapply(dplyr::arrange, t) %>%
                 lapply(dplyr::slice, 1) %>%
                 lapply(dplyr::select, -probs) %>%
                 bind_rows() %>%
                 mutate(
                   name = c(
                     "m_lm"
                   ) %>% as.factor(),
                   a = a
                 ) %>%
                 dplyr::select(a, name, t))
      }
      
      ie_result <- mapply(
        FUN = ie_fun,
        a = F_N,
        MoreArgs = list(dats = list(
          mlm_df
        )),
        SIMPLIFY = FALSE
      )
      ## OK WE"RE DONE
      a.a <- F_N
      hatq <- ie_result %>%
        bind_rows() %>%
        filter(name == "m_lm") %>%
        filter(a == a.a) %>%
        dplyr::select(t) %>%
        unlist() %>%
        unname()
      
      
      hatF_hatT <- (1 / sum(A$weight)) * sum(A$weight * e_ecdf(hatq - lm_pred_A))
      
      
      v1 <- (1 - nA / N) / nA * var(e_ecdf(hatq - lm_pred_A))
      
      # this part has to get the min, fair warning :/
      raw_index <- expandGrid(nThreads = detectCores() -1, return_df = TRUE,i_it = 1:nA, h_it = 1:nA)
      
      Rh_df <- A[raw_index["h_it"] %>% unlist(), ] %>%
        dplyr::select(-c(y, weight, Prob)) %>%
        mutate(Rh = hatq - predict(lm_B, newdata = .)) %>%
        mutate(h_it = raw_index["h_it"]) %>%
        dplyr::select(h_it, Rh)
      
      Ri_df <- A[raw_index["i_it"] %>% unlist(), ] %>%
        dplyr::select(-c(y, weight, Prob)) %>%
        mutate(Ri = hatq - predict(lm_B, newdata = .)) %>%
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
                                        mutate(Rh = hatq - predict(lm_B, newdata = .)) %>%
                                        mutate(Ghat_h = e_ecdf(Rh)) %>%
                                        dplyr::select(Ghat_h) %>%
                                        unlist())^2
      vars_q <- v1 + (v2_1 - v2_2)
      
      
      LL.q_N <- list(mlm_df) %>%
        lapply(dplyr::filter, probs >= a.a - qnorm(1 - alp / 2) * sqrt(vars_q)) %>%
        lapply(function(x) {
          if (nrow(x) == 0) {
            data.frame(
              probs = NA,
              t = NA
            )
          } else {
            x
          }
        }) %>%
        lapply(setNames, c("probs", "t")) %>%
        lapply(dplyr::arrange, t) %>%
        lapply(dplyr::slice, 1) %>%
        bind_rows() %>%
        mutate(a = a.a) %>%
        dplyr::select(t) %>%
        unlist() %>%
        unname()
      
      UL.q_N <- list(mlm_df) %>%
        lapply(dplyr::filter, probs >= a.a + qnorm(1 - alp / 2) * sqrt(vars_q)) %>%
        lapply(function(x) {
          if (nrow(x) == 0) {
            data.frame(
              probs = NA,
              t = NA
            )
          } else {
            x
          }
        }) %>%
        lapply(setNames, c("probs", "t")) %>%
        lapply(dplyr::arrange, t) %>%
        lapply(dplyr::slice, 1) %>%
        bind_rows() %>%
        mutate(a = a.a) %>%
        dplyr::select(t) %>%
        unlist() %>%
        unname()
      
      q_var <- vars_q
      return(list(data.frame(F_N, hatq, q_var, LL.q_N, UL.q_N, hatF_hatT)))
    }
    
    
    # idea here -- don't run asymp variance for boot reps
    if (isboot == TRUE) {
      # Quantile Estimation
      
      # R
      q_mlm <- function(a, b, e_ecdf = e_ecdf) {
        return(e_ecdf(a - b) %>% mean())
      }
      
      vq <- Vectorize(q_mlm,
                      vectorize.args = "a"
      )
      
      mlm_df <- data.frame(
        probs = vq(lm_pred_A, lm_pred_A, e_ecdf = e_ecdf),
        t = lm_pred_A
      ) %>%
        arrange(probs)
      
      ie_fun <- function(dats, a) {
        return(dats %>%
                 lapply(dplyr::filter, probs >= a) %>%
                 lapply(function(x) {
                   if (nrow(x) == 0) {
                     data.frame(
                       probs = NA,
                       t = NA
                     )
                   } else {
                     x
                   }
                 }) %>%
                 lapply(setNames, c("probs", "t")) %>%
                 lapply(dplyr::arrange, t) %>%
                 lapply(dplyr::slice, 1) %>%
                 lapply(dplyr::select, -probs) %>%
                 bind_rows() %>%
                 mutate(
                   name = c(
                     "m_lm"
                   ) %>% as.factor(),
                   a = a
                 ) %>%
                 dplyr::select(a, name, t))
      }
      
      ie_result <- mapply(
        FUN = ie_fun,
        a = F_N,
        MoreArgs = list(dats = list(
          mlm_df
        )),
        SIMPLIFY = FALSE
      )
      
      
      
      final_quant <- ie_result %>%
        bind_rows() %>%
        mutate(F_N = paste0(a * 100, "%")) %>%
        dplyr::select(-a) %>%
        pivot_wider(
          names_from = "F_N",
          values_from = "t"
        ) %>%
        mutate(group = "q_N", val = "mc") %>%
        as.data.frame() %>%
        mutate(
          group = "q_N",
          val = "mc"
        ) %>%
        rbind(q_N %>%
                t() %>%
                as.data.frame() %>%
                mutate(
                  name = "F_N",
                  group = "q_N",
                  val = "mc"
                ))
      
      hatF_hatT <- hatq <- LL.q <- UL.q <- q_var <- c()
      for (i in 1:length(F_N)) {
        a.a <- F_N[i]
        hatq[i] <- ie_result %>%
          bind_rows() %>%
          filter(name == "m_lm") %>%
          filter(a == a.a) %>%
          dplyr::select(t) %>%
          unlist() %>%
          unname()
        
        hatF_hatT[i] <- (1 / sum(A$weight)) * sum(A$weight * e_ecdf(hatq[i] - lm_pred_A))
        
        # since I have to do a for-loop anyway...
        
        LL.q[i] <- UL.q[i] <- q_var[i] <- NA
      }
      
      names(LL.q) <- names(UL.q) <- names(hatq) <- names(hatF_hatT) <- names(q_N)
      
      q_var_df <- rbind(q_N, hatq, q_var, LL.q, UL.q, hatF_hatT) %>%
        as.data.frame() %>%
        rownames_to_column("name") %>%
        mutate(group = "q", val = "var")
      
      final_df <- rbind(
        final_cdf,
        final_quant,
        q_var_df
      ) %>%
        filter(name %in% c(
          "m_lm",
          "hatq",
          "q_var",
          "LL.q",
          "UL.q",
          "hatF_hatT"
        )) %>%
        as.data.frame() %>%
        dplyr::slice(-2) %>%
        rbind(t(vars_to) %>%
                as.data.frame() %>%
                mutate(
                  name = "m_lm_var",
                  group = "q",
                  val = "mc"
                ) %>%
                dplyr::select(colnames(final_cdf))) %>%
        dplyr::select(-c(group, val))
    } else {
      qvar_results <- mcmapply(
        FUN = quant_est,
        F_N = F_N,
        MoreArgs = (
          list(
            "nA" = nA
          )
        ), mc.cores = mc_cores
      ) %>%
        bind_rows()
      
      hatq <- qvar_results %>%
        dplyr::select("hatq") %>%
        unlist()
      
      q_var <- qvar_results %>%
        dplyr::select("q_var") %>%
        unlist()
      
      LL.q_N <- qvar_results %>%
        dplyr::select("LL.q_N") %>%
        unlist()
      
      UL.q_N <- qvar_results %>%
        dplyr::select("UL.q_N") %>%
        unlist()
      
      hatF_hatT <- qvar_results %>%
        dplyr::select("hatF_hatT") %>%
        unlist()
      
      
      names(LL.q_N) <- names(UL.q_N) <- names(hatq) <- names(hatF_hatT) <- names(q_N)
      
      q_var_df <- rbind(q_N, hatq, q_var, LL.q_N, UL.q_N, hatF_hatT) %>%
        as.data.frame() %>%
        rownames_to_column("name") %>%
        mutate(group = "q_N", val = "var")
      
      final_df <- rbind(
        final_cdf,
        q_var_df
      ) %>%
        # filter(name %in% c(
        #   "m_lm",
        #   "hatq",
        #   "q_var",
        #   "LL.q_N",
        #   "UL.q_N"
        # )) %>%
        as.data.frame() %>%
        dplyr::slice(-2) %>%
        rbind(t(vars_to) %>%
                as.data.frame() %>%
                mutate(
                  name = "m_lm_var",
                  group = "q_N",
                  val = "mc"
                ) %>%
                dplyr::select(colnames(final_cdf))) %>%
        dplyr::select(-c(group, val)) %>%
        mutate(nB = nB, miss = da_miss)
      
      # I know this is messy, but this is actually the easiest way to do this
      
      # R
      q_mlm <- function(a, b, e_ecdf = e_ecdf) {
        return(e_ecdf(a - b) %>% mean())
      }
      
      vq <- Vectorize(q_mlm,
                      vectorize.args = "a"
      )
      
      mlm_df <- data.frame(
        probs = vq(lm_pred_A, lm_pred_A, e_ecdf = e_ecdf),
        t = lm_pred_A
      ) %>%
        arrange(probs)
    }
    
    return(list(final_df, mlm_df))
  }
  # selecting A
  
  
  A_perm <- A %>%
    dplyr::select(paste0("x", 1:p), "y") %>%
    mutate(
      Prob = nA / N,
      weight = N / nA,
      fpc = Prob
    )
  
  B_perm_MAR_1nA <- B_MAR_nA
  B_perm_MAR_20nA <- B_MAR_20nA
  B_perm_MNAR_1nA <- B_MNAR_1nA
  B_perm_MNAR_20nA <- B_MNAR_20nA
  
  
  mydesign <- svydesign(
    ids = ~0,
    weights = ~weight,
    fpc = ~fpc,
    data = A_perm
  )
  
  
  bootstrap_rep_design <- as_bootstrap_design(mydesign,
                                              type = "Rao-Wu-Yue-Beaumont",
                                              replicates = L,
                                              samp_method_by_stage = "SRSWOR"
  )
  
  A_repweights_perm <- as_data_frame_with_weights(bootstrap_rep_design,
                                                  full_wgt_name = "FULL_SAMPLE_WGT",
                                                  rep_wgt_prefix = "REP_WGT_"
  ) %>%
    dplyr::select(-weight)
  #
  #
  #
  A_list <- c()
  for (i in 1:L) {
    A_list[[i]] <- A_repweights_perm %>%
      dplyr::select(y, paste0("x", 1:p), paste0("REP_WGT_", i), Prob) %>%
      dplyr::rename(weight = paste0("REP_WGT_", i))
  }
  
  
  # selecting B
  
  ## nB = .05*N
  set.seed(seed)
  B_list.MAR_1nA <- mc_replicate(L, list(B_perm_MAR_1nA[sample(1:(.05*N),
                                                               size = .05*N,
                                                               replace = TRUE
  ), ]), mc.cores = mc_cores)
  
  set.seed(seed)
  B_list.MNAR_1nA <- mc_replicate(L, list(B_perm_MNAR_1nA[sample(1:(.05*N),
                                                                 size = .05*N,
                                                                 replace = TRUE
  ), ]), mc.cores = mc_cores)
  
  ## nB = 20000
  set.seed(seed)
  B_list.MAR_20nA <- mc_replicate(L, list(B_perm_MAR_20nA[sample(1:(20 * .05*N),
                                                                 size = 20 * .05*N,
                                                                 replace = TRUE
  ), ]), mc.cores = mc_cores)
  
  
  set.seed(seed)
  B_list.MNAR_20nA <- mc_replicate(L, list(B_perm_MNAR_20nA[sample(1:(20 * .05*N),
                                                                   size = 20 * .05*N,
                                                                   replace = TRUE
  ), ]), mc.cores = mc_cores)
  
  
  # For each a in A and b in B, compute bootstrap replicates
  #
  results.MAR_1nA <- mcmapply(
    FUN = var_cdf,
    A = A_list,
    B = B_list.MAR_1nA,
    MoreArgs = list(p = p, alp = .10, da_miss = "MAR", isboot = TRUE),
    mc.cores = mc_cores,
    SIMPLIFY = FALSE
  )
  
  results.MAR_20nA <- mcmapply(var_cdf,
                               A = A_list,
                               B = B_list.MAR_20nA,
                               MoreArgs = list("p" = p, "alp" = .10, da_miss = "MAR", isboot = TRUE),
                               mc.cores = mc_cores,
                               SIMPLIFY = FALSE
  )
  
  results.MNAR_1nA <- mcmapply(var_cdf,
                               A = A_list,
                               B = B_list.MNAR_1nA,
                               MoreArgs = list("p" = p, "alp" = .10, da_miss = "MAR", isboot = TRUE),
                               mc.cores = mc_cores,
                               SIMPLIFY = FALSE
  )
  
  results.MNAR_20nA <- mcmapply(var_cdf,
                                A = A_list,
                                B = B_list.MNAR_20nA,
                                MoreArgs = list("p" = p, "alp" = .10, da_miss = "MAR", isboot = TRUE),
                                mc.cores = mc_cores,
                                SIMPLIFY = FALSE
  )
  
  results2 <- list(
    results.MAR_1nA %>% bind_rows() %>% mutate(nB = nA, miss = "MAR"),
    results.MAR_20nA %>% bind_rows() %>% mutate(nB = .20*N, miss = "MAR"),
    results.MNAR_1nA %>% bind_rows() %>% mutate(nB = nA, miss = "MNAR"),
    results.MNAR_20nA %>% bind_rows() %>% mutate(nB = .20*N, miss = "MNAR")
  ) %>%
    lapply(filter, name %in% c("m_lm", "hatq", "hatF_hatT"))
  
  
  # Testing
  #B_perm = Bs_MAR_20nA[[1]]
  var_calcs <- function(A_perm, B_perm, L, alp, da_miss, da_nB) {
    asymp_var_cdf <- var_cdf(
      A = A_perm,
      B = B_perm,
      p = p,
      alp = .10,
      da_miss = da_miss,
      isboot = FALSE
    )
    
    actual <- asymp_var_cdf[[1]]
    
    
    standard_values <-
      actual %>%
      filter(name %in% c(
        "m_lm",
        "hatq",
        "hatF_hatT"
      )) %>%
      dplyr::select(nB, miss, name, everything()) %>%
      pivot_longer(
        cols = 4:ncol(.),
        names_to = "perc",
        values_to = "B_vals"
      )
    
    boot_var <-
      results2 %>%
      bind_rows(.id = "iter") %>%
      dplyr::select(-c(probs, t)) %>%
      filter(miss == da_miss, nB == da_nB) %>%
      filter(name %in% c(
        "m_lm",
        "hatq",
        "hatF_hatT"
      )) %>%
      dplyr::select(iter, name, miss, nB, everything()) %>%
      pivot_longer(
        cols = 5:ncol(.),
        names_to = "perc",
        values_to = "boot_vals"
      ) %>%
      left_join(standard_values,
                by = c("nB", "miss", "name", "perc")
      ) %>%
      group_by(nB, miss, name, perc) %>%
      dplyr::summarize(
        boot_mean_est = mean(boot_vals, na.rm = TRUE),
        B_perm_est = mean(B_vals, na.rm = TRUE),
        boot_var = (1 / L) * sum_((boot_vals - B_vals)^2)
      ) %>%
      dplyr::select(nB, miss, name, perc, boot_var, everything()) %>%
      pivot_longer(cols = 6:ncol(.), names_to = "bootstats") %>%
      pivot_wider(values_from = c("boot_var", "value")) %>%
      filter(bootstats == "boot_mean_est") %>%
      dplyr::select(-bootstats) %>%
      dplyr::select(nB, miss, perc, everything())
    
    pop_quant <- data.frame(F_N = F_N, q_N = q_N, perc = paste0(F_N * 100, "%"))
    
    
    bigguy <-
      actual %>%
      dplyr::select(nB, miss, name, everything()) %>%
      pivot_longer(cols = 4:ncol(.), names_to = "perc") %>%
      pivot_wider() %>%
      left_join(boot_var %>%
                  dplyr::select(-c(value_hatq, value_m_lm, value_hatF_hatT)), by = c("nB", "miss", "perc"))
    
    # reCDF est variance
    mlm_varsum <-
      bigguy %>%
      dplyr::select(
        F_N, nB, miss, perc,
        m_lm,
        m_lm_var,
        boot_var_m_lm
      ) %>%
      setNames(c(
        "pop_quant", "nB", "miss", "perc",
        "est_quant",
        "asymp",
        "boot"
      )) %>%
      dplyr::select(asymp, boot, everything()) %>%
      pivot_longer(
        cols = 1:2,
        names_to = "var_type",
        values_to = "var_val"
      ) %>%
      rowwise() %>%
      mutate(
        LL = est_quant - qnorm(p = 1 - alp / 2) * sqrt(var_val),
        UL = est_quant + qnorm(p = 1 - alp / 2) * sqrt(var_val)
      ) %>%
      mutate(CR = ifelse(pop_quant >= LL & pop_quant <= UL, 1, 0)) %>%
      dplyr::select(
        nB, miss,
        perc,
        pop_quant,
        est_quant,
        var_type,
        CR,
        everything()
      ) %>%
      ungroup() %>%
      mutate(est_type = "cdf")
    
    ## recall, q_var is the estimated variance of F(hat(T))
    
    # quantile
    hatq_varsum <-
      bigguy %>%
      dplyr::select(
        q_N, LL.q_N, UL.q_N, nB, miss, perc,
        hatq,
        boot_var_hatq
      ) %>%
      setNames(c(
        "pop_quant", "LL.q_N", "UL.q_N", "nB", "miss", "perc",
        "est_quant",
        "boot"
      )) %>%
      mutate(asymp = ((UL.q_N - LL.q_N) / (2 * qnorm(1 - (alp / 2))))^2) %>%
      dplyr::select(
        boot, asymp,
        pop_quant, LL.q_N, UL.q_N, nB, miss, perc, est_quant
      ) %>%
      pivot_longer(
        cols = 1:2,
        names_to = "var_type",
        values_to = "var_val"
      ) %>%
      rowwise() %>%
      mutate(
        LL = ifelse(var_type == "boot",
                    est_quant - qnorm(1 - alp / 2) * sqrt(var_val),
                    LL.q_N
        ),
        UL = ifelse(var_type == "boot",
                    est_quant + qnorm(1 - alp / 2) * sqrt(var_val),
                    UL.q_N
        )
      ) %>%
      dplyr::select(-c(LL.q_N, UL.q_N)) %>%
      mutate(CR = ifelse(pop_quant >= LL & pop_quant <= UL, 1, 0)) %>%
      ungroup() %>%
      mutate(est_type = "t")
    
    # new bootstrapped based confidence interval
    crit_vals <-
      boot_var %>%
      dplyr::select(nB, miss, perc, boot_var_hatF_hatT) %>%
      mutate(pop_quant = F_N) %>%
      mutate(
        critval_LL = pop_quant - qnorm(1 - alp / 2) * sqrt(boot_var_hatF_hatT),
        critval_UL = pop_quant + qnorm(1 - alp / 2) * sqrt(boot_var_hatF_hatT)
      ) %>%
      dplyr::select(critval_LL, critval_UL)
    
    mlm_df <- asymp_var_cdf[[2]]
    
    LL_q <- UL_q <- c()
    for (i in 1:nrow(crit_vals)) {
      LL_potential <- try(mlm_df %>%
                            filter(probs >= crit_vals[i, "critval_LL"] %>% unlist()) %>%
                            arrange(t) %>%
                            dplyr::slice(1) %>%
                            dplyr::select(t) %>%
                            unlist(), silent = TRUE)
      UL_potential <- try(mlm_df %>%
                            filter(probs >= crit_vals[i, "critval_UL"] %>% unlist()) %>%
                            arrange(t) %>%
                            dplyr::slice(1) %>%
                            dplyr::select(t) %>%
                            unlist(), silent = TRUE)
      
      LL_q[i] <- ifelse(length(LL_potential) != 0, LL_potential, NA)
      UL_q[i] <- ifelse(length(UL_potential) != 0, UL_potential, NA)
    }
    
    proposed_boot_results <-
      actual %>%
      filter(name == "hatq") %>%
      dplyr::select(name, nB, miss, everything()) %>%
      pivot_longer(cols = 4:ncol(.), names_to = "perc", values_to = "hatq") %>%
      mutate(LL = LL_q, UL = UL_q) %>%
      mutate(pop_quant = q_N) %>%
      mutate(
        CR = ifelse(pop_quant >= LL & pop_quant <= UL, 1, 0),
        var_val = ((UL - LL) / (2 * qnorm(1 - (alp / 2))))^2
      ) %>%
      mutate(
        est_quant = hatq,
        var_type = "boot2",
        est_type = "t"
      ) %>%
      dplyr::select(colnames(hatq_varsum))
    
    hatq_varsum %<>%
      bind_rows(proposed_boot_results) %>%
      arrange(est_quant)
    
    final <- rbind(mlm_varsum, hatq_varsum) %>%
      dplyr::select(est_type, nB, miss, perc, everything())
    return(final)
  }
  
  final_results <- list(
    var_calcs(A_perm, B_perm_MAR_1nA, L, alp, da_miss = "MAR", da_nB = 1 * nA),
    var_calcs(A_perm, B_perm_MAR_20nA, L, alp, da_miss = "MAR", da_nB = .20*N),
    var_calcs(A_perm, B_perm_MNAR_1nA, L, alp, da_miss = "MNAR", da_nB = 1 * nA),
    var_calcs(A_perm, B_perm_MNAR_20nA, L, alp, da_miss = "MNAR", da_nB = .20*N)
  ) %>%
    lapply(bind_rows) %>%
    bind_rows()
  
  return(final_results)
}

results <- c()
for (i in 1:nsim) {
  results[[i]] <- samp.f(
    A = As[[i]],
    B_MAR_nA = Bs_MAR_1nA[[i]],
    B_MAR_20nA = Bs_MAR_20nA[[i]],
    B_MNAR_1nA = Bs_MNAR_1nA[[i]],
    B_MNAR_20nA = Bs_MNAR_20nA[[i]],
    p = p, seed = seed
  )
  # in case computer dies, cache every 100 iters
  if (i %% 100 == 0) {
    setwd('/Users/jeremyflood/Library/CloudStorage/OneDrive-Personal/Documents/Grad School/2024-2025/Fall 2025/reCDF/reCDF/Variance Estimation/Additional Requests/Confirming Variance/smaller N/Data/Cached Iter Files')
    results[(i - 100):100] %>%
      bind_rows() %>%
      mutate(mod = mod) %>%
      openxlsx::write.xlsx(
        paste0("f3_raw_results_iter_", i - 100, "_", i, ".xlsx")
      )
    progress(i, nsim)
  } else {
    progress(i, nsim)
  }
}

cleaned_ddf <- results %>%
  bind_rows() # %>%
# mutate(var_name = var_type) %>%
# dplyr::select(-var_type) #%>%
# mutate(est_name= ifelse(is.na(m_lm) == is.na(F_N) & is.na(F_N) == TRUE, 't', 'cdf'),
#        est_value = ifelse(est_name == 'cdf', m_lm, hatq),
#        actual_value = ifelse(est_name == 'cdf', F_N, q_N)) %>%
# dplyr::select(F_N, q_N, hatq, m_lm, everything())

mc_results <- cleaned_ddf %>%
  # filter(var_name == 'asymp') %>%
  group_by(perc, miss, nB, est_type, var_type) %>%
  dplyr::summarize(MC_var = var(est_quant))


cleaned_results <- cleaned_ddf %>%
  # left_join(mc_results) %>%
  dplyr::select(est_type, est_quant, pop_quant, everything()) %>%
  pivot_longer(cols = 2) %>%
  group_by(est_type, perc, miss, nB, var_type, name) %>%
  dplyr::summarize(
    pop_quant = mean(pop_quant, na.rm = TRUE),
    est_value = mean(value, na.rm= TRUE),
    est_var = mean(var_val, na.rm = TRUE),
    CR = mean(CR, na.rm = TRUE),
    d = mean(UL - LL, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  left_join(mc_results, by = c("est_type", "perc", "miss", "nB", "var_type")) %>%
  mutate(rb = ((MC_var - est_var) / MC_var) * 100) %>%
  mutate(nsim = nsim, L = L)

final_results <- list(
  cleaned_ddf %>% mutate(mod = mod, seed = seed),
  cleaned_results %>% mutate(mod = mod, seed = seed)
)

names(final_results) <- c("raw", "summary")

setwd('/Users/jeremyflood/Library/CloudStorage/OneDrive-Personal/Documents/Grad School/2024-2025/Fall 2025/reCDF/reCDF/Variance Estimation/Additional Requests/Confirming Variance/smaller N/Data')
openxlsx::write.xlsx(final_results, paste0("modf3_results.xlsx"))
