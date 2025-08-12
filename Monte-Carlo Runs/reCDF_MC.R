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
Sample=Vectorize(sample)
Dat = function(x, p, size, seed){
  set.seed(seed)
  return(x[x$part==p,][sample(nrow(x[x$part==p,]), size),])
}
Dat.v = Vectorize(Dat, vectorize.args=c("size", "p"))
require(MASS)
require(mgcv) 
require(parallel)
require(ranger)
require(xgboost)
library(dplyr)
library(mgcv)
library(tidymv)
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

sim.f2 = function(N,
                  nA,
                  nB,
                  nsim,
                  mod, 
                  seed, 
                  r){
  # Superpopulation Model
  if(mod == 'f1'){ 
    set.seed(seed)
    x1 = rnorm(N, mean = 2, sd = 1)
    x2 = rnorm(N, mean = 2, sd = 1)
    x3 = rnorm(N, mean = 4, sd = 1)
    x4 = rnorm(N, mean = 4, sd = 1)
    X=cbind(x1, x2, x3, x4)
    p = ncol(X)
    y = rnorm(N, 4*x1+4*x2 + 2*x3 + 2*x4, 3)
  }  else if(mod == 'f2'){
    set.seed(seed)
    x1= runif(N, min = 0, max = 4)
    x2= runif(N, min = 0, max = 4)
    x3= runif(N, min = 4, max = 8)
    x4 = runif(N, min = 4, max = 8)
    X = cbind(x1, x2, x3, x4) 
    p = ncol(X)
    y = rnorm(N, 4*x1^2 +4*x2^2 + 2*x3^2 + 2*x4^2 + (x1+x2)^2 + (x3+x4)^2, sd =50)
  } else if(mod == 'f3'){
    set.seed(seed)
    x1 = runif(N, min = -1, max = 1)
    x2 = runif(N, min = -1, max = 1)
    x3 = runif(N, min = -1, max = 1)
    x4 = runif(N, min = -1, max = 1)
    X = data.frame(x1, x2, x3, x4)
    p = ncol(X)
    err = rnorm(N, mean =0, sd = sqrt(.5))
    y = -sin(x1)+x2^2+x3-exp(-x4^2)+err
  } else if(mod == 'f4'){
    set.seed(seed)
    x1 = rnorm(N)
    x2 = rnorm(N)
    x3 = rnorm(N)
    x4 = rnorm(N)
    x5 = rnorm(N)
    x6 = rnorm(N)
    X = data.frame(x1, x2, x3, x4, x5, x6)
    p = ncol(X)
    err = rnorm(N)
    y = x1 + .707*x2^2 + 2*ifelse(x3>0, 1, 0) + .873*log(abs(x1))*abs(x3) +
      .894*x2*x4 + 2*ifelse(x5>0, 1, 0) + .46*exp(x6)+err
  } else{
    stop('You must select one of the following models: (f1, f2, f3, f4)')
  }
  
  
  
  pop = as.data.frame(cbind(X, y)) %>%
    mutate(mnar_strat = ifelse(y < median(y), "0", "1") %>% as.factor() 
           
    ) %>%
    arrange(mnar_strat)
  
  
  # population characteristics
  mu = mean(pop$y)
  F_N = c(.01, .10, .25, .50, .75, .90, .99)
  q= quantile(pop$y, probs= F_N, type = 1)
  
  
  # sampling function
  B_generator = function(miss, nB, r){
    sizes_rounded = round(c(r*nB, (1-r)*nB), 
                          digits =0)
    
    sizes = ifelse(sizes_rounded <1, 1, sizes_rounded)
    if(miss == 'MAR'){
      strat = paste0('x', c(abs(cor(X,y)) %>% which.max()))
      pop_restrat= pop %>%
        mutate(restrat = ifelse(get(strat) <= median(get(strat)), 1,
                                2)) %>%
        arrange(restrat)
      s.B = sampling::strata(pop_restrat,
                             stratanames = "restrat",
                             size = sizes,
                             method = 'srswor')
      B = getdata(pop_restrat, s.B) 
    } else if(miss == 'MNAR'){
      s.B = sampling::strata(pop,
                             stratanames = "mnar_strat",
                             size = sizes,
                             method = 'srswor')
      B = getdata(pop, s.B) 
    }
    return(B)
  }
  B = B_generator('MAR', nB, r)
  MI_cdf = function(B, A, p){
    
    # A
    s.A= sample(1:N, nA, replace=FALSE)
    A = pop[s.A, ] %>%
      mutate(prob = nA/N, 
             weight = 1/prob)
    ecdf_A = ecdf(A$y)
    A_plug = ecdf_A(q)
    names(A_plug) = names(q)
    
    # naive estimator
    ecdf_B = ecdf(B$y)
    B_plug = ecdf_B(q)
    names(B_plug) = names(q)
    
    #lmod
    lm_B = lm(y~., 
              data =B[, c(paste0('x', 1:p),
                          'y')])
    lm_pred_B = predict(lm_B, newdata = B)
    lm_pred_A = predict(lm_B, newdata = A)
    
    # test -- lmod but with A
    
    # lm_A = lm(y~., 
    #           data =A[, c(paste0('x', 1:p),
    #                       'y')])
    # lm_Apred_A = predict(lm_A, newdata = A)
    
    ## plug in
    ecdf_lm = ecdf(lm_pred_A)
    lm_plug = ecdf_lm(q)
    names(lm_plug) = names(q)
    
    ## recdf
    
    # reCDF
    v_lm = B$y - lm_pred_B
    #v_lmA = A$y - lm_Apred_A
    e_ecdf= ecdf(v_lm)
    #e_ecdfA = ecdf(v_lmA)
    
    m_lm = vars_to= c()
    
    for(i in 1:length(q)){
      m_lm[i] = e_ecdf(q[i] - lm_pred_A) %>% mean()
      #m_lmA[i] = e_ecdfA(q[i] - lm_Apred_A) %>% mean()
    }
    
    names(m_lm) = names(q)
    
    
    final_cdf =rbind(F_N,
                     A_plug,
                     B_plug,
                     lm_plug,
                     m_lm
                     #m_lmA
    ) %>%
      as.data.frame() %>%
      rownames_to_column('name')  %>%
      mutate(group = 'cdf') 
    
    
    # Quantile Estimation
    
    ## A 
    A_df = data.frame(probs = ecdf_A(A$y),
                      t = A$y) %>% 
      arrange(probs)
    
    # B
    B_df = data.frame(probs = ecdf_B(B$y),
                      t = B$y) %>%
      arrange(probs)
    
    
    # P
    
    plm_df = data.frame(probs = ecdf_lm(lm_pred_A),
                        t = lm_pred_A) %>%
      arrange(probs)
    
    
    # R
    q_mlm = function(a, b, e_ecdf = e_ecdf){
      return(e_ecdf(a - b) %>% mean())
    }
    
    vq = Vectorize(q_mlm, 
                   vectorize.args = 'a')
    
    mlm_df = data.frame(probs = vq(lm_pred_A, lm_pred_A, e_ecdf = e_ecdf),
                        t = lm_pred_A) %>%
      arrange(probs)
    
    ie_fun = function(dats, a){
      
      return(dats %>%
               lapply(dplyr::filter, probs >= a ) %>%
               lapply(function(x) if(nrow(x) == 0) data.frame(probs = NA,
                                                              t = NA) else x) %>%
               lapply(setNames, c('probs', 't')) %>%
               lapply(dplyr::arrange, t) %>%
               lapply(dplyr::slice, 1) %>%
               lapply(dplyr::select, -probs) %>%
               bind_rows() %>%
               mutate(name = c(
                 'A_plug',
                 'B_plug',
                 'lm_plug',
                 'm_lm') %>% as.factor(),
                 a = a ) %>%
               dplyr::select(a, name, t))
      
    }
    
    ie_result = mapply(FUN = ie_fun,
                       a = F_N,
                       MoreArgs = list(dats =list(
                         A_df,
                         B_df,
                         plm_df,
                         mlm_df)),
                       SIMPLIFY=FALSE)
    
    final_quant=  
      ie_result %>%
      bind_rows() %>%
      mutate(F_N = paste0(a*100, '%')) %>%
      dplyr::select(-a) %>%
      pivot_wider(names_from = 'F_N',
                  values_from = 't') %>% 
      mutate(group = 'q') %>%
      as.data.frame() 
    
    # pool results together
    
    final_df =rbind(final_cdf,
                    final_quant
    ) %>%
      as.data.frame() 
    
    return(final_df)
  }
  
  samp.fours = function(samp, nB, r, pop, p, nA){
    
    B_MAR = B_generator(
      miss = 'MAR', nB = nB, r = r)
    B_MNAR = B_generator(
      miss = 'MNAR', nB = nB, r = r)
    
    
    cdf_MAR = MI_cdf(B = B_MAR, A = A, p = p)
    cdf_MNAR = MI_cdf(B = B_MNAR, A = A, p = p) 
    
    
    MC_results = rbind(cdf_MAR,
                       cdf_MNAR) %>%
      mutate(miss = rep(c('MAR', 'MNAR'), each = nrow(.)/2)) %>%
      mutate(nB = nB, nA = nA, nsim = nsim, r=r, mod = mod# to track errors
      )
    
    return(MC_results)
  }
  
  # parallel processing 
  iter = 1:nsim
  numbcores = detectCores()
  RNGkind("L'Ecuyer-CMRG")
  set.seed(seed)
  
  results_r_nB = c()
  for(i in 1:length(nB)){
    results_r_nB[[i]] = pbmclapply(X = iter,  r = r,
                                   nB = nB[i], 
                                   nA = nA,
                                   pop = pop,
                                   p =p, 
                                   FUN = samp.fours, 
                                   mc.cores = 17)
    progress(i, length(nB))
    
  }
  
  results = results_r_nB %>% bind_rows()
  
  
  results_just_A = 
    results %>%
    bind_rows(.id = 'iter') %>%
    mutate(est = name) %>%
    dplyr::select(-name) %>%
    filter(est != 'F_N') %>%
    dplyr::select(iter, nB, nA, nsim, r, mod,# to track errors
                  group, est, miss, everything()) %>%
    pivot_longer(10:ncol(.), names_to = 'perc') %>%
    group_by(r, nB, group, miss, est, perc) %>%
    left_join(data.frame(
      group = rep(c('cdf', 'q'), each = length(F_N)), 
      pop_val = c(F_N, q)
    ) %>% 
      mutate(perc = paste0(F_N*100, '%')),
    by = c('group', 'perc')
    ) %>%
    filter(est == 'A_plug') %>%
    group_by(r, nB, group, miss, perc, est) %>%
    dplyr::summarize(rmse_A = sqrt(mean((value - pop_val)^2, na.rm = TRUE)),
                     bias_A = mean(value, na.rm = TRUE) - mean(pop_val, na.rm=TRUE),
                     var_A = var(value, na.rm = TRUE)) %>%
    dplyr::select(-est) %>%
    ungroup()
  
  results_not_A = 
    results %>%
    bind_rows(.id = 'iter') %>%
    mutate(nB = nB, nA = nA, nsim = nsim, r=r, mod = mod) %>%
    mutate(est = name) %>%
    dplyr::select(-name) %>%
    filter(est != 'F_N') %>%
    dplyr::select(iter, nB, nA, nsim, r, mod,# to track errors
                  group, est, miss, everything())  %>%
    pivot_longer(10:ncol(.), names_to = 'perc') %>%
    group_by(r, nB, group, miss, est, perc) %>%
    left_join(data.frame(
      group = rep(c('cdf', 'q'), each = length(F_N)), 
      pop_val = c(F_N, q)
    ) %>% 
      mutate(perc = paste0(F_N*100, '%')),
    by = c('group', 'perc')
    ) %>%
    dplyr::select(nB, r, perc, group, miss, est, value, pop_val) %>%
    filter(est != 'A_plug') %>%
    group_by(nB, r, group, miss, perc, est) %>%
    dplyr::summarize(rmse = sqrt(mean((value - pop_val)^2, na.rm =TRUE)),
                     bias = mean(value, na.rm = TRUE) - mean(pop_val, na.rm = TRUE),
                     var = var(value, na.rm = TRUE)) %>%
    ungroup()
  
  
  perf_results = results_just_A %>%
    left_join(results_not_A,
              by = c('nB', 'r', 'perc', 'group', 'miss')
    ) %>%
    mutate(rrmse = rmse/rmse_A,
           rbias = bias/bias_A,
           rvar = var/var_A) %>%
    dplyr::select(nB, r, group, miss, perc, est, rrmse, rbias, rvar) %>%
    arrange(group) %>%
    mutate(mod = mod)
  
  
  return(list(results, perf_results))
}

# example run
nsim = 1500
N = 100000
nA = .01*N
nB = c(nA/N, .10, .20)*N
seed=101
r = .15


setwd("/Users/jeremyflood/OneDrive/Documents/Grad School/2024-2025/Spring 2025/reCDF Revisions/MC Code")

# f1
mod = 'f1'
f1_model = sim.f2(N = N,
                  nA = nA,
                  nB = nB,
                  nsim = nsim,
                  mod = 'f1',
                  seed = 101,
                  r = r)

f1_model %>%
  setNames(c('raw', 'perf')) %>%
  openxlsx::write.xlsx('modf1_results.xlsx')

# f2
mod = 'f2'
f2_model = sim.f2(N = N,
                  nA = nA,
                  nB = nB,
                  nsim = nsim,
                  mod = 'f2',
                  seed = 101,
                  r = r)

f2_model %>%
  setNames(c('raw', 'perf')) %>%
  openxlsx::write.xlsx('modf2_results.xlsx')

# f3
mod = 'f3'
f3_model = sim.f2(N = N,
                  nA = nA,
                  nB = nB,
                  nsim = nsim,
                  mod = 'f3',
                  seed = 101,
                  r = r)

f3_model %>%
  setNames(c('raw', 'perf')) %>%
  openxlsx::write.xlsx('modf3_results.xlsx')

# f4
mod = 'f4'
f4_model = sim.f2(N = N,
                  nA = nA,
                  nB = nB,
                  nsim = nsim,
                  mod = 'f4',
                  seed = 101,
                  r = r)

f4_model %>%
  setNames(c('raw', 'perf')) %>%
  openxlsx::write.xlsx('modf4_results.xlsx')

setwd("/Users/jeremyflood/OneDrive/Documents/Grad School/2024-2025/Spring 2025/reCDF Revisions/MC Code")

mod_list = paste0('modf', 1:4, '_results.xlsx')


dat_together = c()

for(i in 1:length(mod_list)){
  dat_together[[i]] = openxlsx::read.xlsx(mod_list[i], sheet = 'perf') %>% mutate(mod = paste0('f', i))
}

est_label =  c('B_plug_cdf' = TeX("$\\hat{F}_{B}$"),
               'lm_plug_cdf' = TeX("$\\hat{F}_{P}"),
               'm_lm_cdf' = TeX("$\\hat{F}_{R}"), 
               'B_plug_q' = TeX("$\\hat{t}_{B}"),
               'lm_plug_q' = TeX("$\\hat{t}_{P}"),
               'm_lm_q' = TeX("$\\hat{t}_{R}")
)


plots_together = c()

for (i in 1:length(nB)){
  plots_together[[i]] = dat_together %>%
    bind_rows() %>%
    filter(nB == nB[i]) %>%
    #filter(est != 'B_plug') %>%
    #filter(miss == 'MAR') %>%
    mutate(new_est = paste0(est, '_', group)) %>%
    dplyr::select(nB, mod, miss, perc, new_est, est, group, rrmse) %>%
    mutate(est = factor(est, labels = c('B_plug' = 'Naive', 'lm_plug' = 'Plug-in', 'm_lm' = 'Residual eCDF'))) %>%
    mutate(mod = factor(mod,
                        labels= c('f1' = TeX("Model $\\xi_{1}"),
                                  'f2' = TeX("Model $\\xi_{2}"),
                                  'f3' = TeX("Model $\\xi_{3}"),
                                  'f4' = TeX("Model $\\xi_{4}")))
    ) %>%
    mutate(cat_est = ifelse(new_est %in% c('B_plug_cdf', 'lm_plug_cdf', 'm_lm_cdf'), 
                            "F(t)", "t(a)")) %>%
    mutate(new_est = factor(new_est, levels = c('B_plug_cdf', 'lm_plug_cdf', 'm_lm_cdf',
                                                'B_plug_q', 'lm_plug_q', 'm_lm_q'),
                            labels = c('B_plug_cdf' = TeX("$\\hat{F}_{B}$"),
                                       'lm_plug_cdf' = TeX("$\\hat{F}_{P}"),
                                       'm_lm_cdf' = TeX("$\\hat{F}_{R}"),
                                       'B_plug_q' = TeX("$\\hat{t}_{B}"),
                                       'lm_plug_q' = TeX("$\\hat{t}_{P}"),
                                       'm_lm_q' = TeX("$\\hat{t}_{R}")))) %>%
    # mutate(nB = factor(nB, 
    #                    labels = c('1000' = TeX("$n_B = n_A$"),
    #                               '10000' = TeX("$n_B = 10n_A$"),
    #                               '20000' = TeX("$n_B = 20n_A$")))) %>%
    mutate(perc_strip = as.numeric(gsub(perc, pattern = '%', replacement = ''))/100) %>%
    mutate(perc_strip =factor(perc_strip)) %>%
    mutate(Missingness = factor(miss)) %>%
    arrange(nB, perc) %>%
    
    ggplot(aes(x = perc_strip, y = rrmse^(1/2), color =est, group = est)) +
    geom_line(linewidth = .5) +
    geom_point(size = 1) +
    facet_grid(scales ='fixed', rows = dplyr::vars(cat_est, Missingness), cols = dplyr::vars(mod),
               labeller = label_parsed) +
    # scale_x_discrete(labels = c('B_plug_cdf' = TeX("$\\hat{F}_{B}$"),
    #                             'lm_plug_cdf' = TeX("$\\hat{F}_{P}"),
    #                             'm_lm_cdf' = TeX("$\\hat{F}_{R}"), 
    #                             'B_plug_q' = TeX("$\\hat{t}_{B}"),
    #                             'lm_plug_q' = TeX("$\\hat{t}_{P}"),
    #                             'm_lm_q' = TeX("$\\hat{t}_{R}")
    # )
    # ) + 
    ylab(TeX("\\sqrt{RMSER}")) +
    xlab('Percentile') +
    theme_bw() +
    theme(text = element_text(size = 15),
          axis.text.x = element_text(angle = 90, hjust = 1)) +
    theme(legend.position = 'top') +
    ggtitle(ifelse(nB[i] == nA, TeX('\\sqrt{RMSER} for $n_{B} = n_{A}$'), 
                   TeX(paste0('\\sqrt{RMSER} for $n_{B} =') %>% paste0(nB[i]/nA) %>% paste0('nA$')))) +
    labs(colour = NULL)
}

setwd("/Users/jeremyflood/OneDrive/Documents/Grad School/2024-2025/Spring 2025/reCDF Revisions/MC Code/Plots/Lower Quality")

ggsave('RRMSE_nB_1.png', plots_together[[1]], dpi=400, width = 10, height = 8)
ggsave('RRMSE_nB_2.png', plots_together[[2]], dpi=400, width = 10, height = 8)
ggsave('RRMSE_nB_3.png', plots_together[[3]], dpi=400, width = 10, height = 8)


# for generating var plots

setwd("/Users/jeremyflood/OneDrive/Documents/Grad School/2024-2025/Spring 2025/reCDF Revisions/Variance Estimation/New Try/Split Jobs/Var Results")


mod_list = c(
  paste0('f1_', 1:3, '.xlsx'),
  paste0('f3_', 1:4, '.xlsx')
)


results = c()

for(i in 1:length(mod_list)){
  results[[i]] = openxlsx::read.xlsx(mod_list[i], sheet = 'raw')
}


cleaned_ddf =  results %>%
  bind_rows() %>%
  mutate(CR = ifelse(pop_quant >= LL & pop_quant <= UL, 1, 0)) 


# mutate(var_name = var_type) %>%
# dplyr::select(-var_type) #%>%
# mutate(est_name= ifelse(is.na(m_lm) == is.na(F_N) & is.na(F_N) == TRUE, 't', 'cdf'),
#        est_value = ifelse(est_name == 'cdf', m_lm, hatq),
#        actual_value = ifelse(est_name == 'cdf', F_N, q_N)) %>%
# dplyr::select(F_N, q_N, hatq, m_lm, everything())

mc_results = cleaned_ddf %>%
  #filter(var_name == 'asymp') %>%
  group_by(mod, perc, miss, nB, est_type, var_type) %>% 
  dplyr::summarize(MC_se = sqrt(var(est_quant)))


cleaned_results = 
  cleaned_ddf %>%
  #left_join(mc_results) %>%
  dplyr::select(est_type, est_quant, pop_quant, bootstats, everything()) %>%
  pivot_longer(cols = 2) %>%
  group_by(mod, est_type, perc, miss, nB, var_type, name, bootstats) %>%
  dplyr::summarize(
    pop_quant = mean(pop_quant),
    est_value = mean(value),
    est_se = mean(sqrt(var_val)),
    CR = mean(CR),
    d = mean(UL - LL),
    boot_mean_val = mean(bootstats_value)) %>%
  ungroup() %>%
  left_join(mc_results, by = c('mod', 'est_type', 'perc', 'miss', 'nB', 'var_type')) %>%
  mutate(rb = (-(MC_se - est_se)/MC_se)*100) %>%
  mutate(nsim = 1500, L = 750) %>%
  filter(bootstats == 'boot_mean_est') 


FN_plots = c()
mods = c(1, 3)

setwd("/Users/jeremyflood/OneDrive/Documents/Grad School/2024-2025/Spring 2025/reCDF Revisions/Variance Estimation/New Try/Split Jobs/Var Results/Lower Quality Plots")

for(i in 1:length(mods)){
  FN_plots[[i]] = cleaned_results %>%
    dplyr::select(-c(name, bootstats, nsim, L, boot_mean_val, pop_quant,MC_se, est_value)) %>%
    filter(mod == paste0('f', mods[i]), est_type == 'cdf' ) %>% # for now
    pivot_longer(cols = 7:ncol(.), names_to = 'perf') %>%
    mutate(perf = factor(perf, levels = c('CR', 'd', 'est_se', 'rb'), labels = c(
      'CR' = TeX("\\% CR"),
      'd' = "AL",
      'est_se' = TeX("\\widehat{SE}"),
      'rb' = TeX("\\% RB")
    ))) %>%
    mutate(est_type = factor(est_type, levels = c('cdf', 't'),
                             labels= c('cdf' = TeX('$F(t)$'),
                                       't' = TeX('$t(\\alpha)$')
                             ))) %>%
    mutate(nB = factor(nB, 
                       levels = c('1000', '20000'),
                       labels = c('1000' = TeX('$n_{B} = n_A$'), '20000' = TeX('$n_{B} = 20n_A$')))) %>%
    mutate(var_name = factor(var_type, labels = c('asymp' = 'V1', 'boot' = 'V2'))) %>%
    ggplot(aes(x = perc, y = value, col = var_name, group = var_name)) +
    geom_line(linewidth = .5) +
    geom_point(size = 1) +
    scale_color_discrete(labels = unname(TeX(c("$V_{1}$ (Asymptotic)", "$V_{2}$ (Bootstrap)")))) +
    facet_grid(scales ='free', col = dplyr::vars(miss, nB), rows = dplyr::vars(perf),
               labeller = label_parsed) +
    theme_bw() +
    ylab("Performance Value") +
    xlab('Percentile') +
    theme(legend.position = 'top') +
    ggtitle(TeX(paste0('$F_{N}(t)$ Variance Estimation (Model $\\xi_', mods[i] , ')'))) +
    labs(colour = "Variance Type:")
  
  ggsave(paste0('var_cdf_f', mods[i], '.png'), FN_plots[[i]], dpi=800, width = 11, height = 8)
}

TN_plots = c()

for(i in 1:length(mods)){
  TN_plots[[i]] = cleaned_results %>%
    dplyr::select(-c(name, bootstats, nsim, L, boot_mean_val, pop_quant,MC_se, est_value)) %>%
    filter(mod == paste0('f', mods[i]), est_type == 't' ) %>% # for now
    pivot_longer(cols = 7:ncol(.), names_to = 'perf') %>%
    mutate(perf = factor(perf, levels = c('CR', 'd', 'est_se', 'rb'), labels = c(
      'CR' = TeX("\\% CR"),
      'd' = "AL",
      'est_se' = TeX("\\widehat{SE}"),
      'rb' = TeX("\\% RB")
    ))) %>%
    mutate(est_type = factor(est_type, levels = c('cdf', 't'),
                             labels= c('cdf' = TeX('$F(t)$'),
                                       't' = TeX('$t(\\alpha)$')
                             ))) %>%
    mutate(nB = factor(nB, 
                       levels = c('1000', '20000'),
                       labels = c('1000' = TeX('$n_{B} = n_A$'), '20000' = TeX('$n_{B} = 20n_A$')))) %>%
    mutate(var_name = factor(var_type, labels = c('asymp' = 'V1', 'boot' = 'V2'))) %>%
    ggplot(aes(x = perc, y = value, col = var_name, group = var_name)) +
    geom_line(linewidth = .5) +
    geom_point(size = 1) +
    scale_color_discrete(labels = unname(TeX(c("$V_{1}$ (Asymptotic)", "$V_{2}$ (Bootstrap)")))) +
    facet_grid(scales ='free', col = dplyr::vars(miss, nB), rows = dplyr::vars(perf),
               labeller = label_parsed) +
    theme_bw() +
    ylab("Performance Value") +
    xlab('Percentile') +
    theme(legend.position = 'top') +
    ggtitle(TeX(paste0('$t_{N}(\\alpha)$ Variance Estimation (Model $\\xi_', mods[i] , ')'))) +
    labs(colour = "Variance Type:")
  
  ggsave(paste0('var_t_f', mods[i], '.png'), TN_plots[[i]], dpi=800, width = 11, height = 8)
}



