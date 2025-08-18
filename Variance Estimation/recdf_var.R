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



L_param = 750
L_parm = 750
L = L_parm
seed = 101
nsim = 300
N = 1e05
nA = 1000
alp = .10
mod = 'f3'
r =.15

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


# MAR strat
MAR_strat = paste0('x', c(cor(X,y) %>% which.max()))


pop = as.data.frame(cbind(X, y)) %>%
  mutate(mnar_strat = ifelse(y < median(y), "0", "1") %>% as.factor() 
         
  ) %>%
  arrange(mnar_strat)


# population characteristics
mu = mean(pop$y)
F_N = c(.01, .10, .25, .50, .75, .90, .99)
q= q_N = quantile(pop$y, probs= F_N, type = 1)

rm(seed)

seed = 101 + as.numeric(format(Sys.time(), "%S")) + as.numeric(format(Sys.time(), "%M")) + as.numeric(format(Sys.time(), "%H"))


As = lapply(1:nsim, function(seed){
  set.seed(seed)
  pop[sample(1:N, nA, replace = FALSE),]
})

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

# B

# sizes2 = round(c(r*nB, (1-r)*nB), 
#                digits =0)
# 
# sizes = ifelse(sizes2 <1, 1, sizes2)

## MAR

Bs_MAR_1nA = mclapply(1:nsim, function(seed){
  B_generator(miss = 'MAR', nB = 1000, r = r)
}, mc.cores = 7)

Bs_MAR_20nA = mclapply(1:nsim, function(seed){
  B_generator(miss = 'MAR', nB = 20*nA, r = r)
}, mc.cores = 7)



## MNAR

Bs_MNAR_1nA = mclapply(1:nsim, function(seed){
  B_generator(miss = 'MNAR', nB = 1000, r = r)
}, mc.cores = 7)

Bs_MNAR_20nA = mclapply(1:nsim, function(seed){
  B_generator(miss = 'MNAR', nB = 20*nA, r = r)
}, mc.cores = 7)


# B_perm_MAR_1nA= Bs_MAR_1nA[[1]]
# B_perm_MAR_20nA = Bs_MAR_20nA[[1]]
# B_perm_MNAR_1nA = Bs_MNAR_1nA[[1]]
# B_perm_MNAR_20nA = Bs_MNAR_20nA[[1]]



samp.f = function(samp, A, B_MAR_nA, B_MAR_20nA, B_MNAR_1nA, B_MNAR_20nA, r, pop, p, miss, L = L_parm, mc_cores = 15, alp = .10, seed){
  options(dplyr.summarise.inform = FALSE)
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
  
  var_cdf = function(B, A, p, alp, da_miss){
    nA = nrow(A)
    nB = nrow(B)
    # defining A
    ecdf_A = ecdf(A$y)
    ecdf_B = ecdf(B$y)
    
    #lmod
    lm_B = lm(y~.,
              data =B[, c(paste0('x', 1:p),
                          'y')])
    lm_pred_B = predict(lm_B, newdata = B)
    lm_pred_A = predict(lm_B, newdata = A)
    
    pA_pred= svydesign(id = ~1,
                       weights = ~weight,
                       fpc=~fpc,
                       data = A %>%
                         mutate(fpc = nrow(A)/sum(A$weight),
                                y = lm_pred_A))
    ecdf_lm = svycdf(~y, pA_pred)
    
    
    # 
    # # plug-in
    
    A_plug = ecdf_A(q)
    B_plug = ecdf_B(q)
    lm_plug = ecdf_lm[[1]](q)
    names(lm_plug) = names(A_plug) = names(lm_plug) =  names(q)
    
    
    # reCDF
    v_lm = B$y - lm_pred_B
    e_ecdf= ecdf(v_lm)
    
    
    
    m_lm = vars_to= c()
    
    for(i in 1:length(q)){
      m_lm[i] = (1/sum(A$weight))*sum(A$weight*e_ecdf(q[i] - lm_pred_A))
      vars_to[i] = (1-nA/N)/nA* var(e_ecdf(q[i] - lm_pred_A))
    }
    
    names(vars_to) = names(m_lm) = names(q)
    
    
    
    final_cdf =rbind(F_N,
                     A_plug,
                     B_plug,
                     lm_plug,
                     m_lm) %>%
      as.data.frame() %>%
      rownames_to_column('name')  %>%
      mutate(group = 'cdf') %>%
      mutate(val = 'mc')
    
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
    
    plm_df = data.frame(probs = ecdf_lm[[1]](lm_pred_A),
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
    
    
    
    final_quant=  ie_result %>%
      bind_rows() %>%
      mutate(F_N = paste0(a*100, '%')) %>%
      dplyr::select(-a) %>%
      pivot_wider(names_from = 'F_N',
                  values_from = 't') %>%
      mutate(group = 'q', val = 'mc') %>%
      as.data.frame() %>%
      mutate(
        group = 'q',
        val = 'mc') %>%
      rbind(q%>%
              t() %>%
              as.data.frame %>%
              mutate(name = 'F_N',
                     group = 'q',
                     val = 'mc'))
    
    
    
    # organizing variances
    
    LL.q= UL.q = q_var = hatq = c()
    for(i in 1:length(F_N)){
      a.a = F_N[i]
      hatq[i] = ie_result %>%
        bind_rows() %>%
        filter(name == 'm_lm') %>%
        filter(a == a.a) %>%
        dplyr::select(t)%>%
        unlist() %>%
        unname()
      
      vars_q = (1-nA/N)/nA* var(e_ecdf(hatq[i]- lm_pred_A))
      
      LL.q[i] =  list(mlm_df) %>%
        lapply(dplyr::filter, probs >= a.a - qnorm(1-alp/2)*sqrt(vars_q)) %>%
        lapply(function(x) if(nrow(x) == 0) data.frame(probs = NA,
                                                       t = NA) else x) %>%
        lapply(setNames, c('probs', 't')) %>%
        lapply(dplyr::arrange, t) %>%
        lapply(dplyr::slice, 1) %>%
        bind_rows() %>%
        mutate(a = a.a) %>%
        dplyr::select(t) %>%
        unlist() %>%
        unname()
      
      UL.q[i]=  list(mlm_df) %>%
        lapply(dplyr::filter, probs >= a.a + qnorm(1-alp/2)*sqrt(vars_q)) %>%
        lapply(function(x) if(nrow(x) == 0) data.frame(probs = NA,
                                                       t = NA) else x) %>%
        lapply(setNames, c('probs', 't')) %>%
        lapply(dplyr::arrange, t) %>%
        lapply(dplyr::slice, 1) %>%
        bind_rows() %>%
        mutate(a = a.a) %>%
        dplyr::select(t) %>%
        unlist() %>%
        unname()
      
      q_var[i] = vars_q
    }
    
    names(LL.q) = names(UL.q) = names(hatq) = names(q)
    
    q_var_df = rbind(q, hatq,  q_var, LL.q, UL.q) %>%
      as.data.frame() %>%
      rownames_to_column('name') %>%
      mutate(group = 'q', val = 'var')
    
    nB = nrow(B)
    
    final_df =rbind(final_cdf,
                    final_quant,
                    q_var_df
    ) %>%
      filter(name %in% c(
        'm_lm',
        'hatq',
        'q_var',
        'LL.q',
        "UL.q"
      )) %>%
      as.data.frame() %>%
      dplyr::slice(-2) %>%
      rbind(t(vars_to) %>%
              as.data.frame() %>%
              mutate(name = 'm_lm_var',
                     group = 'q',
                     val = 'mc')%>%
              dplyr::select(colnames(final_cdf))) %>%
      dplyr::select(-c(group, val)) %>%
      mutate(nB = nB, miss = da_miss)
    
    
    
    
    
    return(final_df)
  }
  # selecting A
  
  
  A_perm = A %>%
    dplyr::select(paste0('x', 1:p),  'y') %>%
    mutate(Prob = nA/N,
           weight = N/nA,
           fpc = Prob)
  
  B_perm_MAR_1nA= B_MAR_nA
  B_perm_MAR_20nA = B_MAR_20nA
  B_perm_MNAR_1nA = B_MNAR_1nA
  B_perm_MNAR_20nA = B_MNAR_20nA
  
  
  mydesign = svydesign(ids = ~0,
                       weights = ~weight,
                       fpc = ~fpc,
                       data = A_perm)
  
  
  bootstrap_rep_design = as_bootstrap_design(mydesign,
                                             type = "Rao-Wu-Yue-Beaumont",
                                             replicates = L,
                                             samp_method_by_stage = 'SRSWOR')
  
  A_repweights_perm = as_data_frame_with_weights(bootstrap_rep_design,
                                                 full_wgt_name = "FULL_SAMPLE_WGT",
                                                 rep_wgt_prefix = "REP_WGT_"  )  %>%
    dplyr::select(-weight)
  #
  #
  #
  A_list =c()
  for(i in 1:L){
    A_list[[i]] = A_repweights_perm %>% dplyr::select(y, paste0("x", 1:p), paste0("REP_WGT_", i), Prob) %>%
      dplyr::rename(weight = paste0("REP_WGT_", i))
  }
  #  
  #  # selecting B
  #  
  #  
  
  set.seed(seed)
  B_list.MAR_1nA =mc_replicate(L, list(B_perm_MAR_1nA[sample(1:1000,
                                                             size = 1000,
                                                             replace = TRUE),]), mc.cores = mc_cores)
  set.seed(seed)
  B_list.MAR_20nA =mc_replicate(L, list(B_perm_MAR_20nA[sample(1:(20*1000),
                                                               size = 20*1000,
                                                               replace = TRUE),]), mc.cores = mc_cores)
  set.seed(seed)
  B_list.MNAR_1nA =mc_replicate(L, list(B_perm_MNAR_1nA[sample(1:1000,
                                                               size = 1000,
                                                               replace = TRUE),]), mc.cores = mc_cores)
  set.seed(seed)
  B_list.MNAR_20nA =mc_replicate(L, list(B_perm_MNAR_20nA[sample(1:(20*1000),
                                                                 size = 20*1000,
                                                                 replace = TRUE),]), mc.cores = mc_cores)
  #  
  #  
  #  
  #  
  #  
  #  # list of A and list of B are selected. Now for each a in A and b in B, repeat technique found in original paper
  #  
  results.MAR_1nA=  mcmapply(FUN = var_cdf,
                             A = A_list,
                             B = B_list.MAR_1nA,
                             MoreArgs = list(p= p, alp= .10, da_miss = 'MAR'),
                             mc.cores = mc_cores,
                             SIMPLIFY = FALSE)
  
  results.MAR_20nA=  mcmapply(var_cdf,
                              A = A_list,
                              B = B_list.MAR_20nA,
                              MoreArgs = list('p' = p, 'alp' = .10, da_miss = 'MAR'),
                              mc.cores = mc_cores,
                              SIMPLIFY = FALSE)
  
  results.MNAR_1nA=  mcmapply(var_cdf,
                              A = A_list,
                              B = B_list.MNAR_1nA,
                              MoreArgs = list('p' = p, 'alp' = .10, da_miss = 'MAR'),
                              mc.cores = mc_cores,
                              SIMPLIFY = FALSE)
  
  results.MNAR_20nA=  mcmapply(var_cdf,
                               A = A_list,
                               B = B_list.MNAR_20nA,
                               MoreArgs = list('p' = p, 'alp' = .10, da_miss = 'MAR'),
                               mc.cores = mc_cores,
                               SIMPLIFY = FALSE)
  
  results2 = list(results.MAR_1nA %>% bind_rows() %>% mutate(nB = nA, miss = 'MAR'),
                  results.MAR_20nA %>% bind_rows() %>% mutate(nB = 20*nA, miss = 'MAR'),
                  results.MNAR_1nA %>% bind_rows() %>% mutate(nB = nA, miss = 'MNAR'),
                  results.MNAR_20nA %>% bind_rows() %>% mutate(nB = 20*nA, miss = 'MNAR')) 
  
  
  var_calcs = function(A_perm, B_perm, L, alp, da_miss, da_nB){
    actual =  var_cdf(
      A = A_perm,
      B = B_perm,
      p = p,
      alp = .10,
      da_miss = da_miss)
    
    standard_values =
      actual %>%
      filter(name %in% c('m_lm',
                         'hatq')) %>%
      dplyr::select(nB, miss, name, everything()) %>%
      pivot_longer(cols =4:ncol(.),
                   names_to = 'perc',
                   values_to = 'standard')
    
    boot_var = results2 %>%
      bind_rows(.id = 'iter') %>%
      filter(miss == da_miss, nB == da_nB) %>%
      filter(name %in% c('m_lm',
                         'hatq')) %>%
      dplyr::select(iter, name,miss, nB, everything()) %>%
      pivot_longer(cols = 5:ncol(.),
                   names_to = 'perc',
                   values_to = 'estimate') %>%
      left_join(standard_values,
                by = c('nB', 'miss', 'name', 'perc')) %>%
      group_by(nB, miss, name, perc) %>%
      dplyr::summarize(boot_mean_est = mean(estimate),
                       B_perm_est = mean(standard),
                       boot_var = (1/L)*sum_((estimate - standard)^2))  %>%
      dplyr::select(nB, miss, name, perc, boot_var, everything()) %>%
      pivot_longer(cols = 6:ncol(.), names_to = 'bootstats') %>%
      pivot_wider(values_from = c('boot_var', 'value')) %>%
      pivot_longer(cols = 7:8, names_to = 'bootstats_est', values_to = 'bootstats_value') %>%
      mutate(bootstats_est = gsub(x = bootstats_est, pattern = 'value_', replacement = '')) %>%
      dplyr::select(nB, miss, perc, boot_var_hatq, boot_var_m_lm, bootstats, bootstats_est, bootstats_value, everything())
    #
    #  # boot_var  =
    #  #   results2 %>%
    #  #   bind_rows(.id = 'iter') %>%
    #  #   filter(name %in% c('m_lm',
    #  #                      'hatq')) %>%
    #  #   dplyr::select(iter, name, miss, nB, everything()) %>%
    #  #   pivot_longer(cols = 5:ncol(.),
    #  #                names_to = 'perc',
    #  #                values_to = 'estimate') %>%
    #  #   filter(nB == da_nB, miss == da_miss) %>%
    #  #   left_join(standard_values,
    #  #             by = c('name', 'perc', 'miss', 'nB')) %>%
    #  #   group_by(name, perc, miss, nB) %>%
    #  #   dplyr::summarize(boot_var = (1/L)*sum_((estimate - standard)^2)) %>%
    #  #   pivot_wider(values_from = 'boot_var') %>%
    #  #   setNames(c('perc',
    #  #              'miss',
    #  #              'nB',
    #  #              'bs_hatq',
    #  #              'bs_mlm'))
    #  #
    #
    
    pop_quant = data.frame(F_N = F_N, q_N = q, perc = paste0(F_N*100, '%'))
    
    
    bigguy = actual %>%
      dplyr::select(nB, miss, name, everything()) %>%
      pivot_longer(cols = 4:ncol(.), names_to = 'perc') %>%
      pivot_wider() %>%
      left_join(boot_var, by = c('nB', 'miss',  'perc')) %>%
      # dplyr::select(nB, miss, name, everything()) %>%
      # pivot_longer(cols = 4:ncol(.),
      #              names_to = 'perc') %>%
      # pivot_wider() %>%
      # left_join(boot_var, by = c('perc', 'miss', 'nB'))  %>%
      left_join(pop_quant, by = 'perc')
    
    
    mlm_varsum =
      bigguy %>%
      dplyr::select(F_N, nB, miss, perc,
                    m_lm,
                    bootstats,
                    bootstats_est,
                    bootstats_value,
                    m_lm_var,
                    boot_var_m_lm) %>%
      setNames(c('pop_quant', 'nB', 'miss', 'perc',
                 'est_quant',
                 'bootstats', 'bootstats_est', "bootstats_value",
                 'asymp',
                 'boot')) %>%
      pivot_longer(cols = 9:ncol(.),
                   names_to = 'var_type',
                   values_to = 'var_val') %>%
      rowwise() %>%
      mutate(LL = est_quant - qnorm(p = 1-alp/2)*sqrt(var_val),
             UL = est_quant + qnorm(p = 1-alp/2)*sqrt(var_val)) %>%
      mutate(CR = ifelse(F_N >= LL & F_N <= UL, 1, 0)) %>%
      dplyr::select(nB, miss,
                    perc,
                    pop_quant,
                    est_quant,
                    var_type,
                    CR,
                    everything()
      ) %>%
      ungroup() %>%
      mutate(est_type = 'cdf') %>%
      filter(bootstats_est == 'm_lm') %>%
      dplyr::select(-bootstats_est)
    
    # q_var is NOT the variance of qhat -- bad naming, Jeremy!! q_var is the estimated variance of F(hat(T))
    
    hatq_varsum =
      bigguy %>%
      dplyr::select(q_N, LL.q, UL.q, nB, miss, perc,
                    hatq,
                    bootstats,
                    bootstats_est,
                    bootstats_value,
                    boot_var_hatq) %>%
      setNames(c('pop_quant', 'LL.q', 'UL.q', 'nB', 'miss', 'perc',
                 'est_quant', 'bootstats', 'bootstats_est', 'bootstats_value',
                 'boot')) %>%
      mutate(asymp = (UL.q - LL.q)^2 / (2*qnorm(1-(alp/2)))) %>%
      dplyr::select(pop_quant, LL.q, UL.q, nB, miss, perc,  est_quant,
                    bootstats, bootstats_est, bootstats_value,
                    boot, asymp) %>%
      pivot_longer(cols = 11:ncol(.),
                   names_to = 'var_type',
                   values_to = 'var_val') %>%
      rowwise() %>%
      mutate(LL = ifelse(var_type == 'boot',
                         est_quant - qnorm(1-alp/2)*sqrt(var_val),
                         LL.q
      ),
      UL = ifelse(var_type == 'boot',
                  est_quant +qnorm(1-alp/2)*sqrt(var_val),
                  UL.q
      )) %>%
      dplyr::select(-c(LL.q, UL.q)) %>%
      # ungroup() %>%
      # left_join(bigguy %>%
      #             dplyr::select(perc, LL.q, UL.q) %>%
      #             mutate(var_type = 'asymp') %>%
      #             dplyr::select(perc, var_type, everything()) %>%
      #             setNames(c('perc', 'var_type', 'LL_other','UL_other')),
      #           by = c('perc', 'var_type')) %>%
      # rowwise() %>%
      # mutate(LL = ifelse(is.na(LL) == FALSE,
      #                    LL, LL_other
      # ),
      # UL =ifelse(is.na(UL) == FALSE,
      #            UL, UL_other)
      # ) %>%
      # dplyr::select(-c(LL_other,
      #                  UL_other)) %>%
      mutate(CR = ifelse(pop_quant >= LL & pop_quant <= UL, 1, 0)) %>%
      ungroup() %>%
      mutate(est_type = 't') %>%
      filter(bootstats_est == 'hatq') %>%
      dplyr::select(-bootstats_est)
    
    final = rbind(mlm_varsum, hatq_varsum) %>%
      dplyr::select(est_type, nB, miss, perc, everything())
    return(final)
  }
  
  final_results = list(var_calcs(A_perm,B_perm_MAR_1nA, L, alp, da_miss = 'MAR', da_nB = 1*nA),
                       var_calcs(A_perm,B_perm_MAR_20nA, L, alp, da_miss = 'MAR', da_nB = 20*nA),
                       var_calcs(A_perm,B_perm_MNAR_1nA, L, alp, da_miss = 'MNAR', da_nB = 1*nA),
                       var_calcs(A_perm,B_perm_MNAR_20nA, L, alp, da_miss = 'MNAR', da_nB = 20*nA)
  ) %>%
    lapply(bind_rows) %>%
    bind_rows() 
  
  return(final_results)
  
}

results = c()
for(i in 1:nsim){
  results[[i]] = samp.f(A = As[[i]],
                        B_MAR_nA = Bs_MAR_1nA[[i]], 
                        B_MAR_20nA = Bs_MAR_20nA[[i]], 
                        B_MNAR_1nA=Bs_MNAR_1nA[[i]], 
                        B_MNAR_20nA = Bs_MNAR_20nA[[i]], 
                        p = p,
                                        L = L, seed = seed)
  progress(i, nsim)
}



cleaned_ddf = results %>%
  bind_rows() #%>%
# mutate(var_name = var_type) %>%
# dplyr::select(-var_type) #%>%
# mutate(est_name= ifelse(is.na(m_lm) == is.na(F_N) & is.na(F_N) == TRUE, 't', 'cdf'),
#        est_value = ifelse(est_name == 'cdf', m_lm, hatq),
#        actual_value = ifelse(est_name == 'cdf', F_N, q_N)) %>%
# dplyr::select(F_N, q_N, hatq, m_lm, everything())

mc_results = cleaned_ddf %>%
  #filter(var_name == 'asymp') %>%
  group_by(perc, miss, nB, est_type, var_type) %>% 
  dplyr::summarize(MC_var = var(est_quant))


cleaned_results = cleaned_ddf %>%
  #left_join(mc_results) %>%
  dplyr::select(est_type, est_quant, pop_quant, bootstats, everything()) %>%
  pivot_longer(cols = 2) %>%
  group_by(est_type, perc, miss, nB, var_type, name, bootstats) %>%
  dplyr::summarize(
    pop_quant = mean(pop_quant),
    est_value = mean(value),
    est_var = mean(var_val),
    CR = mean(CR),
    d = mean(UL - LL),
    boot_mean_val = mean(bootstats_value)) %>%
  ungroup() %>%
  left_join(mc_results, by = c('est_type', 'perc', 'miss', 'nB', 'var_type')) %>%
  mutate(rb = ((MC_var - est_var)/MC_var)*100) %>%
  mutate(nsim = nsim, L = L_parm)

final_results = list(cleaned_ddf %>% mutate(mod = mod, seed = seed), 
                     cleaned_results %>% mutate(mod = mod, seed = seed))

names(final_results) = c('raw', 'summary')

random_val = rgamma(n =1, shape = 500, rate =2)

openxlsx::write.xlsx(final_results, paste0('modf33_results_', seed, '_', round(random_val, digits = 3) *1000, '_.xlsx'))
