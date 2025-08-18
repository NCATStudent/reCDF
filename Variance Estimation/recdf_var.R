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
require(xlsx)
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
require(hablar)
require(pbapply)

sim.var = function(N,
                  nA,
                  nB,
                  nsim,
                  mod, 
                  seed, 
                  r = .15,
                  alp = .10,
                  L = 1000){
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
  } else if(mod == 'f2'){
    set.seed(seed)
    x1 = runif(N, min = -1, max = 1)
    x2 = runif(N, min = -1, max = 1)
    x3 = runif(N, min = -1, max = 1)
    x4 = runif(N, min = -1, max = 1)
    X = data.frame(x1, x2, x3, x4)
    p = ncol(X)
    err = rnorm(N, mean =0, sd = sqrt(.5))
    y = -sin(x1)+x2^2+x3-exp(-x4^2)+err
  } else if(mod == 'f3'){
    set.seed(seed)
    x1= runif(N, min = 0, max = 4) 
    x2= runif(N, min = 0, max = 4) 
    x3= runif(N, min = 4, max = 8)
    x4 = runif(N, min = 4, max = 8)
    X = cbind(x1, x2, x3, x4) 
    p = ncol(X)
    y = rnorm(n = N,
              4*x1 + 4*x2 + 2*x3 + 2*x4 + 
                (x1+x2)^2 + (x3+x4)^2, sd = 17.5) 
    
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
  
  
  pop = as.data.frame(cbind(X, y))  %>%
    mutate(MAR_strat = ifelse(get(MAR_strat) <= median(get(MAR_strat)), 'L',
                              'H') %>% as.factor(),
           MNAR_strat =  ifelse(y <= median(y), 'L',
                                'H') %>% as.factor()
    )
  

  
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
  
  var_cdf = function(B, A, p){
    options(dplyr.summarise.inform = FALSE)
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
    
    
    
    # plug-in
    
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
      dplyr::select(-c(group, val))
    
    
    
    
    
    return(final_df)
  }
  As = lapply(1:nsim, function(seed){
    set.seed(seed)
    pop[sample(1:N, nA, replace = FALSE),]
  })
  
  # B
  
  sizes2 = round(c(r*nB, (1-r)*nB), 
                 digits =0)
  
  sizes = ifelse(sizes2 <1, 1, sizes2)
  
  ## MAR
  
  nl_MAR= table(pop$MAR_strat)[2]
  nh_MAR = table(pop$MAR_strat)[1]
  pop_MAR = pop %>% 
    rowwise() %>%
    mutate(Prob = ifelse(MAR_strat == 'L', sizes2[1]/nl_MAR,  sizes2[2]/nh_MAR 
    )) %>%
    ungroup() %>%
    mutate(index = 1:N)
  
  
  Bs_MAR = pbmclapply(1:nsim, function(seed){
    set.seed(seed)
    getdata(pop_MAR,strata(pop_MAR %>% arrange(MAR_strat),
                           stratanames='MAR_strat',
                           size = sizes,
                           method = 'srswor') )
  }, mc.cores = 17)
  
  
  
  
  ## MNAR
  
  nl_MNAR= table(pop$MNAR_strat)[2]
  nh_MNAR = table(pop$MNAR_strat)[1]
  pop_MNAR = pop %>% 
    rowwise() %>%
    mutate(Prob = ifelse(MNAR_strat == 'L', sizes2[1]/nl_MNAR,  sizes2[2]/nh_MNAR 
    )) %>%
    ungroup() %>%
    mutate(index = 1:N)
  
  Bs_MNAR = pbmclapply(1:nsim, function(seed){
    set.seed(seed)
    getdata(pop_MNAR,strata(pop_MNAR %>% arrange(MNAR_strat),
                            stratanames='MNAR_strat',
                            size = sizes,
                            method = 'srswor') )
  }, mc.cores = 17)
  
  
  
  
  samp.f = function(samp, A, B, r, pop, p, miss){
    options(dplyr.summarise.inform = FALSE)
    # selecting A
    
    
    A_perm = A %>%
      dplyr::select(paste0('x', 1:p),  'y') %>%
      mutate(Prob = nA/N,
             weight = N/nA,
             fpc = Prob)
    B_perm = B
    
    
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
    
    
    
    A_list =c()
    for(i in 1:L){
      A_list[[i]] = A_repweights_perm %>% dplyr::select(y, paste0("x", 1:p), paste0("REP_WGT_", i), Prob) %>%
        dplyr::rename(weight = paste0("REP_WGT_", i))
    }
    
    # selecting B
    
    
    
    B_list =replicate(L, B_perm[sample(1:nB, 
                                       size = nB, 
                                       replace = TRUE),], 
                      simplify = FALSE)
    
    
    
    
    
    # list of A and list of B are selected. Now for each a in A and b in B, repeat technique found in original paper
    
    results2 =  mcmapply(var_cdf,
                         A = A_list,
                         B = B_list,
                         MoreArgs = list('p' = p),
                         mc.cores = 20,
                         SIMPLIFY = FALSE)
    
    actual =  var_cdf(
      A = A_perm,
      B = B_perm,
      p = p)
    
    standard_values =  
      actual %>%
      filter(name %in% c('m_lm',
                         'hatq')) %>%
      dplyr::select(name, everything()) %>%
      pivot_longer(cols =2:ncol(.),
                   names_to = 'perc',
                   values_to = 'standard')
    
    boot_var  = 
      results2 %>%
      bind_rows(.id = 'iter') %>%
      filter(name %in% c('m_lm',
                         'hatq')) %>%
      dplyr::select(iter, name, everything()) %>%
      pivot_longer(cols = 3:ncol(.),
                   names_to = 'perc',
                   values_to = 'estimate') %>%
      left_join(standard_values,
                by = c('name', 'perc')) %>%
      group_by(name, perc) %>%
      dplyr::summarize(boot_var = (1/L)*sum_((estimate - standard)^2)) %>%
      pivot_wider(values_from = 'boot_var') %>%
      setNames(c('perc',
                 'bs_hatq',
                 'bs_mlm'))
    
    bigguy =  actual %>% 
      pivot_longer(cols = 2:ncol(.), 
                   names_to = 'perc') %>%
      pivot_wider() %>%
      left_join(boot_var, by = 'perc') %>%
      dplyr::select(perc,
                    m_lm,
                    hatq,
                    q_var,
                    m_lm_var,
                    bs_hatq,
                    bs_mlm,
                    LL.q,
                    UL.q)
    
    mlm_varsum = 
      bigguy %>%
      dplyr::select(perc, 
                    m_lm,
                    m_lm_var,
                    bs_mlm) %>%
      setNames(c('perc',
                 'm_lm',
                 'asymp',
                 'boot')) %>%
      pivot_longer(cols = 3:ncol(.),
                   names_to = 'var_type',
                   values_to = 'var_val') %>%
      rowwise() %>%
      mutate(LL = m_lm - qnorm(p = 1-alp/2)*sqrt(var_val),
             UL = m_lm + qnorm(p = 1-alp/2)*sqrt(var_val),
             F_N = as.numeric(gsub(perc,
                                   pattern = '%', 
                                   replacement = ''))/100) %>%
      mutate(CR = ifelse(F_N >= LL & F_N <= UL, 1, 0)) %>%
      dplyr::select(
        perc,
        F_N,
        m_lm,
        var_type,
        CR,
        everything()
      )
    
    hatq_varsum =
      bigguy %>%
      dplyr::select(perc,
                    hatq,
                    q_var,
                    bs_hatq) %>%
      setNames(c('perc',
                 'hatq',
                 'asymp',
                 'boot')) %>%
      mutate(q_N = q) %>%
      dplyr::select(perc, hatq, q_N, everything()) %>%
      pivot_longer(cols = 4:ncol(.),
                   names_to = 'var_type',
                   values_to = 'var_val') %>%
      rowwise() %>%
      mutate(LL = ifelse(var_type == 'boot', 
                         hatq - qnorm(1-alp/2)*sqrt(var_val),
                         NA
      ),
      UL = ifelse(var_type == 'boot', 
                  hatq +qnorm(1-alp/2)*sqrt(var_val),
                  NA
      )) %>%
      ungroup() %>%
      left_join(bigguy %>% 
                  dplyr::select(perc, LL.q, UL.q) %>% 
                  mutate(var_type = 'asymp') %>%
                  dplyr::select(perc, var_type, everything()) %>%
                  setNames(c('perc', 'var_type', 'LL_other','UL_other')),
                by = c('perc', 'var_type')) %>%
      rowwise() %>%
      mutate(LL = ifelse(is.na(LL) == FALSE,
                         LL, LL_other
      ),
      UL =ifelse(is.na(UL) == FALSE,
                 UL, UL_other)
      ) %>%
      dplyr::select(-c(LL_other,
                       UL_other)) %>%
      mutate(CR = ifelse(q_N >= LL & q_N <= UL, 1, 0)) %>% 
      dplyr::select(
        perc,
        q_N,
        hatq,
        var_type,
        CR,
        everything()
      ) %>%
      ungroup()
    
    
    return(list(mlm_varsum,
                hatq_varsum))
  }
  
  
  results_MAR =  pbmapply(samp.f,
                          A = As,
                          B=Bs_MAR,
                          MoreArgs = list(p = p,
                                          miss = 'MAR'),
                          SIMPLIFY = FALSE)
  
  results_MNAR =  pbmapply(samp.f,
                           A = As,
                           B=Bs_MNAR,
                           MoreArgs = list(p = p,
                                           miss = 'MNAR'),
                           SIMPLIFY = FALSE)
  
  
  raw_CDF1  = rbind(
    results_MAR %>%
      lapply('[', 1) %>%
      bind_rows() %>%
      mutate(miss = 'MAR'),
    results_MNAR %>%
      lapply('[', 1) %>%
      bind_rows() %>%
      mutate(miss = 'MNAR')
  ) %>%
    mutate(var_name = var_type ) %>%
    dplyr::select(-var_type) 
  
  raw_CDF = raw_CDF1 %>% 
    left_join(raw_CDF1 %>%
                filter(var_name == 'asymp') %>%
                group_by(perc) %>%
                dplyr::summarize(MC_var = var(m_lm)),
              by = 'perc'
    ) %>%
    mutate(mod = mod,
           nB = nB,
           r = r,
           alp = alp,
           L = L)
  
  
  
  
  raw_q1  = rbind(
    results_MAR %>%
      lapply('[', 2) %>%
      bind_rows() %>%
      mutate(miss = 'MAR'),
    results_MNAR %>%
      lapply('[', 2) %>%
      bind_rows() %>%
      mutate(miss = 'MNAR')
  ) %>%
    mutate(var_name = var_type ) %>%
    dplyr::select(-var_type) 
  
  raw_q = raw_q1 %>%
    left_join(raw_q1 %>%
                filter(var_name == 'asymp') %>%
                group_by(perc) %>%
                dplyr::summarize(MC_var = var(hatq)),
              by = 'perc'
    ) %>%
    mutate(mod = mod,
           nB = nB,
           r = r,
           alp = alp,
           L = L)
  
  
  
  
  CDF_results = raw_CDF %>%
    dplyr::select(mod, nB, r, alp, L, miss, perc, var_name, MC_var, everything()) %>%
    pivot_longer(cols = 9:ncol(.)) %>%
    group_by(miss, perc, var_name, name) %>%
    dplyr::summarize(mean = mean(value)) %>%
    pivot_wider(values_from = 'mean') %>%
    arrange(var_name) %>%
    mutate(rb = ((MC_var - var_val)/MC_var)*100) %>%
    mutate(CR = CR * 100)
  
  
  q_results = raw_q %>%
    dplyr::select(mod, nB, r, alp, L, miss, perc, var_name, MC_var, everything()) %>%
    pivot_longer(cols = 9:ncol(.))%>%
    group_by(miss, perc, var_name, name) %>%
    dplyr::summarize(mean = mean(value)) %>%
    pivot_wider(values_from = 'mean') %>%
    arrange(var_name) %>%
    mutate(rb = ((MC_var - var_val)/MC_var)*100) %>%
    mutate(CR = CR * 100)
  
  
  
  final_results = list(raw_CDF, raw_q, 
                       CDF_results, q_results)  %>%
    lapply(ungroup)
  
  names(final_results) = c('raw_CDF',
                           'raw_q',
                           'CDF_results',
                           'q_results')
  
  return(final_results)
}

sim.varV = Vectorize(sim.var, vectorize.args = 'nB')

N = 100000
nsim = 1500
nA = 1000
nB = c(nA, .20*N)
seed=101
F_N = c(.01, .10, .25, .50, .75, .90, .99 )

# model f1
## y1


setwd("/Users/jeremyflood/OneDrive/Documents/Grad School/2024-2025/Spring 2025/reCDF Revisions/Variance Estimation")
f1 =nB %>%
  purrr::map(function(x){sim.varV(#r=r,
    N = N,
    nA = nA,
    nB = x,
    nsim = nsim,
    mod = 'f1',
    seed = seed)})

f1 %>% lapply('[', 1) %>% 
  bind_rows() %>%
  openxlsx::write.xlsx('f1_results_rawCDF.xlsx')

f1 %>% lapply('[', 2) %>% 
  bind_rows() %>%
  openxlsx::write.xlsx('f1_results_rawq.xlsx')

f1 %>% lapply('[', 3) %>% 
  bind_rows() %>%
  openxlsx::write.xlsx('f1_summary_CDF.xlsx')

f1 %>% lapply('[', 4) %>% 
  bind_rows() %>%
  openxlsx::write.xlsx('f1_summary_q.xlsx')


# model f3

f3 =nB %>%
  purrr::map(function(x){sim.varV(#r=r,
    N = N,
    nA = nA,
    nB = x,
    nsim = nsim,
    mod = 'f3',
    seed = seed)})

f3 %>% lapply('[', 1) %>% 
  bind_rows() %>%
  openxlsx::write.xlsx('f3_results_rawCDF.xlsx')

f3 %>% lapply('[', 2) %>% 
  bind_rows() %>%
  openxlsx::write.xlsx('f3_results_rawq.xlsx')

f3 %>% lapply('[', 3) %>% 
  bind_rows() %>%
  openxlsx::write.xlsx('f3_summary_CDF.xlsx')

f3 %>% lapply('[', 4) %>% 
  bind_rows() %>%
  openxlsx::write.xlsx('f3_summary_q.xlsx')


mod_list = 
c(paste0(rep('f1', 2), '_summary_', c('CDF', 'q'), '.xlsx'),
  paste0(rep('f3', 2), '_summary_', c('CDF', 'q'), '.xlsx')
)

n
est_type = rep(c('cdf', 'q'), 2)
mod = rep(c('f1', 'f3'), each = 2)

summary_together = c()
for(i in 1:length(mod_list)){
  summary_together[[i]] = openxlsx::read.xlsx(mod_list[i])
}


plots = c()

for(i in 1:length(nB)){
plots[[i]] = summary_together %>%
  bind_rows() %>%
  mutate(d = UL-LL) %>%
  dplyr::select(-c(MC_var, LL, UL, var_val)) %>% # am taking out MC variance, LL & UL, and variance estimate to clean up the plots
  mutate(mod = rep(c('f1', 'f3'), each = nrow(.)/2)) %>%
  mutate(est_type = rep(rep(c('cdf', 'q'), 2), each = nrow(.)/4) ) %>%
  mutate(nB = rep(nB, each = nrow(.)/8) ) %>%
  mutate(Estimate = ifelse(is.na(m_lm) == TRUE, hatq, m_lm),
         Truth = ifelse(is.na(F_N) == TRUE, q_N, F_N),
         relBias = ((Estimate - Truth)/Truth)*100) %>%
  dplyr::select(-c('hatq', 'm_lm', 'Estimate', 'Truth')) %>%
  dplyr::select(mod, est_type, nB, miss, perc, var_name,  q_N, F_N, relBias, everything()) %>%
  pivot_longer(cols = 9:ncol(.)) %>%
  filter(name != 'relBias') %>%
  #filter(name != 'd') %>%
  filter(nB == nB[i]) %>%
  # mutate(name = ifelse(name == 'relBias', 'rb_est',
  #                      ifelse(name == 'rb', 'rb_var', name))) %>%
  # mutate(est_type = ifelse(name == 'rb_var' & est_type == 'cdf', 'var_cdf',
  #                          ifelse(name == 'rb_var' & est_type == 'q', 'var_q', est_type))) %>%
  
  mutate(name = factor(name, levels = c('CR', 'rb', 'd'), labels = c(
    'CR' = TeX("\\% CR"),
    'rb' = TeX("\\% RB"),
    'd' = "AL"))) %>%
  mutate(est_type = factor(est_type, levels = c('cdf', 'q'),
                           labels= c('cdf' = TeX('$F(t)$'),
                                      'q' = TeX('$t(\\alpha)$')
                                     ))) %>%
  mutate(mod = factor(mod, labels = c('f1' = TeX('Model $\\xi_{1}$'), 'f3' = TeX('Model $\\xi_{3}$')))) %>%
  mutate(var_name = factor(var_name, labels = c('asymp' = 'V1', 'boot' = 'V2'))) %>%
  ggplot(aes(x = perc, y = value, col = var_name, group = var_name)) +
  geom_line(linewidth = .5) +
  geom_point(size = 1) +
  scale_color_discrete(labels = unname(TeX(c("$V_{1}$ (Asymptotic)", "$V_{2}$ (Bootstrap)")))) +
  facet_grid(scales ='free', col = dplyr::vars(mod, miss), rows = dplyr::vars(est_type, name),
             labeller = label_parsed) +
  theme_bw()+
  ylab("Performance Value") +
  xlab('Percentile') +
  theme(legend.position = 'top') +
  ggtitle(TeX(paste0('Variance Estimation for $n_{B} =', ifelse(nB[i]/nA ==1, '', nB[i]/nA), 'n_{A}$')))+
  labs(colour = "Variance Type:")
}

# saving plots in high quality -- don't use Rstudio Export! Very blurry in Overleaf!

setwd("/Users/jeremyflood/OneDrive/Documents/Grad School/2024-2025/Spring 2025/reCDF Revisions/Variance Estimation/Plots")

# 1000 pixels is roughly 10.416666666666666 inches
# 786 pixels is roughly 8.1875 inches

ggsave('var_nA.png', plots[[1]], dpi=600, width = 11, height = 8)
ggsave('var_20nA.png', plots[[2]], dpi=600, width = 11, height = 8)

