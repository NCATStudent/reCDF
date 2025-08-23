# plot generation

setwd("/Users/jeremyflood/Library/CloudStorage/OneDrive-Personal/Documents/Grad School/2024-2025/Fall 2025/reCDF/reCDF/Variance Estimation/Data")

cleaned_results = 
  openxlsx::read.xlsx('final_agg_results.xlsx', sheet = 'summary' ) %>%
  bind_rows() %>%
  mutate(var_type = ifelse(
    est_type == 't' & var_type == 'boot', 'boot_ignore', var_type
  )) %>%
  filter(var_type != 'boot_ignore') %>%
  mutate(var_type = ifelse(
    var_type == 'boot2', 'boot', var_type
  )) %>%
  dplyr::select(-name) # needs to be checked for general case of all data together; my data looks interesting bc for loop failed and needed to depend on cached results


plots = c()

mods = c(1, 3)

FN_plots = TN_plots = c()
for(i in 1:length(mods)){
  FN_plots[[i]] = 
    cleaned_results %>%
    filter(mod == paste0('f', mods[i]), est_type == 'cdf' ) %>% # for now
    mutate(se_val = sqrt(est_var)) %>%
    mutate(rb = (est_var - MC_var) / MC_var * 100) %>%
    dplyr::select(-c(MC_var, pop_quant, est_value, est_var)) %>%
    dplyr::select(CR, se_val, d, rb, everything()) %>%
    pivot_longer(cols = 1:4, names_to = 'perf') %>%
    mutate(perf = factor(perf, levels = c('CR', 'd', 'se_val', 'rb'), labels = c(
      'CR' = TeX("\\% CR"),
      'd' = "AL",
      'se_val' = TeX("\\widehat{SE}"),
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
    scale_color_discrete(labels = unname(TeX(c("Asymptotic", "Bootstrap")))) +
    facet_grid(scales ='free', col = dplyr::vars(miss, nB), rows = dplyr::vars(perf),
               labeller = label_parsed) +
    theme_bw() +
    ylab("Performance Value") +
    xlab('Percentile') +
    theme(legend.position = 'top') +
    ggtitle(TeX(paste0('$F_{N}(t)$ Variance Estimation (Model $\\xi_', mods[i] , ')'))) +
    labs(colour = "Variance Type:")
  
  ggsave(paste0('var_cdf_f', mods[i], '.png'), FN_plots[[i]], dpi=550, width = 11, height = 8)
}

TN_plots = c()

for(i in 1:length(mods)){
  #TN_plots[[i]] = 
  cleaned_results %>%
    filter(mod == paste0('f', mods[i]), est_type == 't' ) %>% # for now
    mutate(se_val = sqrt(est_var)) %>%
    mutate(rb = (est_var - MC_var) / MC_var * 100) %>%
    dplyr::select(-c(MC_var, pop_quant, est_value, est_var)) %>%
    dplyr::select(CR, se_val, d, rb, everything()) %>%
    pivot_longer(cols = 1:4, names_to = 'perf') %>%
    mutate(perf = factor(perf, levels = c('CR', 'd', 'se_val', 'rb'), labels = c(
      'CR' = TeX("\\% CR"),
      'd' = "AL",
      'se_val' = TeX("\\widehat{SE}"),
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
    scale_color_discrete(labels = unname(TeX(c("Asymptotic", "Bootstrap")))) +
    facet_grid(scales ='free', col = dplyr::vars(miss, nB), rows = dplyr::vars(perf),
               labeller = label_parsed) +
    theme_bw() +
    ylab("Performance Value") +
    xlab('Percentile') +
    theme(legend.position = 'top') +
    ggtitle(TeX(paste0('$T_{N}(\\alpha)$ Variance Estimation (Model $\\xi_', mods[i] , ')'))) +
    labs(colour = "Variance Type:")
  
  ggsave(paste0('var_t_f', mods[i], '.png'), TN_plots[[i]], dpi=550, width = 11, height = 8)
}



#for(i in 1:length(nB)){
#plots[[i]] =
summary_together %>%
  bind_rows() %>%
  fil
dplyr::select(-c(MC_var, est_var)) %>% # am taking out MC variance and variance estimate to clean up the plots
  dplyr::select(est_value, CR, d, rb, everything()) %>%
  pivot_longer(cols = 1:4) %>%
  mutate(name = factor(name, 
                       levels = c('est_value', 'CR', 'd', 'rb'),
                       labels = c(
                         'est_value' = TeX("$\\hat{\\theta}$"),
                         'CR' = TeX("\\% CR"),
                         'rb' = TeX("\\% RB"),
                         'd' = "AL"))
  ) %>%
  mutate(est_type = factor(est_type, levels = c('cdf', 'q'),
                           labels= c('cdf' = TeX('$F(t)$'),
                                     't' = TeX('$T(\\alpha)$')
                           ))) %>%
  mutate(mod = factor(mod, labels = c('f1' = TeX('Model $\\xi_{1}$'), 'f3' = TeX('Model $\\xi_{3}$')))) %>%
  rename(c('var_type'= 'var_name')) %>%
  mutate(var_name = factor(var_name, labels = c('asymp' = 'V1', 'boot' = 'V2')))%>%
  ggplot(aes(x = perc, y = value, col = var_name, group = var_name)) +
  geom_line(linewidth = .5) +
  geom_point(size = 1) +
  scale_color_discrete(labels = unname(TeX(c("$V_{1}$ (Asymptotic)", "$V_{2}$ (Bootstrap)")))) +
  facet_grid(scales ='free', col = dplyr::vars(mod, miss), rows = dplyr::vars(est_type, name),
             labeller = label_parsed) +
  theme_bw()#+
ylab("Performance Value") +
  xlab('Percentile') +
  theme(legend.position = 'top') +
  ggtitle(TeX(paste0('Variance Estimation for $n_{B} =', ifelse(nB[i]/nA ==1, '', nB[i]/nA), 'n_{A}$')))+
  labs(colour = "Variance Type:")
#}

# saving plots in high quality -- don't use Rstudio Export! Very blurry in Overleaf!

setwd("/Users/jeremyflood/OneDrive/Documents/Grad School/2024-2025/Spring 2025/reCDF Revisions/Variance Estimation/Plots")

# 1000 pixels is roughly 10.416666666666666 inches
# 786 pixels is roughly 8.1875 inches

ggsave('var_nA.png', plots[[1]], dpi=600, width = 11, height = 8)
ggsave('var_20nA.png', plots[[2]], dpi=600, width = 11, height = 8)

