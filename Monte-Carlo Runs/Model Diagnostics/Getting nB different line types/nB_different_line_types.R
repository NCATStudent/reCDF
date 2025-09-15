setwd("/Users/jeremyflood/Library/CloudStorage/OneDrive-Personal/Documents/Grad School/2024-2025/Fall 2025/reCDF/reCDF/Monte-Carlo Runs/r_15")

mod_list <- paste0("modf", 1:4, "_results.xlsx")


dat_together <- c()

for (i in 1:length(mod_list)) {
  dat_together[[i]] <- openxlsx::read.xlsx(mod_list[i], sheet = "perf") %>% mutate(mod = paste0("f", i))
}

est_label <- c(
  "B_plug_cdf" = TeX("$\\hat{F}_{B}$"),
  "lm_plug_cdf" = TeX("$\\hat{F}_{P}"),
  "m_lm_cdf" = TeX("$\\hat{F}_{R}"),
  "B_plug_q" = TeX("$\\hat{t}_{B}"),
  "lm_plug_q" = TeX("$\\hat{t}_{P}"),
  "m_lm_q" = TeX("$\\hat{t}_{R}")
)


plots_together <- c()

dat_together %>%
  bind_rows() %>%
  # filter(est != 'B_plug') %>%
  # filter(miss == 'MAR') %>%
  mutate(new_est = paste0(est, "_", group)) %>%
  dplyr::select(nB, mod, miss, perc, new_est, est, group, rrmse) %>%
  mutate(est = factor(est, labels = c("B_plug" = "Naive", "lm_plug" = "Plug-in", "m_lm" = "Residual eCDF"))) %>%
  mutate(mod = factor(mod,
    labels = c(
      "f1" = TeX("Model $\\xi_{1}"),
      "f2" = TeX("Model $\\xi_{2}"),
      "f3" = TeX("Model $\\xi_{3}"),
      "f4" = TeX("Model $\\xi_{4}")
    )
  )) %>%
  mutate(cat_est = ifelse(new_est %in% c("B_plug_cdf", "lm_plug_cdf", "m_lm_cdf"),
    "F(t)", "t(a)"
  )) %>%
  mutate(new_est = factor(new_est,
    levels = c(
      "B_plug_cdf", "lm_plug_cdf", "m_lm_cdf",
      "B_plug_q", "lm_plug_q", "m_lm_q"
    ),
    labels = c(
      "B_plug_cdf" = TeX("$\\hat{F}_{B}$"),
      "lm_plug_cdf" = TeX("$\\hat{F}_{P}"),
      "m_lm_cdf" = TeX("$\\hat{F}_{R}"),
      "B_plug_q" = TeX("$\\hat{t}_{B}"),
      "lm_plug_q" = TeX("$\\hat{t}_{P}"),
      "m_lm_q" = TeX("$\\hat{t}_{R}")
    )
  )) %>%
  mutate(nB = factor(nB,
    labels = c(
      "100" = TeX("$n_B = 100$"),
      "1000" = TeX("$n_B = n_A$"),
      "10000" = TeX("$n_B = 10n_A$"),
      "20000" = TeX("$n_B = 20n_A$"),
      "40000" = TeX("$n_B = 40n_A$")
    )
  )) %>%
  mutate(perc_strip = as.numeric(gsub(perc, pattern = "%", replacement = "")) / 100) %>%
  mutate(perc_strip = factor(perc_strip)) %>%
  mutate(Missingness = factor(miss)) %>%
  arrange(nB, perc) %>%
  ggplot(aes(x = perc_strip, y = rrmse^(1 / 2), color = est, group= interaction(est, nB))) +
  geom_line(aes(linetype = nB)) +
  geom_point(size = 1) +
  facet_grid(
    scales = "fixed", rows = dplyr::vars(cat_est, Missingness), cols = dplyr::vars(mod),
    labeller = label_parsed
  ) +
  # scale_x_discrete(labels = c('B_plug_cdf' = TeX("$\\hat{F}_{B}$"),
  #                             'lm_plug_cdf' = TeX("$\\hat{F}_{P}"),
  #                             'm_lm_cdf' = TeX("$\\hat{F}_{R}"),
  #                             'B_plug_q' = TeX("$\\hat{t}_{B}"),
  #                             'lm_plug_q' = TeX("$\\hat{t}_{P}"),
  #                             'm_lm_q' = TeX("$\\hat{t}_{R}")
  # )
  # ) +
  ylab(TeX("\\sqrt{RMSER}")) +
  xlab("Percentile") +
  theme_bw() +
  theme(
    text = element_text(size = 15),
    axis.text.x = element_text(angle = 90, hjust = 1)
  ) +
  theme(legend.position = "top") #+
  for (i in 1:length(nB)) {
    nBs <- nB[i] # needed bc R's dumbass can't do a fucking filter condition correctly.
    plots_together[[i]] <- dat_together %>%
      bind_rows() %>%
      filter(nB == nBs) %>%
      # filter(est != 'B_plug') %>%
      # filter(miss == 'MAR') %>%
      mutate(new_est = paste0(est, "_", group)) %>%
      dplyr::select(nB, mod, miss, perc, new_est, est, group, rrmse) %>%
      mutate(est = factor(est, labels = c("B_plug" = "Naive", "lm_plug" = "Plug-in", "m_lm" = "Residual eCDF"))) %>%
      mutate(mod = factor(mod,
        labels = c(
          "f1" = TeX("Model $\\xi_{1}"),
          "f2" = TeX("Model $\\xi_{2}"),
          "f3" = TeX("Model $\\xi_{3}"),
          "f4" = TeX("Model $\\xi_{4}")
        )
      )) %>%
      mutate(cat_est = ifelse(new_est %in% c("B_plug_cdf", "lm_plug_cdf", "m_lm_cdf"),
        "F(t)", "t(a)"
      )) %>%
      mutate(new_est = factor(new_est,
        levels = c(
          "B_plug_cdf", "lm_plug_cdf", "m_lm_cdf",
          "B_plug_q", "lm_plug_q", "m_lm_q"
        ),
        labels = c(
          "B_plug_cdf" = TeX("$\\hat{F}_{B}$"),
          "lm_plug_cdf" = TeX("$\\hat{F}_{P}"),
          "m_lm_cdf" = TeX("$\\hat{F}_{R}"),
          "B_plug_q" = TeX("$\\hat{t}_{B}"),
          "lm_plug_q" = TeX("$\\hat{t}_{P}"),
          "m_lm_q" = TeX("$\\hat{t}_{R}")
        )
      )) %>%
      # mutate(nB = factor(nB,
      #                    labels = c('1000' = TeX("$n_B = n_A$"),
      #                               '10000' = TeX("$n_B = 10n_A$"),
      #                               '20000' = TeX("$n_B = 20n_A$")))) %>%
      mutate(perc_strip = as.numeric(gsub(perc, pattern = "%", replacement = "")) / 100) %>%
      mutate(perc_strip = factor(perc_strip)) %>%
      mutate(Missingness = factor(miss)) %>%
      arrange(nB, perc) %>%
      ggplot(aes(x = perc_strip, y = rrmse^(1 / 2), color = est, group = est)) +
      geom_line(linewidth = .5) +
      geom_point(size = 1) +
      facet_grid(
        scales = "fixed", rows = dplyr::vars(cat_est, Missingness), cols = dplyr::vars(mod),
        labeller = label_parsed
      ) +
      # scale_x_discrete(labels = c('B_plug_cdf' = TeX("$\\hat{F}_{B}$"),
      #                             'lm_plug_cdf' = TeX("$\\hat{F}_{P}"),
      #                             'm_lm_cdf' = TeX("$\\hat{F}_{R}"),
      #                             'B_plug_q' = TeX("$\\hat{t}_{B}"),
      #                             'lm_plug_q' = TeX("$\\hat{t}_{P}"),
      #                             'm_lm_q' = TeX("$\\hat{t}_{R}")
      # )
      # ) +
      ylab(TeX("\\sqrt{RMSER}")) +
      xlab("Percentile") +
      theme_bw() +
      theme(
        text = element_text(size = 15),
        axis.text.x = element_text(angle = 90, hjust = 1)
      ) +
      theme(legend.position = "top") +
      ggtitle(ifelse(nBs == nA, TeX(paste0("\\sqrt{RMSER} for $n_{B} = n_{A}$, r = ", r)),
        TeX(paste0("\\sqrt{RMSER} for $n_{B} =") %>% paste0(nBs / nA) %>% paste0("nA$") %>% paste0(", r =", r))
      )) +
      labs(colour = NULL)
  }
setwd("/Users/jeremyflood/Library/CloudStorage/OneDrive-Personal/Documents/Grad School/2024-2025/Fall 2025/reCDF/reCDF/Monte-Carlo Runs/Model Diagnostics/Getting nB different line types")

ggsave("RRMSE_nB_1.png", plots_together[[1]], dpi = 400, width = 10, height = 8)
ggsave("RRMSE_nB_2.png", plots_together[[2]], dpi = 400, width = 10, height = 8)
ggsave("RRMSE_nB_3.png", plots_together[[3]], dpi = 400, width = 10, height = 8)
ggsave("RRMSE_nB_4.png", plots_together[[4]], dpi = 400, width = 10, height = 8)
ggsave("RRMSE_nB_5.png", plots_together[[5]], dpi = 400, width = 10, height = 8)