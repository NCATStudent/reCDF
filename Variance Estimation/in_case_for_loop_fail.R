# did your for-loop fail? use this script to help clean up results

setwd("/Users/jeremyflood/Library/CloudStorage/OneDrive-Personal/Documents/Grad School/2024-2025/Fall 2025/reCDF/reCDF/Variance Estimation/Data/Cached Iter Files")


# getting f1 data
f1_path <- "/Users/jeremyflood/Library/CloudStorage/OneDrive-Personal/Documents/Grad School/2024-2025/Fall 2025/reCDF/reCDF/Variance Estimation/Additional Requests/Adding Smaller nB/f1/Cached Iter Files"
setwd(f1_path)
all_files <- list.files(f1_path)

f1_agg <- c()
for (i in 1:length(all_files)) {
  f1_agg[[i]] <- openxlsx::read.xlsx(all_files[i])
  print(paste0("completed ", round(i / length(all_files) * 100, 0)))
}

f3_path <- "/Users/jeremyflood/Library/CloudStorage/OneDrive-Personal/Documents/Grad School/2024-2025/Fall 2025/reCDF/reCDF/Variance Estimation/Additional Requests/Adding Smaller nB/f3/Cached Iter Files"
setwd(f3_path)
all_files_f3 <- list.files(f3_path)


f3_agg <- c()
for (i in 1:length(all_files_f3)) {
  f3_agg[[i]] <- openxlsx::read.xlsx(all_files_f3[i]) %>% mutate(mod = "f3")
  print(paste0("completed ", round(i / length(all_files_f3) * 100, 0)))
}

cleaned_ddf <- rbind(
  f1_agg %>% bind_rows(),
  f3_agg %>% bind_rows()
)


mc_results <- cleaned_ddf %>%
  # filter(var_name == 'asymp') %>%
  group_by(mod, perc, miss, nB, est_type, var_type) %>%
  dplyr::summarize(MC_var = var(est_quant, na.rm = TRUE))


cleaned_results <-
  cleaned_ddf %>%
  # left_join(mc_results) %>%
  dplyr::select(est_type, est_quant, pop_quant, everything()) %>%
  pivot_longer(cols = 2) %>%
  group_by(mod, est_type, perc, miss, nB, var_type, name) %>%
  dplyr::summarize(
    pop_quant = mean(pop_quant, na.rm = TRUE),
    est_value = mean(value, na.rm = TRUE),
    est_var = mean(var_val, na.rm = TRUE),
    CR = mean(CR, na.rm = TRUE),
    d = mean(UL - LL, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  left_join(mc_results, by = c("est_type", "perc", "miss", "nB", "var_type", "mod")) %>%
  mutate(rb = ((MC_var - est_var) / MC_var) * 100) %>%
  mutate(nsim = 1500, L = 1500)

final_results <- list(
  cleaned_ddf,
  cleaned_results
)
names(final_results) <- c("raw", "summary")

setwd("/Users/jeremyflood/Library/CloudStorage/OneDrive-Personal/Documents/Grad School/2024-2025/Fall 2025/reCDF/reCDF/Variance Estimation/Additional Requests/Adding Smaller nB")
openxlsx::write.xlsx(cleaned_results, paste0("final_agg_results.xlsx")) # github doesn't like large files
