# -----------------------------------------------------------
# Figure 1 script
# simulated outcomes in whole pop w/o waning 12+ vs. 18+
# -----------------------------------------------------------

# load packages ---------------------------------------------
library(tidyr)
library(dplyr)
library(ggplot2)
library(cowplot)

# define function for data wrangling ------------------------
wrangle_results <- function(x, times){
  n_sim <- length(x)
  rtn <- list()
  
  for (i in 1:n_sim){
  tmp <- lapply(x[[i]], wide_to_long, times) %>%
    bind_rows() %>%
    mutate(sim = i)
  
  rtn[[i]] <- tmp
  }

  rtn1 <- bind_rows(rtn) %>%
    group_by(time, state, age_group) %>%
    summarise(mean = mean(value),
              lower = quantile(value, probs = 0.025),
              upper = quantile(value, probs = 0.975)) 
  
  return(rtn1)
}

# read in simulation results --------------------------------
file_date <- "2021-08-26"
# 12+
mle_res_12plus <- readRDS(paste0("inst/extdata/results/results_12plus_mle_beta_", file_date, ".rds"))
lower_res_12plus <- readRDS(paste0("inst/extdata/results/results_12plus_lower_beta_", file_date, ".rds"))
upper_res_12plus <- readRDS(paste0("inst/extdata/results/results_12plus_upper_beta_", file_date, ".rds"))

# 18+
mle_res_18plus <- readRDS(paste0("inst/extdata/results/results_18plus_mle_beta_", file_date, ".rds"))
lower_res_18plus <- readRDS(paste0("inst/extdata/results/results_18plus_lower_beta_", file_date, ".rds"))
upper_res_18plus <- readRDS(paste0("inst/extdata/results/results_18plus_upper_beta_", file_date, ".rds"))

# wrangle raw results ----------------------------------------
# 12 +
mle_12plus_all <- wrangle_results(mle_res_12plus$out_all, times)
lower_12plus_all <- wrangle_results(lower_res_12plus$out_all, times)
upper_12plus_all <- wrangle_results(upper_res_12plus$out_all, times)

all_12plus <- bind_rows(mle_12plus_all, lower_12plus_all, upper_12plus_all, .id = "R0") %>%
  mutate(R0 = case_when(
    R0 == 1 ~ "4.6",
    R0 == 2 ~ "3.45",
    R0 == 3 ~ "9"
  ), Scenario = "12+") %>%
  pivot_wider(names_from = state, values_from = mean:upper) %>%
  mutate(h = rep(params$h, 3),
         i1 = rep(params$i1, 3),
         i2 = rep(params$i2, 3),
         d = rep(params$d, 3),
         d_ic = rep(params$d_ic, 3),
         d_hic = rep(params$d_hic, 3),
         cases_mle = params$sigma * (mean_E + mean_Ev_1d + mean_Ev_2d) * params$p_report,
         cases_lower = params$sigma * (lower_E + lower_Ev_1d + lower_Ev_2d) * params$p_report,
         cases_upper = params$sigma * (upper_E + upper_Ev_1d + upper_Ev_2d) * params$p_report,
         hosp_mle = h * (mean_I + mean_Iv_1d + mean_Iv_2d),
         hosp_lower = h * (lower_I + lower_Iv_1d + lower_Iv_2d),
         hosp_upper = h * (upper_I + upper_Iv_1d + upper_Iv_2d),
         ic_mle = i1 * (mean_H + mean_Hv_1d + mean_Hv_2d),
         ic_lower = i1 * (lower_H + lower_Hv_1d + lower_Hv_2d),
         ic_upper = i1 * (upper_H + upper_Hv_1d + upper_Hv_2d),
         deaths_mle = d * (mean_H + mean_Hv_1d + mean_Hv_2d) + d_ic * (mean_IC + mean_ICv_1d + mean_ICv_2d) + d_hic * (mean_H_IC + mean_H_ICv_1d + mean_H_ICv_2d),
         deaths_lower = d * (lower_H + lower_Hv_1d + lower_Hv_2d) + d_ic * (lower_IC + lower_ICv_1d + lower_ICv_2d) + d_hic * (lower_H_IC + lower_H_ICv_1d + lower_H_ICv_2d),
         deaths_upper = d * (upper_H + upper_Hv_1d + upper_Hv_2d) + d_ic * (upper_IC + upper_ICv_1d + upper_ICv_2d) + d_hic * (upper_H_IC + upper_H_ICv_1d + upper_H_ICv_2d)
        ) %>%
  select(Scenario, R0, time, age_group, cases_mle:deaths_upper) %>%
  pivot_longer(cases_mle:deaths_upper, names_to = c("outcome", "estimate"), names_sep = "_", values_to = "value") %>%
  pivot_wider(names_from = "estimate", values_from = "value")

# 18+
mle_18plus_all <- wrangle_results(mle_res_18plus$out_all, times)
lower_18plus_all <- wrangle_results(lower_res_18plus$out_all, times)
upper_18plus_all <- wrangle_results(upper_res_18plus$out_all, times)

all_18plus <- bind_rows(mle_18plus_all, lower_18plus_all, upper_18plus_all, .id = "R0") %>%
  mutate(R0 = case_when(
    R0 == 1 ~ "4.6",
    R0 == 2 ~ "3.45",
    R0 == 3 ~ "9"
  ), Scenario = "18+") %>%
  pivot_wider(names_from = state, values_from = mean:upper) %>%
  mutate(h = rep(params$h, 3),
         i1 = rep(params$i1, 3),
         i2 = rep(params$i2, 3),
         d = rep(params$d, 3),
         d_ic = rep(params$d_ic, 3),
         d_hic = rep(params$d_hic, 3),
         cases_mle = params$sigma * (mean_E + mean_Ev_1d + mean_Ev_2d) * params$p_report,
         cases_lower = params$sigma * (lower_E + lower_Ev_1d + lower_Ev_2d) * params$p_report,
         cases_upper = params$sigma * (upper_E + upper_Ev_1d + upper_Ev_2d) * params$p_report,
         hosp_mle = h * (mean_I + mean_Iv_1d + mean_Iv_2d),
         hosp_lower = h * (lower_I + lower_Iv_1d + lower_Iv_2d),
         hosp_upper = h * (upper_I + upper_Iv_1d + upper_Iv_2d),
         ic_mle = i1 * (mean_H + mean_Hv_1d + mean_Hv_2d),
         ic_lower = i1 * (lower_H + lower_Hv_1d + lower_Hv_2d),
         ic_upper = i1 * (upper_H + upper_Hv_1d + upper_Hv_2d),
         deaths_mle = d * (mean_H + mean_Hv_1d + mean_Hv_2d) + d_ic * (mean_IC + mean_ICv_1d + mean_ICv_2d) + d_hic * (mean_H_IC + mean_H_ICv_1d + mean_H_ICv_2d),
         deaths_lower = d * (lower_H + lower_Hv_1d + lower_Hv_2d) + d_ic * (lower_IC + lower_ICv_1d + lower_ICv_2d) + d_hic * (lower_H_IC + lower_H_ICv_1d + lower_H_ICv_2d),
         deaths_upper = d * (upper_H + upper_Hv_1d + upper_Hv_2d) + d_ic * (upper_IC + upper_ICv_1d + upper_ICv_2d) + d_hic * (upper_H_IC + upper_H_ICv_1d + upper_H_ICv_2d)
  ) %>%
  select(Scenario, R0, time, age_group, cases_mle:deaths_upper) %>%
  pivot_longer(cases_mle:deaths_upper, names_to = c("outcome", "estimate"), names_sep = "_", values_to = "value") %>%
  pivot_wider(names_from = "estimate", values_from = "value")

# make plots --------------------------------------------------
all_res_for_plot <- all_res %>%
  pivot_wider(names_from = "estimate", values_from = "value")

# figure 1 - 12+ vs. 18+, no waning ---------------------------
fig1 <- ggplot(data = all_res_for_plot, aes(x = date, y = mle, fill = R0, 
                                            linetype = Scenario)) +
  geom_line() +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = R0), alpha = 0.3) +
  labs(y = "Value", x = "Time (days)") +
  ylim(0,NA) +
  scale_x_date(date_breaks = "2 weeks", date_labels = "%d %b %Y") +
  theme(legend.position = "bottom",
        panel.background = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
        axis.text.y = element_text(size = 14),
        strip.text.x = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        axis.title=element_text(size=14,face="bold")) +
  facet_wrap(~outcome, scales = "free_y", nrow = 4)
fig1

ggsave(filename = "inst/extdata/results/figure 1 alt.jpg", plot = fig1,
       units = "in", height = 10, width = 8, dpi = 300)



# old code 
