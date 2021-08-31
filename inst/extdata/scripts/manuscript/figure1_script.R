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
wrangle_results <- function(my_list){
  n_sim <- length(my_list)
  rtn <- list()
  for (s in 1:n_sim){
    # convert output from each simulation into long data frame
    tmp <- lapply(my_list[[s]], wide_to_long) %>%
      bind_rows() %>%
      mutate(sim = s)
    
    rtn[[s]] <- tmp
    
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
mle_12plus_all <- wrangle_results(mle_res_12plus$out_all)
lower_12plus_all <- wrangle_results(lower_res_12plus$out_all)
upper_12plus_all <- wrangle_results(upper_res_12plus$out_all)

mle_18plus_all <- wrangle_results(mle_res_18plus$out_all)
lower_18plus_all <- wrangle_results(lower_res_18plus$out_all)
upper_18plus_all <- wrangle_results(upper_res_18plus$out_all)

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
