# -----------------------------------------------------------
# Figure 1 script
# simulated outcomes in whole pop w/o waning 12+ vs. 18+
# -----------------------------------------------------------

# load packages ---------------------------------------------
library(dplyr)
library(ggplot2)
library(cowplot)

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

# create a joint dataframe ----------------------------------
all_res_12plus <- bind_rows(mle_res_12plus$df_total, lower_res_12plus$df_total, upper_res_12plus$df_total, .id = "R0") %>%
  pivot_longer(cases_mle:deaths_upper, names_to = c("outcome", "estimate"), names_sep = "_", values_to = "value") %>%
  mutate(R0 = case_when(
    R0 == 1 ~ "4.6",
    R0 == 2 ~ "3.45",
    R0 == 3 ~ "9"
  ),
  date = time + as.Date("2020-01-01"),
  outcome = factor(case_when(
    outcome == "cases" ~ "Cases",
    outcome == "hosp" ~ "Hospital Admissions",
    outcome == "ic" ~ "IC Admissions",
    outcome == "deaths" ~ "Deaths"
  ), levels = c("Cases", "Hospital Admissions", "IC Admissions", "Deaths")),
  Scenario = "12+"
  )

all_res_18plus <- bind_rows(mle_res_18plus$df_total, lower_res_18plus$df_total, upper_res_18plus$df_total, .id = "R0") %>%
  pivot_longer(cases_mle:deaths_upper, names_to = c("outcome", "estimate"), names_sep = "_", values_to = "value") %>%
  mutate(R0 = case_when(
    R0 == 1 ~ "4.6",
    R0 == 2 ~ "3.45",
    R0 == 3 ~ "9"
  ),
  date = time + as.Date("2020-01-01"),
  outcome = factor(case_when(
    outcome == "cases" ~ "Cases",
    outcome == "hosp" ~ "Hospital Admissions",
    outcome == "ic" ~ "IC Admissions",
    outcome == "deaths" ~ "Deaths"
  ), levels = c("Cases", "Hospital Admissions", "IC Admissions", "Deaths")),
  Scenario = "18+"
  )

all_res <- bind_rows(all_res_12plus, all_res_18plus)
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
  facet_wrap(~outcome, scales = "free")
fig1

ggsave(filename = "inst/extdata/results/figure 1.jpg", plot = fig1,
       units = "in", height = 10, width = 12, dpi = 300)