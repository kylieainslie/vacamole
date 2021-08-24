# manuscript plots
library(dplyr)
library(ggplot2)

# read in simulation results --------------------------------
# 12+
mle_res_12plus <- readRDS("inst/extdata/results/df_basis_12plus_mle_beta_23aug.rds")
lower_res_12plus <- readRDS("inst/extdata/results/df_basis_12plus_lower_beta_23aug.rds")
upper_res_12plus <- readRDS("inst/extdata/results/df_basis_12plus_upper_beta_23aug.rds")

# 18+
mle_res_18plus <- readRDS("inst/extdata/results/df_basis_18plus_mle_beta_23aug.rds")
lower_res_18plus <- readRDS("inst/extdata/results/df_basis_18plus_lower_beta_23aug.rds")
upper_res_18plus <- readRDS("inst/extdata/results/df_basis_18plus_upper_beta_23aug.rds")

# create a joint dataframe2 ----------------------------------
all_res_12plus <- bind_rows(mle_res_12plus, lower_res_12plus, upper_res_12plus, .id = "R0") %>%
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

all_res_18plus <- bind_rows(mle_res_18plus, lower_res_18plus, upper_res_18plus, .id = "R0") %>%
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
# figure 1 - 12+ vs. 18+, no waning
fig1 <- ggplot(data = all_res, aes(x = date, y = mle, fill = R0, 
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

# figure 2 - 12+ vs. 18+, waning
fig1 <- ggplot(data = all_res, aes(x = date, y = mle, fill = R0)) +
  geom_line() +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = R0), alpha = 0.3) +
  labs(y = "Value", x = "Time (days)") +
  ylim(0,NA) +
  scale_x_date(date_breaks = "2 weeks", date_labels = "%d %b %Y") +
  theme(legend.position = "bottom",
        panel.background = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_wrap(~outcome, scales = "free")
fig1

ggsave(filename = "inst/extdata/results/figure 1.jpg", plot = fig1)

# figure S1 - model fits to Osiris cases with confidence bounds
