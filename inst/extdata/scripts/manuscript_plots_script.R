# manuscript plots
library(dplyr)
library(ggplot2)
library(cowplot)
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

# figure SX - vaccination coverage
# first data wrangle for plot
age_dist_10 <- c(0.10319920, 0.11620856, 0.12740219, 0.12198707, 0.13083463, 
                 0.14514332, 0.12092904, 0.08807406, 0.03976755, 0.007398671)
n <- 17407585 # Dutch population size
n_vec_10 <- n * age_dist_10

# 12 plus --------------------------------------------------------
basis_12plus_long <- basis_12plus %>%
  mutate(pf_d1_9 = (pf_d1_9 * n_vec_10[9] + pf_d1_10 * n_vec_10[10])/sum(n_vec_10[9:10]),
         pf_d2_9 = (pf_d2_9 * n_vec_10[9] + pf_d2_10 * n_vec_10[10])/sum(n_vec_10[9:10]),
         mo_d1_9 = (mo_d1_9 * n_vec_10[9] + mo_d1_10 * n_vec_10[10])/sum(n_vec_10[9:10]),
         mo_d2_9 = (mo_d2_9 * n_vec_10[9] + mo_d2_10 * n_vec_10[10])/sum(n_vec_10[9:10]),
         az_d1_9 = (az_d1_9 * n_vec_10[9] + az_d1_10 * n_vec_10[10])/sum(n_vec_10[9:10]),
         az_d2_9 = (az_d2_9 * n_vec_10[9] + az_d2_10 * n_vec_10[10])/sum(n_vec_10[9:10]),
         ja_d1_9 = (ja_d1_9 * n_vec_10[9] + ja_d1_10 * n_vec_10[10])/sum(n_vec_10[9:10]),
         ja_d2_9 = (ja_d2_9 * n_vec_10[9] + ja_d2_10 * n_vec_10[10])/sum(n_vec_10[9:10])
  ) %>%
  select(date, pf_d1_1:ja_d2_9, -pf_d1_10, -pf_d2_10, -mo_d1_10, -mo_d2_10, -az_d1_10, 
         -az_d2_10, -ja_d1_10, -ja_d2_10) %>%
  pivot_longer(!date, 
               names_to = c("vaccine", "dose", "age_group"),
               names_sep = "_",
               values_to = "coverage") %>%
  mutate(Dose = as.factor(case_when(
    dose == "d1" ~ 1,
    dose == "d2" ~ 2
  )),
  date = as.Date(date, format = "%m/%d/%Y"),
  coverage = ifelse(coverage == 0, NA, coverage),
  Vaccine = case_when(
    vaccine == "az" ~ "AstraZeneca",
    vaccine == "pf" ~ "Pfizer/BioNTech",
    vaccine == "mo" ~ "Moderna",
    vaccine == "ja" ~ "Janssen"
  ),
  age_group = case_when(
    age_group == 1 ~ "0-9",
    age_group == 2 ~ "10-19",
    age_group == 3 ~ "20-29",
    age_group == 4 ~ "30-39",
    age_group == 5 ~ "40-49",
    age_group == 6 ~ "50-59",
    age_group == 7 ~ "60-69",
    age_group == 8 ~ "70-79",
    age_group == 9 ~ "80+"
  )) %>%
  filter(date > as.Date("2020-12-31"))

# plot
figsxa <- ggplot(data = basis_12plus_long, 
                aes(x = date, y = coverage, color = Vaccine, linetype = Dose)) +
  geom_line() +
  labs(y = "Vaccine Coverage", x = "Date") +
  ylim(0,1) +
  scale_x_date(date_breaks = "2 months", date_labels = "%d %b %Y") +
  facet_grid(.~age_group) +
  theme(#legend.position = "none",
        panel.background = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
        axis.text.y = element_text(size = 14),
        strip.text.x = element_text(size = 14),
        #legend.text = element_text(size = 14),
        #legend.title = element_text(size = 14),
        axis.title=element_text(size=14,face="bold"))
figsxa

# 18 plus --------------------------------------------------------
basis_18plus_long <- basis_18plus %>%
  mutate(pf_d1_9 = (pf_d1_9 * n_vec_10[9] + pf_d1_10 * n_vec_10[10])/sum(n_vec_10[9:10]),
         pf_d2_9 = (pf_d2_9 * n_vec_10[9] + pf_d2_10 * n_vec_10[10])/sum(n_vec_10[9:10]),
         mo_d1_9 = (mo_d1_9 * n_vec_10[9] + mo_d1_10 * n_vec_10[10])/sum(n_vec_10[9:10]),
         mo_d2_9 = (mo_d2_9 * n_vec_10[9] + mo_d2_10 * n_vec_10[10])/sum(n_vec_10[9:10]),
         az_d1_9 = (az_d1_9 * n_vec_10[9] + az_d1_10 * n_vec_10[10])/sum(n_vec_10[9:10]),
         az_d2_9 = (az_d2_9 * n_vec_10[9] + az_d2_10 * n_vec_10[10])/sum(n_vec_10[9:10]),
         ja_d1_9 = (ja_d1_9 * n_vec_10[9] + ja_d1_10 * n_vec_10[10])/sum(n_vec_10[9:10]),
         ja_d2_9 = (ja_d2_9 * n_vec_10[9] + ja_d2_10 * n_vec_10[10])/sum(n_vec_10[9:10])
  ) %>%
  select(date, pf_d1_1:ja_d2_9, -pf_d1_10, -pf_d2_10, -mo_d1_10, -mo_d2_10, -az_d1_10, 
         -az_d2_10, -ja_d1_10, -ja_d2_10) %>%
  pivot_longer(!date, 
               names_to = c("vaccine", "dose", "age_group"),
               names_sep = "_",
               values_to = "coverage") %>%
  mutate(Dose = as.factor(case_when(
    dose == "d1" ~ 1,
    dose == "d2" ~ 2
  )),
  date = as.Date(date, format = "%m/%d/%Y"),
  coverage = ifelse(coverage == 0, NA, coverage),
  Vaccine = case_when(
    vaccine == "az" ~ "AstraZeneca",
    vaccine == "pf" ~ "Pfizer/BioNTech",
    vaccine == "mo" ~ "Moderna",
    vaccine == "ja" ~ "Janssen"
  ),
  age_group = case_when(
    age_group == 1 ~ "0-9",
    age_group == 2 ~ "10-19",
    age_group == 3 ~ "20-29",
    age_group == 4 ~ "30-39",
    age_group == 5 ~ "40-49",
    age_group == 6 ~ "50-59",
    age_group == 7 ~ "60-69",
    age_group == 8 ~ "70-79",
    age_group == 9 ~ "80+"
  )) %>%
  filter(date > as.Date("2020-12-31"))

# plot
figsxb <- ggplot(data = basis_18plus_long, 
                 aes(x = date, y = coverage, color = Vaccine, linetype = Dose)) +
  geom_line() +
  labs(y = "Vaccine Coverage", x = "Date") +
  ylim(0,1) +
  scale_x_date(date_breaks = "2 months", date_labels = "%d %b %Y") +
  facet_grid(.~age_group) +
  theme(legend.position = "bottom",
        panel.background = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
        axis.text.y = element_text(size = 14),
        strip.text.x = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        axis.title=element_text(size=14,face="bold"))
figsxb

figsx_no_legend <- plot_grid(figsxa + theme(legend.position = "none"), 
                             figsxb + theme(legend.position = "none"), 
                             labels = "AUTO", nrow = 2)
legend <- get_legend(
  figsxb + theme(legend.box.margin = margin(12, 0, 0, 0))
)

figsx <- plot_grid(figsx_no_legend, legend, rel_heights = c(3, .4), nrow = 2)
figsx

ggsave(filename = "inst/extdata/results/figure SX.jpg", plot = figsx,
       units = "in", height = 6, width = 10, dpi = 300)

