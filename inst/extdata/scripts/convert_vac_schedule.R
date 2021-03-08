# convert cumulative vaccination schedule to non-cumulative ---------------------------------------
library(readr)
library(tidyr)
library(dplyr)
library(ggplot2)
library(cowplot)
# read in vac schedule file -----------------------------------------------------------------------
cum_vac_schedule_orig <- read_csv("inst/extdata/data/Cum_upt_B_parallel_20210225.csv")
cum_vac_schedule_no2dose<- read_csv("inst/extdata/data/Cum_upt_B_parallel_20210225 scno2nddose.csv")
cum_vac_schedule_delay3mo <- read_csv("inst/extdata/data/Cum_upt_B_parallel_20210225 sc2nddose3mo.csv")

#cum_vac_schedule_orig_w_jansen <- read_csv("inst/extdata/data/Cum_upt_B_parallel_20210225 + Jansen.csv")
# to combine age groups 9 and 10
age_dist_10 <- c(0.10319920, 0.11620856, 0.12740219, 0.12198707, 0.13083463, 
              0.14514332, 0.12092904, 0.08807406, 0.03976755, 0.007398671)
n <- 17407585 # Dutch population size
n_vec_10 <- n * age_dist_10

# take the difference for each row
# original
vac_schedule_orig <- data.frame(diff(as.matrix(cum_vac_schedule_orig[-1,-1]))) %>%
  add_row(cum_vac_schedule_orig[1,-1],.before = 1) %>%
  mutate(date = seq.Date(from = as.Date("2021-01-04"), to = as.Date("2021-12-30"), by = 1),
         pf_d1_9 = (pf_d1_9 * n_vec_10[9] + pf_d1_10 * n_vec_10[10])/sum(n_vec_10[9:10]),
         pf_d2_9 = (pf_d2_9 * n_vec_10[9] + pf_d2_10 * n_vec_10[10])/sum(n_vec_10[9:10]),
         mo_d1_9 = (mo_d1_9 * n_vec_10[9] + mo_d1_10 * n_vec_10[10])/sum(n_vec_10[9:10]),
         mo_d2_9 = (mo_d2_9 * n_vec_10[9] + mo_d2_10 * n_vec_10[10])/sum(n_vec_10[9:10]),
         az_d1_9 = (az_d1_9 * n_vec_10[9] + az_d1_10 * n_vec_10[10])/sum(n_vec_10[9:10]),
         az_d2_9 = (az_d2_9 * n_vec_10[9] + az_d2_10 * n_vec_10[10])/sum(n_vec_10[9:10]),
         # ja_d1_9 = (Ja_d1_9 * n_vec_10[9] + Ja_d1_10 * n_vec_10[10])/sum(n_vec_10[9:10]),
         # ja_d2_9 = (Ja_d2_9 * n_vec_10[9] + Ja_d2_10 * n_vec_10[10])/sum(n_vec_10[9:10])
         ) %>%
  select(date, pf_d1_1:az_d2_9, -pf_d1_10, -pf_d2_10, -mo_d1_10, -mo_d2_10, -az_d1_10 #, 
         #-Ja_d1_10, -Ja_d2_10
         ) 

# no second dose
vac_schedule_no2dose <- data.frame(diff(as.matrix(cum_vac_schedule_no2dose[-1,-1]))) %>%
  add_row(cum_vac_schedule_no2dose[1,-1],.before = 1) %>%
  mutate(date = seq.Date(from = as.Date("2021-01-04"), to = as.Date("2021-12-30"), by = 1),
         pf_d1_9 = (pf_d1_9 * n_vec_10[9] + pf_d1_10 * n_vec_10[10])/sum(n_vec_10[9:10]),
         pf_d2_9 = (pf_d2_9 * n_vec_10[9] + pf_d2_10 * n_vec_10[10])/sum(n_vec_10[9:10]),
         mo_d1_9 = (mo_d1_9 * n_vec_10[9] + mo_d1_10 * n_vec_10[10])/sum(n_vec_10[9:10]),
         mo_d2_9 = (mo_d2_9 * n_vec_10[9] + mo_d2_10 * n_vec_10[10])/sum(n_vec_10[9:10]),
         az_d1_9 = (az_d1_9 * n_vec_10[9] + az_d1_10 * n_vec_10[10])/sum(n_vec_10[9:10]),
         az_d2_9 = (az_d2_9 * n_vec_10[9] + az_d2_10 * n_vec_10[10])/sum(n_vec_10[9:10])
  ) %>%
  select(date, pf_d1_1:az_d2_9, -pf_d1_10, -pf_d2_10, -mo_d1_10, -mo_d2_10, -az_d1_10) %>%
  filter(date > as.Date("2021-03-08"))

# delay second dose for 3 months
vac_schedule_delay3mo <- data.frame(diff(as.matrix(cum_vac_schedule_delay3mo[-1,-1]))) %>%
  add_row(cum_vac_schedule_delay3mo[1,-1],.before = 1) %>%
  mutate(date = seq.Date(from = as.Date("2021-01-04"), to = as.Date("2021-12-30"), by = 1),
         pf_d1_9 = (pf_d1_9 * n_vec_10[9] + pf_d1_10 * n_vec_10[10])/sum(n_vec_10[9:10]),
         pf_d2_9 = (pf_d2_9 * n_vec_10[9] + pf_d2_10 * n_vec_10[10])/sum(n_vec_10[9:10]),
         mo_d1_9 = (mo_d1_9 * n_vec_10[9] + mo_d1_10 * n_vec_10[10])/sum(n_vec_10[9:10]),
         mo_d2_9 = (mo_d2_9 * n_vec_10[9] + mo_d2_10 * n_vec_10[10])/sum(n_vec_10[9:10]),
         az_d1_9 = (az_d1_9 * n_vec_10[9] + az_d1_10 * n_vec_10[10])/sum(n_vec_10[9:10]),
         az_d2_9 = (az_d2_9 * n_vec_10[9] + az_d2_10 * n_vec_10[10])/sum(n_vec_10[9:10])
  ) %>%
  select(date, pf_d1_1:az_d2_9, -pf_d1_10, -pf_d2_10, -mo_d1_10, -mo_d2_10, -az_d1_10) %>%
  filter(date > as.Date("2021-03-08"))


# put the allocation up to 8 March from original schedule into new delay schedules
vac_schedule_orig_8March <- vac_schedule_orig %>%
  filter(date <= as.Date("2021-03-08"))

# filter for dates before 1 February 2021
before_feb <- vac_schedule_orig %>%
  filter(date < as.Date("2021-02-01")) %>%
  select(-date) %>%
  summarise_all(sum) %>%
  mutate(date = as.Date("2021-01-31")) %>%
  select(date, pf_d1_1:az_d2_9)

# filter so start date is 1 Feb 2021 with row for Jan 31 to reflect people who have already been vaccinated
vac_schedule_orig_new <- vac_schedule_orig %>%
  filter(date > as.Date("2021-01-31")) %>%
  add_row(before_feb, .before = 1)
vac_schedule_no2dose_new <- bind_rows(vac_schedule_orig_8March, vac_schedule_no2dose) %>%
  filter(date > as.Date("2021-01-31")) %>%
  add_row(before_feb, .before = 1)
vac_schedule_delay3mo_new <- bind_rows(vac_schedule_orig_8March, vac_schedule_delay3mo) %>%
  filter(date > as.Date("2021-01-31")) %>%
  add_row(before_feb, .before = 1)

# calculate number of doses administered
end_date <- as.Date("2021-07-01")
# convert proportion vaccinated into number vaccinated
num_vac_orig <- sweep(vac_schedule_orig_new[,-1], 2, n_vec, "*") %>%
  mutate(date = as.Date(vac_schedule_orig_new$date, "%m/%d/%Y")) %>%
  select(date, pf_d1_1:az_d2_9) %>%
  filter(date <= end_date) %>%
  summarise_at(vars(pf_d1_1:az_d2_9), sum)
  
num_vac_no2dose <- sweep(vac_schedule_no2dose_new[,-1], 2, n_vec, "*") %>%
  mutate(date = as.Date(vac_schedule_no2dose_new$date, "%m/%d/%Y")) %>%
  select(date, pf_d1_1:az_d2_9) %>%
  filter(date <= end_date) %>%
  summarise_at(vars(pf_d1_1:az_d2_9), sum)

num_vac_delay3mo <- sweep(vac_schedule_delay3mo_new[,-1], 2, n_vec, "*") %>%
  mutate(date = as.Date(vac_schedule_delay3mo_new$date, "%m/%d/%Y")) %>%
  select(date, pf_d1_1:az_d2_9) %>%
  filter(date <= end_date) %>%
  summarise_at(vars(pf_d1_1:az_d2_9), sum)

num_vac_all <- bind_rows(num_vac_orig,
                         num_vac_delay3mo,
                         num_vac_no2dose,
                         .id = "scenario") %>%
  mutate(scenario = case_when(
    scenario == 1 ~ "original",
    scenario == 2 ~ "delay 3 months",
    scenario == 3 ~ "no 2nd doses"
  )) %>%
  # select(-date) %>%
  pivot_longer(cols = pf_d1_1:az_d2_9,names_to = c("vaccine", "dose", "age_group"),
               names_pattern = "(.*)_(.*)_(.*)", values_to = "n") %>%
  mutate(vaccine = case_when(
    vaccine == "pf" ~ "pfizer",
    vaccine == "mo" ~ "moderna",
    vaccine == "az" ~ "astrazeneca"),
    dose = case_when(
      dose == "d1" ~ 1,
      dose == "d2" ~ 2),
    dose = as.factor(dose),
    age_group = case_when(
      age_group == 1 ~ "0-9",
      age_group == 2 ~ "10-19",
      age_group == 3 ~ "20-29",
      age_group == 4 ~ "30-39",
      age_group == 5 ~ "40-49",
      age_group == 6 ~ "50-59",
      age_group == 7 ~ "60-69",
      age_group == 8 ~ "70-79",
      age_group == 9 ~ "80+"#,
      #age_group == 10 ~ "90+"
    )#,
    #n_divide_100000 = n / 100000
  )

ymax = max(num_vac_all$n)

p_orig <- ggplot(data = num_vac_all %>% filter(scenario == "original"), 
                 aes(x = age_group, y = n, fill = dose)) +
  geom_bar(stat = "identity", position = 'dodge') +
  scale_y_continuous(limits = c(0,ymax)) +
  #geom_text(aes(label = n), position = position_dodge(width=0.9), vjust=-0.25) +
  facet_grid(vars(vaccine), scales = "free") 

p_delay3mo <- ggplot(data = num_vac_all %>% filter(scenario == "delay 3 months"), 
                 aes(x = age_group, y = n, fill = dose)) +
  geom_bar(stat = "identity", position = 'dodge') +
  scale_y_continuous(limits = c(0,ymax)) +
  facet_grid(vars(vaccine), scales = "free")

p_no2dose <- ggplot(data = num_vac_all %>% filter(scenario == "no 2nd doses"), 
                 aes(x = age_group, y = n, fill = dose)) +
  geom_bar(stat = "identity", position = 'dodge') +
  scale_y_continuous(limits = c(0,ymax)) +
  facet_grid(vars(vaccine), scales = "free")


# arrange the three plots in a single row
prow <- plot_grid(
  p_orig     + theme(legend.position="none"),
  p_delay3mo + theme(legend.position="none"),
  p_no2dose  + theme(legend.position="none"),
  align = 'vh',
  labels = c("A", "B", "C"),
  hjust = -1,
  nrow = 1
)

# extract the legend from one of the plots
legend <- get_legend(
  # create some space to the left of the legend
  p_orig + theme(legend.box.margin = margin(0, 0, 0, 12))
)

# add the legend to the row we made earlier. Give it one-third of 
# the width of one plot (via rel_widths).
p <- plot_grid(prow, legend, rel_widths = c(3, .4))
ggsave("inst/extdata/results/vac_allocation_by_age_group.jpg",
       plot = p,
       height = 6,
       width = 15,
       dpi = 300)


# summarise over age groups
totals <- num_vac_all %>%
  group_by(scenario, vaccine, dose) %>%
  summarise_at(.vars = "n", sum)
#write_csv(totals, "")

ymax2 <- max(totals$n)

p_totals_orig <- ggplot(data = totals %>% filter(scenario == "original"), 
                        aes(x = vaccine, y = n, fill = dose)) +
  geom_bar(stat = "identity", position = 'dodge') +
  scale_y_continuous(limits = c(0,ymax2)) +
  geom_text(aes(label = round(n,0)), position = position_dodge(width=0.9), vjust=-0.25)
  # facet_grid(vars(vaccine), scales = "free")

p_totals_delay3mo <- ggplot(data = totals %>% filter(scenario == "delay 3 months"), 
                        aes(x = vaccine, y = n, fill = dose)) +
  geom_bar(stat = "identity", position = 'dodge') +
  scale_y_continuous(limits = c(0,ymax2)) +
  geom_text(aes(label = round(n,0)), position = position_dodge(width=0.9), vjust=-0.25)

p_totals_no2dose <- ggplot(data = totals %>% filter(scenario == "no 2nd doses"), 
                        aes(x = vaccine, y = n, fill = dose)) +
  geom_bar(stat = "identity", position = 'dodge') +
  scale_y_continuous(limits = c(0,ymax2)) +
  geom_text(aes(label = round(n,0)), position = position_dodge(width=0.9), vjust=-0.25)

# arrange the three plots in a single row
prow2 <- plot_grid(
  p_totals_orig     + theme(legend.position="none"),
  p_totals_delay3mo + theme(legend.position="none"),
  p_totals_no2dose  + theme(legend.position="none"),
  align = 'vh',
  labels = c("A", "B", "C"),
  hjust = -1,
  nrow = 1
)

# extract the legend from one of the plots
legend2 <- get_legend(
  # create some space to the left of the legend
  p_totals_orig + theme(legend.box.margin = margin(0, 0, 0, 12))
)

# add the legend to the row we made earlier. Give it one-third of 
# the width of one plot (via rel_widths).
p2 <- plot_grid(prow2, legend2, rel_widths = c(3, .4))
ggsave("inst/extdata/results/vac_allocation_totals.jpg",
       plot = p2,
       height = 6,
       width = 15,
       dpi = 300)
