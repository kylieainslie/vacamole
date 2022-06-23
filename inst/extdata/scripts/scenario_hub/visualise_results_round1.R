# ------------------------------------------------------------------------------
# Visualise scenario hub round 1 results
# ------------------------------------------------------------------------------

# load required packages -------------------------------------------------------
library(tidyr)
library(dplyr)
library(ggplot2)

# read in results df -----------------------------------------------------------
df_round1 <- readRDS("inst/extdata/results/scenario_hub/2022-05-22-rivm-vacamole.rds")

# summarise over all age groups
df_all_age_groups <- df_round1 %>%
  group_by(scenario_id, target_variable, date, sample) %>%
  summarise(sum = sum(value)) %>%
  ungroup() %>%
  group_by(scenario_id, target_variable, date) %>%
  summarise(mean  = mean(sum),
            q025 = quantile(sum, probs = 0.025),
            q975 = quantile(sum, probs = 0.975)) %>%
  select(date, scenario_id, target_variable, mean:q975_max)

# summarise over 60+
df_60_plus <- df_round1 %>%
  filter(age_group %in% c("7", "8", "9")) %>%
  group_by(scenario_id, target_variable, date, sample) %>%
  summarise(sum = sum(value)) %>%
  ungroup() %>%
  group_by(scenario_id, target_variable, date) %>%
  summarise(mean  = mean(sum),
            q025 = quantile(sum, probs = 0.025),
            q975 = quantile(sum, probs = 0.975)) %>%
  select(date, scenario_id, target_variable, mean:q975)

# plot -------------------------------------------------------------------------
p_all <- ggplot(data = df_all_age_groups %>%
                  filter(target_variable %in% c("inc hosp", "inc icu", "inc death")) %>%
                  mutate(scenario_id = factor(scenario_id, 
                                              levels = c("C-2022-05-22", "D-2022-05-22",
                                                         "A-2022-05-22", "B-2022-05-22")),
                         target_variable = factor(target_variable,
                                                  levels = c("inc infection", "inc case",
                                                             "inc hosp", "inc icu", "inc death"))),
                aes(x = date, y = mean, color = scenario_id, fill = scenario_id)) +
  geom_line() +
  geom_ribbon(aes(ymin = q025, ymax = q975), alpha = 0.1, color = NA) +
  # scale_color_discrete(limits = c("A-2022-05-22", "B-2022-05-22","C-2022-05-22", "D-2022-05-22")) +
  # scale_fill_discrete(limits = c("A-2022-05-22", "B-2022-05-22","C-2022-05-22", "D-2022-05-22")) +
  xlab("Date") + 
  ylab("Mean value") +
  #ylim(0,NA) +
  scale_x_date(date_breaks = "1 month", date_labels = "%d %b %Y") +
  theme(legend.position = "bottom",
        panel.background = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
        axis.text.y = element_text(size = 14),
        strip.text.x = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.box="vertical",
        axis.title=element_text(size=14,face="bold")) +
  guides(fill=guide_legend("Scenario ID"), colour = guide_legend("Scenario ID")) +
  facet_wrap(~target_variable, nrow = 1)
p_all

# plot by epiweek
# p <- ggplot(data = df_round1 %>%
#               filter(target_variable == "inc case"), aes(x = epiweek, y = med, color = scenario_id)) +
#   geom_point(position = position_dodge(width = .75)) +
#   geom_pointrange(aes(ymin=q025, ymax=q975), position = position_dodge(width = .75)) +
#   scale_color_viridis(discrete = TRUE)
# p + facet_wrap(~target_variable, scales = "free")

# percent reduction ------------------------------------------------------------
# total population
peak_inc <- df_round1 %>%
  group_by(scenario_id, target_variable, date, sample) %>%
  summarise(sum = sum(value)) %>%
  ungroup() %>%
  group_by(scenario_id, target_variable, sample) %>%
  summarise(max = max(sum)) %>%
  ungroup() %>%
  group_by(scenario_id, target_variable) %>%
  summarise(max_mean = mean(max),
            max_q025 = quantile(max, probs = 0.025),
            max_q975 = quantile(max, probs = 0.975))

peak_incAB <- peak_inc %>%
  filter(scenario_id %in% c("A-2022-05-22", "B-2022-05-22")) %>%
  group_by(target_variable) %>%
  mutate(perc_change_mean = scales::percent((max_mean - max_mean[scenario_id == 'A-2022-05-22'])/max_mean[scenario_id == 'A-2022-05-22']),
         perc_change_q025 = scales::percent((max_q025 - max_q025[scenario_id == 'A-2022-05-22'])/max_q025[scenario_id == 'A-2022-05-22']),
         perc_change_q975 = scales::percent((max_q975 - max_q975[scenario_id == 'A-2022-05-22'])/max_q975[scenario_id == 'A-2022-05-22']))

peak_incCD <- peak_inc %>%
  filter(scenario_id %in% c("C-2022-05-22", "D-2022-05-22")) %>%
  group_by(target_variable) %>%
  mutate(perc_change_mean = scales::percent((max_mean - max_mean[scenario_id == 'C-2022-05-22'])/max_mean[scenario_id == 'C-2022-05-22']),
         perc_change_q025 = scales::percent((max_q025 - max_q025[scenario_id == 'C-2022-05-22'])/max_q025[scenario_id == 'C-2022-05-22']),
         perc_change_q975 = scales::percent((max_q975 - max_q975[scenario_id == 'C-2022-05-22'])/max_q975[scenario_id == 'C-2022-05-22']))
