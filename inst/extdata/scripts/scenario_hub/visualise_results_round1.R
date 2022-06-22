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
  summarise(med  = median(sum),
            q025 = quantile(sum, probs = 0.025),
            q975 = quantile(sum, probs = 0.975)) %>%
  select(date, scenario_id, target_variable, med:q975)

# summarise over 60+
df_60_plus <- df_round1 %>%
  filter(age_group %in% c("7", "8", "9")) %>%
  group_by(scenario_id, target_variable, date) %>%
  summarise(sum = sum(value)) %>%
  ungroup() %>%
  group_by(scenario_id, target_variable, date) %>%
  summarise(med  = median(sum),
            q025 = quantile(sum, probs = 0.025),
            q975 = quantile(sum, probs = 0.975)) %>%
  select(date, scenario_id, target_variable, med:q975)

# plot -------------------------------------------------------------------------
p_all <- ggplot(data = df_all_age_groups %>%
                  filter(target_variable == "inc case"),
                aes(x = date, y = med, color = scenario_id, fill = scenario_id)) +
  geom_line() +
  geom_ribbon(aes(ymin = q025, ymax = q975), alpha = 0.2)
p_all

# plot by epiweek
# p <- ggplot(data = df_round1 %>%
#               filter(target_variable == "inc case"), aes(x = epiweek, y = med, color = scenario_id)) +
#   geom_point(position = position_dodge(width = .75)) +
#   geom_pointrange(aes(ymin=q025, ymax=q975), position = position_dodge(width = .75)) +
#   scale_color_viridis(discrete = TRUE)
# p + facet_wrap(~target_variable, scales = "free")
