# ------------------------------------------------------------------------------
# Visualise scenario hub round 2 results
# ------------------------------------------------------------------------------

# load required packages -------------------------------------------------------
library(tidyr)
library(dplyr)
library(ggplot2)

# read in results df -----------------------------------------------------------
df_round2 <- readRDS("inst/extdata/results/scenario_hub/2022-07-24-rivm-vacamole.rds")
#df_round2_no_boost <- readRDS("inst/extdata/results/scenario_hub/2022-07-24-rivm-vacamole_no_boost.rds")
# summarise over all age groups
df_all <- df_round2 %>%
  #bind_rows(., df_round2_no_boost) %>%
  group_by(scenario_id, target_variable, date, sample) %>%
  summarise(sum = sum(value)) %>%
  ungroup() %>%
  mutate(target_variable = factor(target_variable,
                           levels = c("inc infection", "inc case",
                                      "inc hosp", "inc icu", "inc death")),
         sample = factor(sample))

df_summary <- df_all %>%
  group_by(scenario_id, target_variable, date) %>%
  summarise(mean  = median(sum),
            q025 = quantile(sum, probs = 0.025),
            q25  = quantile(sum, probs = 0.25),
            q75  = quantile(sum, probs = 0.75),
            q975 = quantile(sum, probs = 0.975)
            ) %>%
  select(date, scenario_id, target_variable, mean:q975) #%>%
  # mutate(VE = case_when(
  #   scenario_id == "A-2022-07-24" ~ "Optimistic",
  #   scenario_id == "B-2022-07-24" ~ "Optimistic",
  #   scenario_id == "C-2022-07-24" ~ "Pessimistic",
  #   scenario_id == "D-2022-07-24" ~ "Pessimistic"
  # ))

# plot -------------------------------------------------------------------------
# mean line with ribbon 
p_ribbon <- ggplot(data = df_summary %>%
                  filter(target_variable %in% c("inc infection")#,
                         #date < as.Date("2022-10-01")
                         ), # "inc hosp", "inc icu", "inc death"
                 aes(x = date, y = mean, color = scenario_id, fill = scenario_id)) +
  geom_line() +
  geom_ribbon(aes(ymin = q025, ymax = q975), alpha = 0.1, color = NA) +
  #geom_ribbon(aes(ymin = q25, ymax = q75), alpha = 0.15, color = NA) +
  xlab("Date") + 
  ylab("Mean value") +
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
  facet_grid(.~target_variable) +
  # annotate("rect", xmin = as.Date("2022-09-22"), xmax = as.Date("2022-12-15"), ymin = 0, ymax = 200000, 
  #          alpha = .5)
  geom_vline(xintercept = as.Date("2022-09-15"), linetype = "dashed", color = "grey70")
p_ribbon

ggsave(filename = "/rivm/s/ainsliek/results/scenario_hub/round2/case_plot_round2.jpg", 
       plot = p_ribbon,
       units = "in", height = 8, width = 13, dpi = 300)

# individual lines
p_lines <- ggplot(data = df_all %>%
         filter(target_variable %in% c("inc case"),
                sample %in% sample(.data$sample, 9)
                #scenario_id %in% c("A-2022-07-24", "B-2022-07-24")
                ), #"inc hosp", "inc icu", "inc death"
       aes(x = date, 
           y = sum, 
           #group = sample, 
           color = scenario_id)) +
  geom_line() +
  facet_wrap(~sample)
p_lines

  xlab("Date") + 
  ylab("Mean value") +
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
# annotate("rect", xmin = as.Date("2022-09-22"), xmax = as.Date("2022-12-15"), ymin = 0, ymax = 200000, 
#          alpha = .5)
p_lines +
  geom_vline(xintercept = as.Date("2022-09-15"), linetype = "dashed", color = "grey70")

# percent reduction ------------------------------------------------------------
# total population
peak_inc <- df_all %>%
  group_by(scenario_id, target_variable, sample) %>%
  summarise(max = max(sum)) %>%
  ungroup() %>%
  group_by(scenario_id, target_variable) %>%
  summarise(max_mean = mean(max),
            max_q025 = quantile(max, probs = 0.025),
            max_q975 = quantile(max, probs = 0.975))

peak_incAB <- #peak_inc %>%
  df_all %>%
  group_by(scenario_id, target_variable, sample) %>%
  filter(scenario_id %in% c("A-2022-07-24", "B-2022-07-24")) %>%
  summarise(max = max(sum)) %>%
  mutate(perc_change = scales::percent((max - max[scenario_id == 'A-2022-07-24'])/max[scenario_id == 'A-2022-07-24']))

peak_incCD <- peak_inc %>%
  filter(scenario_id %in% c("C-2022-07-24", "D-2022-07-24")) %>%
  group_by(target_variable) %>%
  mutate(perc_change_mean = scales::percent((max_mean - max_mean[scenario_id == 'C-2022-07-24'])/max_mean[scenario_id == 'C-2022-07-24']),
         perc_change_q025 = scales::percent((max_q025 - max_q025[scenario_id == 'C-2022-07-24'])/max_q025[scenario_id == 'C-2022-07-24']),
         perc_change_q975 = scales::percent((max_q975 - max_q975[scenario_id == 'C-2022-07-24'])/max_q975[scenario_id == 'C-2022-07-24']))
