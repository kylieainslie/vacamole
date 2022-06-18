# ------------------------------------------------------------------------------
# Visualise scenario hub round 1 results
# ------------------------------------------------------------------------------

# load required packages -------------------------------------------------------
library(tidyr)
library(dplyr)
library(ggplot2)

# read in results df -----------------------------------------------------------
df_round1 <- readRDS("inst/extdata/results/scenario_hub/2021-05-22-rivm-vacamole.rds")

# make epiweek a factor, so correct order is maintained for plotting
df_plot <- df_round1 %>%
  mutate(epiweek = factor(epiweek, levels = c(21:52,1:20))) %>%
  # get mean and CI of samples
  group_by(scenario_id, epiweek) %>%
  summarise(med = median(value), 
            q025 = quantile(value, probs = 0.025),
            q975 = quantile(value, probs = 0.975))

# plot -------------------------------------------------------------------------
p <- ggplot(data = df_plot, aes(x = epiweek, y = med, color = scenario_id)) +
  geom_point() #+
  #geom_errorbar(aes(ymin=q025, ymax=q975), width=.1)
  #geom_ribbon(aes(ymin = q025, ymax = q975))
p
