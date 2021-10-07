# -----------------------------------------------------------
# Figure 2c script
# bar chart of cases by vac status
# -----------------------------------------------------------

totals <- all_res %>%
  filter(outcome == "Daily Cases",
         R0 == "4.6") %>%
  group_by(Scenario, R0, Immunity) %>%
  summarise_at(.vars = c("mle", "lower", "upper"), .funs = "sum")

data_for_bar_plot <- all_res %>%
  filter(outcome == "Daily Cases",
         R0 == "4.6") %>%
  group_by(Scenario, R0, Immunity, vac_status) %>%
  summarise_at(.vars = c("mle", "lower", "upper"), .funs = "sum") %>%
  left_join(., totals, by = c("Scenario", "R0", "Immunity")) %>%
  mutate(mle_prop = mle.x/mle.y,
         lower_prop = lower.x/lower.y,
         upper_prop = upper.x/upper.y,
         vac_status = factor(case_when(
           vac_status == "unvac" ~ "Unvaccinated",
           vac_status == "partvac" ~ "Partially Vaccinated",
           vac_status == "fullvac" ~ "Fully Vaccinated"
         ), levels = c("Unvaccinated", "Partially Vaccinated", "Fully Vaccinated")))

fig2c <- ggplot(data = data_for_bar_plot, 
                aes(x=vac_status, y=mle_prop, fill=Immunity)) +
  geom_bar(stat="identity", position=position_dodge()) +
  geom_errorbar(aes(ymin=lower_prop, ymax=upper_prop), width=.2,
                position=position_dodge(.9)) +
  ylim(0,1) +
  labs(y = "Proportion of Daily Cases", x = "Vaccination Status") +
  facet_wrap(~Scenario, nrow = 2) +
  theme(legend.position = "bottom",
        panel.background = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
        axis.text.y = element_text(size = 14),
        strip.text.x = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        axis.title=element_text(size=14,face="bold"))
fig2c

fig2 <- plot_grid(fig2ab, fig2c, labels = c("", "C"), rel_widths = c(1, 0.7), ncol = 2)
fig2

ggsave(filename = "inst/extdata/results/figure 2.jpg", plot = fig2,
       units = "in", height = 10, width = 12, dpi = 300)

