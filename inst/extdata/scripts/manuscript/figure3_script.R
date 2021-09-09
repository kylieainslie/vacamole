# -------------------------------------
# plot of outcomes by age group
# -------------------------------------
# need to run figure1_script.R and figure2_script.R first
# -------------------------------------------------------------
# a) no waning
fig3a <- ggplot(data = all_res %>%
                   filter(outcome != "Daily Deaths",
                          Immunity == "No Waning") %>%
                   group_by(Scenario, Immunity, date, age_group, outcome) %>%
                   summarise_at(.vars = c("mle", "lower", "upper"), .funs = "sum"), 
                 aes(x = date, y = mle, fill = age_group, linetype = Scenario)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = age_group), alpha = 0.3) +
  geom_line(aes(color = age_group)) +
  labs(y = "Value", x = "Date") +
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
  facet_wrap(~outcome, scales = "free_y", nrow = 3)
fig3a

ggsave(filename = "inst/extdata/results/figure 3a.jpg", plot = fig3a,
       units = "in", height = 8, width = 8, dpi = 300)
# -------------------------------------------------------------
# b) waning
fig3b <- ggplot(data = all_res %>%
                   filter(outcome != "Daily Deaths",
                          Immunity == "Waning") %>%
                   group_by(Scenario, Immunity, date, age_group, outcome) %>%
                   summarise_at(.vars = c("mle", "lower", "upper"), .funs = "sum"), 
                 aes(x = date, y = mle, fill = age_group, linetype = Scenario)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = age_group), alpha = 0.3) +
  geom_line(aes(color = age_group)) +
  labs(y = "Value", x = "Date") +
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
  facet_wrap(~outcome, scales = "free_y", nrow = 3)
fig3b

# -------------------------------------------------------------
# combine plots
fig3_no_legend <- plot_grid(fig3a + theme(legend.position = "none"), 
                             fig3b + theme(legend.position = "none"), 
                             labels = "AUTO", nrow = 1)

legend3 <- get_legend(
  figs3a + theme(legend.box.margin = margin(0, 0, 0, 12))
)

fig3 <- plot_grid(fig3_no_legend, legend3, rel_heights = c(3, .4), nrow = 2)
fig3

# save output -------------------------------------------------
ggsave(filename = "inst/extdata/results/figure 3.jpg", plot = fig3,
       units = "in", height = 12, width = 10, dpi = 300)

