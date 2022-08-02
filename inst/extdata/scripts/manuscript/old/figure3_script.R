# -----------------------------------------------------------
# Figure 3 script
# simulated outcomes w/ and w/o waning 5+ vs. 12+ vs. 18+
# -----------------------------------------------------------

# source files ----------------------------------------------
source("inst/extdata/scripts/manuscript/data_wrangling_for_figures.R")
source("inst/extdata/scripts/manuscript/tables_script.R")
# -----------------------------------------------------------
save_path <- "/rivm/s/ainsliek/results/impact_vac/"

# Figure 3 --------------------------------------------------
# 5+ vs. 12+ vs. 18+, no waning vs. waning ------------------

dat_fig3 <- all_res_for_plot %>%
  filter(outcome == "Daily Cases",
         date >= as.Date("2021-11-01")) %>%
  group_by(Immunity, Scenario, age_group2, date, outcome) %>%
  summarise_at(.vars = c("mle", "lower", "upper"), .funs = "sum")

fig3 <- ggplot(data = dat_fig3 %>%
                  filter(Immunity == "No Waning"), 
                aes(x = date, y = mle, fill = age_group2,linetype = Scenario)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = age_group2), alpha = 0.3) +
  geom_line(aes(color = age_group2), size = 1, alpha = 0.1) +
  geom_line(data = dat_fig3 %>%
              filter(Immunity == "Waning") %>%
              mutate(age_group3 = age_group2),
            aes(x = date, y = mle, linetype = Scenario, color = age_group2), size = 1) +
  scale_linetype_manual(values = c("dotted", "dashed", "solid")) +
  labs(y = "Cases per day", x = "Date of infection") +
  ylim(0,NA) +
  scale_x_date(date_breaks = "2 weeks", date_labels = "%d %b %Y") +
  theme(legend.position = "bottom",
        panel.background = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
        axis.text.y = element_text(size = 14),
        strip.text.x = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.box="vertical",
        axis.title=element_text(size=14,face="bold")) +
  guides(fill=guide_legend("Age Group"), colour = guide_legend("Age Group"),
         linetype = guide_legend("Strategy")) #+
#facet_wrap(~age_group2, scales = "free_y", nrow = 1)
fig3

# ggsave(filename = paste0(save_path, "figure 3.jpg"), plot = fig3,
#        units = "in", height = 10, width = 12, dpi = 300)

# Supplemental figure --------------------------------------
# cases by age group with waning ---------------------------

dat_figS_waning <- all_res_for_plot %>%
  filter(outcome == "Daily Cases",
         date >= as.Date("2021-11-01")) %>%
  group_by(Immunity, Scenario, age_group, date, outcome) %>%
  summarise_at(.vars = c("mle", "lower", "upper"), .funs = "sum")

aaas_cols <- pal_aaas("default")(3)
my_colors <- c(aaas_cols[3], aaas_cols[2], gray.colors(7, start = 0.2, end = 0.8))

figS_waning <- ggplot(data = dat_figS_waning %>%
                  filter(Immunity == "Waning"), 
                aes(x = date, y = mle, fill = age_group)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = age_group), alpha = 0.3) +
  geom_line(aes(color = age_group), size = 1) +
  scale_fill_manual(values = my_colors) +
  scale_color_manual(values = my_colors) +
  labs(y = "Daily Cases", x = "Date") +
  ylim(0,NA) +
  scale_x_date(date_breaks = "2 weeks", date_labels = "%d %b %Y") +
  guides(fill = guide_legend("Age Group"), color = guide_legend("Age Group")) +
  theme(legend.position = "bottom",
        panel.background = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
        axis.text.y = element_text(size = 14),
        strip.text.y = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        axis.title=element_text(size=14,face="bold")) +
  facet_grid(Scenario~., scales = "free_y")
figS_waning

ggsave(filename = paste0(save_path, "figure_waning_by_age_group.jpg"), plot = figS_waning,
       units = "in", height = 10, width = 12, dpi = 300)


