# -------------------------------------
# Figure 2 script
# plot of outcomes by age group
# -------------------------------------
library(ggsci)
library(gridExtra)

# source files ------------------------------------------------
source("inst/extdata/scripts/manuscript/data_wrangling_for_figures.R")
# -------------------------------------------------------------
save_path <- "/rivm/s/ainsliek/results/impact_vac/"

# Figure 2 
dat_fig2 <- all_res_for_plot %>%
  filter(outcome == "Daily Cases",
         Immunity == "No Waning",
         date >= as.Date("2021-11-01")) %>%
  group_by(Immunity, Scenario, date, age_group, outcome) %>%
  summarise_at(.vars = c("mle", "lower", "upper"), .funs = "sum")

aaas_cols <- pal_aaas("default")(3)
my_colors <- c(aaas_cols[3], aaas_cols[2], gray.colors(7, start = 0.2, end = 0.8))

fig2 <- ggplot(data = dat_fig2, 
                 aes(x = date, y = mle, fill = age_group)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = age_group), alpha = 0.3) +
  geom_line(aes(color = age_group)) +
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
        strip.text.x = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        axis.title=element_text(size=14,face="bold")) +
  facet_grid(Scenario~., scales = "free_y")
fig2

ggsave(filename = paste0(save_path, "figure2_new_color_palette.jpg"), plot = fig3,
       units = "in", height = 8, width = 12, dpi = 300)

# -------------------------------------------------------------
# b) hosp and IC admissions is supplemental figure
dat_figS4 <- all_res_for_plot %>%
  filter(outcome %in% c("Hospital Admissions", "IC Admissions"),
         Immunity == "No Waning",
         date >= as.Date("2021-11-01")) %>%
  group_by(Immunity, Scenario, date, age_group, outcome) %>%
  summarise_at(.vars = c("mle", "lower", "upper"), .funs = "sum")

figS4 <- ggplot(data = dat_figS4, 
                 aes(x = date, y = mle, fill = age_group)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = age_group), alpha = 0.3) +
  geom_line(aes(color = age_group)) +
  scale_fill_manual(values = my_colors) +
  scale_color_manual(values = my_colors) +
  labs(y = "Value", x = "Date") +
  ylim(0,NA) +
  scale_x_date(date_breaks = "2 weeks", date_labels = "%d %b %Y") +
  guides(fill = guide_legend("Age Group"), color = guide_legend("Age Group")) +
  theme(legend.position = "bottom",
        panel.background = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
        axis.text.y = element_text(size = 14),
        strip.text.x = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        axis.title=element_text(size=14,face="bold")) +
  facet_grid(outcome~Scenario, scales = "free_y")
figS4

# add text annotation
# figS4 + 
#   annotate(
#     geom = "curve", x = as.Date("2022-03-07"), y = 500, xend = as.Date("2022-02-01"), yend = 400, 
#     curvature = .3, arrow = arrow(length = unit(2, "mm"))
#   ) +
#   annotate(geom = "text", x = as.Date("2022-03-08"), y = 500, label = "", hjust = "left")

# ggsave(filename = paste0(save_path, "figureS4.jpg"), plot = figS4,
#        units = "in", height = 8, width = 12, dpi = 300)
