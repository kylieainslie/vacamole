# -------------------------------------
# plot of outcomes by age group
# -------------------------------------
# need to run figure1_script.R and figure2_script.R first
source("inst/extdata/scripts/manuscript/figure1_script.R")
source("inst/extdata/scripts/manuscript/figure2_script.R")
library(patchwork)
# -------------------------------------------------------------
# a) no waning
dat_fig3 <- all_res %>%
  filter(outcome != "Daily Deaths",
         Immunity == "No Waning") %>%
  group_by(Variant, Scenario, date, age_group, outcome) %>%
  summarise_at(.vars = c("mle", "lower", "upper"), .funs = "sum")

fig3 <- ggplot(data = dat_fig3, 
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
  facet_grid(outcome~Variant, scales = "free_y")
fig3

ggsave(filename = "inst/extdata/results/figure 3 no inserts.jpg", plot = fig3a,
       units = "in", height = 8, width = 8, dpi = 300)

# add inserts
# helper function
annotation_custom2 <- function (grob, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, data) 
{
  layer(data = data, stat = StatIdentity, position = PositionIdentity, 
        geom = ggplot2:::GeomCustomAnn,
        inherit.aes = TRUE, params = list(grob = grob, 
                                          xmin = xmin, xmax = xmax, 
                                          ymin = ymin, ymax = ymax))
}

# plot to be inset
inset_plot <- ggplot(data = dat_fig3 %>% 
                       filter(Variant == "Delta",
                              outcome == "Daily Cases"), 
                     aes(x = date, y = mle, fill = age_group, linetype = Scenario)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = age_group), alpha = 0.3) +
  geom_line(aes(color = age_group)) +
  labs(y = "Value", x = "Date") +
  ylim(0,NA) +
  scale_x_date(limits = c(as.Date("2021-08-01"), as.Date("2021-11-01")), 
               date_breaks = "2 weeks", date_labels = "%d %b %Y") +
  theme(legend.position = "bottom",
        panel.background = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
        axis.text.y = element_text(size = 14),
        strip.text.x = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        axis.title=element_text(size=14,face="bold"))

# add inset to main plot
fig3 +
  annotation_custom2(grob=ggplotGrob(inset_plot), 
                     data = data.frame(outcome="Daily Cases", Variant = "Delta"),
                     ymin = 5000, ymax=15000, xmin=as.Date("2021-11-08"), xmax=as.Date("2022-03-25"))  

# save plot
ggsave(filename = "inst/extdata/results/figure 3.jpg", plot = fig3a,
       units = "in", height = 8, width = 8, dpi = 300)

# -------------------------------------------------------------
# b) waning - make supplemental figure
dat_figS4 <- all_res %>%
  filter(outcome != "Daily Deaths",
         Immunity == "Waning") %>%
  group_by(Variant, Scenario, date, age_group, outcome) %>%
  summarise_at(.vars = c("mle", "lower", "upper"), .funs = "sum")

fig3b <- ggplot(data = dat_fig3b, 
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
  facet_grid(outcome~Variant, scales = "free_y")
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

