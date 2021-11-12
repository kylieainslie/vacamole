# -----------------------------------------------------------
# Figure 1 script
# simulated outcomes in whole pop w/o waning 12+ vs. 18+
# -----------------------------------------------------------

# source files ----------------------------------------------
source("inst/extdata/scripts/manuscript/data_wrangling_for_figures.R")
source("inst/extdata/scripts/manuscript/tables_script.R")

# Figure 1 (main) -------------------------------------------
fig1a <- ggplot(data = all_res_for_plot %>%
                        filter(outcome == "Daily Cases",
                               date >= as.Date("2021-11-01"),
                               Immunity = "No Waning") %>%
                        group_by(Scenario, age_group2, date, outcome) %>%
                        summarise_at(.vars = c("mle", "lower", "upper"), .funs = "sum"), 
                aes(x = date, y = mle, fill = age_group2,linetype = Scenario)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = age_group2), alpha = 0.3) +
  geom_line(aes(color = age_group2), size = 1) +
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
         linetype = guide_legend("Strategy"))
fig1a

# Figure 1 (inset) ------------------------------------------
# bar plot of percent difference ----------------------------
fig1_inset <- ggplot(data = table1a %>%
                     filter(Scenario != "Vaccination of 18+"), 
                     aes(x = outcome, y = abs(perc_diff), fill = age_group2)) +
  geom_bar(stat = "Identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = abs(perc_diff_upper), ymax = abs(perc_diff_lower), width = 0.2),
                position = position_dodge(0.9)) +
  labs(x = "", y = "Percent Reduction (%)", fill = "Age Group") +
  scale_x_discrete(labels = c("Daily\nCases", "Hospital\nAdmissions", "IC\n Admissions")) +
  facet_wrap(~Scenario, nrow = 2) +
  theme(legend.position = "bottom",
        panel.background = element_blank(),
        axis.text.x = element_text(size = 12), #angle = 45, hjust = 1, 
        axis.text.y = element_text(size = 14),
        strip.text.x = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        axis.title=element_text(size=14,face = "bold")) +
  guides(fill=guide_legend(nrow=2,byrow=TRUE))
fig1_inset

fig1 <- fig1a + annotation_custom(ggplotGrob(fig1_inset), 
                                  xmin = as.Date("2022-02-01"), xmax = as.Date("2022-04-04"),
                                  ymin = 15000, ymax = 100000)
fig1

# save output -------------------------------------------------
# ggsave(filename = "/rivm/s/ainsliek/results/impact_vac/figure_1_w_inset.pdf", plot = fig1,
#        units = "in", height = 10, width = 12, dpi = 300)
