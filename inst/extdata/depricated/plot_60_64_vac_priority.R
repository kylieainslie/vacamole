# Plot results from SEIR compartmental model with 2 vaccine doses

# Load packages ----------------------------------------------------
library(ggplot2)
library(gridExtra)
library(dplyr)
library(cowplot)

# Combine results from model runs ----------------------------------
# non-constant FOI -------------------------------------------------
# read in individual results files
pfizer_60_64 <- readRDS("inst/extdata/results/new_pfizer_60_64_output.rds") %>%
  mutate(vac_type = "Pfizer",
         age_group = "60-64",
         hosp_admissions = c(rep(NA, 11), head(hosp_admissions,-11)))

pfizer_65_69 <- readRDS("inst/extdata/results/new_pfizer_65_69_output.rds") %>%
  mutate(vac_type = "Pfizer",
         age_group = "65-69",
         hosp_admissions = c(rep(NA, 11), head(hosp_admissions,-11)))

az_10_60_64 <- readRDS("inst/extdata/results/new_AZ_10_60_64_output.rds") %>%
  mutate(vac_type = "AstraZeneca (10% VE)",
         age_group = "60-64",
         hosp_admissions = c(rep(NA, 11), head(hosp_admissions,-11)))

az_30_60_64 <- readRDS("inst/extdata/results/new_AZ_30_60_64_output.rds") %>%
  mutate(vac_type = "AstraZeneca (30% VE)",
         age_group = "60-64",
         hosp_admissions = c(rep(NA, 11), head(hosp_admissions,-11)))

az_60_60_64 <- readRDS("inst/extdata/results/new_AZ_60_60_64_output.rds") %>%
  mutate(vac_type = "AstraZeneca (62% VE)",
         age_group = "60-64",
         hosp_admissions = c(rep(NA, 11), head(hosp_admissions,-11)))

# combine files based on strategy
pfizer_only <- pfizer_65_69 %>%
  mutate(inc_sum = incidence + pfizer_60_64$incidence,
         hosp_sum = hosp_admissions + pfizer_60_64$hosp_admissions,
         vac_type_new = "Pfizer (all)") %>%
  dplyr::select(time, vac_type_new, inc_sum, hosp_sum)

pfizer_az_10 <- pfizer_65_69 %>%
  mutate(inc_sum = incidence + az_10_60_64$incidence,
         hosp_sum = hosp_admissions + az_10_60_64$hosp_admissions,
         vac_type_new = "Astrazeneca (10% VE, 60-64), Pfizer (65-69)") %>%
  select(time, vac_type_new, inc_sum, hosp_sum) 

pfizer_az_30 <- pfizer_65_69 %>%
  mutate(inc_sum = incidence + c(az_30_60_64$incidence, 
                                 rep(az_30_60_64$incidence[dim(az_30_60_64)[1]], 
                                     dim(pfizer_65_69)[1] - dim(az_30_60_64)[1])),
         hosp_sum = hosp_admissions + c(az_30_60_64$hosp_admissions, 
                                        rep(az_30_60_64$hosp_admissions[dim(az_30_60_64)[1]], 
                                            dim(pfizer_65_69)[1] - dim(az_30_60_64)[1])),
         vac_type_new = "Astrazeneca (30% VE, 60-64), Pfizer (65-69)") %>%
  select(time, vac_type_new, inc_sum, hosp_sum)

pfizer_az_60 <- pfizer_65_69 %>%
  mutate(inc_sum = incidence + az_60_60_64$incidence,
         hosp_sum = hosp_admissions + az_60_60_64$hosp_admissions,
         vac_type_new = "Astrazeneca (62% VE, 60-64), Pfizer (65-69)") %>%
  select(time, vac_type_new, inc_sum, hosp_sum) 

# create data set for plotting
data_for_plot <- rbind(pfizer_only, pfizer_az_10, pfizer_az_30, pfizer_az_60) %>%
  mutate(date = time + as.Date("2021-01-21"))

# Make summaries ---------------------------------------------------
summary_data <- data_for_plot %>%
  group_by(vac_type_new) %>%
  summarise_at(.vars = c("inc_sum", "hosp_sum"), .funs = sum, na.rm = TRUE) %>%
  mutate(perc_diff_inc = ((inc_sum * 100)/inc_sum[3]) - 100,
         perc_diff_hosp = ((hosp_sum * 100)/hosp_sum[3]) - 100)

# Make plot --------------------------------------------------------

# colors <- c("Incidence" = "red", "Hospital Admissions" = "Green")
# dates <- seq.Date(from = as.Date("2021-01-21"), by = "day", length.out = 201)

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

cols <- gg_color_hue(4)
# non-constant FOI (figure 1) --------------------------------------
az_dose1_start <- vac_start_day$AstraZeneca_2[1] + as.Date("2021-01-21")
az_dose2_start <- vac_start_day$AstraZeneca_2[2] + as.Date("2021-01-21")
p_dose1_start <- vac_start_day$Pfizer_2[1] + as.Date("2021-01-21")
p_dose2_start <- vac_start_day$Pfizer_2[2] + as.Date("2021-01-21")


# incidence plot
g <- ggplot(data_for_plot, aes(x = date, y = inc_sum, color = vac_type_new)) +
  geom_line() +
  geom_vline(xintercept = az_dose1_start,linetype="dashed", color = "grey70") +
  geom_vline(xintercept = az_dose2_start,linetype="dashed", color = "grey70") +
  geom_vline(xintercept = p_dose1_start,linetype="dotted", color = "grey70") +
  geom_vline(xintercept = p_dose2_start,linetype="dotted", color = "grey70") +
  labs(y = "Incidence of Infections", x = "Date", color = "Vaccine Type") +
  scale_x_date(date_breaks = "2 weeks", date_labels = "%d %b") +
  scale_color_manual(values = cols) +
  #annotation_custom(tableGrob(Value), xmin=125, xmax=200, ymin=25, ymax=100) +
  theme(legend.position = "bottom",
        panel.background = element_blank(),
        axis.text.x = element_text(angle = 45))

# hospitalisation plot
g2 <- ggplot(data_for_plot, aes(x = date, y = hosp_sum, color = vac_type_new)) +
  geom_line() +
  geom_vline(xintercept = az_dose1_start,linetype="dashed", color = "grey70") +
  geom_vline(xintercept = az_dose2_start,linetype="dashed", color = "grey70") +
  geom_vline(xintercept = p_dose1_start,linetype="dotted", color = "grey70") +
  geom_vline(xintercept = p_dose2_start,linetype="dotted", color = "grey70") +
  labs(y = "Hospital Admissions", x = "Date", color = "Vaccine Type") +
  scale_x_date(date_breaks = "2 weeks", date_labels = "%d %b") +
  scale_color_manual(values=cols) +
  # scale_y_continuous(limits = c(0,5)) +
  #annotation_custom(tableGrob(Value), xmin=125, xmax=200, ymin=25, ymax=100) +
  theme(legend.position = "bottom",
        panel.background = element_blank(),
        axis.text.x = element_text(angle = 45))

g_both <- plot_grid(g + theme(legend.position="none"), 
                    g2 + theme(legend.position="none"),
                    nrow = 1)
# extract a legend that is laid out horizontally
legend_b <- get_legend(
  g + 
    guides(color = guide_legend(nrow = 1)) +
    theme(legend.position = "bottom")
)

# add the legend underneath the row we made earlier. Give it 10%
# of the height of one plot (via rel_heights).
fig3 <- plot_grid(g_both, legend_b, ncol = 1, rel_heights = c(1, .1))

ggsave(filename = "inst/extdata/results/figure3.png",
       plot = fig3,
       width = 12,
       height = 6,
       dpi = 300)
