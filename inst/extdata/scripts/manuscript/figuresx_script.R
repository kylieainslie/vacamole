# ----------------------------------------------------------------
# figure SX - vaccination coverage
# ----------------------------------------------------------------

# data wrangling -------------------------------------------------

age_dist_10 <- c(0.10319920, 0.11620856, 0.12740219, 0.12198707, 0.13083463, 
                 0.14514332, 0.12092904, 0.08807406, 0.03976755, 0.007398671)
n <- 17407585 # Dutch population size
n_vec_10 <- n * age_dist_10

# 12 plus --------------------------------------------------------
basis_12plus_long <- basis_12plus %>%
  mutate(pf_d1_9 = (pf_d1_9 * n_vec_10[9] + pf_d1_10 * n_vec_10[10])/sum(n_vec_10[9:10]),
         pf_d2_9 = (pf_d2_9 * n_vec_10[9] + pf_d2_10 * n_vec_10[10])/sum(n_vec_10[9:10]),
         mo_d1_9 = (mo_d1_9 * n_vec_10[9] + mo_d1_10 * n_vec_10[10])/sum(n_vec_10[9:10]),
         mo_d2_9 = (mo_d2_9 * n_vec_10[9] + mo_d2_10 * n_vec_10[10])/sum(n_vec_10[9:10]),
         az_d1_9 = (az_d1_9 * n_vec_10[9] + az_d1_10 * n_vec_10[10])/sum(n_vec_10[9:10]),
         az_d2_9 = (az_d2_9 * n_vec_10[9] + az_d2_10 * n_vec_10[10])/sum(n_vec_10[9:10]),
         ja_d1_9 = (ja_d1_9 * n_vec_10[9] + ja_d1_10 * n_vec_10[10])/sum(n_vec_10[9:10]),
         ja_d2_9 = (ja_d2_9 * n_vec_10[9] + ja_d2_10 * n_vec_10[10])/sum(n_vec_10[9:10])
  ) %>%
  select(date, pf_d1_1:ja_d2_9, -pf_d1_10, -pf_d2_10, -mo_d1_10, -mo_d2_10, -az_d1_10, 
         -az_d2_10, -ja_d1_10, -ja_d2_10) %>%
  pivot_longer(!date, 
               names_to = c("vaccine", "dose", "age_group"),
               names_sep = "_",
               values_to = "coverage") %>%
  mutate(Dose = as.factor(case_when(
    dose == "d1" ~ 1,
    dose == "d2" ~ 2
  )),
  date = as.Date(date, format = "%m/%d/%Y"),
  coverage = ifelse(coverage == 0, NA, coverage),
  Vaccine = case_when(
    vaccine == "az" ~ "AstraZeneca",
    vaccine == "pf" ~ "Pfizer/BioNTech",
    vaccine == "mo" ~ "Moderna",
    vaccine == "ja" ~ "Janssen"
  ),
  age_group = case_when(
    age_group == 1 ~ "0-9",
    age_group == 2 ~ "10-19",
    age_group == 3 ~ "20-29",
    age_group == 4 ~ "30-39",
    age_group == 5 ~ "40-49",
    age_group == 6 ~ "50-59",
    age_group == 7 ~ "60-69",
    age_group == 8 ~ "70-79",
    age_group == 9 ~ "80+"
  )) %>%
  filter(date > as.Date("2020-12-31"))

# plot
figsxa <- ggplot(data = basis_12plus_long, 
                 aes(x = date, y = coverage, color = Vaccine, linetype = Dose)) +
  geom_line() +
  labs(y = "Vaccine Coverage", x = "Date") +
  ylim(0,1) +
  scale_x_date(date_breaks = "3 months", date_labels = "%d %b %Y") +
  facet_grid(.~age_group) +
  theme(#legend.position = "none",
    panel.background = element_blank(),
    axis.text.x = element_blank(), #element_text(angle = 45, hjust = 1, size = 14),
    axis.title.x = element_blank(),
    axis.text.y = element_text(size = 14),
    strip.text.x = element_text(size = 14),
    #legend.text = element_text(size = 14),
    #legend.title = element_text(size = 14),
    axis.title=element_text(size=14,face="bold"))
figsxa

# 18 plus --------------------------------------------------------
basis_18plus_long <- basis_18plus %>%
  mutate(pf_d1_9 = (pf_d1_9 * n_vec_10[9] + pf_d1_10 * n_vec_10[10])/sum(n_vec_10[9:10]),
         pf_d2_9 = (pf_d2_9 * n_vec_10[9] + pf_d2_10 * n_vec_10[10])/sum(n_vec_10[9:10]),
         mo_d1_9 = (mo_d1_9 * n_vec_10[9] + mo_d1_10 * n_vec_10[10])/sum(n_vec_10[9:10]),
         mo_d2_9 = (mo_d2_9 * n_vec_10[9] + mo_d2_10 * n_vec_10[10])/sum(n_vec_10[9:10]),
         az_d1_9 = (az_d1_9 * n_vec_10[9] + az_d1_10 * n_vec_10[10])/sum(n_vec_10[9:10]),
         az_d2_9 = (az_d2_9 * n_vec_10[9] + az_d2_10 * n_vec_10[10])/sum(n_vec_10[9:10]),
         ja_d1_9 = (ja_d1_9 * n_vec_10[9] + ja_d1_10 * n_vec_10[10])/sum(n_vec_10[9:10]),
         ja_d2_9 = (ja_d2_9 * n_vec_10[9] + ja_d2_10 * n_vec_10[10])/sum(n_vec_10[9:10])
  ) %>%
  select(date, pf_d1_1:ja_d2_9, -pf_d1_10, -pf_d2_10, -mo_d1_10, -mo_d2_10, -az_d1_10, 
         -az_d2_10, -ja_d1_10, -ja_d2_10) %>%
  pivot_longer(!date, 
               names_to = c("vaccine", "dose", "age_group"),
               names_sep = "_",
               values_to = "coverage") %>%
  mutate(Dose = as.factor(case_when(
    dose == "d1" ~ 1,
    dose == "d2" ~ 2
  )),
  date = as.Date(date, format = "%m/%d/%Y"),
  coverage = ifelse(coverage == 0, NA, coverage),
  Vaccine = case_when(
    vaccine == "az" ~ "AstraZeneca",
    vaccine == "pf" ~ "Pfizer/BioNTech",
    vaccine == "mo" ~ "Moderna",
    vaccine == "ja" ~ "Janssen"
  ),
  age_group = case_when(
    age_group == 1 ~ "0-9",
    age_group == 2 ~ "10-19",
    age_group == 3 ~ "20-29",
    age_group == 4 ~ "30-39",
    age_group == 5 ~ "40-49",
    age_group == 6 ~ "50-59",
    age_group == 7 ~ "60-69",
    age_group == 8 ~ "70-79",
    age_group == 9 ~ "80+"
  )) %>%
  filter(date > as.Date("2020-12-31"))

# plot -----------------------------------------------------------
figsxb <- ggplot(data = basis_18plus_long, 
                 aes(x = date, y = coverage, color = Vaccine, linetype = Dose)) +
  geom_line() +
  labs(y = "Vaccine Coverage", x = "Date") +
  ylim(0,1) +
  scale_x_date(date_breaks = "3 months", date_labels = "%d %b") +
  facet_wrap(~age_group, nrow = 1) +
  theme(legend.position = "bottom",
        panel.background = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        #axis.title.x = element_blank(),
        axis.text.y = element_text(size = 14),
        strip.text.x = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        axis.title=element_text(size=14,face="bold"))
figsxb

figsx_no_legend <- plot_grid(figsxa + theme(legend.position = "none"), 
                             figsxb + theme(legend.position = "none"),
                             rel_heights = c(0.75,1),
                             labels = "AUTO", nrow = 2)
legend <- get_legend(
  figsxb + theme(legend.box.margin = margin(12, 0, 0, 0))
)

figsx <- plot_grid(figsx_no_legend, legend, rel_heights = c(3, .4), nrow = 2)
figsx

ggsave(filename = "inst/extdata/results/figure SX.jpg", plot = figsx,
       units = "in", height = 8, width = 12, dpi = 300)
