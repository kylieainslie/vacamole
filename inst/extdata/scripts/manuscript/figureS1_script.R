# -----------------------------------------------------------
# Figure S1 script
# Model fit to data with CIs
# -----------------------------------------------------------

# load packages ---------------------------------------------
library(dplyr)
library(ggplot2)
#library(cowplot)

# read in model fit data set ---------------------------------
file_date <- "2021-10-01"

model_fit <- readRDS(paste0("inst/extdata/results/model_fits/model_fit_df_", file_date, ".rds"))
# subset for period before childhood vaccination started
model_fit_sub <- model_fit %>%
  filter(date < as.Date("2021-06-23"))


p <- ggplot(data = model_fit_sub, aes(x = date, y = mle, linetype="solid")) +
  geom_point(data = model_fit_sub, aes(x = date, y = real, color = "Osiris notifications")) +
  geom_line() +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = "95% Confidence bounds"), alpha = 0.3) +
  scale_color_manual(values = c("red"),
                     labels = c("Osiris notifications")) +
  scale_fill_manual(values = c("grey70")) +
  scale_linetype_manual(values=c(1), labels = c("Model Fit")) +
  #scale_shape_manual(values=c(NA,20)) +
  labs(y = "Daily Cases", x = "Date") +
  ylim(0,NA) +
  scale_x_date(date_breaks = "1 month", date_labels = "%d %b %Y") +
  theme(legend.position = "bottom",
        panel.background = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
        axis.text.y = element_text(size = 14),
        strip.text.x = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_blank(),
        axis.title=element_text(size=14))
p

ggsave(filename = "inst/extdata/results/figure S1.jpg", plot = p,
       units = "in", height = 8, width = 10, dpi = 300)
