# plot results -----------------------------------------------------------------------
library(ggplot2)
library(dplyr)
library(cubature)
library(EpiEstim)
library(scales)

# load results -----------------------------------------------------------------------
setwd("./../..")
results_y2o <- readRDS("results/model_results_beta03_AZ_young_to_old.rds")
results_o2y <- readRDS("results/model_results_beta03_AZ_old_to_young.rds")
results_60_plus <- readRDS("results/model_results_beta03_AZ_60_plus_first.rds")

# plot Rt ----------------------------------------------------------------------------
# determine generation interval by determining the convolution
# of the latent and infectious periods
pdf_latent <- function(x) dweibull(x, shape = 1.97, scale = 9.09)
#pdf_infectious_a <- function(y) dgamma(y, shape = 2.1, scale = 4)
pdf_infectious <- function(z) dgamma(z, shape = 2.4, scale = 4)

pdf_generation <- function(t){
  cubintegrate(function(x){pdf_latent(x = x) * pdf_infectious(z = t - x)},
               lower = c(0),
               upper = c(Inf)
  )$integral
}

vpdf_generation <- Vectorize(pdf_generation, vectorize.args = "t")
t = seq(0, 112, 1)
plot(t, vpdf_generation(t), type = "l")

# estimate Rt using incidence data and generation interval
R_t <- list()
inc_dat_list <- list()
tag <- c("old_to_young", "young_to_old", "60_plus_first")
inc_y2o <- results_y2o$new_infections
inc_o2y <- results_o2y$new_infections
inc_60p <- results_60_plus$new_infections

for (i in 1:3){
  if (i == 1){ inc_mat <- inc_o2y
  } else if (i == 2) { inc_mat <- inc_y2o
  } else { inc_mat <- inc_60p}
  
  
  dates <- as.Date(vaccines_per_day$date)
  incidence <- apply(inc_mat, 2, sum)
  inc_dat <- data.frame(dates = dates, I = incidence)
  inc_dat_list[[i]] <- data.frame(inc_dat, scenario = tag[i])
  
  res_non_parametric_si <- estimate_R(inc_dat, 
                                    method="non_parametric_si",
                                    config = make_config(list(
                                      si_distr = vpdf_generation(t)))
                                    )
  
  R_t[[i]] <- res_non_parametric_si$R %>%
    mutate(scenario = tag[i],
           date_start = inc_dat$dates[.data$t_start],
           date_end = inc_dat$dates[.data$t_start][.data$t_end]) %>%
    rename(r_mean = .data$`Mean(R)`, r_q2.5 = .data$`Quantile.0.025(R)`,
           r_q97.5 = .data$`Quantile.0.975(R)`,
           r_median = .data$`Median(R)`) %>%
    select(.data$date_start, .data$date_end, .data$scenario, .data$r_mean, 
           .data$r_q2.5, .data$r_q97.5)
}

# bind rows to create single data frame 
r_dat <- bind_rows(R_t)
inc_dat_all <- bind_rows(inc_dat_list)

# plot Rt
p <- ggplot(data = r_dat, aes(x = date_start, y = r_mean, group = scenario)) +
  geom_line(aes(color = scenario)) +
  geom_ribbon(aes(ymin = r_q2.5, ymax = r_q97.5), alpha = 0.1) +
  xlab("Date") + ylab("") +
  ylim(0, 3) +
  scale_x_date(labels = date_format("%Y-%m-%d")) +
  geom_hline(yintercept = 1, linetype="dashed", colour = "black") +
  #geom_hline(yintercept = 0, colour = "black") +
  guides(color=guide_legend(nrow=1,byrow=TRUE)) +
  theme(legend.position = "bottom",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave("results/Rt_plot.pdf", p)

# plot incidence --------------------------------------------------------------------
inc_dat_all <- inc_dat_all %>%
  group_by(scenario) %>%
  mutate(perc_final_size = I / sum(I)) %>%
  ungroup()

p2 <- ggplot(data = inc_dat_all, aes(x = dates, y = I, fill = scenario)) +
  geom_bar(stat="identity", position=position_dodge()) +
  xlab("Date") + ylab("Incidence") +
  scale_x_date(labels = date_format("%Y-%m-%d")) +
  theme(legend.position = "bottom",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave("results/incidence_plot.pdf", p2)

p3 <- ggplot(data = inc_dat_all, aes(x = dates, y = perc_final_size, fill = scenario)) +
  geom_bar(stat="identity", position=position_dodge()) +
  xlab("Date") + ylab("Percentage of Final Size") +
  scale_x_date(labels = date_format("%Y-%m-%d")) +
  theme(legend.position = "bottom",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave("results/perc_final_size_plot.pdf", p2)
