# R0 plot ---------------------------------------------
# -----------------------------------------------------
# days <- seq(yday(as.Date("2021-01-31")), yday(as.Date("2021-12-31")), by = 1)
# breakpoints <- yday(c(
#   as.Date("2021-02-28"),
#   as.Date("2021-03-21"),
#   as.Date("2021-04-01"),
#   as.Date("2021-04-08"),
#   as.Date("2021-04-15"),
#   as.Date("2021-04-22"),
#   as.Date("2021-04-30"),
#   as.Date("2021-05-25")
# ))
# beta_mle_vals <- c(0.00045, 0.00038, 0.00041, 0.00052, 0.00057, 0.00061, 0.00062, 0.00061)
# betas <- c(
#   rep(beta_mle_vals[1], length(days[1]:breakpoints[1])),
#   rep(beta_mle_vals[2], length(breakpoints[1]:breakpoints[2]) - 1),
#   rep(beta_mle_vals[3], length(breakpoints[2]:breakpoints[3]) - 1),
#   rep(beta_mle_vals[4], length(breakpoints[3]:breakpoints[4]) - 1),
#   rep(beta_mle_vals[5], length(breakpoints[4]:breakpoints[5]) - 1),
#   rep(beta_mle_vals[6], length(breakpoints[5]:breakpoints[6]) - 1),
#   rep(beta_mle_vals[7], length(breakpoints[6]:breakpoints[7]) - 1),
#   rep(beta_mle_vals[8], length(breakpoints[7]:days[length(days)]) - 1)
# )
# beta_t_ <- betas * (1 + 0.15 * cos(2 * pi * days / 365.24))
# R0 <- (beta_t_ / g) * rho
# R0_dat <- date <- data.frame(
#   date = seq.Date(as.Date("2021-01-31"), as.Date("2021-12-31"), by = 1),
#   R0 = R0,
#   roll_mean_R0 = zoo::rollmean(R0, k = 7, fill = NA)
# )
# r0_plot <- ggplot(data = R0_dat, aes(x = date, y = roll_mean_R0)) +
#   geom_line() +
#   labs(y = "Basic reproduction number (R0)", x = "Date") +
#   theme( # legend.position = "bottom",
#     panel.background = element_blank(),
#     axis.text = element_text(size = 14),
#     axis.text.x = element_text(size = 14),
#     axis.text.y = element_text(size = 14)
#   )
# r0_plot

# ggsave("inst/extdata/results/r0_plot_16june.jpg",
#   plot = r0_plot,
#   # height = 8,
#   # width = 12,
#   dpi = 300
# )