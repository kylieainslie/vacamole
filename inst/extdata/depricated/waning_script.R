# looking at waning rates

t_vec <- seq(0, 365, by = 0.1)
waning <- 1 - (1 / (1 + exp(-0.05 * (t_vec - 180))))
plot(t_vec,waning,type = "l")
