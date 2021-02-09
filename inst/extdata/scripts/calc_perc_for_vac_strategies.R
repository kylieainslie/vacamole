
n <- 17407585 
age_dist <- c(0.10319920, 0.11620856, 0.12740219, 0.12198707, 0.13083463, 
              0.14514332, 0.12092904, 0.08807406, 0.04622194)
n_vec <- n * age_dist
n_vec[2] <- n_vec[2] * 0.2 # only 18 and 19 year olds will be vaccinated
uptake <- 0.85
perc_already_vaccinated <- c(0, 0.0737, 0.0737, 0.0736, 0.0822, 0.0773, 1, 1, 1)
vac_per_week <- 800000

num_weeks <- (n_vec*(uptake - perc_already_vaccinated)) / 800000
perc_per_week <- 1/num_weeks
num_days <- num_weeks*7
perc_per_day <- perc_per_week/7
perc_per_day[c(1,7:9)] <- NA
num_days[c(1,7:9)] <- NA

df <- data.frame(age_group = 1:9,
                 n = n_vec,
                 perc_per_day = perc_per_day,
                 num_days = num_days)

((vac_per_week - (0.02268998 * df$n[5]))/800000) * df$perc_per_day[4]
(df$num_days[5] * df$perc_per_day[5]) - (0.06484154 + 13*df$perc_per_day[4])

# when 18-30 treated as single group
num_weeks_18_30 <- sum((n_vec[c(2:3)]*(uptake - perc_already_vaccinated[c(2:3)])) / 800000)
perc_per_week_18_30 <- 1/num_weeks_18_30
num_days_18_30 <- num_weeks_18_30*7
perc_per_day_18_30 <- perc_per_week_18_30/7
