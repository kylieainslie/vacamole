# find mle of beta from SEIR model fit to real data

# --------------------------------------------------
# First perform a breakpoint analysis to find time 
# points over which to separately maximise likelihood
library(segmented)
library(ggplot2)
library(tidyverse)
library(lubridate)

source("inst/extdata/scrips/model_run_helper.R")
# read in OSIRIS data
osiris <- readRDS("inst/extdata/data/Osiris_Data_20210602_1041_agg.rds")
# osiris1 <- osiris %>%
#   group_by(date) %>%
#   summarise_at(.vars = "n", .funs = "sum") %>%
#   filter(date >= as.Date("2021-01-31")) %>%
#   rename(inc = n)

# remove last week of points
osiris1 <- osiris %>% 
  filter(date <= tail(date,1) - 7)

# plot data
p <- ggplot(osiris1, aes(x = date, y = inc)) +
  geom_point()
p

my_glm <- glm(inc ~ date, data = osiris1) 
my_coef <- coef(my_glm)

p1 <- p + geom_abline(intercept = my_coef[1],
                      slope = my_coef[2])
p1

# now for the actual breakpoint analysis
my_seg <- segmented(my_glm,
                    seg.Z = ~ date,
                    psi = list(date = c(as.Date("2021-03-07"), 
                                        as.Date("2021-03-21"),
                                        as.Date("2021-04-01"),
                                        as.Date("2021-04-17"))))
summary(my_seg)
# the breakpoints
my_seg$psi

# the slopes
my_slopes <- slope(my_seg)
#beta_mult <- (my_slopes$date[,1]/my_slopes$date[1,1])

#  get the fitted data
my_fitted <- fitted(my_seg)
my_model <- data.frame(Date = osiris1$date, Daily_Cases = my_fitted)
p + geom_line(data = my_model, aes(x = Date, y = Daily_Cases), color = "blue") 
# --------------------------------------------------

# --------------------------------------------------
# re-run mle for each set of time points
beta_range <- seq(0.00001, 0.001, by = 0.00001)

breakpoints <- yday(c(as.Date("2021-02-28"),
                      as.Date("2021-03-21"),
                      as.Date("2021-04-01"),
                      as.Date("2021-04-08"),
                      as.Date("2021-04-15"),
                      as.Date("2021-04-22"),
                      as.Date("2021-04-30"),
                      as.Date("2021-05-25")
                      )
                    ) - yday(as.Date("2021-01-31"))

# breakpoints <- yday(seq(osiris1$date[1], tail(osiris1$date,1), by = 14)) - yday(as.Date("2021-01-31"))
# breakpoints <- breakpoints[-1]
lik_out <- matrix(, nrow = length(beta_range), ncol = length(breakpoints))
mles <- c(rep(NA, length(breakpoints)))
out_mle <- list()
daily_cases_mle <- list()

for (j in 1:length(breakpoints)) {
  if (j == 1){
    times <- seq(0, breakpoints[j], by = 1)
    init <- c(t = 0,
              S = init_states_dat$init_S,
              Shold_1d = empty_state,
              Sv_1d = empty_state,
              Shold_2d = empty_state,
              Sv_2d = empty_state,
              E = init_states_dat$init_E,
              Ev_1d = empty_state,
              Ev_2d = empty_state,
              I = init_states_dat$init_I,
              Iv_1d = empty_state,
              Iv_2d = empty_state,
              H = init_states_dat$n_hosp,
              Hv_1d = empty_state,
              Hv_2d = empty_state,
              H_IC = init_states_dat$n_hosp_after_ic,
              H_ICv_1d = empty_state,
              H_ICv_2d = empty_state,
              IC = init_states_dat$n_ic,
              ICv_1d = empty_state,
              ICv_2d = empty_state,
              D = empty_state,
              R = init_states_dat$n_recovered,
              Rv_1d = empty_state,
              Rv_2d = empty_state
    )
    params$c_start <- params$c_lockdown
    params$keep_cm_fixed <- TRUE
  } else {
    times <- seq(breakpoints[j-1], breakpoints[j], by = 1)
    init <- c(t = times[1],
            S = as.numeric(tail(out_mle[[j-1]]$S,1)),
            Shold_1d = as.numeric(tail(out_mle[[j-1]]$Shold_1d,1)),
            Sv_1d = as.numeric(tail(out_mle[[j-1]]$Sv_1d,1)),
            Shold_2d = as.numeric(tail(out_mle[[j-1]]$Shold_2d,1)),
            Sv_2d = as.numeric(tail(out_mle[[j-1]]$Sv_2d,1)),
            E = as.numeric(tail(out_mle[[j-1]]$E,1)),
            Ev_1d = as.numeric(tail(out_mle[[j-1]]$Ev_1d,1)),
            Ev_2d = as.numeric(tail(out_mle[[j-1]]$Ev_2d,1)),
            I = as.numeric(tail(out_mle[[j-1]]$I,1)),
            Iv_1d = as.numeric(tail(out_mle[[j-1]]$Iv_1d,1)),
            Iv_2d = as.numeric(tail(out_mle[[j-1]]$Iv_2d,1)),
            H = as.numeric(tail(out_mle[[j-1]]$H,1)),
            Hv_1d = as.numeric(tail(out_mle[[j-1]]$Hv_1d,1)),
            Hv_2d = as.numeric(tail(out_mle[[j-1]]$Hv_2d,1)),
            H_IC = as.numeric(tail(out_mle[[j-1]]$H_IC,1)),
            H_ICv_1d = as.numeric(tail(out_mle[[j-1]]$H_ICv_1d,1)),
            H_ICv_2d = as.numeric(tail(out_mle[[j-1]]$H_ICv_2d,1)),
            IC = as.numeric(tail(out_mle[[j-1]]$IC,1)),
            ICv_1d = as.numeric(tail(out_mle[[j-1]]$ICv_1d,1)),
            ICv_2d = as.numeric(tail(out_mle[[j-1]]$ICv_2d,1)),
            D = as.numeric(tail(out_mle[[j-1]]$D,1)),
            R = as.numeric(tail(out_mle[[j-1]]$R,1)),
            Rv_1d = as.numeric(tail(out_mle[[j-1]]$Rv_1d,1)),
            Rv_2d = as.numeric(tail(out_mle[[j-1]]$Rv_2d,1))
            )
    params$c_start <- params$c_very_relaxed
    params$keep_cm_fixed <- FALSE
    }

  osiris2 <- osiris1[times+1,]
  for (i in 1:length(beta_range)){
    print(i)
    params$beta <- beta_range[i]

    seir_out <- lsoda(init,times,age_struct_seir_ode,params) #
    seir_out <- as.data.frame(seir_out)
    test_out <- postprocess_age_struct_model_output(seir_out)

    daily_cases <- params$sigma * rowSums(test_out$E + test_out$Ev_1d + test_out$Ev_2d) * params$p_report

    lik <- sum(dpois(osiris2$inc,daily_cases,log=TRUE))
    print(lik)
    lik_out[i,j] <- lik
  }
  
  mles[j] <- beta_range[which(lik_out[,j] == max(lik_out[,j], na.rm = TRUE))] 
              # 45, 38, 41, 47* (52), 51* (57), 55* (61), 56* (62), 62* (61) 
  params$beta <- mles[j] 
  seir_out <- lsoda(init,times,age_struct_seir_ode,params) #
  seir_out <- as.data.frame(seir_out)
  out_mle[[j]] <- postprocess_age_struct_model_output(seir_out)
  daily_cases_mle[[j]] <- params$sigma * rowSums(out_mle[[j]]$E + out_mle[[j]]$Ev_1d + out_mle[[j]]$Ev_2d) * params$p_report
  
  plot(daily_cases_mle[[j]]~times, type = "l", ylim = c(0,8000))
  points(osiris2$inc~times, col = "red", pch = 16)
  
} # end of for loop over breakpoints

last_date_in_osiris <- "2020-05-25"
saveRDS(out_mle, file = paste0("output_from_fits_", last_date_in_osiris, ".rds"))
saveRDS(daily_cases_mle, file = paste0("cases_from_fits_", last_date_in_osiris, ".rds"))

# --------------------------------------------------
#  combine all piecewise results to plot together
cases_all <- unlist(daily_cases_mle)
cases_all1 <- unique(cases_all) # remove duplicate time points
times_all <- 1:length(cases_all1)
model_fit <- data.frame(time = times_all, cases = cases_all1)
saveRDS(model_fit, file = paste0("model_fit_df_", last_date_in_osiris,".rds"))

plot(osiris1$inc~seq(1, dim(osiris1)[1], by = 1),col="red",pch=16, 
     xlab="Time (days)",ylab="Incidence", ylim = c(0,7600)) 
lines(times_all, cases_all1, col="blue")
legend("bottomright",c("Osiris Data","Model Fit"),
       col=c("red","blue"),lty=c(0,1),pch=c(16,NA))
# --------------------------------------------------

# --------------------------------------------------
# save initial conditions for forward simulation
init_forward <- c(t = tail(times_all,1) - 1,
             S = as.numeric(tail(out_mle[[8]]$S,1)),
             Shold_1d = as.numeric(tail(out_mle[[8]]$Shold_1d,1)),
             Sv_1d = as.numeric(tail(out_mle[[8]]$Sv_1d,1)),
             Shold_2d = as.numeric(tail(out_mle[[8]]$Shold_2d,1)),
             Sv_2d = as.numeric(tail(out_mle[[8]]$Sv_2d,1)),
             E = as.numeric(tail(out_mle[[8]]$E,1)),
             Ev_1d = as.numeric(tail(out_mle[[8]]$Ev_1d,1)),
             Ev_2d = as.numeric(tail(out_mle[[8]]$Ev_2d,1)),
             I = as.numeric(tail(out_mle[[8]]$I,1)),
             Iv_1d = as.numeric(tail(out_mle[[8]]$Iv_1d,1)),
             Iv_2d = as.numeric(tail(out_mle[[8]]$Iv_2d,1)),
             H = as.numeric(tail(out_mle[[8]]$H,1)),
             Hv_1d = as.numeric(tail(out_mle[[8]]$Hv_1d,1)),
             Hv_2d = as.numeric(tail(out_mle[[8]]$Hv_2d,1)),
             H_IC = as.numeric(tail(out_mle[[8]]$H_IC,1)),
             H_ICv_1d = as.numeric(tail(out_mle[[8]]$H_ICv_1d,1)),
             H_ICv_2d = as.numeric(tail(out_mle[[8]]$H_ICv_2d,1)),
             IC = as.numeric(tail(out_mle[[8]]$IC,1)),
             ICv_1d = as.numeric(tail(out_mle[[8]]$ICv_1d,1)),
             ICv_2d = as.numeric(tail(out_mle[[8]]$ICv_2d,1)),
             D = as.numeric(tail(out_mle[[8]]$D,1)),
             R = as.numeric(tail(out_mle[[8]]$R,1)),
             Rv_1d = as.numeric(tail(out_mle[[8]]$Rv_1d,1)),
             Rv_2d = as.numeric(tail(out_mle[[8]]$Rv_2d,1))
  )  
saveRDS(init_forward, file = paste0("init_conditions_", last_date_in_osiris,".rds"))

# --------------------------------------------------

