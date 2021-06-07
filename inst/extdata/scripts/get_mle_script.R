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
beta_range <- seq(0.00001, 0.002, by = 0.00001)

breakpoints <- yday(c(as.Date("2021-03-07"), 
                      as.Date("2021-03-21"),
                      as.Date("2021-04-01"),
                      as.Date("2021-04-17"),
                      as.Date("2021-05-25")
                      )
                    ) - yday(as.Date("2021-01-31"))
lik_out <- data.frame(beta = beta_range, 
                      lnlike_t1 = c(rep(NA, length(beta_range))),
                      lnlike_t2 = c(rep(NA, length(beta_range))),
                      lnlike_t3 = c(rep(NA, length(beta_range))),
                      lnlike_t4 = c(rep(NA, length(beta_range))),
                      lnlike_t5 = c(rep(NA, length(beta_range))))


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
    params$keep_cm_fixed <- FALSE
  } else {
    times <- seq(breakpoints[j-1], breakpoints[j], by = 1)
    init <- c(t = times[1],
            S = as.numeric(tail(out$S,1)),
            Shold_1d = as.numeric(tail(out$Shold_1d,1)),
            Sv_1d = as.numeric(tail(out$Sv_1d,1)),
            Shold_2d = as.numeric(tail(out$Shold_2d,1)),
            Sv_2d = as.numeric(tail(out$Sv_2d,1)),
            E = as.numeric(tail(out$E,1)),
            Ev_1d = as.numeric(tail(out$Ev_1d,1)),
            Ev_2d = as.numeric(tail(out$Ev_2d,1)),
            I = as.numeric(tail(out$I,1)),
            Iv_1d = as.numeric(tail(out$Iv_1d,1)),
            Iv_2d = as.numeric(tail(out$Iv_2d,1)),
            H = as.numeric(tail(out$H,1)),
            Hv_1d = as.numeric(tail(out$Hv_1d,1)),
            Hv_2d = as.numeric(tail(out$Hv_2d,1)),
            H_IC = as.numeric(tail(out$H_IC,1)),
            H_ICv_1d = as.numeric(tail(out$H_ICv_1d,1)),
            H_ICv_2d = as.numeric(tail(out$H_ICv_2d,1)),
            IC = as.numeric(tail(out$IC,1)),
            ICv_1d = as.numeric(tail(out$ICv_1d,1)),
            ICv_2d = as.numeric(tail(out$ICv_2d,1)),
            D = as.numeric(tail(out$D,1)),
            R = as.numeric(tail(out$R,1)),
            Rv_1d = as.numeric(tail(out$Rv_1d,1)),
            Rv_2d = as.numeric(tail(out$Rv_2d,1))
            )
    params$c_start <- params$c_relaxed
    params$keep_cm_fixed <- FALSE
  }
  osiris2 <- osiris1[times+1,]
  for (i in 1:length(beta_range)){
    print(i)
    params$beta <- beta_range[i]

    seir_out <- lsoda(init,times,age_struct_seir_ode,params) #
    seir_out <- as.data.frame(seir_out)
    out <- postprocess_age_struct_model_output(seir_out)

    daily_cases <- params$sigma * rowSums(out$E + out$Ev_1d + out$Ev_2d) * params$p_report

    lik <- sum(dpois(osiris2$inc,daily_cases,log=TRUE))
    print(lik)
    lik_out[i,j+1] <- lik
  }
}

params$beta <- lik_out$beta[45] 
seir_out <- lsoda(init,times,age_struct_seir_ode,params) #
seir_out <- as.data.frame(seir_out)
out <- postprocess_age_struct_model_output(seir_out)
daily_cases <- params$sigma * rowSums(out$E + out$Ev_1d + out$Ev_2d) * params$p_report

plot(daily_cases~times, type = "l")
points(osiris2$inc~times, col = "red", pch = 16)


last_date_in_osiris <- tail(osiris2)
write.csv(lik_out, file = paste0("mle_betas_", last_date_in_osiris, ".csv"))

mle_pos <- apply(lik_out[,-1], 2, function(x) which(x == max(x)))
# # 40, 47, 41, 46
mle_betas <- lik_out$beta[c(mle_pos)] 
mle_betas <- beta_range[c(40, 47, 41, 46)]
plot(lik_out$beta,lik_out$lnlike_t1, type = "l")
abline(v = mle_betas[1], col = "blue")
# --------------------------------------------------

# --------------------------------------------------
# run model for mles
# must do it piecewise
rtn <- list()
# t = 0 to 35
times1 <- seq(0, breakpoints[1], by = 1)
init1 <- c(t = 0,                  
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

    params$beta <- lik_out$beta[40] 
    params$keep_cm_fixed <- FALSE
    params$c_start <-params$c_lockdown
    seir_out <- lsoda(init1,times1,age_struct_seir_ode,params) #
    seir_out <- as.data.frame(seir_out)
    out1 <- postprocess_age_struct_model_output(seir_out)
    
    daily_cases1 <- params$sigma * rowSums(out1$E + out1$Ev_1d + out1$Ev_2d) * params$p_report
    rtn[[1]] <- daily_cases1
    

# t = 35 to 49
  times2 <- seq(breakpoints[1], breakpoints[2], by = 1)
  init2 <- c(t = times2[1],
             S = as.numeric(tail(out1$S,1)),
             Shold_1d = as.numeric(tail(out1$Shold_1d,1)),
             Sv_1d = as.numeric(tail(out1$Sv_1d,1)),
             Shold_2d = as.numeric(tail(out1$Shold_2d,1)),
             Sv_2d = as.numeric(tail(out1$Sv_2d,1)),
             E = as.numeric(tail(out1$E,1)),
             Ev_1d = as.numeric(tail(out1$Ev_1d,1)),
             Ev_2d = as.numeric(tail(out1$Ev_2d,1)),
             I = as.numeric(tail(out1$I,1)),
             Iv_1d = as.numeric(tail(out1$Iv_1d,1)),
             Iv_2d = as.numeric(tail(out1$Iv_2d,1)),
             H = as.numeric(tail(out1$H,1)),
             Hv_1d = as.numeric(tail(out1$Hv_1d,1)),
             Hv_2d = as.numeric(tail(out1$Hv_2d,1)),
             H_IC = as.numeric(tail(out1$H_IC,1)),
             H_ICv_1d = as.numeric(tail(out1$H_ICv_1d,1)),
             H_ICv_2d = as.numeric(tail(out1$H_ICv_2d,1)),
             IC = as.numeric(tail(out1$IC,1)),
             ICv_1d = as.numeric(tail(out1$ICv_1d,1)),
             ICv_2d = as.numeric(tail(out1$ICv_2d,1)),
             D = as.numeric(tail(out1$D,1)),
             R = as.numeric(tail(out1$R,1)),
             Rv_1d = as.numeric(tail(out1$Rv_1d,1)),
             Rv_2d = as.numeric(tail(out1$Rv_2d,1))
  )
  
  params$beta <- lik_out$beta[47]
  params$c_start <- params$c_relaxed
  params$keep_cm_fixed <- TRUE
  
  seir_out <- lsoda(init2,times2,age_struct_seir_ode,params) #
  seir_out <- as.data.frame(seir_out)
  out2 <- postprocess_age_struct_model_output(seir_out)
  
  daily_cases2 <- params$sigma * rowSums(out2$E + out2$Ev_1d + out2$Ev_2d) * params$p_report
  rtn[[2]] <- daily_cases2
  
# t = 49 to 60
  times3 <- seq(breakpoints[2], breakpoints[3], by = 1)
  init3 <- c(t = times3[1],
             S = as.numeric(tail(out2$S,1)),
             Shold_1d = as.numeric(tail(out2$Shold_1d,1)),
             Sv_1d = as.numeric(tail(out2$Sv_1d,1)),
             Shold_2d = as.numeric(tail(out2$Shold_2d,1)),
             Sv_2d = as.numeric(tail(out2$Sv_2d,1)),
             E = as.numeric(tail(out2$E,1)),
             Ev_1d = as.numeric(tail(out2$Ev_1d,1)),
             Ev_2d = as.numeric(tail(out2$Ev_2d,1)),
             I = as.numeric(tail(out2$I,1)),
             Iv_1d = as.numeric(tail(out2$Iv_1d,1)),
             Iv_2d = as.numeric(tail(out2$Iv_2d,1)),
             H = as.numeric(tail(out2$H,1)),
             Hv_1d = as.numeric(tail(out2$Hv_1d,1)),
             Hv_2d = as.numeric(tail(out2$Hv_2d,1)),
             H_IC = as.numeric(tail(out2$H_IC,1)),
             H_ICv_1d = as.numeric(tail(out2$H_ICv_1d,1)),
             H_ICv_2d = as.numeric(tail(out2$H_ICv_2d,1)),
             IC = as.numeric(tail(out2$IC,1)),
             ICv_1d = as.numeric(tail(out2$ICv_1d,1)),
             ICv_2d = as.numeric(tail(out2$ICv_2d,1)),
             D = as.numeric(tail(out2$D,1)),
             R = as.numeric(tail(out2$R,1)),
             Rv_1d = as.numeric(tail(out2$Rv_1d,1)),
             Rv_2d = as.numeric(tail(out2$Rv_2d,1))
  )
  
  params$beta <- lik_out$beta[41]
  
  seir_out <- lsoda(init3,times3,age_struct_seir_ode,params) #
  seir_out <- as.data.frame(seir_out)
  out3 <- postprocess_age_struct_model_output(seir_out)
  
  daily_cases3 <- params$sigma * rowSums(out3$E + out3$Ev_1d + out3$Ev_2d) * params$p_report
  rtn[[3]] <- daily_cases3
  
# t = 60 to 79
  times4 <- seq(breakpoints[3], breakpoints[4], by = 1)
  init4 <- c(t = times4[1],
             S = as.numeric(tail(out3$S,1)),
             Shold_1d = as.numeric(tail(out3$Shold_1d,1)),
             Sv_1d = as.numeric(tail(out3$Sv_1d,1)),
             Shold_2d = as.numeric(tail(out3$Shold_2d,1)),
             Sv_2d = as.numeric(tail(out3$Sv_2d,1)),
             E = as.numeric(tail(out3$E,1)),
             Ev_1d = as.numeric(tail(out3$Ev_1d,1)),
             Ev_2d = as.numeric(tail(out3$Ev_2d,1)),
             I = as.numeric(tail(out3$I,1)),
             Iv_1d = as.numeric(tail(out3$Iv_1d,1)),
             Iv_2d = as.numeric(tail(out3$Iv_2d,1)),
             H = as.numeric(tail(out3$H,1)),
             Hv_1d = as.numeric(tail(out3$Hv_1d,1)),
             Hv_2d = as.numeric(tail(out3$Hv_2d,1)),
             H_IC = as.numeric(tail(out3$H_IC,1)),
             H_ICv_1d = as.numeric(tail(out3$H_ICv_1d,1)),
             H_ICv_2d = as.numeric(tail(out3$H_ICv_2d,1)),
             IC = as.numeric(tail(out3$IC,1)),
             ICv_1d = as.numeric(tail(out3$ICv_1d,1)),
             ICv_2d = as.numeric(tail(out3$ICv_2d,1)),
             D = as.numeric(tail(out3$D,1)),
             R = as.numeric(tail(out3$R,1)),
             Rv_1d = as.numeric(tail(out3$Rv_1d,1)),
             Rv_2d = as.numeric(tail(out3$Rv_2d,1))
  )
  
  params$beta <- lik_out$beta[46]
  
  seir_out <- lsoda(init4,times4,age_struct_seir_ode,params) #
  seir_out <- as.data.frame(seir_out)
  out4 <- postprocess_age_struct_model_output(seir_out)
  
  daily_cases4 <- params$sigma * rowSums(out4$E + out4$Ev_1d + out4$Ev_2d) * params$p_report
  rtn[[4]] <- daily_cases4
# --------------------------------------------------

# --------------------------------------------------
# save initial conditions for forward simulation
  init_forward <- c(t = times4[length(times4)],
             S = as.numeric(tail(out4$S,1)),
             Shold_1d = as.numeric(tail(out4$Shold_1d,1)),
             Sv_1d = as.numeric(tail(out4$Sv_1d,1)),
             Shold_2d = as.numeric(tail(out4$Shold_2d,1)),
             Sv_2d = as.numeric(tail(out4$Sv_2d,1)),
             E = as.numeric(tail(out4$E,1)),
             Ev_1d = as.numeric(tail(out4$Ev_1d,1)),
             Ev_2d = as.numeric(tail(out4$Ev_2d,1)),
             I = as.numeric(tail(out4$I,1)),
             Iv_1d = as.numeric(tail(out4$Iv_1d,1)),
             Iv_2d = as.numeric(tail(out4$Iv_2d,1)),
             H = as.numeric(tail(out4$H,1)),
             Hv_1d = as.numeric(tail(out4$Hv_1d,1)),
             Hv_2d = as.numeric(tail(out4$Hv_2d,1)),
             H_IC = as.numeric(tail(out4$H_IC,1)),
             H_ICv_1d = as.numeric(tail(out4$H_ICv_1d,1)),
             H_ICv_2d = as.numeric(tail(out4$H_ICv_2d,1)),
             IC = as.numeric(tail(out4$IC,1)),
             ICv_1d = as.numeric(tail(out4$ICv_1d,1)),
             ICv_2d = as.numeric(tail(out4$ICv_2d,1)),
             D = as.numeric(tail(out4$D,1)),
             R = as.numeric(tail(out4$R,1)),
             Rv_1d = as.numeric(tail(out4$Rv_1d,1)),
             Rv_2d = as.numeric(tail(out4$Rv_2d,1))
  )  
saveRDS(init_forward, file = "init_conditions_2021-04-20.rds")

#  combine all piecewise results to plot together
times_all <- c(times1, times2, times3, times4)
cases_all <- c(rtn[[1]], rtn[[2]], rtn[[3]], rtn[[4]])
model_fit <- data.frame(time = times_all, cases = cases_all)
saveRDS(model_fit, file = "daily_cases_from_model_fit_2021-04-20.rds")

plot(osiris$inc~seq(1, dim(osiris)[1], by = 1),col="red",pch=16, 
     xlab="Time (days)",ylab="Incidence", ylim = c(0,7000)) 
lines(times_all, cases_all, col="blue")
legend("bottomright",c("Osiris Data","Model Fit"),
       col=c("red","blue"),lty=c(0,1),pch=c(16,NA))
# --------------------------------------------------

