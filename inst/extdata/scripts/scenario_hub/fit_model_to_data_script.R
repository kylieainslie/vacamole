# ------------------------------------------------------------------
# Fit model to data script
# ------------------------------------------------------------------
# Fit model with vaccination (4 doses) and waning to data
# ------------------------------------------------------------------

# Load required packages -------------------------------------------
library(deSolve)
library(lubridate)
library(dplyr)
library(tidyr)
library(readxl)
library(readr)
library(lubridate)
library(ggplot2)
# ------------------------------------------------------------------
# Define model -----------------------------------------------------
source("R/age_struct_seir_ode2.R")
# -------------------------------------------------------------------

# Specify initial conditions ----------------------------------------
# define population size (by age group)
n_vec <- c(0.10319920, 0.11620856, 0.12740219, 0.12198707, 0.13083463,
           0.14514332, 0.12092904, 0.08807406, 0.04622194) * 17407585 # Dutch population size
empty_state <- c(rep(0, 9)) # vector of zeros
seed_age_group <- sample(1:9,1)
inf_seed_vec <- empty_state
inf_seed_vec[seed_age_group] <- 1

s_vec     <- n_vec - inf_seed_vec; sv1_vec <- empty_state; sv2_vec <- empty_state; sv3_vec <- empty_state; sv4_vec <- empty_state; sv5_vec <- empty_state;
shold1_vec<- empty_state; shold2_vec <- empty_state; shold3_vec <- empty_state; shold4_vec <- empty_state; shold5_vec <- empty_state;
e_vec     <- empty_state; ev1_vec <- empty_state; ev2_vec <- empty_state; ev3_vec <- empty_state; ev4_vec <- empty_state; ev5_vec <- empty_state;
i_vec     <- inf_seed_vec; iv1_vec <- empty_state; iv2_vec <- empty_state; iv3_vec <- empty_state; iv4_vec <- empty_state; iv5_vec <- empty_state;
h_vec     <- empty_state; hv1_vec <- empty_state; hv2_vec <- empty_state; hv3_vec <- empty_state; hv4_vec <- empty_state; hv5_vec <- empty_state
ic_vec    <- empty_state; icv1_vec <- empty_state; icv2_vec <- empty_state; icv3_vec <- empty_state; icv4_vec <- empty_state; icv5_vec <- empty_state;
hic_vec   <- empty_state; hicv1_vec <- empty_state; hicv2_vec <- empty_state; hicv3_vec <- empty_state; hicv4_vec <- empty_state; hicv5_vec <- empty_state;
d_vec     <- empty_state
r_vec     <- empty_state; rv1_vec <- empty_state; rv2_vec <- empty_state; rv3_vec <- empty_state; rv4_vec <- empty_state; rv5_vec <- empty_state;
r1_vec    <- empty_state; r1v1_vec <- empty_state; r1v2_vec <- empty_state; r1v3_vec <- empty_state; r1v4_vec <- empty_state; r1v5_vec <- empty_state;
r2_vec    <- empty_state; r2v1_vec <- empty_state; r2v2_vec <- empty_state; r2v3_vec <- empty_state; r2v4_vec <- empty_state; r2v5_vec <- empty_state;
r3v1_vec  <- empty_state; r3v2_vec <- empty_state; r3v3_vec <- empty_state; r3v4_vec <- empty_state; r3v5_vec <- empty_state;
r3_vec    <- n_vec - s_vec - sv1_vec - sv2_vec - sv3_vec - sv4_vec - sv5_vec -
  shold1_vec - shold2_vec - shold3_vec - shold4_vec - shold5_vec -
  e_vec - ev1_vec - ev2_vec - ev3_vec - ev4_vec - ev5_vec -
  i_vec - iv1_vec - iv2_vec - iv3_vec - iv4_vec - iv5_vec -
  h_vec - hv1_vec - hv2_vec - hv3_vec - hv4_vec - hv5_vec -
  ic_vec - icv1_vec - icv2_vec - icv3_vec - icv4_vec - icv5_vec -
  hic_vec - hicv1_vec - hicv2_vec - hicv3_vec - hicv4_vec - hicv5_vec -
  d_vec - dv1_vec - dv2_vec - dv3_vec - dv4_vec - dv5_vec -
  r_vec - rv1_vec - rv2_vec - rv3_vec - rv4_vec - rv5_vec -
  r1_vec - r1v1_vec - r1v2_vec - r1v3_vec - r1v4_vec - r1v5_vec -
  r2_vec - r2v1_vec - r2v2_vec - r2v3_vec - r2v4_vec - r2v5_vec -
  r3v1_vec - r3v2_vec - r3v3_vec - r3v4_vec - r3v5_vec

init_t0 <- c(t        = 0,
             S        = s_vec, Sv_1d = sv1_vec, Sv_2d = sv2_vec, Sv_3d = sv3_vec, Sv_4d = sv4_vec, Sv_5d = sv5_vec,
             Shold_1d = shold1_vec, Shold_2d = shold2_vec, Shold_3d = shold3_vec, Shold_4d = shold4_vec, Shold_5d = shold5_vec,
             E        = e_vec, Ev_1d = ev1_vec, Ev_2d = ev2_vec, Ev_3d = ev3_vec, Ev_4d = ev4_vec, Ev_5d = ev5_vec,
             I        = i_vec, Iv_1d = iv1_vec, Iv_2d = iv2_vec, Iv_3d = iv3_vec, Iv_4d = iv4_vec, Iv_5d = iv5_vec,
             H        = h_vec, Hv_1d = hv1_vec, Hv_2d = hv2_vec, Hv_3d = hv3_vec, Hv_4d = hv4_vec, Hv_5d = hv5_vec,
             IC       = ic_vec, ICv_1d = icv1_vec, ICv_2d = icv2_vec, ICv_3d = icv3_vec, ICv_4d = icv4_vec, ICv_5d = icv5_vec,
             H_IC     = hic_vec, H_ICv_1d = hicv1_vec, H_ICv_2d = hicv2_vec, H_ICv_3d = hicv3_vec, H_ICv_4d = hicv4_vec, H_ICv_5d = hicv5_vec,
             D        = d_vec,
             R        = r_vec, Rv_1d = rv1_vec, Rv_2d = rv2_vec, Rv_3d = rv3_vec, Rv_4d = rv4_vec, Rv_5d = rv5_vec,
             R_1w     = r1_vec, Rv_1d_1w = r1v1_vec, Rv_2d_1w = r1v2_vec, Rv_3d_1w = r1v3_vec, Rv_4d_1w = r1v4_vec, Rv_5d_1w = r1v5_vec,
             R_2w     = r2_vec, Rv_1d_2w = r2v1_vec, Rv_2d_2w = r2v2_vec, Rv_3d_2w = r2v3_vec, Rv_4d_2w = r2v4_vec, Rv_5d_2w = r2v5_vec,
             R_3w     = r3_vec, Rv_1d_3w = r3v1_vec, Rv_2d_3w = r3v2_vec, Rv_3d_3w = r3v3_vec, Rv_4d_3w = r3v4_vec, Rv_5d_3w = r3v5_vec
)

# Specify model parameters ------------------------------------------
# define contact/transmission matrix --------------------------------
# path <- "inst/extdata/inputs/contact_matrices/converted/"
path <- "/rivm/s/ainsliek/data/contact_matrices/converted/" 
april_2017     <- readRDS(paste0(path,"transmission_matrix_april_2017.rds"))
april_2020     <- readRDS(paste0(path,"transmission_matrix_april_2020.rds"))
june_2020      <- readRDS(paste0(path,"transmission_matrix_june_2020.rds"))
september_2020 <- readRDS(paste0(path,"transmission_matrix_september_2020.rds"))
february_2021  <- readRDS(paste0(path,"transmission_matrix_february_2021.rds"))
june_2021      <- readRDS(paste0(path,"transmission_matrix_june_2021.rds"))
november_2021  <- readRDS(paste0(path,"transmission_matrix_november_2021.rds"))

# put contact matrices into a list for input into fit_to_data_func()
contact_matrices <- list(
  april_2017 = april_2017,
  april_2020 = april_2020,
  june_2020 = june_2020,
  september_2020 = september_2020,
  february_2021 = february_2021,
  june_2021 = june_2021,
  november_2021 = november_2021
)

# probabilities -------------------------------------------------------
dons_probs <- read_xlsx("inst/extdata/inputs/ProbabilitiesDelays_20210107.xlsx")
p_infection2admission <- dons_probs$P_infection2admission
p_admission2death <- dons_probs$P_admission2death
p_admission2IC <- dons_probs$P_admission2IC
p_IC2hospital <- dons_probs$P_IC2hospital
p_hospital2death <- c(rep(0, 5), 0.01, 0.04, 0.12, 0.29) # (after ICU)
p_reported_by_age <- c(0.29, 0.363, 0.381, 0.545, 0.645, 0.564, 0.365, 0.33, 0.409) # from Jantien

# delays --------------------------------------------------------------
time_symptom2admission <- c(2.29, 5.51, 5.05, 5.66, 6.55, 5.88, 5.69, 5.09, 4.33) # assume same as infectious2admission
time_admission2discharge <- 7.9
time_admission2IC <- 2.28
time_IC2hospital <- 15.6
time_hospital2discharge <- 10.1 # (after ICU)
time_admission2death <- 7
time_IC2death <- 19
time_hospital2death <- 10 # (after ICU)

# define transition rates ---------------------------------------------
i2r    <- (1-p_infection2admission) / 2                    # I -> R
i2h    <- p_infection2admission / time_symptom2admission   # I -> H

h2ic   <- p_admission2IC / time_admission2IC               # H -> IC
h2d    <- p_admission2death / time_admission2death         # H -> D
h2r    <- (1 - (p_admission2IC + p_admission2death)) / time_admission2discharge
                                                           # H -> R

ic2hic <- p_IC2hospital / time_IC2hospital                 # IC -> H_IC
ic2d   <- (1 - p_IC2hospital) / time_IC2death              # IC -> D

hic2d  <- p_hospital2death / time_hospital2death           # H_IC -> D
hic2r  <- (1 - p_hospital2death) / time_hospital2discharge # H_IC -> R

dt <- 1

# ve estimates ------------------------------------------------------
# list containing the following named lists:
# delays, ve_inf, ve_hosp, ve_trans
# each named list has the following named elements:
# pfizer, moderna, astrazeneca, jansen
ve_params <- readRDS("inst/extdata/inputs/ve_params.rds")
ve_params$delay_rate <- lapply(ve_params$delays, function(x) 1/x)
# vaccination schedule ----------------------------------------------
# read in vaccination schedule
vac_schedule <- read_csv("inst/extdata/inputs/vac_schedule_real_w_4th_and_5th_dose.csv") %>%
  select(-X1)

# source function to convert vac schedule for model input -----------
source("R/convert_vac_schedule2.R")
source("R/calc_ve_w_waning.R")
# convert vaccination schedule for input into model
vac_rates_wt <- convert_vac_schedule2(
  vac_schedule = vac_schedule,
  delay = ve_params$delay_rate,
  ve = ve_params$ve_inf$wildtype,
  hosp_multiplier = ve_params$ve_hosp$wildtype,
  ve_trans = ve_params$ve_trans$wildtype,
  wane = TRUE)

vac_rates_alpha <- convert_vac_schedule2(
  vac_schedule = vac_schedule,
  delay = ve_params$delays,
  ve = ve_params$ve_inf$alpha,
  hosp_multiplier = ve_params$ve_hosp$alpha,
  ve_trans = ve_params$ve_trans$alpha,
  wane = TRUE)

vac_rates_delta <- convert_vac_schedule2(
  vac_schedule = vac_schedule,
  delay = ve_params$delays,
  ve = ve_params$ve_inf$delta,
  hosp_multiplier = ve_params$ve_hosp$delta,
  ve_trans = ve_params$ve_trans$delta,
  wane = TRUE)

vac_rates_omicron <- convert_vac_schedule2(
  vac_schedule = vac_schedule,
  delay = ve_params$delays,
  ve = ve_params$ve_inf$omicron,
  hosp_multiplier = ve_params$ve_hosp$omicron,
  ve_trans = ve_params$ve_trans$omicron,
  wane = TRUE)

# make into a list for input into fit_to_data_func()
vac_rates_list <- list(
  wildtype = vac_rates_wt,
  alpha = vac_rates_alpha,
  delta = vac_rates_delta,
  omicron = vac_rates_omicron
)

lapply(vac_rates_wt, )
# model input parameters ---------------------------------------------
# parameters must be in a named list
params <- list(dt = dt,
               N = n_vec,
               beta = 0.0004734092/dt,
               beta1 = 0.14/dt,
               sigma = 0.5/dt,
               gamma = i2r/dt,
               h = i2h/dt,
               i1 = h2ic/dt,
               d = h2d/dt,
               r = h2r/dt,
               i2 = ic2hic/dt,
               d_ic = ic2d/dt,
               d_hic = hic2d/dt,
               r_ic = hic2r/dt,
               epsilon = 0.00/dt,
               omega = 0.0038/dt,
               p_report = p_reported_by_age,
               contact_mat = april_2017,
               calendar_start_date = as.Date("2020-01-01"),
               no_vac = FALSE
)

# -------------------------------------------------------------------
# Define likelihood function ----------------------------------------
# we're fitting the transmission probability (beta) and an
# over-dispersion parameter (alpha)
likelihood_func_test <- function(x, t, data, params, init) {
  # parameters to be estimated
  beta <- x[1]
  alpha <- x[2]
  print(x)
  # observed daily cases
  inc_obs <- data$inc
  
  # run model with current parameter values
  params$beta <- x[1]/10000
  rk45 <- rkMethod("rk45dp7")
  seir_out <- ode(init, t, age_struct_seir_ode2, params, method = rk45) # , rtol = 1e-08, hmax = 0.02
  out <- as.data.frame(seir_out)
  
  # modeled cases
  daily_cases <- rowSums(params$sigma * (out[c(paste0("E",1:9))] + out[c(paste0("Ev_1d",1:9))] +
                                         out[c(paste0("Ev_2d",1:9))] + out[c(paste0("Ev_3d",1:9))] +
                                         out[c(paste0("Ev_4d",1:9))] + out[c(paste0("Ev_5d",1:9))]) * params$p_report)
  daily_cases <- ifelse(daily_cases == 0, 0.0001, daily_cases) # prevent likelihood function function from being Inf
  
  # log-likelihood function
  # lik <- sum(dpois(x = inc_obs,lambda = incidence,log=TRUE))
  lik <- -sum(stats::dnbinom(x = inc_obs, mu = daily_cases, size = alpha, log = TRUE))
  
  # print(lik)
  lik
}
# -------------------------------------------------------------------
# -------------------------------------------------------------------

# define fit windows ------------------------------------------------
df_breakpoints <- read_csv2("inst/extdata/inputs/breakpoints_for_model_fit_v3.csv") %>%
  mutate(date = as.Date(date, format = "%d-%m-%Y"),
         time = as.numeric(date - date[1])) %>%
  select(date, time, variant, contact_matrix)

bp_for_fit <- df_breakpoints#[1:5,]
n_bp <- dim(bp_for_fit)[1] - 1

# specify initial values and bounds for fitted parameters -----------
fit_params <- list(
  init_value = c(4, 1),
  lower_bound = c(0.5, 0.0001),
  upper_bound = c(Inf, Inf)
)

# create empty containers to store outputs from fit -----------------
init_cond <- list()     # store initial conditions for each window
out <- list()           # model output for each time point
times <- list()         # time points
cases <- list()         # daily cases to plot against real data
mles <- list()          # MLEs for each time window
beta_draws <- list()    # store 200 parameter draws
ci_out <- list()        # store model outputs for confidence bounds
ci_cases <- list()

susceptibles <- list()
exposed <- list()
infected <- list()
hospitalised <- list()
ic <- list()
hosp_after_ic <- list()
recovered <- list()
deaths <- list()
recovered1 <- list()
recovered2 <- list()
recovered3 <- list()
# load case data ----------------------------------------------------
data_date <- "2022-05-21"
case_data <- readRDS(paste0("inst/extdata/data/case_data_upto_", data_date, ".rds"))

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Fit model to data -------------------------------------------------
init_cond[[1]] <- init_t0

# loop over time windows --------------------------------------------
for (j in 1:n_bp) {
  
  print(paste(paste0(j,")"),"Fitting from", bp_for_fit$date[j], "to", bp_for_fit$date[j+1]))
  
  # set time sequence -----------------------------------------------
  times[[j]] <- seq(bp_for_fit$time[j], bp_for_fit$time[j+1], by = 1)
  
  # set contact matrix for time window ------------------------------
  if (bp_for_fit$contact_matrix[j+1] == "april_2017"){contact_matrix <- contact_matrices$april_2017 #; print("april_2017")
  } else if (bp_for_fit$contact_matrix[j+1] == "april_2020"){contact_matrix <- contact_matrices$april_2020 #; print("april_2020")
  } else if (bp_for_fit$contact_matrix[j+1] == "june_2020"){contact_matrix <- contact_matrices$june_2020 #; print("june_2020")
  } else if (bp_for_fit$contact_matrix[j+1] == "september_2020"){contact_matrix <- contact_matrices$september_2020 #; print("septemeber_2020")
  } else if (bp_for_fit$contact_matrix[j+1] == "february_2021"){contact_matrix <- contact_matrices$february_2021 #; print("february_2021")
  } else if (bp_for_fit$contact_matrix[j+1] == "june_2021"){contact_matrix <- contact_matrices$june_2021 #; print("june_2021")
  } else {contact_matrix <- contact_matrices$november_2021} 
  
  # change contact matrix in params list ----------------------------
  params$contact_mat <- contact_matrix$mean
  
  # set no_vac = TRUE before vaccination program starts -------------
  if(bp_for_fit$date[j+1] <= as.Date("2021-01-04")){ params$no_vac <- TRUE
  } else {params$no_vac <- FALSE}
  
  # set vaccination characteristics depending on variant ------------
  if(!params$no_vac){
    # set VE for time window depending on which variant was dominant
    if (bp_for_fit$variant[j+1] == "wildtype"){params$vac_inputs <- vac_rates_list$wildtype
    } else if (bp_for_fit$variant[j+1] == "alpha"){params$vac_inputs <- vac_rates_list$alpha
    } else if (bp_for_fit$variant[j+1] == "delta"){params$vac_inputs <- vac_rates_list$delta
    } else {params$vac_inputs <- vac_rates_list$omicron}
  }
  
  # subset data for time window -------------------------------------
  case_data_sub <- case_data[times[[j]] + 1, ]
  
  # run optimization procedure --------------------------------------
  res <- optim(par = fit_params$init_value, 
               fn = likelihood_func_test,
               method = "L-BFGS-B",
               lower = fit_params$lower_bound,
               upper = fit_params$upper_bound,
               t = times[[j]],
               data = case_data_sub,
               params = params,
               init = init_cond[[j]],
               hessian = TRUE
  )
  
  # store MLE --------------------------------------------------------
  mles[[j]] <- c(beta = res$par[1]/10000, alpha = res$par[2])
  print(mles[[j]])
  
  # Run model --------------------------------------------------------
  params$beta <- res$par[1]/10000
  rk45 <- rkMethod("rk45dp7")
  seir_out <- ode(init_cond[[j]], times[[j]], age_struct_seir_ode2,  
                  params, method = rk45) # , rtol = 1e-08, hmax = 0.02
  
  # store outputs ----------------------------------------------------
  out[[j]] <- as.data.frame(seir_out) 
  cases[[j]] <-  rowSums(params$sigma * (out[[j]][c(paste0("E",1:9))] + out[[j]][c(paste0("Ev_1d",1:9))] +
                                         out[[j]][c(paste0("Ev_2d",1:9))] + out[[j]][c(paste0("Ev_3d",1:9))] +
                                         out[[j]][c(paste0("Ev_4d",1:9))] + out[[j]][c(paste0("Ev_5d",1:9))]) * params$p_report)
  
  # plot for quick check of fit --------------------------------------
  plot(case_data_sub$inc ~ times[[j]], pch = 16, col = "red", 
       ylim = c(0, max(case_data_sub$inc,cases[[j]])))
  lines(cases[[j]] ~ times[[j]]) 
  
  #------------------------------------------------------------------
  # get confidence bounds -------------------------------------------
  # draw 200 parameter values
  parameter_draws <- mvtnorm::rmvnorm(200, res$par, solve(res$hessian))
  beta_draws[[j]] <- data.frame(beta = (parameter_draws[,1]/10000)) %>%
    mutate(index = 1:200)
  
  # run model for each beta draw (with different contact matrix) ----
  print("Getting confidence bounds...")
  #ci_out[[j]] <- list()
  ci_cases[[j]] <- list()
  for(i in 1:200){
    params$beta <- beta_draws[[j]][i,1]
    params$contact_mat <- contact_matrix[[i]]
    seir_out_ci <- ode(init_cond[[j]], times[[j]], age_struct_seir_ode2,  
                       params, method = rk45) #, rtol = 1e-08, hmax = 0.02
    seir_out_ci1 <- as.data.frame(seir_out_ci) 
    ci_cases[[j]][[i]] <-  rowSums(params$sigma * seir_out_ci1[c(paste0("E",1:9))] * params$p_report)
  }
  ci_out[[j]] <- do.call("rbind", ci_cases[[j]])
  # -----------------------------------------------------------------
  
  # update initial conditions for next time window
  init_cond[[j+1]] <- unlist(tail(out[[j]],1)[names(out[[j]]) != "time"])
  
  # output error message if negative compartment values
  if(any(init_cond[[j+1]] < 0)){
    stop("Negative compartment values")
  }
  
  # output error message if sum of initial conditions != N
  test <- all.equal(sum(init_cond[[j+1]][-1]), sum(params$N))
  if(!test){
    stop("Sum of compartments not equal to total population size")
  }
  
  # ------------------------------------------------------------------  
  # get number of people in each compartment (to check for unusual 
  # behaviour)
  susceptibles[[j]]  <- rowSums(out[[j]][,c(paste0("S",1:9))])
  exposed[[j]]       <- rowSums(out[[j]][,c(paste0("E",1:9))])
  infected[[j]]      <- rowSums(out[[j]][,c(paste0("I",1:9))])
  hospitalised[[j]]  <- rowSums(out[[j]][,c(paste0("H",1:9))])
  ic[[j]]            <- rowSums(out[[j]][,c(paste0("IC",1:9))])
  hosp_after_ic[[j]] <- rowSums(out[[j]][,c(paste0("H_IC",1:9))])
  deaths[[j]]        <- rowSums(out[[j]][,c(paste0("D",1:9))])
  recovered[[j]]     <- rowSums(out[[j]][,c(paste0("R",1:9))]) 
  recovered1[[j]]    <- rowSums(out[[j]][,c(paste0("R_1w",1:9))]) 
  recovered2[[j]]    <- rowSums(out[[j]][,c(paste0("R_2w",1:9))]) 
  recovered3[[j]]    <- rowSums(out[[j]][,c(paste0("R_3w",1:9))]) 
  
} # end of for loop over breakpoints

# save outputs ------------------------------------------------------
path_out <- "inst/extdata/results/model_fits/"
# get confidence bounds
ci_out_wide <- do.call("cbind", ci_out)
bounds <- apply(ci_out_wide, 2, quantile, probs = c(0.025, 0.975)) # get quantiles

x_axis <- unlist(times)
df_model_fit <- data.frame(time = x_axis, 
                           date = params$calendar_start_date + x_axis,
                           obs = case_data$inc[x_axis + 1], 
                           mle = unlist(cases), 
                           lower = bounds[1,], 
                           upper = bounds[2,])

saveRDS(df_model_fit,
        file = paste0(path_out, "model_fit_df_from_", df_model_fit$date[1],"_to_",
                      tail(df_model_fit$date,1), ".rds"))

rtn <- list(out_mle = out,
            cases_mle = cases,
            cases_ci = ci_out)
saveRDS(rtn, file = paste0(path_out,"model_fit_outputs_from_",df_model_fit$date[1],"_to_",
                           tail(df_model_fit$date,1), ".rds"))

# -------------------------------------------------------------------
# -------------------------------------------------------------------
# Plot output -------------------------------------------------------
# plot SEIR compartments
plot(unlist(susceptibles) ~ x_axis, type = "l", ylim = c(0, sum(params$N)))
abline(h = sum(params$N), lty = "dashed")
lines(unlist(recovered) ~ x_axis, type = "l", col = "blue") #, ylim = c(0,max(recovered))
lines(unlist(recovered1) ~ x_axis, col = "blue", lty = "dashed")
lines(unlist(recovered2) ~ x_axis, col = "blue", lty = "dotted")
lines(unlist(recovered3) ~ x_axis, col = "blue", lty = "twodash")
lines(unlist(exposed) ~ x_axis, col = "green")
lines(unlist(infected) ~ x_axis, col = "red")

# plot severe disease compartments
plot(unlist(hospitalised) ~ x_axis, type = "l", col = "orange", 
     ylim = c(min(unlist(ic)),max(unlist(hospitalised),unlist(deaths))))
lines(unlist(ic) ~ x_axis, col = "pink", type = "l")
lines(unlist(hosp_after_ic) ~ x_axis, col = "purple")
lines(unlist(deaths) ~ x_axis, col = "grey")


p <- ggplot(data = df_model_fit, aes(x = date, y = mle, linetype="solid")) +
  geom_point(data = df_model_fit, aes(x = date, y = obs, color = "Osiris notifications")) +
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
# --------------------------------------------------------------------