# Fit model to data to get appropriate initial conditions ---------------

# Load packages ----------------------------------------------------
library(deSolve)
library(reshape2)
library(ggplot2)
library(dplyr)
library(stringr)
library(tidyr)
library(readxl)
library(rARPACK)
library(readr)
library(lubridate)
library(FSA)

# Source functions ------------------------------------------------------
source("R/age_struct_seir_ode.R")
source("R/postprocess_age_struct_model_output.R")
source("R/get_foi.R")
source("R/get_beta.R")
source("R/get_vac_rate.R")
source("R/get_vac_rate_2.R")

# read in hospital admission counts -------------------------------------
nice <- readRDS("inst/extdata/data/nice_admissions_20210212.rds") #%>%
  #filter(AdmissionDate >= "2020-01-01")
# make dates consecutive
nice1 <- complete(nice, AdmissionDate = full_seq(AdmissionDate, 1)) %>%
  mutate(n = ifelse(is.na(n), 0, n))

# probabilities
dons_probs <- read_xlsx("inst/extdata/data/ProbabilitiesDelays_20210107.xlsx")
p_infection2admission <- dons_probs$P_infection2admission
p_admission2death <- dons_probs$P_admission2death
p_admission2IC <- dons_probs$P_admission2IC
p_IC2hospital <- dons_probs$P_IC2hospital
p_IC2death <- 1-p_IC2hospital
p_hospital2death <- c(rep(0,5), 0.01, 0.04, 0.12, 0.29) #(after ICU)
p_reported_by_age <- c(0.29, 0.363, 0.381, 0.545, 0.645, 0.564, 0.365, 0.33, 0.409) # from Jantien
p_inf_by_age <- c(0.018, 0.115, 0.156, 0.118, 0.142, 0.199, 0.114, 0.062, 0.054 + 0.023)          

# delays
time_symptom2admission <- c(2.29, 5.51, 5.05, 5.66, 6.55, 5.88, 5.69, 5.09, 4.33) # assume same as infectious2admission
time_admission2discharge <- 7.9
time_admission2IC <- 2.28
time_IC2hospital <- 15.6
time_hospital2discharge <- 10.1 #(after ICU)
time_admission2death <- 7 
time_IC2death <- 19 
time_hospital2death <- 10 #(after ICU)

# age distribution and pop size 
age_dist <- c(0.10319920, 0.11620856, 0.12740219, 0.12198707, 0.13083463, 
              0.14514332, 0.12092904, 0.08807406, 0.04622194)
n <- 17407585                                 # Dutch population size
n_vec <- n * age_dist

# contact matrices
contact_matrices_all <- readRDS("inst/extdata/data/contact_matrices_for_model_input.rds")
c1 <- as.matrix(contact_matrices_all$baseline[,-1])
c2 <- as.matrix(contact_matrices_all$april2020[,-1])
c3 <- as.matrix(contact_matrices_all$june2020[,-1])
c4 <- as.matrix(contact_matrices_all$september2020[,-1])

# relative susceptibility/infectiousness
rel_trans <- c(1.000, 3.051, 5.751, 3.538, 3.705, 4.365, 5.688, 5.324, 7.211)
get_transmission_matrix <- function(x, contact_mat){
  tmp <- sweep(contact_mat, 1, x, "*")
  rtn <- sweep(tmp, 2, x, "*")
  return(rtn)
}

t1 <- get_transmission_matrix(rel_trans, c1)
t2 <- get_transmission_matrix(rel_trans, c2)
t3 <- get_transmission_matrix(rel_trans, c3)
t4 <- get_transmission_matrix(rel_trans, c4)

# parameter inputs
s <- 0.5
g <- 0.5
r0 <- 3.33966     
init_i <- 0.544645
tmp <- get_beta(R0 = r0, contact_matrix = t1, N = n_vec, sigma = s, 
                gamma = s) 
beta <- tmp$beta
h <- p_infection2admission/time_symptom2admission
i1 <- p_admission2IC/time_admission2IC
i2 <- p_IC2hospital/time_IC2hospital
d <- p_admission2death/time_admission2death
d_ic <- p_IC2death/time_IC2death
d_hic <- p_hospital2death/time_hospital2death
r <- (1 - p_admission2death)/time_admission2discharge
r_ic <- (1 - p_IC2death)/time_hospital2discharge

# data
nice_sub <- nice1 %>% 
  filter(AdmissionDate > as.Date("2019-12-30") & AdmissionDate <= as.Date("2020-03-25")) %>%
  mutate(moving_avg3 = zoo::rollmean(n, k = 3, fill = NA)) %>%
  filter(!is.na(moving_avg3))

# write wrapper function for seir model code that outputs only hospital
# admission counts ------------------------------------------------------
fit_to_data_wrapper <- function(x, params){
 
  r0 <- params[1]
  init_i <- params[2]
  #delta1 <- params[1]
  #t01 <- params[2]
  #t02 <- params[3]
  #delta2 <- params[2]
  # delta3 <- params[3]
  # delta4 <- params[4]
  
  # list of dates at different intervention points
  t_vec <- lubridate::yday(as.Date(c("2020-01-01",
                                     "2020-03-24"#,
                                     #"2020-06-01"#,
                                     #"2020-09-01"#,
                                     #"2021-01-01",
                                     #"2021-02-01"
                                     )))
  #t_vec[2] <- t_vec[2] + t01
  #t_vec[length(t_vec)] <- t_vec[length(t_vec)] + 366 #2020 was a leap year
  
  s <- 0.5
  g <- 0.5
  # r0 <- 2.95550 
  # init_i <- 0.489163
  # # r0 <- 3.33966     
  # init_i <- 0.544645
  # r0 <- 3.48164 
  # init_i <- 0.299867
  
  tmp <- get_beta(R0 = r0, contact_matrix = t1, N = n_vec, sigma = s, 
                  gamma = s) 
  beta <- tmp$beta
  slope <- 0
  # loop over time periods -------------------------------------------
  rtn <- list() # store output for each iteration here
  for (i in 1:(length(t_vec)-1)){
    delta <- (i == 1) * 1 #+ (i == 2) * delta1 #+ (i == 3) * delta2 #+ (i == 4) * delta3 + (i == 5) * delta4
    #cm <- (i == 1) * c1 + (i == 2) * c2 + (i == 3) * c3 + (i == 4) * c4

    params <- list(beta = beta,                    # transmission rate
                   gamma = g,                      # 1/gamma = infectious period
                   sigma = s,                      # 1/sigma = latent period
                   N = n_vec,                      # Population (no need to change)
                   h = h,                          # Rate from infection to hospital admission/ time from infection to hosp admission
                   i1 = i1,
                   i2 = i2,
                   d = d, 
                   d_ic = d_ic,
                   d_hic = d_hic,
                   r = r,
                   r_ic = r_ic,
                   p_report = p_reported_by_age,
                   c_lockdown = t1,
                   c_relaxed = t4,
                   c_very_relaxed = t3,
                   c_normal = t1,
                   force_relax = NULL,
                   #vac_schedule = old_to_young,
                   #ve = ve,
                   #delay = delays,
                   use_cases = FALSE,              # use cases as criteria to change contact matrices. If FALSE, IC admissions used.
                   thresh_l = 0,#14.3/100000 * sum(n_vec),        # 10 for IC admissions
                   thresh_m = sum(n_vec),
                   thresh_u = sum(n_vec),#35.7/100000 * sum(n_vec),        # 20 for IC admissions
                   no_vac = TRUE
    )
    
  # Specify initial values -------------------------------------------
    #t0 <- (i == 1) * t01 + (i == 2) * t02 #+ (i == 3) * t03 #+ (i == 4) * t04 + (i == 5) * t05

    if (i == length(t_vec) - 1){
      #out_length <- length(seq(floor(t_vec[i]), floor(t_vec[i+1]), by = 1))
      times <- seq(floor(t_vec[i]), floor(t_vec[i+1]), by = 1)
    } else {
      #out_length <- length(seq(floor(t_vec[i]), floor(t_vec[i+1]) - 1, by = 1))
      times <- seq(floor(t_vec[i]), floor(t_vec[i+1]) - 1, by = 1)
    }
  timeInt <- times[2]-times[1] 
  empty_state <- c(rep(0,9))
  
  if (i == 1){
    i_start <- init_i * p_inf_by_age
    init <- c(t = times[1],                  
              S = params$N - (init_i * p_inf_by_age),
              Shold_1d = empty_state,
              Sv_1d = empty_state,
              Shold_2d = empty_state,
              Sv_2d = empty_state,
              E = empty_state,
              Ev_1d = empty_state,
              Ev_2d = empty_state,
              I = (init_i * p_inf_by_age),
              Iv_1d = empty_state,
              Iv_2d = empty_state,
              H = empty_state,
              Hv_1d = empty_state,
              Hv_2d = empty_state,
              H_IC = empty_state,
              H_ICv_1d = empty_state,
              H_ICv_2d = empty_state,
              IC = empty_state,
              ICv_1d = empty_state,
              ICv_2d = empty_state,
              D = empty_state,
              R = empty_state,
              Rv_1d = empty_state,
              Rv_2d = empty_state
    )                      
    
  } else {init <- init_new}
  # Solve model ------------------------------------------------------
  seir_out <- lsoda(init,times,age_struct_seir_ode,params)
  seir_out <- as.data.frame(seir_out)
  out <- postprocess_age_struct_model_output(seir_out)
  init_new <- c(t = dim(out$S)[1],
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
  #print(init_i)
  hosp_by_age_group <- sweep((out$I + out$Iv_1d + out$Iv_2d), 2, h, "*")
  out_vec <- rowSums(hosp_by_age_group)
  
  # if(length(out_vec) != out_length){
  #   out_vec <- c(rep(0,out_length-length(times)), out_vec)
  # }
  rtn[[i]] <- out_vec
  }
  
  rtn2 <- unlist(rtn)
  
  return(rtn2)
}

# use nls to fit non-linear least squares ----------------------------
delta_start <- c(rep(0.69, 4))
delta_lower <- c(rep(0.001, 4))
delta_upper <- c(rep(1, 4))
t0_start <- c(1.1, 88, 153)
# t0_lower <- c(0, 70)
# t0_upper <- c(30, 95)

# fit model
model_fit <- nls(moving_avg3 ~ fit_to_data_wrapper(AdmissionDate, params), 
                 data = nice_sub, # where x and y are
                 start = list(params = c(2.5, 0.03)), #initial conditions
                 lower = c(2, 0), # lower bound of parameter values
                 upper = c(4, Inf), # upper bound of parameter values
                 trace = T,
                 algorithm = "port", # must be used to specify bounds
                 control = list(maxiter = 100000, 
                                warnOnly=T)
                 )
# plot ---------------------------------------------------------------
# loess_fit <- loess(n ~ time, data=nice1, span=0.5)
# smoothed <- predict(loess_fit) 
ymax <- max(max(model_fit$m$fitted()),nice_sub$moving_avg3)
plot(moving_avg3~AdmissionDate,data=nice_sub,type = "l", ylim = c(0,ymax),
     xlab = "Date", ylab = "Hospital Admissions") 
lines(x = nice_sub$AdmissionDate, y = model_fit$m$fitted(), 
      type = "l", col = "red")
# lines(smoothed, x = nice1$AdmissionDate, col = "blue")

