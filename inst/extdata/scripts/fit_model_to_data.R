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
nice <- readRDS("inst/extdata/data/nice_admissions_20210212.rds") %>%
  filter(AdmissionDate >= "2020-01-01")
# make dates consecutive
nice1 <- complete(nice, AdmissionDate = full_seq(AdmissionDate, 1)) %>%
  mutate(n = ifelse(is.na(n), 0, n))

# probabilities ---------------------------------------------------------
p_infection2admission <- c(0.003470, 0.000377, 0.000949, 0.003880, 0.008420, 
                           0.016500,0.025100, 0.049400, 0.046300)
p_admission2death <- c(0.00191, 0.00433, 0.00976, 0.02190, 0.02500, 0.04010,
                       0.10600, 0.22900, 0.31100)

# contact matrix --------------------------------------------------------
contact_matrices_all <- readRDS("inst/extdata/data/contact_matrices_for_model_input.rds")
c1 <- as.matrix(contact_matrices_all$baseline[,-1])
c2 <- as.matrix(contact_matrices_all$april2020[,-1])
c3 <- as.matrix(contact_matrices_all$june2020[,-1])
c4 <- as.matrix(contact_matrices_all$september2020[,-1])
# age distribution and pop size -----------------------------------------
age_dist <- c(0.10319920, 0.11620856, 0.12740219, 0.12198707, 0.13083463, 
              0.14514332, 0.12092904, 0.08807406, 0.04622194)
n <- 17407585                                 # Dutch population size
n_vec <- n * age_dist

empty_state <- c(rep(0,9))

# write wrapper function for seir model code that outputs only hospital
# admission counts ------------------------------------------------------
fit_to_data_wrapper <- function(x, params){ #, t0, delta, s, g
 
  r0 <- params[1]
  init_i <- params[2]
  #s <- params[3]
  #g <- params[4]
  s <- 0.5
  g <- 0.5
  delta1 <- params[3]
  delta2 <- params[4]
  delta3 <- params[5]
  delta4 <- params[6]
  
  # list of dates at different intervention points
  t_vec <- lubridate::yday(as.Date(c("2020-01-01", 
                                     "2020-03-15", 
                                     #"2020-05-11",
                                     "2020-06-01", 
                                     #"2020-07-01", 
                                     #"2020-08-16",
                                     "2020-09-29", 
                                     #"2020-10-14",
                                     "2020-12-15",
                                     "2021-02-11"))) 
  t_vec[length(t_vec)] <- t_vec[length(t_vec)] + 366 #2020 was a leap year
  rtn <- list()
  
  tmp <- get_beta(R0 = r0, contact_matrix = c1, N = n_vec, sigma = s, 
                  gamma = s) 
  beta <- tmp$beta
  
  # loop over time periods -------------------------------------------
  for (i in 1:(length(t_vec)-1)){
    delta <- (i == 1) * 1 + (i == 2) * delta1 + (i == 3) * delta2 + (i == 4) * delta3 + (i == 5) * delta4
    params <- list(beta = beta * delta,           # transmission rate
                 gamma = g,                      # 1/gamma = infectious period
                 sigma = s,                      # 1/sigma = latent period
                 N = n_vec,                      # Population (no need to change)
                 h = p_infection2admission,      # Rate from infection to hospital admission
                 d = p_admission2death,          # Rate from admission to death
                 r = 0.0206,                     # Rate from admission to recovery
                 c_start = c1,
                 no_vac = TRUE
    )
    
  # Specify initial values -------------------------------------------
    if (i == length(t_vec) - 1){
      out_length <- length(seq(t_vec[i], t_vec[i+1], by = 1))
      times <- seq(t_vec[i], t_vec[i+1], by = 1) 
    } else {
      out_length <- length(seq(t_vec[i], t_vec[i+1] - 1, by = 1))
      times <- seq(t_vec[i], t_vec[i+1] - 1, by = 1) 
    }

  timeInt <- times[2]-times[1]             
  init <- c(t = times[1],                  
            S = params$N - c(rep(init_i/9, 9)),
            Shold_1d = empty_state,
            Sv_1d = empty_state,
            Shold_2d = empty_state,
            Sv_2d = empty_state,
            E = empty_state,
            Ev_1d = empty_state,
            Ev_2d = empty_state,
            I = c(rep(init_i/9, 9)),
            Iv_1d = empty_state,
            Iv_2d = empty_state,
            H = empty_state,
            Hv_1d = empty_state,
            Hv_2d = empty_state,
            D = empty_state,
            R = empty_state,
            Rv_1d = empty_state,
            Rv_2d = empty_state
  )                      

  # Solve model ------------------------------------------------------
  seir_out <- lsoda(init,times,age_struct_seir_ode,params)

  seir_out <- as.data.frame(seir_out)
  out <- postprocess_age_struct_model_output(seir_out)
  out_vec <- rowSums(out$H + out$Hv_1d + out$Hv_2d)
  if(length(out_vec) != out_length){
    out_vec <- c(rep(0,out_length-length(times)), out_vec)
  }
  rtn[[i]] <- out_vec
  }
  
  rtn2 <- unlist(rtn)
  
  return(rtn2)
}

# use nls to fit non-linear least squares ----------------------------
delta_start <- c(mean(c2/c1), mean(c3/c1), mean(c4/c1), 1)
delta_lower <- c(min(c2/c1), min(c3/c1), min(c4/c1), 0.05)
delta_upper <- c(max(c2/c1), max(c3/c1), max(c4/c1), 1)

model_fit <- nls(n ~ fit_to_data_wrapper(AdmissionDate, params), #, t0, delta, s, g
                 data = nice1, # where x and y are
                 start = list(params = c(3, 0.03, delta_start)), #initial conditions s = 0.5, g = 0.5
                 lower = c(0, 0, delta_lower), # lower bound of parameter values s = 0.1, g = 0.1
                 upper = c(10, 10e6, delta_upper), # upper bound of parameter values s = 1, g = 1
                 trace = T,
                 algorithm = "port", # must be used to specify bounds
                 control = list(maxiter = 50000, 
                                warnOnly=T)
                 )
summary(model_fit)
# plot ---------------------------------------------------------------
plot(n~AdmissionDate,data=nice1,type = "l")
abline(v = t_vec)
lines(x = nice1$AdmissionDate, y = model_fit$m$fitted(), 
      type = "l", col = "red")

