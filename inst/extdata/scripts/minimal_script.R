# Minimal script for age-structured SEIR model
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
#source("R/get_vac_rate_2.R")

# Model inputs ----------------------------------------------------------
# probabilities 
p_infection2admission <- c(0.003470, 0.000377, 0.000949, 0.003880, 0.008420, 
                           0.016500,0.025100, 0.049400, 0.046300)
p_admission2death <- c(0.00191, 0.00433, 0.00976, 0.02190, 0.02500, 0.04010,
                       0.10600, 0.22900, 0.31100)

# contact matrix 
contact_matrices_all <- readRDS("inst/extdata/data/contact_matrices_for_model_input.rds")
c1 <- as.matrix(contact_matrices_all$baseline[,-1])
c2 <- as.matrix(contact_matrices_all$april2020[,-1])
c3 <- as.matrix(contact_matrices_all$june2020[,-1])
c4 <- as.matrix(contact_matrices_all$september2020[,-1])
# age distribution and pop size 
age_dist <- c(0.10319920, 0.11620856, 0.12740219, 0.12198707, 0.13083463, 
              0.14514332, 0.12092904, 0.08807406, 0.04622194)
n <- 17407585                                 # Dutch population size
n_vec <- n * age_dist

# time vector
t_vec <- lubridate::yday(as.Date(c("2020-01-01",
                                   "2020-03-15",
                                   "2020-06-01"#,
                                   # "2020-09-01",
                                   # "2021-12-15",
                                   # "2021-02-11"
                         )))
t_vec[length(t_vec)] <- t_vec[length(t_vec)] + 366 #2020 was a leap year
t_index <- length(t_vec)
# parameter values
s <- 0.2
g <- 0.125
r0 <- 3.33966     
init_i <- 0.544645

tmp <- get_beta(R0 = r0, contact_matrix = c1, N = n_vec, sigma = s, 
                gamma = s) 

params <- list(beta = tmp$beta,                # transmission rate
               gamma = g,                      # 1/gamma = infectious period
               sigma = s,                      # 1/sigma = latent period
               N = n_vec,                      # Population (no need to change)
               h = p_infection2admission,      # Rate from infection to hospital admission
               time_inf2hosp = 21,
               d = p_admission2death,          # Rate from admission to death
               r = 1-p_admission2death,                     # Rate from admission to recovery
               time_hosp2rec = 10.8,
               c1 = c1,
               c2 = c2,
               c3 = c3,
               c4 = c4,
               t_vec = t_vec,
               no_vac = TRUE
)

# initial conditions
empty_state <- c(rep(0,9))
times <- seq(t_vec[1],t_vec[t_index], by = 1)
timeInt <- times[2]-times[1]             
init <- c(t = times[1],                  
          S = params$N - (init_i * prob_inf_by_age),
          Shold_1d = empty_state,
          Sv_1d = empty_state,
          Shold_2d = empty_state,
          Sv_2d = empty_state,
          E = empty_state,
          Ev_1d = empty_state,
          Ev_2d = empty_state,
          I = (init_i * prob_inf_by_age),
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

# Solve model --------------------------------------------------------
seir_out <- lsoda(init,times,age_struct_seir_ode,params)
seir_out <- as.data.frame(seir_out)
out <- postprocess_age_struct_model_output(seir_out)

# estimates
# lambda <- get_foi(out, beta = params$beta, c_main = c1, N = n_vec)
# inc <- rowSums(out$S * lambda)
hosp_by_age_group <- sweep((out$I + out$Iv_1d + out$Iv_2d), 2, params$h/params$time_inf2hosp, "*")
hosp <- rowSums(hosp_by_age_group)

# quick plots
plot(times, rowSums(out$S), type = "l")
plot(times, rowSums(out$E), type = "l")
plot(times, rowSums(out$I), type = "l")
plot(times, rowSums(out$R), type = "l")

plot(times, hosp, col = "red", type = "l")
