# Simulate age-structured SEIR compartmental model with 2 vaccine doses

# Load packages ----------------------------------------------------
library(deSolve)
library(reshape2)
library(ggplot2)
#library(readxl)
source("R/age_struct_seir_ode.R")

# Input parameters:
params <- list(beta = 0.61, 
               gamma = 0.5,                   # R0 = beta/gamma
               sigma = 0.5,                   # 1/sigma = latent period
               N = c(500000, 500000),                  # Population (no need to change)
               vac_per_day = 25000,           # Number of vaccines per day (dose 1)
               vac_per_day2 = 25000,          # Number of vaccines per day (dose 2)
               tv = 10,                       # Time vaccination starts (dose 1)
               tv2 = 50,                      # Time vaccination starts (dose 2)
               delay = 14,                    # Delay from vaccination to protection (days)
               delay2 = 14,                   # Delay for dose 2
               eta = 1- 0.58,                 # 1 - VE (dose 1)
               eta2 = 1- 0.62,                # 1 - VE (dose 2)
               uptake = 0.85,                 # Proportion of population able and willing to be vaccinated
               h = 0.0251,                    # Rate from infection to hospital admission
               d = 0.106,                     # Rate from admission to death
               r = 0.0206,                    # Rate from admission to recovery
               C = matrix(c(3,2,2,3), nrow = 2)
) 

# Specify initial values -------------------------------------------
times <- seq(0,200,length.out = 201)     # Vector of times
timeInt <- times[2]-times[1]             # Time interval (for technical reasons)
init <- c(t = times[1],                  # Initial conditions
          S = params$N - c(20,50),
          Shold = c(0,0),
          Sv = c(0,0),
          Shold2 = c(0,0),
          Sv2 = c(0,0),
          E = c(0,0),
          Ev = c(0,0),
          Ev2 = c(0,0),
          I = c(20,50),
          Iv = c(0,0),
          Iv2 = c(0,0),
          H = c(0,0),
          D = c(0,0),
          R = c(0,0),
          Rv = c(0,0),
          Rv2 = c(0,0))                      

# Solve model ------------------------------------------------------
seir_out <- lsoda(init,times,age_struct_seir_ode,params)
seir_out <- as.data.frame(seir_out)
# Summarise results ------------------------------------------------
beta <- params$beta * timeInt
eta <- params$eta
eta2 <- params$eta2
N <- params$N
h <- params$h
gamma <- params$gamma
C <- params$C

lambda <- beta * C %*% rowSums(seir_out[,c(11:13)])/N
time <- seir_out$time
S <- seir_out$S1 + seir_out$S2
Shold <- seir_out$Shold1 + seir_out$Shold2
Sv <- seir_out$Sv1 + seir_out$Sv2
Shold2 <- seir_out$Shold21 + seir_out$Shold22
Sv2 <- seir_out$Sv21 + seir_out$Sv22
inc <- (S + Shold + (eta * (Sv + Shold2)) + (eta2 * Sv2)) * lambda
I <- seir_out[,11]
Iv <- seir_out[,12]
Iv2 <- seir_out[,13]
hosp <- h * (I + Iv + Iv2)

# Create object for plotting:
df <- data.frame(time = time, 
                 incidence = inc, 
                 hosp_admissions = hosp,
                 hosp_admissions_from_inc = hosp2)
