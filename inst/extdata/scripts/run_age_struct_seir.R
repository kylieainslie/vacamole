# Simulate age-structured SEIR compartmental model with 2 vaccine doses

# Load packages ----------------------------------------------------
library(deSolve)
library(reshape2)
library(ggplot2)
library(dplyr)
library(stringr)
library(tidyr)
#library(readxl)
source("R/age_struct_seir_ode.R")
source("R/postprocess_age_struct_model_output.R")
source("R/get_foi.R")

# Input parameters:
params <- list(beta = 0.61, 
               gamma = 0.5,                   # R0 = beta/gamma
               sigma = 0.5,                   # 1/sigma = latent period
               N = c(500000, 500000),         # Population (no need to change)
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
          Shold_1d = c(0,0),
          Sv_1d = c(0,0),
          Shold_2d = c(0,0),
          Sv_2d = c(0,0),
          E = c(0,0),
          Ev_1d = c(0,0),
          Ev_2d = c(0,0),
          I = c(20,50),
          Iv_1d = c(0,0),
          Iv_2d = c(0,0),
          H = c(0,0),
          Hv_1d = c(0,0),
          Hv_2d = c(0,0),
          D = c(0,0),
          R = c(0,0),
          Rv_1d = c(0,0),
          Rv_2d = c(0,0)
          )                      

# Solve model ------------------------------------------------------
seir_out <- lsoda(init,times,age_struct_seir_ode,params)
seir_out <- as.data.frame(seir_out)
out <- postprocess_age_struct_model_output(seir_out)

# Summarise results ------------------------------------------------
beta <- params$beta * timeInt
eta <- params$eta
eta2 <- params$eta2
N <- params$N
h <- params$h
gamma <- params$gamma
C <- params$C

lambda <- get_foi(dat = out, beta = beta, contact_matrix = C, N = N)
time <- seir_out$time
inc <- (out$S + out$Shold_1d + (eta * (out$Sv_1d + out$Shold_2d)) + (eta2 * out$Sv_2d)) * lambda
hosp <- h * (out$I + out$Iv_1d + out$Iv_2d)

# Create object for plotting ---------------------------------------
# convert from wide to long format
inc_long <- inc %>% 
  mutate(time = time) %>%
  pivot_longer(cols = starts_with("S"), 
               names_to = "age_group", 
               names_prefix = "S",
               values_to = "incidence")

hosp_long <- hosp %>%
  mutate(time = time) %>%
  pivot_longer(cols = starts_with("I"), 
               names_to = "age_group", 
               names_prefix = "I",
               values_to = "hosp_admissions")

df <- left_join(inc_long, hosp_long, by = c("time", "age_group")) %>%
  pivot_longer(cols = c("incidence", "hosp_admissions"),
               names_to = "outcome",
               values_to = "value") # %>%
  # mutate(outcome = factor(outcome, levels = c("Incidence", "Hospital Admissions")))

# Make plot ---------------------------------------------------------
g <- ggplot(df, aes(x = time, y = value, color = age_group)) +
  geom_line() +
  labs(y = "Value", x = "Time (days)", color = "Age Group") +
  theme(legend.position = "bottom",
        panel.background = element_blank()) +
  facet_wrap(~outcome, scales = "free")
plot(g)

