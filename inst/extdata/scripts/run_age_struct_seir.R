# Simulate age-structured SEIR compartmental model with 2 vaccine doses

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

# Source functions -------------------------------------------------
source("R/age_struct_seir_ode.R")
source("R/postprocess_age_struct_model_output.R")
source("R/get_foi.R")
source("R/get_beta.R")
source("R/get_vac_rate.R")
source("R/get_vac_rate_2.R")

# load data ---------------------------------------------------------
# probabilities
#dons_probs <- read_xlsx("inst/extdata/data/ProbabilitiesDelays_20210107.xlsx")
p_infection2admission <- c(0.003470, 0.000377, 0.000949, 0.003880, 0.008420, 0.016500,
                           0.025100, 0.049400, 0.046300)
p_admission2death <- c(0.00191, 0.00433, 0.00976, 0.02190, 0.02500, 0.04010,
                       0.10600, 0.22900, 0.31100)

# read in contact matrices
contact_matrices_all <- readRDS("inst/extdata/data/contact_matrices_for_model_input.rds")

# read in vac schedules
old_to_young <- read_xlsx("inst/extdata/data/old_to_young_az_pf_only.xlsx", sheet = 1)
no_vac <- data.frame(date = old_to_young$date, old_to_young[,-1] * 0)
young_to_old <- read_xlsx("inst/extdata/data/old_to_young_az_pf_only.xlsx", sheet = 2)
alternative <- read_xlsx("inst/extdata/data/old_to_young_az_pf_only.xlsx", sheet = 3)

# age distribution and pop size
age_dist <- c(0.10319920, 0.11620856, 0.12740219, 0.12198707, 0.13083463, 
              0.14514332, 0.12092904, 0.08807406, 0.04622194)
n_age_groups <- length(age_dist)
n <- 17407585                                 # Dutch population size

# proportion initially infected (from coronadashboard.nl)
infectious_total <- 103000
prob_inf_by_age <- c(0.018, 0.115, 0.156, 0.118, 0.142, 0.199, 0.114, 0.062, 0.054 + 0.023)          

# proportion initially recovered (seroprevalence - from Don)
# Method: extrapolated seroprev = [measured seroprev September] * 
#                                 [cumulative hospitalisations today] / [cumulative hospitalisations September]
prop_rec = c(0.01120993, 0.09663659, 0.24141186, 0.11004723, 0.10677859, 
             0.11977255, 0.11904044, 0.11714503, 0.11347191)

# vaccinations per day
ve <- list(pfizer = c(0.926, 0.948), 
           moderna = c(0.896, 0.941), 
           astrazeneca = c(0.583, 0.621))

h <- p_infection2admission
init_states <- list(E = c(3044 * 3 * prob_inf_by_age), # current number of cases * 8
                    I = c(infectious_total * prob_inf_by_age), # number of infectious on Feb 1, 2021 (from corona dashboard) 
                    H = c(h/sum(h) * 270), # from NICE data for Feb 1, 2021
                    R = c(n * age_dist * prop_rec))

test_state <- c(rep(100,n_age_groups)) # for testing
empty_state <- c(rep(0,n_age_groups))
# lsqfit in matlab

# Input parameters:
n_vec <- n * age_dist
c1 <- as.matrix(contact_matrices_all$baseline[,-1])
c2 <- as.matrix(contact_matrices_all$april2020[,-1])
c3 <- as.matrix(contact_matrices_all$september2020[,-1])
s <- 0.5 #0.2174  # from David
g <- 0.5 #0.4762  # from David
susceptibles <- n_vec-init_states$E-init_states$I-init_states$R
tmp <- get_beta(R0 = 3, 
                contact_matrix = c1, 
                N = n_vec, 
                sigma = s, 
                gamma = g,
                Reff = 1.13,
                contact_matrix2 = c2,
                init_s = susceptibles,
                eta = 1-ve$pfizer[1],
                eta2 = 1-ve$pfizer[2]
                ) 
b <- tmp$beta 
b2 <- tmp$beta2

tag <- "ccm_no_vac"

params <- list(beta = b2,                      # transmission rate
               gamma = g,                      # 1/gamma = infectious period
               sigma = s,                      # 1/sigma = latent period
               N = n_vec,                      # Population (no need to change)
               vac_per_day = NULL,             # Number of vaccines per day (dose 1)
               vac_per_day2 = NULL,            # Number of vaccines per day (dose 2)
               t_ve_pf = 1,                       # Time vaccination starts (dose 1)
               t_ve_az = 126,                      # Time vaccination starts (dose 2)
               delay = 14,                     # Delay from vaccination to protection (days)
               delay2 = 14,                    # Delay for dose 2
               #eta = 1- 0.70,                 # 1 - VE (dose 1)
               #eta2 = 1- 0.90,                # 1 - VE (dose 2)
               #uptake = 0.85,                 # Proportion of population able and willing to be vaccinated
               h = p_infection2admission,      # Rate from infection to hospital admission
               d = p_admission2death,          # Rate from admission to death
               r = 0.0206,                     # Rate from admission to recovery
               c_start = c2,
               c_lockdown = c2,
               c_relaxed = c3,
               vac_schedule = no_vac,
               ve = ve,
               no_vac = FALSE,
               constant_contact_matrix = FALSE,
               one_cp = TRUE
)

# Specify initial values -------------------------------------------
t_max <- dim(old_to_young)[1]
times <- seq(0,t_max,length.out = t_max+1)     # Vector of times
timeInt <- times[2]-times[1]             # Time interval (for technical reasons)
init <- c(t = times[1],                  # Initial conditions
          S = params$N - (init_states$E + init_states$I + init_states$R),
          Shold_1d = empty_state,
          Sv_1d = empty_state,
          Shold_2d = empty_state,
          Sv_2d = empty_state,
          E = init_states$E,
          Ev_1d = empty_state,
          Ev_2d = empty_state,
          I = init_states$I,
          Iv_1d = empty_state,
          Iv_2d = empty_state,
          H = init_states$H,
          Hv_1d = empty_state,
          Hv_2d = empty_state,
          D = empty_state,
          R = init_states$R,
          Rv_1d = empty_state,
          Rv_2d = empty_state
          )                      

# Solve model ------------------------------------------------------
seir_out <- lsoda(init,times,age_struct_seir_ode,params)
seir_out <- as.data.frame(seir_out)
out <- postprocess_age_struct_model_output(seir_out)

# Summarise results ------------------------------------------------
beta <- params$beta * timeInt
res <- get_vac_rate_2(times, params$vac_schedule, params$ve, 
                      t_ve_pf = params$t_ve_pf, t_ve_az = params$t_ve_az)
eta_dat <- res %>%
  select(time, age_group, eta) %>%
  pivot_wider(names_from = age_group, names_prefix = "eta",
              values_from = eta)
eta2_dat <- res %>%
  select(time, age_group, eta2) %>%
  pivot_wider(names_from = age_group, names_prefix = "eta",
              values_from = eta2)
N <- params$N
h <- params$h
gamma <- params$gamma
contact_mat <- params$C
#time_inf_to_hosp <- 11

lambda_est <- get_foi(dat = out, 
                      beta = beta, 
                      c_main = c2, 
                      N = N,
                      constant_contact_matrix = TRUE,
                      c_lockdown = c2,
                      c_relaxed = c3,
                      times = times,
                      eta = eta_dat[,-1],
                      eta2 = eta2_dat[,-1])
time <- seir_out$time
inc <- (out$S + out$Shold_1d + (eta_dat[,-1] * (out$Sv_1d + out$Shold_2d)) + (eta2_dat[,-1] * out$Sv_2d)) * lambda_est
hosp <- sweep(inc, 2, h, "*")

# quick check
plot(times, rowSums(inc), type = "l")
plot(times, rowSums(hosp), type = "l")

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
  pivot_longer(cols = starts_with("S"), 
               names_to = "age_group", 
               names_prefix = "S",
               values_to = "hosp_admissions")

df <- left_join(inc_long, hosp_long, by = c("time", "age_group")) %>%
  pivot_longer(cols = c("incidence", "hosp_admissions"),
               names_to = "outcome",
               values_to = "value") %>%
  mutate(date = time + as.Date("2021-02-01")) %>%
  select(time, date, age_group, outcome, value)


df_summary <- df %>%
  group_by(time, date, outcome) %>%
  summarise_at(.vars = "value", .funs = sum) %>%
  mutate(scenario = tag)

saveRDS(df_summary, paste0("inst/extdata/results/res_",tag,".rds"))

# Make summary table ------------------------------------------------
summary_tab <- df_summary %>%
  group_by(outcome) %>%
  summarise_at(.vars = "value", .funs = sum) 

# Make plot ---------------------------------------------------------
# summary over all age groups
p <- ggplot(df_summary %>% filter(time > 3), aes(x = date, y = value, color = outcome)) +
  geom_line() +
  labs(y = "Value", x = "Time (days)", color = "Outcome") +
  ylim(0,NA) +
  scale_x_date(date_breaks = "2 weeks", date_labels = "%d %b") +
  geom_vline(xintercept = as.Date("2021-06-07"),linetype="dashed", color = "grey70") +
  theme(legend.position = "bottom",
        panel.background = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_wrap(~outcome, scales = "free")
p

# by age group
p2 <- ggplot(df %>% filter(time > 3,
                           age_group == 4), aes(x = date, y = value, color = age_group)) +
  geom_line() +
  labs(y = "Value", x = "Time (days)", color = "Outcome") +
  ylim(0,1000) +
  scale_x_date(date_breaks = "2 weeks", date_labels = "%d %b") +
  geom_vline(xintercept = as.Date("2021-06-07"),linetype="dashed", color = "grey70") +
  theme(legend.position = "bottom",
        panel.background = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_wrap(~outcome, scales = "free")
p2

