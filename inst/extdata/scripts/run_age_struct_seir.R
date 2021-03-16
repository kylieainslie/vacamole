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
source("R/choose_contact_matrix.R")
source("R/get_foi.R")
source("R/get_vac_rate.R")
source("R/get_vac_rate_2.R")
source("R/summarise_results.R")
source("inst/extdata/scripts/convert_vac_schedule.R")

# load data ---------------------------------------------------------
# probabilities -----------------------------------------------------
dons_probs <- read_xlsx("inst/extdata/data/ProbabilitiesDelays_20210107.xlsx")
p_infection2admission <- dons_probs$P_infection2admission
p_admission2death <- dons_probs$P_admission2death
p_admission2IC <- dons_probs$P_admission2IC
p_IC2hospital <- dons_probs$P_IC2hospital
p_IC2death <- 1-p_IC2hospital
p_hospital2death <- c(rep(0,5), 0.01, 0.04, 0.12, 0.29) #(after ICU)
p_reported_by_age <- c(0.29, 0.363, 0.381, 0.545, 0.645, 0.564, 0.365, 0.33, 0.409) # from Jantien
p_reported_all <- 0.428 # from Jantien
p_inf_by_age <- c(0.018, 0.115, 0.156, 0.118, 0.142, 0.199, 0.114, 0.062, 0.054 + 0.023)          
p_recovered <- c(0.01120993, 0.09663659, 0.24141186, 0.11004723, 0.10677859, 0.11977255, 
                0.11904044, 0.11714503, 0.11347191)

# delays ------------------------------------------------------------
time_symptom2admission <- c(2.29, 5.51, 5.05, 5.66, 6.55, 5.88, 5.69, 5.09, 4.33) # assume same as infectious2admission
time_admission2discharge <- 7.9
time_admission2IC <- 2.28
time_IC2hospital <- 15.6
time_hospital2discharge <- 10.1 #(after ICU)
time_admission2death <- 7 
time_IC2death <- 19 
time_hospital2death <- 10 #(after ICU)

# age distribution and pop size -------------------------------------
age_dist <- c(0.10319920, 0.11620856, 0.12740219, 0.12198707, 0.13083463, 
              0.14514332, 0.12092904, 0.08807406, 0.04622194)
n <- 17407585                                 # Dutch population size
n_vec <- n * age_dist

# contact matrices --------------------------------------------------
contact_matrices_all <- readRDS("inst/extdata/data/contact_matrices_for_model_input.rds")
c1 <- as.matrix(contact_matrices_all$baseline_sym[,-1])
c2 <- as.matrix(contact_matrices_all$april2020_sym[,-1])
c3 <- as.matrix(contact_matrices_all$june2020_sym[,-1])
c4 <- as.matrix(contact_matrices_all$september2020_sym[,-1])
c5 <- as.matrix(contact_matrices_all$february2021_sym[,-1])

# Convert from number of contacts to rate of contacts ---------------
N_diag <- diag(1/n_vec)
m1 <- c1 %*% N_diag
m2 <- c2 %*% N_diag
m3 <- c3 %*% N_diag
m4 <- c4 %*% N_diag
m5 <- c5 %*% N_diag
  
# relative susceptibility/infectiousness ----------------------------
rel_trans <- c(1.000, 3.051, 5.751, 3.538, 3.705, 4.365, 5.688, 5.324, 7.211)
get_transmission_matrix <- function(x, contact_mat){
  # multiply by relative susc/inf
  tmp <- sweep(contact_mat, 1, x, "*")  # rows
  rtn <- sweep(tmp, 2, x, "*") # columns
  
  # output
  return(rtn)
}

t1 <- get_transmission_matrix(rel_trans, m1)
t2 <- get_transmission_matrix(rel_trans, m2)
t3 <- get_transmission_matrix(rel_trans, m3)
t4 <- get_transmission_matrix(rel_trans, m4)
t5 <- get_transmission_matrix(rel_trans, m5)

# parameter inputs -------------------------------------------------
s <- 0.5
g <- 0.5
r0 <- 2.3 

# determine transmission rate (beta) for r0 ------------------------
S = diag(n_vec - 1)
rho = as.numeric(eigs(S %*% t1,1)$values)
beta = (r0/rho) * g
# check
K = (1/g) * beta * S %*% t1
as.numeric(eigs(K,1)$values) # this should be r0

# define state transition rates ------------------------------------
h <- p_infection2admission/time_symptom2admission
i1 <- p_admission2IC/time_admission2IC
i2 <- p_IC2hospital/time_IC2hospital
d <- p_admission2death/time_admission2death
d_ic <- p_IC2death/time_IC2death
d_hic <- p_hospital2death/time_hospital2death
r <- (1 - p_admission2death)/time_admission2discharge
r_ic <- (1 - p_IC2death)/time_hospital2discharge

# read in vac schedules --------------------------------------------
# this now comes from convert_vac_schedule.R
# prop_o2y
# prop_y2o
# prop_alt

# vaccinations params ----------------------------------------------
ve <- list(pfizer = c(0.926, 0.948), 
           moderna = c(0.896, 0.941), 
           astrazeneca = c(0.583, 0.621),
           jansen = c(0.661))

delays <- list(pfizer = c(14, 7), 
               moderna = c(14, 14), 
               astrazeneca = c(21,14),
               jansen = c(14))

ve_sa20 <- list(pfizer = c(0.741, 0.758), 
           moderna = c(0.717, 0.753), 
           astrazeneca = c(0.466, 0.497),
           jansen = c(0.529))

ve_sa50 <- list(pfizer = c(0.463, 0.474), 
                moderna = c(0.448, 0.471), 
                astrazeneca = c(0.269, 0.311),
                jansen = c(0.331))

# initial states ---------------------------------------------------
# Jacco's suggested way to determine initial conditions
init_states_dat <- data.frame(age_group = c("0-9", "10-19", "20-29", "30-39", "40-49", 
                                            "50-59", "60-69", "70-79","80+"),
                              n = n_vec,
                              # from Scott
                              n_recovered = c(30720, 397100, 642600, 419000, 412200, 505900,
                                              349100, 206800, 115900 + 33200),
                              # from sitrep for 26 januari tot 2 februari: 
                              # https://www.rivm.nl/coronavirus-covid-19/actueel/wekelijkse-update-epidemiologische-situatie-covid-19-in-nederland)
                              n_cases = c(835, 2851, 4591, 3854, 3925, 5191, 3216, 1819, 
                                          1376 + 485),
                              # from NICE data (n_hosp/n_ic refers to occupancy on 1 Feb 2021)
                              n_hosp = c(2, 1, 8, 19, 29, 76, 142, 159, 186),
                              n_ic = c(0, 2, 6, 9, 25, 83, 181, 150, 12),
                              # from NICE data: people with length of stay >= 9 days
                              n_hosp_after_ic = c(2, 2, 9, 15, 45, 158, 321, 392, 266) 
                              ) %>%
  mutate(n_infections = n_cases * 3, #p_reported_by_age,
         init_E = n_infections * (2/7),
         init_I = n_infections * (2/7),
         init_S = n - n_recovered - init_E - init_I  - n_hosp - n_ic - n_hosp_after_ic)

# determine transmission rate for reff ------------------------------
reff <- 1.04 # from RIVM open data for 1 Feb 2021 (midpoint between 
#              0.94 (wt) and 1.13 (UK variant))
S2 = diag(init_states_dat$init_S)
rho2 = as.numeric(eigs(S2 %*% t5,1)$values)
beta2 = reff/rho2 * g
# check
B <- t5
K2 = beta2 * (1/g) * S2 %*% B
as.numeric(eigs(K2,1)$values) # this should be r0

# callibrate distribution of cases across age groups ----------------
w <- eigen(K2)$vectors[,1]/sum(eigen(K2)$vectors[,1]) # should match dist_cases (I think!)
x <- init_states_dat$n_cases / sum(init_states_dat$n_cases)
A <- diag(x/w)
B_prime <- A %*% B
rho2_prime = as.numeric(eigs(S2 %*% B_prime,1)$values)
beta2_prime = reff/rho2_prime * g
K2_prime <- beta2_prime * (1/g) * S2 %*% B_prime
dom_eig_vec <- eigen(K2_prime)$vectors[,1]
w_prime <- dom_eig_vec/sum(dom_eig_vec)
as.numeric(eigs(K2_prime,1)$values)

# Specify initial values -------------------------------------------
empty_state <- c(rep(0,9))
t_max <- dim(vac_schedule)[1] - 1
times <- seq(0,t_max, by = 1)     # Vector of times 242 = Sept 30, 2021
timeInt <- times[2]-times[1]      # Time interval (for technical reasons)
init <- c(t = times[1],                  
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

# Input parameters -------------------------------------------------
params <- list(beta = beta2_prime,           # transmission rate
               gamma = g,                      # 1/gamma = infectious period
               sigma = s,                      # 1/sigma = latent period
               delta = 1,                      # scaling constant 
               N = n_vec,                      # Population (no need to change)
               h = h,                          # Rate from infection to hospital admission/ time from infection to hosp admission
               i1 = i1,
               i2 = i2,
               d = d, 
               d_ic = d_ic,
               d_hic = d_hic,
               r = r,
               r_ic = r_ic,
               p_report = 1/3, #p_reported_by_age,
               c_start = B_prime,
               c_lockdown = B_prime,
               c_relaxed = t4,
               c_very_relaxed = t3,
               c_normal = t1,
               vac_schedule = prop_no_vac_healthy,
               ve = ve,
               delay = delays,
               use_cases = TRUE,              # use cases as criteria to change contact matrices. If FALSE, IC admissions used.
               thresh_n = 0.5/100000 * sum(n_vec),
               thresh_l = 5/100000 * sum(n_vec),           # 3 for IC admissions
               thresh_m = 14.3/100000 * sum(n_vec),        # 10 for IC admissions
               thresh_u = 35.7/100000 * sum(n_vec),        # 20 for IC admissions
               no_vac = TRUE,
               #force_relax = NULL,                          # time step when measures are forced to relax regardless of criteria
               t_start_end = 0#,                           # time step when starting contact matrix ends and criteria are used to decide contact matrix
               #init_lambda = beta2_prime * B_prime %*% init_states_dat$init_I
)

# Solve model ------------------------------------------------------
seir_out <- lsoda(init,times,age_struct_seir_ode,params) #
seir_out <- as.data.frame(seir_out)
out <- postprocess_age_struct_model_output(seir_out)

# quick check ------------------------------------------------------
cases <- (params$sigma * (out$E + out$Ev_1d + out$Ev_2d)) / 3
plot(times, rowSums(cases), type = "l")
hosp_admin <- sweep((out$I + out$Iv_1d + out$Iv_2d), 2, h, "*")
plot(times, rowSums(hosp_admin), type = "l")
ic_admin <- sweep(out$H + out$Hv_1d + out$Hv_2d, 2, i1, "*")
plot(times, rowSums(ic_admin), type = "l")

# Summarise results ------------------------------------------------
tag <- "no_vac_re-lockdown_15March"
results <- summarise_results(out, params, start_date = "2021-01-31", times = times)
saveRDS(results, paste0("inst/extdata/results/res_",tag,".rds"))

# Make summary table ------------------------------------------------
summary_tab <- results$df_summary %>%
  group_by(outcome) %>%
  summarise_at(.vars = "value", .funs = sum) 
summary_tab

# Make plot ---------------------------------------------------------
# summary over all age groups
p <- ggplot(results$df_summary %>%
              filter(outcome == "new_infections"), aes(x = date, y = value, color = outcome)) +
  geom_line() +
  labs(y = "Value", x = "Time (days)", color = "Outcome") +
  ylim(0,NA) +
  scale_x_date(date_breaks = "2 weeks", date_labels = "%d %b") +
  #geom_vline(xintercept = as.Date("2021-06-07"),linetype="dashed", color = "grey70") +
  theme(legend.position = "bottom",
        panel.background = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_wrap(~outcome, scales = "free")
p

# by age group
p2 <- ggplot(results$df, aes(x = date, y = value, color = age_group)) +
  geom_line() +
  labs(y = "Value", x = "Time (days)", color = "Outcome") +
  ylim(0,NA) +
  scale_x_date(date_breaks = "2 weeks", date_labels = "%d %b") +
  geom_vline(xintercept = as.Date("2021-06-07"),linetype="dashed", color = "grey70") +
  theme(legend.position = "bottom",
        panel.background = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_wrap(~outcome, scales = "free")
p2

# plot vaccination
ggplot(results$df_vac, aes(x = date, y = value, color = dose)) +
  geom_line() +
  labs(y = "Value", x = "Time (days)", color = "Dose") +
  ylim(0,NA) +
  scale_x_date(date_breaks = "2 weeks", date_labels = "%d %b") +
  #geom_vline(xintercept = as.Date("2021-06-07"),linetype="dashed", color = "grey70") +
  theme(legend.position = "bottom",
        panel.background = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)) #+
  #facet_wrap(~outcome, scales = "free")

