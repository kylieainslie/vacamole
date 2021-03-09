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
source("R/summarise_results.R")

# load data ---------------------------------------------------------
# probabilities
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
c5 <- readRDS("inst/extdata/data/Contactpatterns_PICO4_10y.rds")
# population distributions
demo <- data.frame(age_group = c("[0,10)", "[10,20)", "[20,30)", "[30,40)","[40,50)", "[50,60)",
                                 "[60,70)", "[70,80)", "[80,Inf]"), 
                   frac_2017 = c(0.105, 0.118, 0.126, 0.120, 0.138, 0.145, 0.122, 0.0808, 0.0447), 
                   frac_2019 = c(0.103, 0.116, 0.127, 0.122, 0.131, 0.145, 0.121, 0.0881, 0.0462),
                   frac_2020 = c(0.102, 0.115, 0.128, 0.123, 0.127, 0.145, 0.121, 0.0904, 0.0472))

# relative susceptibility/infectiousness
rel_trans <- c(1.000, 3.051, 5.751, 3.538, 3.705, 4.365, 5.688, 5.324, 7.211)
get_transmission_matrix <- function(x, contact_mat, year){
  # first make symmetric
  if (year == 2017){
    cm <- contact_mat / demo$frac_2017
  } else if (year == 2019){
    cm <- contact_mat / demo$frac_2019
  } else if (year == 2020) {
    cm <- contact_mat / demo$frac_2020
  }
  # multiply by relative susc/inf
  tmp <- sweep(cm, 1, x, "*")
  rtn <- sweep(tmp, 2, x, "*")
  # output
  return(rtn)
}

t1 <- get_transmission_matrix(rel_trans, c1, year = 2017)
t2 <- get_transmission_matrix(rel_trans, c2, year = 2019)
t3 <- get_transmission_matrix(rel_trans, c3, year = 2019)
t4 <- get_transmission_matrix(rel_trans, c4, year = 2019)
t5 <- get_transmission_matrix(rel_trans, c5, year = 2020)

# parameter inputs
s <- 0.5
g <- 0.5
r0 <- 2.3 
reff <- 1.13
# determine transmission rate for r0
S = diag(n_vec - 1)
rho = as.numeric(eigs(S %*% t1,1)$values)
beta = r0/rho
# check
K = beta * S %*% t1
as.numeric(eigs(K,1)$values) # this should be r0

# determine transmission rate fro reff
S = diag(n_vec - 1)
rho = as.numeric(eigs(S %*% t1,1)$values)
beta = r0/rho
# check
K = beta * S %*% t1
as.numeric(eigs(K,1)$values) # this should be r0



# tmp <- get_beta(R0 = r0, contact_matrix = t1, N = n_vec, sigma = s, 
#                 gamma = g) 
# beta <- tmp$beta
h <- p_infection2admission/time_symptom2admission
i1 <- p_admission2IC/time_admission2IC
i2 <- p_IC2hospital/time_IC2hospital
d <- p_admission2death/time_admission2death
d_ic <- p_IC2death/time_IC2death
d_hic <- p_hospital2death/time_hospital2death
r <- (1 - p_admission2death)/time_admission2discharge
r_ic <- (1 - p_IC2death)/time_hospital2discharge

# read in vac schedules
# this now comes from convert_vac_schedule.R
vac_schedule <- vac_schedule_orig_new
# old_to_young <- read_xlsx("inst/extdata/data/old_to_young_az_pf_only.xlsx", sheet = 1)
# no_vac <- data.frame(date = old_to_young$date, old_to_young[,-1] * 0)
# young_to_old <- read_xlsx("inst/extdata/data/old_to_young_az_pf_only.xlsx", sheet = 2)
# alternative <- read_xlsx("inst/extdata/data/old_to_young_az_pf_only.xlsx", sheet = 3)

# vaccinations params
ve <- list(pfizer = c(0.926, 0.948), 
           moderna = c(0.896, 0.941), 
           astrazeneca = c(0.583, 0.621))

delays <- list(pfizer = c(14, 7), 
               moderna = c(14, 14), 
               astrazeneca = c(21,14))

# ve_sa <- list(pfizer = c(0.51, 0.948), 
#               moderna = c(0.51, 0.941), 
#               astrazeneca = c(0.676, 0.495))

# initial states
# for vaccinated
# before_feb object is created in convert_vac_schedule.R
# init_sv_1d <- unlist(before_feb %>% select(pf_d1_1:pf_d1_9) + before_feb %>% select(mo_d1_1:mo_d1_9)) * n_vec
# init_sv_2d <- unlist(before_feb %>% select(pf_d2_1:pf_d2_9) + before_feb %>% select(mo_d2_1:mo_d2_9)) * n_vec

# init_states <- list(Sv_1d = init_sv_1d,
#                     Sv_2d = init_sv_2d,
#                     E = c((3244 / (s * p_reported_all)) * p_inf_by_age), # number of infections that will lead to same number 
#                         # of new cases as of Feb 1, 2021
#                     I = c(77000 * p_inf_by_age), # number of infectious on Feb 1, 2021 (from corona dashboard) so that
#                         # at the initial time step there are 215 hospital admissions
#                     H = c(p_infection2admission/sum(p_infection2admission) * 500), # number of hospitalisations so that at initial
#                         # time step there are 31 IC admissions (from coronadashboard data for Feb 1, 2021)
#                     H_IC = c(p_IC2hospital/sum(p_IC2hospital) * 1131), # so that total number of hospital occupancy (exc IC) is 1631 (from coronadashboard for 1 Feb 2021)
#                     IC = c(p_admission2IC/sum(p_admission2IC) * 639), # from coronadashboard Feb 1, 2021
#                     R = c(3000000 * p_recovered))

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
                              n_hosp = c(4, 3, 17, 34, 74, 234, 463, 551, 452),
                              n_ic = c(0, 2, 6, 9, 25, 83, 181, 150, 12)
                              ) %>%
  mutate(n_infections = n_cases * 3,
         init_E = n_infections * (2/7),
         init_I = n_infections * (2/7),
         init_S = n - n_recovered - init_E - init_I - n_hosp - n_ic)

empty_state <- c(rep(0,9))

# Specify initial values -------------------------------------------
t_max <- dim(vac_schedule)[1] - 1
times <- seq(0,t_max, by = 1)     # Vector of times
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
          H = empty_state, #init_states_dat$n_hosp,
          Hv_1d = empty_state,
          Hv_2d = empty_state,
          H_IC = empty_state,
          H_ICv_1d = empty_state,
          H_ICv_2d = empty_state,
          IC = empty_state, #init_states_dat$n_ic,
          ICv_1d = empty_state,
          ICv_2d = empty_state,
          D = empty_state,
          R = init_states_dat$n_recovered,
          Rv_1d = empty_state,
          Rv_2d = empty_state
)                      

# Input parameters -------------------------------------------------
params <- list(beta = beta,                   # transmission rate
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
               c_lockdown = t2,
               c_relaxed = t4,
               c_very_relaxed = t3,
               c_normal = t1,
               vac_schedule = vac_schedule,
               ve = ve,
               delay = delays,
               use_cases = FALSE,              # use cases as criteria to change contact matrices. If FALSE, IC admissions used.
               thresh_l = 3, #5/100000 * sum(n_vec),            # 3 for IC admissions
               thresh_m = 10, #14.3/100000 * sum(n_vec),        # 10 for IC admissions
               thresh_u = 20, #35.7/100000 * sum(n_vec),        # 20 for IC admissions
               #thresh_cushion = 1/100000 * sum(n_vec),      # cushion so integrator doesn't get stuck at change point (0 for IC)
               no_vac = FALSE,
               force_relax = NULL                              # time step when measures are forced to relax regardless of criteria
)

# Solve model ------------------------------------------------------
seir_out <- lsoda(init,times,age_struct_seir_ode,params) #
seir_out <- as.data.frame(seir_out)
out <- postprocess_age_struct_model_output(seir_out)

# quick check ------------------------------------------------------
infs <- (out$I + out$Iv_1d + out$Iv_2d)
plot(times, rowSums(infs), type = "l")
ic_admin <- sweep(out$H + out$Hv_1d + out$Hv_2d, 2, i1, "*")
plot(times, rowSums(ic_admin), type = "l")

# Summarise results ------------------------------------------------
tag <- "original_sa"

results <- summarise_results(out, params, start_date = "2021-02-01")
saveRDS(results$df_summary, paste0("inst/extdata/results/res_",tag,".rds"))

# Make summary table ------------------------------------------------
summary_tab <- results$df_summary %>%
  group_by(outcome) %>%
  summarise_at(.vars = "value", .funs = sum) 

# Make plot ---------------------------------------------------------
# summary over all age groups
p <- ggplot(results$df_summary, aes(x = date, y = value, color = outcome)) +
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

