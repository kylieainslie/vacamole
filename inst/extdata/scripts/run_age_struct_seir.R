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

# parameter inputs
s <- 0.2
g <- 0.125
r0 <- 2.95550 
init_i <- 0.489163
tmp <- get_beta(R0 = r0, contact_matrix = c1, N = n_vec, sigma = s, 
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

# contact matrices
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

# read in vac schedules
old_to_young <- read_xlsx("inst/extdata/data/old_to_young_az_pf_only.xlsx", sheet = 1)
no_vac <- data.frame(date = old_to_young$date, old_to_young[,-1] * 0)
young_to_old <- read_xlsx("inst/extdata/data/old_to_young_az_pf_only.xlsx", sheet = 2)
alternative <- read_xlsx("inst/extdata/data/old_to_young_az_pf_only.xlsx", sheet = 3)

# vaccinations params
ve <- list(pfizer = c(0.926, 0.948), 
           moderna = c(0.896, 0.941), 
           astrazeneca = c(0.583, 0.621))

delays <- list(pfizer = c(14, 7), 
               moderna = c(14, 14), 
               astrazeneca = c(21,141))

# initial states
init_states <- list(E = c((3245 / p_reported_all) * p_inf_by_age * 8.5), # cases as of Feb 1, 2021
                    I = c(96617 * p_inf_by_age), # number of infectious on Feb 1, 2021 (from corona dashboard) 
                    H = c(p_infection2admission/sum(p_infection2admission) * 1631), # from NICE data for Feb 1, 2021
                    IC = c(p_admission2IC/sum(p_admission2IC) * 639), # from coronadashboard Feb 1, 2021
                    R = c(3000000 * p_recovered))

empty_state <- c(rep(0,n_age_groups))

tag <- "1cp_old_to_young"

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
               contact_mat = c2,
               c_relaxed = c4,
               vac_schedule = old_to_young,
               delay = delays,
               ic_thresh_l = 10,
               ic_thresh_u = 25,
               no_vac = FALSE
)

# Specify initial values -------------------------------------------
t_max <- dim(old_to_young)[1]
times <- seq(0,t_max,length.out = t_max+1)     # Vector of times
timeInt <- times[2]-times[1]             # Time interval (for technical reasons)
init <- c(t = times[1],                  
          S = params$N - init_states$E - init_states$I - init_states$H - init_states$IC - init_states$R,
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
          H_IC = empty_state,
          H_ICv_1d = empty_state,
          H_ICv_2d = empty_state,
          IC = init_states$IC,
          ICv_1d = empty_state,
          ICv_2d = empty_state,
          D = empty_state,
          R = init_states$R,
          Rv_1d = empty_state,
          Rv_2d = empty_state
)                      


# Solve model ------------------------------------------------------
seir_out <- lsoda(init,times,age_struct_seir_ode,params)
seir_out <- as.data.frame(seir_out)
out <- postprocess_age_struct_model_output(seir_out)

# quick plots
plot(times, rowSums(out$I), type = "l")
plot(times, rowSums(out$H), type = "l")
plot(times, rowSums(out$IC), type = "l")

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
                      eta2 = eta2_dat[,-1],
                      one_cp = params$one_cp)
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
p2 <- ggplot(df %>% filter(time > 3), aes(x = date, y = value, color = age_group)) +
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

