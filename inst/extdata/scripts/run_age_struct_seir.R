# Simulate age-structured SEIR compartmental model with 2 vaccine doses

# Load packages ----------------------------------------------------
library(deSolve)
library(reshape2)
library(ggplot2)
library(dplyr)
library(stringr)
library(tidyr)
library(readxl)

# Source functions -------------------------------------------------
source("R/age_struct_seir_ode.R")
source("R/postprocess_age_struct_model_output.R")
source("R/get_foi.R")
source("R/get_beta.R")

# load data ---------------------------------------------------------
# probabilities
dons_probs <- read_xlsx("inst/extdata/data/ProbabilitiesDelays_20210107.xlsx")
# contact matrix
contact_matrices_all <- read.delim("inst/extdata/data/S2_contact_matrices_withPico3.tsv", header=TRUE, allowEscapes=FALSE, sep="\t")
contact_matrix_april2020 <- contact_matrices_all %>%
  filter(survey == "April 2020") %>%
  filter(contact_type == "community") %>%
  select(-survey, -contact_type) %>%
  pivot_wider(., names_from = cont_age, values_from = m_est)

contact_matrix_baseline <- contact_matrices_all %>%
  filter(survey == "baseline") %>%
  filter(contact_type == "community") %>%
  select(-survey, -contact_type) %>%
  pivot_wider(., names_from = cont_age, values_from = m_est)

contact_matrix_input <- as.matrix(contact_matrix_april2020[-1,-c(1,2)])
#rownames(contact_matrix_input) <- contact_matrix_april2020$part_age[-1]

# age distribution and pop size
age_dist <- c(0.10319920, 0.11620856, 0.12740219, 0.12198707, 0.13083463, 
              0.14514332, 0.12092904, 0.08807406, 0.04622194)
n_age_groups <- length(age_dist)
n <- 17407585                                 # Dutch population size

# proportion initially infected (from coronadashboard.nl)
infectious_total <- 113264
prob_inf_by_age <- c(0.018, 0.115, 0.156, 0.118, 0.142, 0.199, 0.114, 0.062, 0.054 + 0.023)          

# proportion initially recovered (seroprevalence - from Don)
# Method: extrapolated seroprev = [measured seroprev September] * [cumulative hospitalisations today] / [cumulative hospitalisations September]
prop_rec = c(0.01120993, 0.09663659, 0.24141186, 0.11004723, 0.10677859, 0.11977255, 0.11904044, 0.11714503, 0.11347191)

# vaccinations per day
vac_scheduleA <- read_csv("inst/extdata/data/cum_upt_A.csv")
ve <- list(pfizer = c(0.926, 0.948), 
           moderna = c(0.896, 0.941), 
           astrazeneca = c(0.583, 0.621))

init_states <- list(E = c(4060 * 8 * prob_inf_by_age),
                    I = c(infectious_total * prob_inf_by_age),
                    R = c(n * age_dist * prop_rec))
test_state <- c(rep(100,n_age_groups))
empty_state <- c(rep(0,n_age_groups))

# Input parameters:
c <- as.matrix(contact_matrix_baseline[-1,-c(1,2)])
s <- 0.2174  # from David
g <- 0.4762  # from David
b <- get_beta(R0 = 3,contact_matrix = c, N = n * age_dist,sigma = s,gamma = g) 

params <- list(beta = b, 
               gamma = g,                   # R0 = beta/gamma
               sigma = s,                   # 1/sigma = latent period
               N = n * age_dist,              # Population (no need to change)
               vac_per_day = 0,           # Number of vaccines per day (dose 1)
               vac_per_day2 = 0,              # Number of vaccines per day (dose 2)
               tv = 14,                       # Time vaccination starts (dose 1)
               tv2 = 56,                      # Time vaccination starts (dose 2)
               delay = 14,                    # Delay from vaccination to protection (days)
               delay2 = 14,                   # Delay for dose 2
               eta = 1- 0.70,                 # 1 - VE (dose 1)
               eta2 = 1- 0.90,                # 1 - VE (dose 2)
               uptake = 0.85,                 # Proportion of population able and willing to be vaccinated
               h = dons_probs$P_infection2admission, # Rate from infection to hospital admission
               d = dons_probs$P_admission2death,     # Rate from admission to death
               r = 0.0206,                    # Rate from admission to recovery
               C = c,
               vac_schedule = vac_scheduleA,
               ve = ve,
               constant_foi = FALSE,
               init_inf = init_states$I
)

# Specify initial values -------------------------------------------
times <- seq(0,200,length.out = 201)     # Vector of times
timeInt <- times[2]-times[1]             # Time interval (for technical reasons)
init <- c(t = times[1],                  # Initial conditions
          S = params$N - test_state, #(init_states$E + init_states$I + init_states$R),
          Shold_1d = empty_state,
          Sv_1d = empty_state,
          Shold_2d = empty_state,
          Sv_2d = empty_state,
          E = empty_state, #init_states$E,
          Ev_1d = empty_state,
          Ev_2d = empty_state,
          I = test_state, #init_states$I,
          Iv_1d = empty_state,
          Iv_2d = empty_state,
          H = empty_state,
          Hv_1d = empty_state,
          Hv_2d = empty_state,
          D = empty_state,
          R = empty_state, #init_states$R,
          Rv_1d = empty_state,
          Rv_2d = empty_state
          )                      

# Solve model ------------------------------------------------------
seir_out <- lsoda(init,times,age_struct_seir_ode,params)
seir_out <- as.data.frame(seir_out)
out <- postprocess_age_struct_model_output(seir_out)

# quick check
plot(times, out$I[,1], type = "l")

# Summarise results ------------------------------------------------
beta <- params$beta * timeInt
eta <- params$eta
eta2 <- params$eta2
N <- params$N
h <- params$h
gamma <- params$gamma
C <- params$C
time_inf_to_hosp <- 11

lambda <- get_foi(dat = out, beta = beta, contact_matrix = C, N = N)
time <- seir_out$time
inc <- (out$S + out$Shold_1d + (eta * (out$Sv_1d + out$Shold_2d)) + (eta2 * out$Sv_2d)) * lambda
hosp <- sweep(inc, 2, h, "*")

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
               values_to = "value") # %>%
  # mutate(outcome = factor(outcome, levels = c("Incidence", "Hospital Admissions")))

df_summary <- df %>%
  group_by(time, outcome) %>%
  summarise_at(.vars = "value", .funs = sum)

# Make summary table ------------------------------------------------
summary_tab <- df_summary %>%
  group_by(outcome) %>%
  summarise_at(.vars = "value", .funs = sum)

# Make plot ---------------------------------------------------------
# summary over all age groups
g_sum <- ggplot(df_summary, aes(x = time, y = value, color = outcome)) +
  geom_line() +
  labs(y = "Value", x = "Time (days)", color = "Age Group") +
  theme(legend.position = "bottom",
        panel.background = element_blank()) +
  facet_wrap(~outcome, scales = "free")
g_sum

# stratify by age group
g_age <- ggplot(df, aes(x = time, y = value, color = age_group)) +
  geom_line() +
  labs(y = "Value", x = "Time (days)", color = "Age Group") +
  theme(legend.position = "bottom",
        panel.background = element_blank()) +
  facet_wrap(~outcome, scales = "free")
g_age
