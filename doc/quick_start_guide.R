## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- install_package, echo=TRUE, eval=TRUE, message=FALSE, warning=FALSE-----
devtools::install_github("kylieainslie/vacamole")
library(vacamole)

# Load other packages
library(lubridate)
library(deSolve)
library(dplyr)
library(tidyr)
library(readxl)
library(kableExtra)

## ----synth_vac_sched_preview, echo=FALSE, eval=TRUE, results='markup'---------
vac_schedule <- read_xlsx("../inst/extdata/inputs/vac_schedule_18plus_synth.xlsx") %>%
  select(-starts_with("X")) %>%
  mutate(date = seq.Date(as.Date("2020-01-01"), as.Date("2020-12-31"), by = 1)) %>%
  filter(date < as.Date("2020-01-15")) %>%
  select(date:pf_d1_9)

vac_schedule %>%
  kbl() %>%
  kable_styling()

## ----convert_vac_sched, echo=TRUE, eval=FALSE, message = FALSE, warning=FALSE----
#  # read in vaccination schedule
#  vac_schedule <- read_xlsx("../inst/extdata/inputs/vac_schedule_18plus_synth.xlsx", col_types = c("date", rep("numeric", 80))) %>%
#    select(-starts_with("X"))
#  
#  # read in vaccine effectiveness parameters
#  # these will be used when converting the vaccine schedule into weighted vaccine effectiveness, delay to protection, and vaccination rate
#  ve_info <- readRDS("../inst/extdata/inputs/ve_params.rds")
#  
#  # convert vaccination schedule
#  # the vaccination schedule is assumed to be cumulative over time
#  vac_schedule1 <- convert_vac_schedule(
#    vac_schedule = vac_schedule,
#    ve = ve_info[[1]],
#    delay = ve_info[[2]],
#    hosp_multiplier = ve_info[[3]],
#    ve_trans = ve_info[[4]]
#  )

## ----params, echo=TRUE, eval=FALSE--------------------------------------------
#  # read in transition rates
#  transition_rates <- readRDS("../inst/extdata/inputs/transition_rates.rds")
#  
#  # read in transmission matrix
#  april_2017  <- readRDS("../inst/extdata/inputs/contact_matrices/converted/transmission_matrix_april_2017.rds")
#  
#  # define population size (by age group)
#  age_dist <- c(0.10319920, 0.11620856, 0.12740219, 0.12198707, 0.13083463,0.14514332, 0.12092904, 0.08807406, 0.04622194)
#  n <- 17407585 # Dutch population size
#  n_vec <- n * age_dist
#  
#  # create named list of parameter
#  params <- list(beta = 0.0004,
#                 beta1 = 0.14,
#                 gamma = 0.5,
#                 sigma = 0.5,
#                 epsilon = 0.01,
#                 N = n_vec,
#                 h = transition_rates$h,
#                 i1 = transition_rates$i1,
#                 i2 = transition_rates$i2,
#                 d = transition_rates$d,
#                 d_ic = transition_rates$d_ic,
#                 d_hic = transition_rates$d_hic,
#                 r = transition_rates$r,
#                 r_ic = transition_rates$r_ic,
#                 p_report = 1/3,
#                 c_start = april_2017,
#                 keep_cm_fixed = TRUE,
#                 vac_inputs = vac_schedule1,
#                 use_cases = TRUE,
#                 no_vac = FALSE,
#                 t_calendar_start = yday(as.Date("2020-01-01")),
#                 beta_change = NULL
#  )
#  

## ----initial_conditions, echo=TRUE, eval=FALSE, message=FALSE-----------------
#  # Specify initial conditions ---------------------------------
#  empty_state <- c(rep(0, 9)) # vector of zeros
#  
#  init <- c(
#    t = 0,
#    S = c(n_vec[1:4], n_vec[5]-1, n_vec[6:9]),
#    Shold_1d = empty_state,
#    Sv_1d = empty_state,
#    Shold_2d = empty_state,
#    Sv_2d = empty_state,
#    E = empty_state,
#    Ev_1d = empty_state,
#    Ev_2d = empty_state,
#    I = c(rep(0,4),1,rep(0,4)),
#    Iv_1d = empty_state,
#    Iv_2d = empty_state,
#    H = empty_state,
#    Hv_1d = empty_state,
#    Hv_2d = empty_state,
#    H_IC = empty_state,
#    H_ICv_1d = empty_state,
#    H_ICv_2d = empty_state,
#    IC = empty_state,
#    ICv_1d = empty_state,
#    ICv_2d = empty_state,
#    D = empty_state,
#    R = empty_state,
#    Rv_1d = empty_state,
#    Rv_2d = empty_state
#  )
#  
#  # create vector of time points
#  times <- seq(0,nrow(vac_schedule)-1, by = 1)

## ----run_model, echo=TRUE, eval=FALSE-----------------------------------------
#  # Solve model ------------------------------------------------------
#  seir_out <- lsoda(init,                 # initial conditions
#                    times,                # time vector
#                    age_struct_seir_ode,  # model function
#                    params)               # parameters
#  
#  # post-process the model output ------------------------------------
#  seir_out <- as.data.frame(seir_out)
#  out <- postprocess_age_struct_model_output(seir_out)

