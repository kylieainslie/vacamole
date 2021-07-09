# Run vacamole simulation
library(readxl)
library(dplyr)
library(tidyr)
devtools::load_all()
# source('/s-schijf/ainsliek/vacamole/R/run_model.R')

# Read in data probabilities
dons_probs <- read_xlsx("inst/extdata/data/ProbabilitiesDelays_20210107.xlsx")
contact_matrices_all <- read.delim("inst/extdata/data/S2_contact_matrices_withPico3.tsv", header=TRUE, allowEscapes=FALSE, sep="\t")
contact_matrix_april2020 <- contact_matrices_all %>%
  filter(survey == "April 2020") %>%
  filter(contact_type == "community") %>%
  select(-survey, -contact_type) %>%
  pivot_wider(., names_from = cont_age, values_from = m_est)

contact_matrix_input <- as.matrix(contact_matrix_april2020[-1,-1])
rownames(contact_matrix_input) <- contact_matrix_april2020$part_age[-1]

age_dist <- c(0.10319920, 0.11620856, 0.12740219, 0.12198707, 0.13083463, 0.14514332, 0.12092904,
              0.08807406, 0.04622194)
n_age_groups <- length(age_dist)

# define input values
n <- 50000                                  # population size
N <- 17407585                               # Dutch population size
scale_val <- n/N                            # scale for percentage of Dutch population (for vaccine distribution)

# vaccines per day
vaccines_per_day <- read_xlsx("inst/extdata/data/vaccines_per_day_by_type.xlsx", sheet = 1)
vaccines_per_day_scaled <- round(vaccines_per_day[,-1]*scale_val,0)
vac_per_day_mat <- as.matrix(vaccines_per_day_scaled)

n_postcodes <- n/100                        # number of postcodes
t_max <- dim(vaccines_per_day)[1]           # number of time points
# proportion initially infected (from coronadashboard.nl)
infectious_prop_total <- 140833/17407585
prop_inf_by_age <- c(0.018, 0.115, 0.156, 0.118, 0.142, 0.199, 0.114, 0.062, 0.054 + 0.023) * infectious_prop_total           
# proportion initially recovered (seroprevalence - from Don)
# Method: 
#   extrapolated seroprev = [measured seroprev September] * [cumulative hospitalisations today] / [cumulative hospitalisations September]
prop_rec = c(0.01120993, 0.09663659, 0.24141186, 0.11004723, 0.10677859, 0.11977255, 0.11904044, 0.11714503, 0.11347191)
prop_vh = 0.3                               # proportion vaccine hesitant

# susceptibility matrix
sigma_mat = matrix(rep(1,n_age_groups*n_age_groups), nrow = length(age_dist))

# transmission probability matrix
beta_mat = matrix(rep(0.3,length(age_dist)*length(age_dist)), nrow = length(age_dist))

# household information (needs to be improved based on Jantien's data)
hh_info_dat <- data.frame(hh_size = c(1,2,2,3,3,4,4,5,10,4), 
                          hh_type = c("Single", 
                                      "Couple","Single, 1 Child",
                                      "Couple, 1 Child", "Single, 2 Children",
                                      "Couple, 2 Children", "Single, 3 Children",
                                      "Couple, 3 Children",
                                      "Multi-person (care home)", "Multi-person (student)"),
                          hh_prob = c(0.385, 0.283, 0.045, 0.095, 0.022, 0.114, 0.007, 0.045, 0.001, 0.004))

# order of age groups for vaccination
vac_order <- matrix(c(rep(0, 3*7)), nrow = 3)
vac_order[1,] <- c(9,8,7,6,5,4,3)
vac_order[2,] <- c(9,8,7,6,5,4,3)
vac_order[3,] <- c(9,8,7,6,5,4,3) # c(3,4,5,6,7,8,9)

ve_mat <- matrix(c(0.926, 0.524, 0.5269,
                   0.948, 0.924, 0.9), nrow = 3)
rownames(ve_mat) <- names(vaccines_per_day)[-1]


# probabilities
probs_by_age = data.frame(age_group = 1:n_age_groups,
                          prob_high_risk = c(rep(0,9)), # not sure of prop high risk by age, so leaving out for now
                          prob_asympt = c(0.796, 0.796, 0.762, 0.731, 0.695, 0.655, 0.612, 0.551, 0.551), # from Scott and Fumi's paper'
                          prob_severe = dons_probs$`Pr(Se|I)`,
                          prob_hosp = dons_probs$P_infection2admission,
                          prob_non_hosp_death = dons_probs$P_non_hosp_death,
                          prob_ICU = dons_probs$P_admission2IC,
                          prob_hosp_death = dons_probs$P_admission2death,
                          prob_ICU_death = 1-dons_probs$P_IC2hospital)


# tag for output file
tag <- "_beta03_AZ_old_to_young"

# run simulation
# profvis({
system.time(
results <- run_model(pop_size = n,
                           age_dist = age_dist,
                           hh_info = hh_info_dat,
                           n_postcodes = n_postcodes,
                           prop_inf_by_age = prop_inf_by_age,
                           prop_rec_by_age = prop_rec,
                           t_max = t_max,
                           vac_mech = 1,
                           contact_mat = contact_matrix_input,
                           prop_in_pc = 0.85,
                           sigma = sigma_mat,
                           beta = beta_mat,
                           tau = 1,
                           ve = ve_mat,
                           prop_asympt_vec = probs_by_age$prob_asympt,
                           prop_severe = probs_by_age$prob_severe,
                           prop_hosp = probs_by_age$prob_hosp,
                           prop_non_hosp_death = probs_by_age$prob_non_hosp_death,
                           prop_ICU = probs_by_age$prob_ICU,
                           prop_hosp_death = probs_by_age$prob_hosp_death,
                           prop_ICU_death = probs_by_age$prob_ICU_death,
                           vaccines_per_day = vac_per_day_mat,
                           perc_vac_hesitant = prop_vh,
                           perc_high_risk = probs_by_age$prob_high_risk,
                           time_in_hh_isolation = 14,
                           time_btw_doses = 21,
                           time_to_protection = 14,
                           vac_age_group_order = vac_order
  )
)
# })
# save results
saveRDS(results, paste0("model_results", tag, ".rds"))
