### Vacamole 
#' Dynamic agent-based model of infection and vaccination
#'
#' This function initializes the population before running the model.
#' @param pop_size Number of individuals in the population.
#' @param age_dist Age distribution: proportion of individuals in each age group
#' @param hh_info_mat data frame with household distribution information
#' @param n_postcodes number of distict geographic regions to divide population into
#' if 0 then there will be no household structure
#' @param prop_inf_by_age vector of proportion of each age group infected at the start of the simulation
#' @param t_max integer of final time step (days)
#' @param vaccine_mech integer denoting vaccine mechanism:
#'          1 - vaccine reduces susceptibility
#'          2 - vaccine reduces the probability of becoming symptomatic
#'          3 - vaccine protects against severe symptoms
#' @param contact_matrix matrix of daily contacts for each age group (should be of dimension length(age_dist) x length(age_dist))
#' @param prop_in_pc proportion of contacts that are within an individual's postcode
#' @param sigma matrix of susceptibilities for each age group
#' @param beta matrix of age-stratified transmission rates
#' @param tau relative transmissibility of asymptomatic infection
#' @param ve vaccine efficacy
#' @param prop_asymp_vec vector of proportions of asymptomatic infections by age group
#' @param prop_severe
#' @param prop_hosp
#' @param prop_non_hosp_death
#' @param prop_ICU
#' @param prop_hosp_death
#' @param prop_ICU_death
#' @param vaccines_per_day
#' @param perc_vac_hesitant
#' @param perc_high_risk
#' @return List of summary results
#' @keywords vacamole
#' @export
run_model <- function(pop_size,
                      age_dist,
                      hh_info,
                      n_postcodes,
                      prop_inf_by_age,
                      prop_rec_by_age,
                      t_max,
                      vac_mech,
                      contact_mat,
                      prop_in_pc,
                      sigma,
                      beta,
                      tau,
                      ve,
                      prop_asympt_vec,
                      prop_severe,
                      prop_hosp,
                      prop_non_hosp_death,
                      prop_ICU,
                      prop_hosp_death,
                      prop_ICU_death,
                      vaccines_per_day,
                      perc_vac_hesitant,
                      perc_high_risk,
                      time_in_hh_isolation,
                      time_btw_doses,
                      time_to_protection,
                      vac_age_group_order
) {

# create a population with households stratified by age
  # columns of pop_mat:
  # 1. Individual ID
  # 2. Household ID
  # 3. Household size
  # 4. Postal code
  # 5. Age group
  # 6. High-risk status
  # 7. Vaccine Hesitancy
  # 8. State
  
pop_mat = create_pop_w_households(pop_size, 
                                  age_dist, 
                                  hh_info, 
                                  n_postcodes, 
                                  perc_high_risk, 
                                  perc_vac_hesitant)

n = dim(pop_mat)[1]
n_ages = length(age_dist)

# define some empty vectors
time_dose1 <- c(rep(0,n))
time_dose2 <- c(rep(0,n))
vac_status <- c(rep(0,n))
vac_type <- c(rep(0,n))
protection_status <- c(rep(0,n))

# vectors to keep track of states
infection_to_symptoms <- c(rep(0,n))
infectious_period <- c(rep(0,n))
time_of_infection <- c(rep(0,n))
time_of_death <- c(rep(0,n))
household_isolation <- c(rep(0,n))
time_sympt_to_hosp <- c(rep(0,n))
latent_period = round(rweibull(n, 1.97, 9.09),0);
time_hosp_to_discharge = round(rnbinom(n, size = 1.73, mu = 7.90),0);  # from Don
time_hosp_to_death = round(rnbinom(n, size = 1.73, mu = 7.90),0);      # from Don
time_hosp_to_ICU = round(rnbinom(n, size = 0.41, mu = 2.28),0);        # from Don
time_ICU_to_hosp = round(rnbinom(n, size = 1.13, mu = 15.6),0);        # from Don
time_ICU_step_down_care = round(rnbinom(n, 2.73, mu = 10.1),0); # from Don

# output matrices (stratified by age)
empty <- matrix(, nrow = n_ages, ncol = t_max)
rtn <- list(total_infections = empty,
            total_asymptomatic = empty,
            total_symptomatic = empty,
            total_severe = empty,
            total_hospitalisations = empty,
            total_ICU = empty,
            total_recovered = empty,
            total_dead = empty,
            total_vaccinations = empty,
            total_protected = empty,
            new_infections = empty
            )
# change state of initially infected and recovered (previously infected)
state_col <- pop_mat[,"State"]
state_col = determine_intitial_states(pop_mat[,"Age_Group"], 
                                      prop_inf_by_age, 
                                      prop_rec_by_age)

contact_info <- get_total_community_contacts(contact_mat)

# start time loop
for (t in 1:t_max){
cat("Time step:", t, "\n")

# determine who is vaccinated in that time step and change vac status
tmp <- vaccinate(num_vaccines = vaccines_per_day[t,], 
                 ages = pop_mat[,"Age_Group"],
                 high_risk_status = pop_mat[,"High-Risk_Status"], 
                 vaccine_hesitancy = pop_mat[,"Vac_Hesitancy"],
                 vac_status = vac_status,
                 vac_type = vac_type,
                 time_step = t,
                 time_dose1 = time_dose1,
                 time_dose2 = time_dose2,
                 protection_status = protection_status,
                 time_btw_doses = time_btw_doses,
                 time_to_protection = time_to_protection,
                 age_group_vac_order = vac_order)

# update vac status and time of dose 1
vac_status <- tmp$vac_status
vac_type <- tmp$vac_type
protection_status <- tmp$protection_status
time_dose1 <- tmp$time_dose1
time_dose2 <- tmp$time_dose2

# get contacts for only infectious individuals for time step t
infectious_ids <- which(state_col %in% c(2:4))
contacts_list <- get_contacts_2(contact_info, 
                           infectious_ids,
                           pop_mat[,"ID"], 
                           pop_mat[,"Household_ID"], 
                           pop_mat[,"Age_Group"],
                           pop_mat[,"Postcode"], 
                           prop_in_pc, 
                           household_isolation)
# cat("Nodes:", contacts_list$nodes, "\n")
# cat("Edges:", contacts_list$edges, "\n")

contacts <- cbind(contacts_list$nodes, contacts_list$edges)

if(length(which(contacts[,1] == 0)) > 0){
  contacts <- contacts[-which(contacts[,1] == 0),]
}
# start individual loop for transmission process
for (i in 1:n){

if (state_col[i] == 8) {next} # if person i is dead, skip to next person

# get some info about person i
age_i = pop_mat[i, "Age_Group"]
state_i = state_col[i]
vac_status_i = protection_status[i]
vac_type_i = vac_type[i]
hh_id_i = pop_mat[i,"Household_ID"]
time_sympt_to_hosp_i = time_sympt_to_hosp[i]
time_hosp_to_discharge_i = time_hosp_to_discharge[i]
time_hosp_to_death_i = time_hosp_to_death[i]
time_hosp_to_ICU_i = time_hosp_to_ICU[i]
time_ICU_to_hosp_i = time_ICU_to_hosp[i]
time_ICU_step_down_care_i = time_ICU_step_down_care[i]
latent_period_i = latent_period[i]
infectious_period_i = infectious_period[i]

# ve for vaccine type used on person i
ve_i <- c(0,0)
if(vac_type_i > 0){ ve_i = ve[vac_type_i,]}

# determine contacts: household members + edges in contacts
neighbors <- which(pop_mat[seq((i-10):(i+10)),"Household_ID"] == hh_id_i)
hh_i_members <- pop_mat[neighbors, "ID"]
hh_i_members_minus_i <- hh_i_members[!hh_i_members == i]
nodes_i <- which(contacts[,1] == i)
edges_i <- which(contacts[,2] == i)
non_hh_contacts_i = c(contacts[nodes_i, 2], contacts[edges_i, 1])
contacts_i = c(hh_i_members_minus_i, non_hh_contacts_i)
states_c = state_col[contacts_i]

# only make contacts for people with infectious contacts
if (!all(states_c %in% c(0,1,5,6,7,8))){
  ages_c = pop_mat[contacts_i,"Age_Group"]


  # make contacts for indiv
  tmp1 = make_contacts(state_i,
                       age_i,
                       states_c,
                       contacts_i,
                       ages_c,
                       sigma,
                       beta,
                       vac_mech,
                       vac_status_i,
                       ve_i,
                       tau,
                       t)

  # update state and time_of infection
  state_col[i] <- tmp1$state_i
  time_of_infection[i] <- tmp1$time_of_infection_i
  time_of_infection_i <- time_of_infection[i]
  state_i_updated <- state_col[i]
} else {
  state_i_updated <- state_col[i]
  time_of_infection_i <- time_of_infection[i]
}
# change states
if(state_col[i] %in% c(0,7,8)){
  next
} else {
  tmp2 = change_states(i-1,
                     state_i_updated,
                     age_i,
                     time_of_infection_i,
                     latent_period_i,
                     infectious_period_i,
                     vac_mech,
                     vac_status_i,
                     ve_i,
                     t,
                     prop_asympt_vec,
                     prop_severe,
                     prop_hosp,
                     prop_non_hosp_death,
                     prop_hosp_death,
                     prop_ICU,
                     prop_ICU_death,
                     hh_i_members,
                     household_isolation,
                     time_in_hh_isolation,
                     time_sympt_to_hosp_i, 
                     time_hosp_to_discharge_i,
                     time_hosp_to_death_i,
                     time_hosp_to_ICU_i,
                     time_ICU_to_hosp_i,
                     time_ICU_step_down_care_i)

  # update state and hh_isolation and infectious period
  state_col[i] <- tmp2$state_i
  household_isolation <- tmp2$household_isolation
  infectious_period[i] <- tmp2$infectious_period_i
  time_sympt_to_hosp[i] <- tmp2$time_sympt_to_hosp_i
  }
} # end loop over individuals

# output  
num_infections_t <- length(which(state_col == 1))
num_infectious_t <- length(which(state_col == 2)) +  length(which(state_col == 3))

# loop through age groups and report total counts by age group
for (a in unique(pop_mat[,"Age_Group"])){
  age_a_pos <- which(pop_mat[,"Age_Group"] == a)
  state_col_age_a = state_col[age_a_pos]                 # get states for people in age group a
  time_of_infection_age_a = time_of_infection[age_a_pos] # get new infections for people in age group a
  vac_status_age_a = vac_status[age_a_pos]               # get states for people in age group a
  protection_status_age_a = protection_status[age_a_pos] # get states for people in age group a

  rtn$total_infections[a, t] = length(which(state_col_age_a == 1)) 
  rtn$total_asymptomatic[a, t] = length(which(state_col_age_a == 2))
  rtn$total_symptomatic[a, t] = length(which(state_col_age_a == 3))
  rtn$total_severe[a, t] = length(which(state_col_age_a == 4)) 
  rtn$total_hospitalisations[a, t] = length(which(state_col_age_a == 5)) 
  rtn$total_ICU[a, t] = length(which(state_col_age_a == 6))
  rtn$total_recovered[a, t] = length(which(state_col_age_a == 7)) 
  rtn$total_dead[a, t] = length(which(state_col_age_a == 8))
  rtn$total_vaccinations[a, t] = length(which(vac_status_age_a != 0))
  rtn$total_protected[a, t] = length(which(protection_status_age_a != 0)) 
  rtn$new_infections[a, t] =length(which(time_of_infection_age_a == t)) 
}

} # end loop over time steps  

return(rtn)
}
