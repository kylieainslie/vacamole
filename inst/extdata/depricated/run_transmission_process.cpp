#include <Rcpp.h>
#include "create_pop_w_households.h"
#include "determine_initial_states.h"
#include "get_total_community_contacts.h"
#include "get_contacts_2.h"
#include "vaccinate.h"
#include "helpers.h"
#include "make_contacts.h"
#include "change_states.h"

using namespace Rcpp;

// Main transmission function (currently just an outline now)
//
// @param pop_size Number of individuals in the population.
// @param age_dist Age distribution: proportion of individuals in each age group
// @param hh_info_mat data frame with household distribution information
// @param n_postcodes number of distict geographic regions to divide population into
//        if 0 then there will be no household structure
// @param prop_inf_by_age vector of proportion of each age group infected at the start of the simulation
// @param t_max integer of final time step (days)
// @param vaccine_mech integer denoting vaccine mechanism:
//          1 - vaccine reduces susceptibility
//          2 - vaccine reduces the probability of becoming symptomatic
//          3 - vaccine protects against severe symptoms
// @param contact_matrix matrix of daily contacts for each age group (should be of dimension length(age_dist) x length(age_dist))
// @param prop_in_pc proportion of contacts that are within an individual's postcode
// @param sigma matrix of susceptibilities for each age group
// @param beta matrix of age-stratified transmission rates
// @param tau relative transmissibility of asymptomatic infection
// @param ve vaccine efficacy
// @param prop_asymp_vec vector of proportions of asymptomatic infections by age group
// @param prop_severe
// @param prop_hosp
// @param prop_non_hosp_death
// @param prop_ICU
// @param prop_hosp_death
// @param prop_ICU_death
// @param vaccines_per_day
// @param perc_vac_hesitant
// @param perc_high_risk
// @return vectors of XX
// @keywords vac-a-mole
// @export

// [[Rcpp::export]]
List run_transmission_process(
         int pop_size,
         NumericVector age_dist,
         DataFrame hh_info,
         int n_postcodes,
         NumericVector prop_inf_by_age,
         NumericVector prop_rec_by_age,
         int t_max,
         int vac_mech,
         NumericMatrix contact_mat,
         double prop_in_pc,
         NumericMatrix sigma,
         NumericMatrix beta,
         double tau,
         NumericMatrix ve,
         NumericVector prop_asympt_vec,
         NumericVector prop_severe,
         NumericVector prop_hosp,
         NumericVector prop_non_hosp_death,
         NumericVector prop_ICU,
         NumericVector prop_hosp_death,
         NumericVector prop_ICU_death,
         NumericMatrix vaccines_per_day,
         double perc_vac_hesitant,
         NumericVector perc_high_risk,
         int time_in_hh_isolation,
         int time_btw_doses,
         int time_to_protection,
         NumericMatrix vac_age_group_order
         ) {
  
  // create a population with households stratified by age
  NumericMatrix pop_mat = create_pop_w_households(pop_size, 
                                                  age_dist, 
                                                  hh_info, 
                                                  n_postcodes, 
                                                  perc_high_risk, 
                                                  perc_vac_hesitant);
  // columns of pop_mat:
  //  0. Individual ID
  //  1. Household ID
  //  2. Household size
  //  3. Postal code
  //  4. Age group
  //  5. High-risk status
  //  6. Vaccine Hesitancy
  //  7. State
  
  int n = pop_mat.nrow();
  int n_ages = age_dist.length();
  //Rcpp::Rcout << "N: "<< n << std::endl;
  // define some vectors
  NumericVector prop_infs = prop_inf_by_age;
  NumericVector prop_recs = prop_rec_by_age;
  NumericVector age_groups = pop_mat(_, 4);
  NumericVector hh_ids = pop_mat(_,1);
  NumericVector pop_ids = pop_mat(_,0);
  NumericVector high_risk_status = pop_mat(_,5);
  NumericVector vac_hesitant = pop_mat(_,6);
  NumericVector state_col = pop_mat(_,7);
  NumericVector time_dose1(n);
  NumericVector time_dose2(n);
  NumericMatrix vac_order = vac_age_group_order;
  NumericVector vac_status(n);
  NumericVector vac_type(n);
  NumericVector protection_status(n);
  
  // vectors to keep track of states
  NumericVector latent_period = round(rweibull(n, 1.97, 9.09),0);
  NumericVector infection_to_symptoms(n);
  NumericVector infectious_period(n); // empty for now, but will be filled in later
  NumericVector time_of_infection(n);
  NumericVector time_of_death(n);
  NumericVector household_isolation(n);
  NumericVector time_sympt_to_hosp(n);
  NumericVector time_hosp_to_discharge = round(rnbinom_mu(n, 1.73, 7.90),0);  // from Don
  NumericVector time_hosp_to_death = round(rnbinom_mu(n, 1.73, 7.90),0);      // from Don
  NumericVector time_hosp_to_ICU = round(rnbinom_mu(n, 0.41, 2.28),0);        // from Don
  NumericVector time_ICU_to_hosp = round(rnbinom_mu(n, 1.13, 15.6),0);        // from Don
  NumericVector time_ICU_step_down_care = round(rnbinom_mu(n, 2.73, 10.1),0); // from Don
    
  // output matrices (stratified by age)
  NumericMatrix empty(n_ages, t_max);
  NumericVector R_vec(t_max);
  
  // create list of output matrices (all are stratified by age)
  List rtn = List::create(Named("total_infections") = clone(empty), 
                          _["total_asymptomatic"] = clone(empty),
                          _["total_symptomatic"] = clone(empty),
                          _["total_severe"] = clone(empty),
                          _["total_hospitalisations"] = clone(empty),
                          _["total_ICU"] = clone(empty),
                          _["total_recovered"] = clone(empty),
                          _["total_dead"] = clone(empty),
                          _["total_vaccinations"] = clone(empty),
                          _["total_protected"] = clone(empty),
                          _["R"] = R_vec,
                          _["new_infections"] = clone(empty) //, 
                          // _["new_asymptomatic"] = clone(empty),
                          // _["new_symptomatic"] = clone(empty),
                          // _["new_severe"] = clone(empty),
                          // _["new_hospitalisations"] = clone(empty),
                          // _["new_ICU"] = clone(empty),
                          // _["new_recovered"] = clone(empty),
                          // _["new_dead"] = clone(empty),
                          // _["new_vaccinations"] = clone(empty)
                          );
  
  // change state of initially infected and recovered (previously infected)
  state_col = determine_intitial_states(age_groups, 
                                        prop_infs, 
                                        prop_recs);
  //Rcpp::Rcout << "Number initially infected: "<< which_cpp(state_col, 3, "in").length() << std::endl;
  List contact_info = get_total_community_contacts(as<NumericMatrix>(contact_mat));
  
  // start time loop
  for (int t = 0; t < t_max; ++t){
    Rcpp::Rcout << "Time step: " << t << std::endl;
    
    // determine who is vaccinated in that time step and change vac status
    //Rcpp::Rcout << "Vaccinate " << std::endl;
    List tmp = vaccinate(vaccines_per_day(t,_), 
                         age_groups, 
                         high_risk_status, 
                         vac_hesitant,
                         vac_status,
                         vac_type,
                         t,
                         time_dose1,
                         time_dose2,
                         protection_status,
                         time_btw_doses,
                         time_to_protection,
                         vac_order);
    
    // update vac status and time of dose 1
    vac_status = tmp["vac_status"];
    vac_type = tmp["vac_type"];
    protection_status = tmp["protection_status"];
    time_dose1 = tmp["time_dose1"];
    time_dose2 = tmp["time_dose2"];
    
    //Rcpp::Rcout << "Vaccinated w/ 1st dose at time t: " << which_cpp(vac_status, 1, "in").length() << std::endl;
    //Rcpp::Rcout << "Vaccinated w/ 2nd dose at time t: " << which_cpp(vac_status, 2, "in").length() << std::endl;
    
    // get contacts for all individuals for time step t
    //Rcpp::Rcout << "Get contacts" << std::endl;
    NumericMatrix contacts = get_contacts_2(contact_info, 
                                            pop_mat(_,0), 
                                            pop_mat(_,1), 
                                            age_groups,
                                            pop_mat(_,3), 
                                            prop_in_pc, 
                                            household_isolation);
    
    //start individual loop for transmission process
    for (int i = 0; i < n; ++i){
      
      //Rcpp::Rcout << "Individual: " << i+1 << std::endl;
      
      if (state_col[i] == 8) {continue;} // if person i is dead, skip to next person
      
      // get some info about person i
      int age_i = pop_mat(i, 4);
      int state_i = state_col[i];
      int vac_status_i = protection_status[i];
      int vac_type_i = vac_type[i];
      int hh_id_i = hh_ids[i];
      int time_sympt_to_hosp_i = time_sympt_to_hosp[i];
      int time_hosp_to_discharge_i = time_hosp_to_discharge[i];
      int time_hosp_to_death_i = time_hosp_to_death[i];
      int time_hosp_to_ICU_i = time_hosp_to_ICU[i];
      int time_ICU_to_hosp_i = time_ICU_to_hosp[i];
      int time_ICU_step_down_care_i = time_ICU_step_down_care[i];
      int latent_period_i = latent_period[i];
      int infectious_period_i = infectious_period[i];
      
      NumericVector ve_i(2);
      if(vac_type_i > 0){ ve_i = ve(vac_type_i - 1,_);}
      // determine contacts: household members + edges in contacts
      NumericVector hh_i_members = in_group(hh_ids, pop_ids, hh_id_i);
      IntegerVector pos_contacts_i = which_cpp(contacts(_,1), i+1, "in");
      NumericVector non_hh_contacts_i = subset_vec(pos_contacts_i, contacts(_,1));
      IntegerVector contacts_i = union_(as<IntegerVector>(hh_i_members), as<IntegerVector>(non_hh_contacts_i)) - 1;
      //Rcpp::Rcout << "Contacts i: " << contacts_i << std::endl;
      NumericVector states_c = subset_vec(contacts_i, state_col);
      
      
      //Rcpp::Rcout << "Make contacts " << std::endl;
      
     // make contacts for indiv
      List tmp1 = make_contacts(state_i,
                               age_i,
                               states_c,
                               contacts_i,
                               age_groups,
                               sigma,
                               beta,
                               vac_mech,
                               vac_status_i,
                               ve_i,
                               tau,
                               t);
      
      //Rcpp::Rcout << "End make contacts " << std::endl;
      
     // update state and time_of infection
      state_col[i] = tmp1["state_i"];
      time_of_infection[i] = tmp1["time_of_infection_i"];
      int time_of_infection_i = time_of_infection[i];
      int state_i_updated = state_col[i];

      
      //Rcpp::Rcout << "Change states " << std::endl;
      
     // change states
       List tmp2 = change_states(i,
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
                                 time_ICU_step_down_care_i);
        
        //Rcpp::Rcout << "End change states " << std::endl;
        
        // update state and hh_isolation and infectious period
        state_col[i] = tmp2["state_i"];
        household_isolation = tmp2["household_isolation"];
        infectious_period[i] = tmp2["infectious_period_i"];
        time_sympt_to_hosp[i] = tmp2["time_sympt_to_hosp_i"]; 
       
        // output counts of each state
        // if (vac_status[i] == 1){num_vaccinations(age_i-1, t) += 1;}
        // if (state_col[i] == 1) {num_infections(age_i-1, t) += 1;}
        // if (state_col[i] == 2) {num_asymptomatic(age_i-1, t) += 1;}
        // if (state_col[i] == 3 | state_col[i] == 4) {num_symptomatic(age_i-1, t) += 1;}
        // if (state_col[i] == 4) {num_severe(age_i-1, t) += 1;}
        // if (state_col[i] == 5) {num_hospitalised(age_i-1, t) += 1;}
        // if (state_col[i] == 6) {num_ICU(age_i-1, t) += 1;}
        // if (state_col[i] == 7) {num_recovered(age_i-1, t) += 1;}
        // if (state_col[i] == 8) {num_dead(age_i-1, t) += 1;}
        // 
    } // end loop over individuals
    
    //Rcpp::Rcout << "time of infection " << time_of_infection << std::endl;
    
  // output  
      int num_infections_t = 0;
    
      if (t == 0) {
        num_infections_t = which_cpp(state_col, 1, "in").length();
      } else { num_infections_t = which_cpp(time_of_infection, t, "in").length(); }
      
      int num_infectious_t = which_cpp(state_col, 2, "in").length() +  which_cpp(state_col, 3, "in").length();
      //Rcpp::Rcout << "n infections at time t: " << num_infections_t << std::endl;
      //Rcpp::Rcout << "n infectious at time t: " << num_infectious_t << std::endl;
      
      as<NumericVector>(rtn["R"])[t] = num_infections_t / num_infectious_t; 
      
    // loop through age groups and report total counts by age group
    for (int a = 0; a < n_ages; ++a){
      IntegerVector age_a_pos = which_cpp(age_groups, a + 1, "in"); // get people in age group a
      NumericVector state_col_age_a = subset_vec(age_a_pos, state_col); // get states for people in age group a
      NumericVector time_of_infection_age_a = subset_vec(age_a_pos, time_of_infection); // get new infections for people in age group a
      NumericVector vac_status_age_a = subset_vec(age_a_pos, vac_status); // get states for people in age group a
      NumericVector protection_status_age_a = subset_vec(age_a_pos, protection_status); // get states for people in age group a
      
      as<NumericMatrix>(rtn["total_infections"])(a, t) = which_cpp(state_col_age_a, 1, "in").length(); 
      as<NumericMatrix>(rtn["total_asymptomatic"])(a, t) = which_cpp(state_col_age_a, 2, "in").length();
      as<NumericMatrix>(rtn["total_symptomatic"])(a, t) = which_cpp(state_col_age_a, 3, "in").length();
      as<NumericMatrix>(rtn["total_severe"])(a, t) = which_cpp(state_col_age_a, 4, "in").length(); 
      as<NumericMatrix>(rtn["total_hospitalisations"])(a, t) = which_cpp(state_col_age_a, 5, "in").length(); 
      as<NumericMatrix>(rtn["total_ICU"])(a, t) = which_cpp(state_col_age_a, 6, "in").length(); 
      as<NumericMatrix>(rtn["total_recovered"])(a, t) = which_cpp(state_col_age_a, 7, "in").length(); 
      as<NumericMatrix>(rtn["total_dead"])(a, t) = which_cpp(state_col_age_a, 8, "in").length(); 
      as<NumericMatrix>(rtn["total_vaccinations"])(a, t) = which_cpp(vac_status_age_a, 0, "notin").length();
      as<NumericMatrix>(rtn["total_protected"])(a, t) = which_cpp(protection_status_age_a, 0, "notin").length(); 
      as<NumericMatrix>(rtn["new_infections"])(a, t) = which_cpp(time_of_infection_age_a, t, "in").length(); 
    }
     
  } // end loop over time steps  

  
  
  return rtn;
}


/*** R
n <- 500
age_dist <- c(0.10319920, 0.11620856, 0.12740219, 0.12198707, 0.13083463, 0.14514332, 0.12092904,
              0.08807406, 0.04622194)
hh_info_dat <- data.frame(hh_size = c(1,2,2,3,3,4,4,5,10,4),
                          hh_type = c("Single",
                                      "Couple","Single, 1 Child",
                                      "Couple, 1 Child", "Single, 2 Children",
                                      "Couple, 2 Children", "Single, 3 Children",
                                      "Couple, 3 Children",
                                      "Multi-person (care home)", "Multi-person (student)"),
                          hh_prob = c(0.385, 0.283, 0.045, 0.095, 0.022, 0.114, 0.007, 0.045, 0.001, 0.004))

contact_matrices_all <- read.delim("../inst/extdata/data/S2_contact_matrices_withPico3.tsv")
contact_matrix_april2020 <- contact_matrices_all %>%
  filter(survey == "April 2020") %>%
  filter(contact_type == "community") %>%
  select(-survey, -contact_type) %>%
  pivot_wider(., names_from = cont_age, values_from = m_est)

contact_matrix_input <- as.matrix(contact_matrix_april2020[-1,-1])
rownames(contact_matrix_input) <- contact_matrix_april2020$part_age[-1]

prop_inf_vec <- c(rep(0.2/9, 9))
n_postcodes <- floor(0.02 * n)
t_max <- 5
sigma_mat = matrix(rep(1,length(age_dist)*length(age_dist)), nrow = length(age_dist))
beta_mat = matrix(rep(0.2,length(age_dist)*length(age_dist)), nrow = length(age_dist))
probs_by_age = data.frame(age_group = 1:9,
                          prob_high_risk = c(0.01, 0, 0, 0.1, 0.2, 0.2, 0.3, 0.4, 0.5),
                          prob_asympt = runif(9,0.45, 0.65),
                          prob_severe = c(0.01, 0, 0, 0.1, 0.2, 0.2, 0.3, 0.4, 0.5),
                          prob_hosp = c(0.01, 0, 0, 0.1, 0.2, 0.2, 0.3, 0.4, 0.5),
                          prob_non_hosp_death = c(rep(0,5), rep(0.1, 3), 0.2),
                          prob_ICU = c(rep(0,5), 0.05, 0.1, 0.15, 0.2),
                          prob_hosp_death = c(rep(0.01, 5), 0.1, 0.1, 0.2, 0.3),
                          prob_ICU_death = c(rep(0.0, 5), 0.1, 0.2, 0.3, 0.5),
                          prob_rec = c(rep(0.1/9, 9))
                          )

vaccines_per_day = matrix(c(rep(2, 3*t_max)), nrow = t_max)
vac_order <- matrix(c(rep(0, 3*7)), nrow = 3)
vac_order[1,] <- c(9,8,7,6,5,4,3)
vac_order[2,] <- c(9,8,7,6,5,4,3)
vac_order[3,] <- c(3,4,5,6,7,8,9)

ve_mat <- matrix(c(0.524, 0.524, 0,
                   0.924, 0.924, 0.9), nrow = 3)
rownames(ve_mat) <- names(vaccines_per_day)[-1]


results <-
run_transmission_process(pop_size = n,
     age_dist = age_dist,
     hh_info = hh_info_dat,
     n_postcodes = n_postcodes,
     prop_inf_by_age = prop_inf_vec,
     prop_rec_by_age = probs_by_age$prob_rec,
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
     vaccines_per_day = vaccines_per_day,
     perc_vac_hesitant = 0.4,
     perc_high_risk = probs_by_age$prob_high_risk,
     time_in_hh_isolation = 14,
     time_btw_doses = 21,
     time_to_protection = 14,
     vac_age_group_order = vac_order
     )
*/
