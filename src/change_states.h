#include <Rcpp.h>
using namespace Rcpp;

#ifndef CHANGE_STATES_H
#define CHANGE_STATES_H

List change_states(int id,
                   int state_i,
                   int age_i,
                   int time_of_infection_i,
                   int latent_period_i,
                   int infectious_period_i,
                   int vac_mech,
                   int vac_status_i,
                   NumericVector ve,
                   int time_step,
                   NumericVector prop_asympt_vec,
                   NumericVector prop_severe,
                   NumericVector prop_hosp,
                   NumericVector prop_non_hosp_death,
                   NumericVector prop_hosp_death,
                   NumericVector prop_ICU,
                   NumericVector prop_ICU_death,
                   NumericVector hh_i_members,
                   NumericVector household_isolation,
                   int time_in_hh_isolation,
                   int time_sympt_to_hosp_i, 
                   int time_hosp_to_discharge_i,
                   int time_hosp_to_death_i,
                   int time_hosp_to_ICU_i,
                   int time_ICU_to_hosp_i,
                   int time_ICU_step_down_care_i
                  );
  

#endif