#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
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
                  ) {
  
  //Rcpp::Rcout << "state at beginning of change_states(): " << state_i << std::endl;
  
  NumericVector sympt_to_hosp_dist_mu = NumericVector::create(2.29, 5.51, 5.06, 5.66,  // from Don
                                                              6.55, 5.88, 5.69, 5.09,
                                                              4.33);
  
  // generate random numbers
  NumericVector rdm = runif(7);
  
  NumericVector hh_i_member_vec = hh_i_members;
  
  // from exposed to infectious
  // there is probability of symptomatic or asymptomatic
  if(state_i == 1){
    if(time_step - latent_period_i == time_of_infection_i){
      // if person i is vaccinated and vaccination prevents symptoms
      if ((vac_mech == 2) & (vac_status_i > 0)){
        if ((vac_status_i == 1) & (rdm[0] < 1 - (prop_asympt_vec[age_i - 1] * (1 - ve[0])))){
          state_i = 2; // asymptomatic
          infectious_period_i = R::rgamma(2.1,4);
        } else if ((vac_status_i == 2) & (rdm[0] < 1 - (prop_asympt_vec[age_i - 1] * (1 - ve[1])))){
          state_i = 2; // asymptomatic
          infectious_period_i = R::rgamma(2.1,4);
        } else {
          state_i = 3; // symptomatic
          infectious_period_i = R::rgamma(2.9,4);
          
          // person i and their household members are in isolation
          household_isolation[id] = 1;
          for (int h = 0; h < hh_i_member_vec.length(); ++h){
            household_isolation[hh_i_member_vec[h] - 1] = 1;
          }
        }
      // if person i is not vaccinated or vaccination does not prevent symptoms  
    } else {
        if( rdm[1] < prop_asympt_vec[age_i-1]){
          state_i = 2; // asymptomatic
          infectious_period_i = R::fround(R::rgamma(2.1,4),0);
          //Rcpp::Rcout << "infectious period: " << infectious_period_i << std::endl;
        } else {
          state_i = 3; // symptomatic
          infectious_period_i = R::fround(R::rgamma(2.9,4),0);
          time_sympt_to_hosp_i = R::fround(Rf_rnbinom_mu(1.77, sympt_to_hosp_dist_mu[age_i - 1]),0);
          //Rcpp::Rcout << "infectious period: " << infectious_period_i << std::endl;
          // person i and their household members are in isolation
          household_isolation[id] = 1;
          for (int h = 0; h < hh_i_member_vec.length(); ++h){
            household_isolation[hh_i_member_vec[h] - 1] = 1;
          }
        }
      }
    } //else{state_i = 1;}
  }
  // remove household isolation after time_in_hh_isolation days
  if (household_isolation[id] == 1){
    if(time_step == time_of_infection_i + latent_period_i + time_in_hh_isolation){
      household_isolation[id] = 0;
        for (int h = 0; h < hh_i_member_vec.length(); ++h){
          household_isolation[hh_i_member_vec[h] - 1] = 0;
        }
    } else {household_isolation[id] = 1;}
  }
  
  // from asymptomatic to recovered
  if (state_i == 2){
    if(time_step == time_of_infection_i + latent_period_i + infectious_period_i){
    state_i = 7; // recovered
    } else {state_i = 2;}
  }
  
  // from symptomatic to severe
  if (state_i == 3){
    // if person i is not vaccinated or vaccination does not prevent severe disease
    if((vac_mech != 3) | (vac_status_i == 0)){
      // person may have severe disease
      if ( rdm[2] < prop_severe[age_i-1] ){
        state_i = 4; // severe disease
      } else if (time_step == time_of_infection_i + latent_period_i + infectious_period_i){
        state_i = 7; // if not severe and end of infectious period -> recovered
      }
    } else{
      // if person i is vaccinated and vaccination does prevent severe disease
      // person may have severe disease
      if((vac_status_i == 1) & (rdm[2] < (1 - ve[0]) * prop_severe[age_i-1])){
        state_i = 4; // severe disease
      } else if ((vac_status_i == 2) & (rdm[2] < (1 - ve[1]) * prop_severe[age_i-1])){
        state_i = 4; // severe disease
      } else if (time_step == time_of_infection_i + latent_period_i + infectious_period_i){
        state_i = 7; // if not severe and end of infectious period -> recovered
      }
    }
  }
  
  // from severe to hospital OR recovery OR death
  if ((state_i == 3) | (state_i == 4)){
    if( (rdm[3] < prop_hosp[age_i - 1]) & (time_step == time_of_infection_i + latent_period_i + time_sympt_to_hosp_i)){
      state_i = 5; // hospital
    } else if ( (rdm[3] < prop_non_hosp_death[age_i - 1]) & (time_step == time_of_infection_i + latent_period_i + infectious_period_i)){
      state_i = 8; // dead
    } else if ( time_step == time_of_infection_i + latent_period_i + infectious_period_i){
      state_i = 7; // recovered
    }
  }
  
  // from hospital to ICU OR recovery OR death
  if (state_i == 5){
    if( (rdm[4] < prop_ICU[age_i - 1]) & (time_step == time_of_infection_i + latent_period_i + time_sympt_to_hosp_i +
        time_hosp_to_ICU_i)){
      state_i = 6; // ICU
    } else if ( (rdm[4] < prop_hosp_death[age_i - 1]) & (time_step == time_of_infection_i + latent_period_i + infectious_period_i +
      time_sympt_to_hosp_i + time_hosp_to_death_i)){
      state_i = 8; // dead
    } else if ( time_step == time_of_infection_i + latent_period_i + infectious_period_i + time_sympt_to_hosp_i +
      time_hosp_to_discharge_i){
      state_i = 7; // recovered
    }
  }
  
  // from ICU to recovery OR death
  if (state_i == 6){
    if ( (rdm[5] < prop_ICU_death[age_i - 1]) ){
      state_i = 8; // dead
    } else if ( time_step == time_of_infection_i + latent_period_i + infectious_period_i + time_sympt_to_hosp_i +
              time_hosp_to_ICU_i + time_ICU_to_hosp_i + time_ICU_step_down_care_i){
      state_i = 7; // recovered
    }
  }
  
  //Rcpp::Rcout << "state at end of change_states(): " << state_i << std::endl;
  
  // output
  List rtn = List::create(Named("state_i") = state_i,
                          Named("infectious_period_i") = infectious_period_i,
                          Named("time_sympt_to_hosp_i") = time_sympt_to_hosp_i,
                          Named("household_isolation") = household_isolation);
  
  return rtn;
}

/*** R
n <- 1
i <- 1
ages <- sample(1:9, n, replace = TRUE)
state_i <- 1 
vac_status <- sample(0:1, n, replace = TRUE)
household_isolation <- c(rep(0,4))
time_of_infection_i <- 1
latent_period_i <- 5
infectious_period_i <- 0
time_sympt_to_hosp <- round(rweibull(n, 0.845, 5.506),0)
time_hosp_to_discharge <- round(rlnorm(n, 2.724, 0.720),0)
time_hosp_to_death <- round(rlnorm(n, 2.573, 0.842))
time_hosp_to_ICU <- round(rlnorm(n, 2.724, 0.720)/2,0)
time_ICU_to_hosp <- round(rlnorm(n, 2.079, 0.919),0)
time_ICU_step_down_care <- round(rlnorm(n, 2.724, 0.720)/2,0)

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

for (t in 1:20){
cat("Time step: ",t, "\n")

out <-
change_states(id = 1,
              state_i = state_i,
              age_i = ages[i],
              time_of_infection_i = time_of_infection_i,
              latent_period_i = latent_period_i,
              infectious_period_i = infectious_period_i,
              vac_mech = 1,
              vac_status_i = vac_status[i],
              ve = c(0.524, 0.928),
              time_step = t,
              prop_asympt_vec = probs_by_age$prob_asympt,
              prop_severe = probs_by_age$prob_severe,
              prop_hosp = probs_by_age$prob_hosp,
              prop_non_hosp_death = probs_by_age$prob_non_hosp_death,
              prop_hosp_death = probs_by_age$prob_hosp_death,
              prop_ICU = probs_by_age$prob_ICU,
              prop_ICU_death = probs_by_age$prob_ICU_death,
              hh_i_members = c(1,2),
              household_isolation = household_isolation,
              time_in_hh_isolation = 10,
              time_sympt_to_hosp_i = time_sympt_to_hosp[i], 
              time_hosp_to_discharge_i = time_hosp_to_discharge[i],
              time_hosp_to_death_i = time_hosp_to_death[i],
              time_hosp_to_ICU_i = time_hosp_to_ICU[i],
              time_ICU_to_hosp_i = time_ICU_to_hosp[i],
              time_ICU_step_down_care_i = time_ICU_step_down_care[i]
              )

state_i <- out$state_i
cat("State: ",state_i, "\n")

infectious_period_i <- out$infectious_period_i
cat("Infectious period: ", infectious_period_i, "\n")

time_sympt_to_hosp_i <- out$time_sympt_to_hosp_i
cat("Time from symptoms to hospital: ", time_sympt_to_hosp_i, "\n")

household_isolation <- out$household_isolation
cat("Household isolation: ", household_isolation, "\n")

}

*/
