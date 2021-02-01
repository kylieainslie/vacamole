#include <Rcpp.h>
using namespace Rcpp;

// Make contacts
// @param state_i state of person i
// @param age_i age of person i
// @param states_c states of contacts of person i
// @param age_vec vector of ages
// @param sigma susceptibility matrix
// @param beta transmission matrix
// @param vac_mech vaccine mechanism
// @param vac_status_i vaccination status of person i
// @param ve vaccine effectiveness
// @param tau relative infectiousness of asymptomatic person compared to symptomatic
// @return List containing the state of person i and the time of infection 
// @keywords vac-a-mole
// @export

// [[Rcpp::export]]
List make_contacts(int state_i,
                  int age_i,
                  NumericVector states_c,
                  IntegerVector contacts_i,
                  NumericVector age_vec,
                  NumericVector sigma,
                  NumericVector beta,
                  int vac_mech,
                  int vac_status_i,
                  NumericVector ve,
                  double tau,
                  int time_step
                 ) {
  
  // initialise time_of_infection
  int time_of_infection = 0;
  
  // loop over contacts
  for (int c = 0; c < contacts_i.length(); ++c){
    
    // if an individual is in hospital, we assume no contacts 
    if ((state_i == 5) | (state_i == 6)){break;}
    
    // generate random numbers
    NumericVector rdm3 = runif(4);
    
    // info about contact
    int age_c = age_vec[c];
    int contact_state = states_c[c];
    
    // if contact is not infectious, continue to next contact
    if((contact_state == 0) | (contact_state == 1) | (contact_state == 7) | (contact_state == 8)){continue;}
    
    // if susceptible: probability of getting infected from a SYMPTOMATIC infectious contact
    if((state_i == 0) & (vac_status_i == 0) & (contact_state == 3)){
      //Rcpp::Rcout << "Random number "<< rdm3[0] << std::endl;
      //Rcpp::Rcout << "Prob of infection "<< sigma(age_i-1, age_c-1) * beta(age_i-1, age_c-1) << std::endl;
      
      if(rdm3[0] < sigma(age_i-1, age_c-1) * beta(age_i-1, age_c-1)){
        state_i = 1;
        time_of_infection = time_step;
        break;
      }
    }
    // if vaccinated susceptible probability of getting infected from a SYMPTOMATIC infectious contact
    if((state_i == 0) & (vac_status_i > 0) & (contact_state == 3)){
      if(vac_mech == 1){
        //Rcpp::Rcout << "Random number "<< rdm3[1] << std::endl;
        //Rcpp::Rcout << "Prob of infection "<< (1 - ve) * sigma(age_i-1, age_c-1) * beta(age_i-1, age_c-1) << std::endl;
        
        if((vac_status_i == 1) & (rdm3[1] < (1 - ve[0]) * sigma(age_i-1, age_c-1) * beta(age_i-1, age_c-1))){
          state_i = 1;
          time_of_infection = time_step;
          break;
        } else if ((vac_status_i == 2) & (rdm3[1] < (1 - ve[1]) * sigma(age_i-1, age_c-1) * beta(age_i-1, age_c-1))) {
          state_i = 1;
          time_of_infection = time_step;
          break;
        }
      } else{
        //Rcpp::Rcout << "Random number "<< rdm3[1] << std::endl;
        //Rcpp::Rcout << "Prob of infection "<< sigma(age_i-1, age_c-1) * beta(age_i-1, age_c-1) << std::endl;
        
        if(rdm3[1] < sigma(age_i-1, age_c-1) * beta(age_i-1, age_c-1)){
          state_i = 1;
          time_of_infection = time_step;
          break;
        }
      }
    }
    // if susceptible: probability of getting infected from a ASYMPTOMATIC infectious contact
    if((state_i == 0) & (vac_status_i == 0) & (contact_state == 2)){
      //Rcpp::Rcout << "Random number "<< rdm3[2] << std::endl;
      //Rcpp::Rcout << "Prob of infection "<< sigma(age_i-1, age_c-1) * tau * beta(age_i-1, age_c-1) << std::endl;
      
      if(rdm3[2] < sigma(age_i-1, age_c-1) * tau * beta(age_i-1, age_c-1)){
        state_i = 1;
        time_of_infection = time_step;
        break;
      }
    }
    // if vaccinated susceptible probability of getting infected from a ASYMPTOMATIC infectious contact
    if((state_i == 0) & (vac_status_i > 0) & (contact_state == 2)){
      if(vac_mech == 1){
        //Rcpp::Rcout << "Random number: "<< rdm3[3] << std::endl;
        //Rcpp::Rcout << "Prob of infection: "<< (1 - ve) * sigma(age_i-1, age_c-1) * tau * beta(age_i-1, age_c-1) << std::endl;
        if((vac_status_i == 1) & (rdm3[3] < (1 - ve[0]) * sigma(age_i-1, age_c-1) * tau * beta(age_i-1, age_c-1))){
          state_i = 1;
          time_of_infection = time_step;
          break;
        } else if((vac_status_i == 2) & (rdm3[3] < (1 - ve[1]) * sigma(age_i-1, age_c-1) * tau * beta(age_i-1, age_c-1))){
          state_i = 1;
          time_of_infection = time_step;
          break;
        }
      } else{
        //Rcpp::Rcout << "Random number: "<< rdm3[3] << std::endl;
        //Rcpp::Rcout << "Prob of infection: "<< sigma(age_i-1, age_c-1) * tau * beta(age_i-1, age_c-1) << std::endl;
        if(rdm3[3] < sigma(age_i-1, age_c-1) * tau * beta(age_i-1, age_c-1)){
          state_i = 1;
          time_of_infection = time_step;
          break;
        }
      }
    }
    continue; // continue to next contact if non of the conditions above are met
  } // end loop over individal i's contacts  
  // output 
  List rtn = List::create(Named("state_i") = state_i,
                          Named("time_of_infection_i") = time_of_infection);
  
  return rtn;
  
}


/*** R
# n <- 100
# age_dist <- c(0.10319920, 0.11620856, 0.12740219, 0.12198707, 0.13083463, 0.14514332, 0.12092904,
#               0.08807406, 0.04622194)
# hh_info_dat <- data.frame(hh_size = c(1,2,2,3,3,4,4,5,10,4), 
#                           hh_type = c("Single", 
#                                       "Couple","Single, 1 Child",
#                                       "Couple, 1 Child", "Single, 2 Children",
#                                       "Couple, 2 Children", "Single, 3 Children",
#                                       "Couple, 3 Children",
#                                       "Multi-person (care home)", "Multi-person (student)"),
#                           hh_prob = c(0.385, 0.283, 0.045, 0.095, 0.022, 0.114, 0.007, 0.045, 0.001, 0.004))
# prop_inf_vec <- c(rep(0.2/9, 9))
# high_risk_vec <- c(0.01, 0.05, 0.05, 0.07, 0.1, 0.1, 0.15, 0.25, 0.5)
# n_postcodes <- 5
# contact_matrix_wo_hh <- readRDS("/s-schijf/ainsliek/vac-a-mole/inst/extdata/data/Contactmatrix_withoutHH_2020-12-04.rds")
# sigma_mat <- matrix(rep(1,length(age_dist)*length(age_dist)), nrow = length(age_dist))
# beta_mat <- matrix(rep(0.2,length(age_dist)*length(age_dist)), nrow = length(age_dist))
# 
# pop_mat <- create_pop_w_households(n, age_dist, hh_info_dat, n_postcodes, high_risk_vec, perc_vac_hesitant = 0.3)
# ages <- pop_mat[,"Age_Group"]
# contacts <- get_contacts(contact_matrix_wo_hh, pop_mat[,"ID"], pop_mat[,"Household_ID"], pop_mat[,"Age_Group"],
#                         pop_mat[,"Postcode"], 0.9, household_isolation = c(rep(0,n)))
# state_col <- determine_intitial_states(pop_mat[,"Age_Group"], prop_inf = c(rep(0.2/9, 9)), prop_rec = c(rep(0.1/9, 9)));
# state_col[which(state_col == 1)] <- 3
# 
# out <- list()
#   for (i in 1:dim(pop_mat)[1]){
#     print(i)
#     contacts_i <- which(contacts[i,] == 1)
# 
#     out[[i]] <- make_contacts(state_i = state_col[i],
#                               age_i = ages[i],
#                               states_c = state_col[contacts_i],
#                               contacts_i = contacts_i,
#                               age_vec = ages,
#                               sigma = sigma_mat,
#                               beta = beta_mat,
#                               vac_mech = 1,
#                               vac_status_i = 0,
#                               ve = c(0.534, 0.928),
#                               tau = 1,
#                               time_step = 1
#                               )
#     state_col[i] <- out[[i]]$state_i
#   }
# which(state_col == 1)

*/
