#include <Rcpp.h>
#include "helpers.h"

using namespace Rcpp;

// Determine inititial infected individuals
// @param ages vector of ages for each individual in the population
// @param prop_inf vector of proportions of each age group to be infected
// @param prop_rec vector of proportions of each age group recovered
// @return vector of states
// @keywords vac-a-mole
// @export
// [[Rcpp::export]]
List vaccinate(NumericVector num_vaccines,
               NumericVector ages,
               NumericVector high_risk_status,
               NumericVector vaccine_hesitancy,
               NumericVector vac_status,
               NumericVector vac_type,
               int time_step,
               NumericVector time_dose1,
               NumericVector time_dose2,
               NumericVector protection_status,
               int time_btw_doses,
               int time_to_protection,
               NumericMatrix age_group_vac_order) {
  
  //int n = ages.length();
  NumericVector vac_hesitant = vaccine_hesitancy;
  NumericVector age_vec = ages;
  NumericVector high_risk = high_risk_status;

  // get people who are not vaccine hesitant
  IntegerVector which_not_vac_hesitant = which_cpp(vac_hesitant, 0, "in");
  // loop through vaccine types
  for (int type = 0; type < num_vaccines.length(); ++type){
    int vac_counter = 0;
    NumericVector vac_order = age_group_vac_order(type,_);
    
    // Rcpp::Rcout << "type: " << type + 1 << std::endl;
    
    for (int a = 0; a < vac_order.length(); ++a){
      int age_a = vac_order[a];
      // Rcpp::Rcout << "Age group: " << age_a << std::endl;
      
      IntegerVector which_a = which_cpp(age_vec, age_a, "in");
      IntegerVector to_vac_a = intersect(which_a, which_not_vac_hesitant);
      IntegerVector which_high_risk = which_cpp(high_risk, 1, "in");
      IntegerVector to_vac_a_high_risk = intersect(which_a, which_high_risk);
      int n_to_vac_a = to_vac_a.length();
      int n_to_vac_a_hr = to_vac_a_high_risk.length();
      int vac_counter_a = 0;
      
      for (int i = 0; i < n_to_vac_a; ++i){
        
        if(vac_counter == num_vaccines[type]){break;} // if all the vaccines are used up for that time point, stop
        
        int id = to_vac_a[i]; // this is the position
        // Rcpp::Rcout << "ID: " << id << std::endl;
        // Rcpp::Rcout << "vac_status_i: " << vac_status[id] << std::endl;
        // Rcpp::Rcout << "vac_type_i " << vac_type[id] << std::endl;
        // 
        // change protection status
        if ((vac_status[id] == 1) & (time_step == time_dose1[id] + time_to_protection)){
          protection_status[id] = 1;
        } else if ((vac_status[id] == 2) & (time_step == time_dose2[id] + time_to_protection)){
          protection_status[id] = 2;
        }
        
        // high risk individuals first
        if (vac_counter_a < n_to_vac_a_hr){
          if (high_risk_status[id] == 1){
            // dose 1
            if (vac_status[id] == 0){
              vac_status[id] = 1;
              vac_type[id] = type + 1;
              vac_counter += 1;
              vac_counter_a +=1;
              time_dose1[id] = time_step;
            }
            // dose 2
            if ((vac_status[id] == 1) & (time_step > time_dose1[id] + time_btw_doses) & (vac_type[id] == (type + 1))){
              vac_status[id] = 2;
              time_dose2[id] = time_step;
              vac_counter += 1;
              vac_counter_a +=1;
            }
          } 
        } else {
          //Rcpp::Rcout << "flag" << std::endl;
          
          // dose 1
          if (vac_status[id] == 0){
            vac_status[id] = 1;
            vac_type[id] = type + 1;
          
            //Rcpp::Rcout << "vac_type_i " << vac_type[id] << std::endl;
            
            time_dose1[id] = time_step;
            vac_counter += 1;
            vac_counter_a +=1;
          }
          // dose 2
          if ((vac_status[id] == 1) & (time_step > time_dose1[id] + time_btw_doses) & (vac_type[id] == (type + 1))){
            vac_status[id] = 2;
            time_dose2[id] = time_step;
            vac_counter += 1;
            vac_counter_a += 1;
          }
        }
      } // end loop over individuals in age group a
    } // end loop over age groups
  } // end loop over vaccine types
  
  // output
  List rtn = List::create(Named("vac_status") = vac_status,
                          Named("vac_type") = vac_type,
                          Named("time_dose1") = time_dose1,
                          Named("time_dose2") = time_dose2,
                          Named("protection_status") = protection_status
  );
  
  return rtn;
}


/*** R
# n <- 100
# num_vac <- c(10, 10, 20)
# age_dist <- c(0.10319920, 0.11620856, 0.12740219, 0.12198707, 0.13083463, 0.14514332,
#               0.12092904, 0.08807406, 0.04622194)
# ages <- sample(1:9, n, prob = age_dist, replace = TRUE)
# high_risk_vec <- sample(c(0,1), n, prob = c(0.9, 0.1), replace = TRUE)
# vac_hesitant_vec <- sample(c(0,1), n, prob = c(0.95, 0.05), replace = TRUE)
# vac_status <- c(rep(0,n))
# vac_type <- c(rep(0,n))
# vac_order <- matrix(c(rep(0, 3*7)), nrow = 3)
# vac_order[1,] <- c(9,8,7,6,5,4,3)
# vac_order[2,] <- c(9,8,7,6,5,4,3)
# vac_order[3,] <- c(3,4,5,6,7,8,9)
# 
#   results <- vaccinate(num_vac,
#                        ages,
#                        high_risk_vec,
#                        vac_hesitant_vec,
#                        vac_status,
#                        vac_type,
#                        time_step = 22,
#                        time_dose1 = c(rep(0,n)),
#                        time_dose2 = c(rep(0,n)),
#                        protection_status = c(rep(0,n)),
#                        time_btw_doses = 1,
#                        time_to_protection = 14,
#                        age_group_vac_order = vac_order)
# 
#   vac_dat <- data.frame(id = 1:n,
#                         age = ages,
#                         high_risk_status = high_risk_vec,
#                         vac_hesitant = vac_hesitant_vec,
#                         vac_status = results$vac_status,
#                         vac_type = results$vac_type,
#                         time_dose1 = results$time_dose1,
#                         time_dose2 = results$time_dose2
#   )
  */
