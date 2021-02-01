#include <Rcpp.h>
#include "assign_hh.h"
#include "helpers.h"

using namespace Rcpp;

// A function that creates a population with households and is stratified by age groups
//
// @param pop_size Number of individuals in the population.
// @param age_dist Age distribution: proportion of individuals in each age group
// @param hh_info_mat data frame with household distribution information
// @param n_postcodes number of distict geographic regions to divide population into
//        if 0 then there will be no household structure
// @param perc_high_risk Vector of percent of population that's high risk by age group
// @param perc_vac_hesitant percentage of population that's vaccine hesitant
// @return A matrix of population with covariates with columns (0-indexed):
//         0. ID
//         1. Household ID
//         2. Household size
//         3. Postal code
//         4. Age group
//         5. High risk status
//         6. Vaccine Hesitancy
//         7. State
// @keywords vac-a-mole
// @export
// [[Rcpp::export]]
NumericMatrix create_pop_w_households(const int pop_size,
                                      const NumericVector age_dist,
                                      const DataFrame hh_info,
                                      const int n_postcodes,
                                      const NumericVector perc_high_risk,
                                      const double perc_vac_hesitant) {
  
// create sequences for sampling
  IntegerVector age_groups = seq(1,9); //Rcpp::IntegerVector::create(1, 2, 3, 4, 5, 6, 7, 8, 9);
  IntegerVector postcodes = seq(1001, 1000 + n_postcodes); //Rcpp::IntegerVector::create(1, 2, 3, 4, 5, 6, 7, 8, 9, 10);
  IntegerVector age_groups_adults_all = seq(3,9);
  IntegerVector age_groups_students = seq(2,3);
  IntegerVector age_groups_parents = seq(4,6);
  IntegerVector age_groups_elderly = seq(7,9);

// create households
    DataFrame hh_dat = assign_hh(pop_size, hh_info);
    int n_hh = as<NumericVector>(hh_dat[0]).length();
    //Rcpp::Rcout << n_hh << std::endl;
    NumericVector hh_postcodes = as<NumericVector>(sample(postcodes, n_hh, true)); // assign each hh a postcode
    
// determine actual population size based on households (not always = pop_size)
   double n = sum(as<NumericVector>(hh_dat[2]));

// create population matrix
   NumericMatrix pop_mat(n, 8);
   colnames(pop_mat) = CharacterVector::create("ID", "Household_ID", "HH_Size", 
            "Postcode", "Age_Group", "High-Risk_Status", "Vac_Hesitancy", "State");
   
// create individual IDs
   IntegerVector ids = seq(1,n); // individual IDs
  
// expand household IDs to have same number of obs as n
   NumericVector hh_ids = rep_cpp(as<NumericVector>(hh_dat[0]),as<NumericVector>(hh_dat[2]));
   NumericVector hh_sizes = rep_cpp(as<NumericVector>(hh_dat[2]),as<NumericVector>(hh_dat[2]));
   NumericVector hh_postcodes_long = rep_cpp(hh_postcodes,as<NumericVector>(hh_dat[2]));
  
// determine vaccine hesitant households
  // pick a single postcode (should only account for approx 1.5% of population)
  // to be vaccine hesitant - this accounts for Orthodox Protestant Denominations
   int vh_pc = as<int>(sample(postcodes, 1));
   Rcpp::Rcout << "Vaccine hesitant postcode: " << vh_pc << std::endl;
   double n_vh_pc = which_cpp(hh_postcodes_long, vh_pc, "in").length();
   Rcpp::Rcout << "Number of people in vaccine hesitant postcode: " << n_vh_pc << std::endl;
  // determine other people to be vaccine hesitant
   double vh_other = perc_vac_hesitant-(n_vh_pc/n);
   Rcpp::Rcout << "Vaccine hesitant prop (other): " << vh_other << std::endl;
// counter for people in the same household
   int counter = 1;
// household types for each hh_id
   CharacterVector hh_type_vec = as<CharacterVector>(hh_dat[1]);
   
// assign age groups to people in households
    for (int i = 0; i < n; ++i){
      
      NumericVector rdm = runif(2); // random numbers
      
      // assign vaccine hesitant individuals
      if (hh_postcodes_long[i] == vh_pc){pop_mat(i,6) = 1;
      } else if (rdm[0] < vh_other) {pop_mat(i,6) = 1;}
      
      // get HH id and size for individual i
        int hh_id = hh_ids[i];
        //int hh_size = hh_sizes[i];
        std::string hh_type = as<std::string>(hh_type_vec[hh_id - 1]);
       // Rcpp::Rcout << hh_type << std::endl;
        
      // advance counter if not first person in household
      if (hh_id == hh_ids[i-1]){counter += 1;
        }
      // Rcpp::Rcout << counter << std::endl;
      // Single
      if ((hh_type == "Single")) {
        pop_mat(i,4) = as<int>(sample(age_groups_adults_all, 1, age_dist[(age_groups_adults_all-1)]));
        continue;
      } 
      
      // Couple, no children
      if (hh_type == "Couple"){
        if (counter == 1){ // first adult
          pop_mat(i,4) = as<int>(sample(age_groups_adults_all, 1, age_dist[(age_groups_adults_all-1)]));
        } else if (counter == 2){ // second adult (same age group as first)
          pop_mat(i,4) = pop_mat(i-1, 4);
        }
      }
      
      // Couple with children
      if ((hh_type == "Couple, 1 Child") | (hh_type == "Couple, 2 Children") | (hh_type == "Couple, 3 Children")){
        if (counter == 1){ // first adult
          pop_mat(i,4) = as<int>(sample(age_groups_parents, 1, age_dist[(age_groups_parents-1)]));
        } else if (counter == 2){ // second adult (same age group as first)
          pop_mat(i,4) = pop_mat(i-1, 4);
        } else if (counter == 3){
          if (pop_mat(i-1,4) == 4){pop_mat(i,4) = 1;
          } else if (pop_mat(i-1,4) > 4) {pop_mat(i,4) = 2;}
        } else if (counter == 4){
          pop_mat(i,4) = pop_mat(i-1, 4); // 2nd child is in same age group as 1st
        } else if (counter == 5){
          pop_mat(i,4) = pop_mat(i-1, 4); // 3rd child is in same age group as 1st and 2nd
        }
      }
      
      // Single with children 
      if ((hh_type == "Single, 1 Child") | (hh_type == "Single, 2 Children") | (hh_type == "Single, 3 Children")){
        if (counter == 1){ // first adult
          pop_mat(i,4) = as<int>(sample(age_groups_parents, 1, age_dist[(age_groups_parents-1)]));
        } else if (counter == 2){ // first child
          if (pop_mat(i-1,4) == 4){pop_mat(i,4) = 1;
          } else if (pop_mat(i-1,4) > 4) {pop_mat(i,4) = 2;}
        } else if (counter == 3){
          pop_mat(i,4) = pop_mat(i-1, 4); // 2nd child is in same age group as 1st
        } else if (counter == 4){
          pop_mat(i,4) = pop_mat(i-1, 4); // 3rd child is in same age group as 1st and 2nd
        }
      }
      
      // Multi-person (care home)
      if (hh_type == "Multi-person (care home)"){
        pop_mat(i,4) = as<int>(sample(age_groups_elderly, 1, age_dist[(age_groups_elderly-1)]));
        pop_mat(i,5) = 1; //assign high-risk status
      }
      
      // Multi-person (student/young professional)
      if (hh_type == "Multi-person (student)"){
        pop_mat(i,4) = as<int>(sample(age_groups_students, 1, age_dist[(age_groups_students-1)]));
      }
      
    // return counter to 1 if last person in household
    if (hh_id != hh_ids[i+1]){counter = 1;}
    
    // assign high-risk status
    if(pop_mat(i,5) == 0){
      if (rdm[1] < perc_high_risk[pop_mat(i,4)-1]){
        pop_mat(i,5) = 1; //assign high-risk status
      }
    }
    
    }
    
    
    // put vectors into matrix for output
    pop_mat(_,0) = ids;
    pop_mat(_,1) = as<IntegerVector>(hh_ids);
    pop_mat(_,2) = as<IntegerVector>(hh_sizes);
    pop_mat(_,3) = as<IntegerVector>(hh_postcodes_long);
    
    return(pop_mat);
}

/*** R
n <- 1000
age_dist <- c(0.10319920, 0.11620856, 0.12740219, 0.12198707, 0.13083463, 0.14514332, 
              0.12092904, 0.08807406, 0.04622194)
hh_info_dat <- data.frame(hh_size = c(1,2,2,3,3,4,4,5,10,4), 
                          hh_type = c("Single", 
                                      "Couple","Single, 1 Child",
                                      "Couple, 1 Child", "Single, 2 Children",
                                      "Couple, 2 Children", "Single, 3 Children",
                                      "Couple, 3 Children",
                                      "Multi-person (care home)", "Multi-person (student)"),
                                      hh_prob = c(0.385, 0.283, 0.045, 0.095, 0.022, 0.114, 0.007, 0.045, 0.001, 0.004))
high_risk_vec <- c(0.01, 0.05, 0.05, 0.07, 0.1, 0.1, 0.15, 0.25, 0.5)

pop <- create_pop_w_households(pop_size = n, 
                               age_dist = age_dist, 
                               hh_info = hh_info_dat, 
                               n_postcodes = 50,
                               perc_high_risk = high_risk_vec,
                               perc_vac_hesitant = 0.4)
*/
