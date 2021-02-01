#include <Rcpp.h>
using namespace Rcpp;

// Determine inititial infected individuals
// @param ages vector of ages for each individual in the population
// @param prop_inf vector of proportions of each age group to be infected
// @param prop_rec vector of proportions of each age group recovered

// [[Rcpp::export]]
NumericVector determine_intitial_states(NumericVector ages, 
                                        NumericVector prop_inf,
                                        NumericVector prop_rec) {
  
  int n = ages.length();
  NumericVector state(n); // initialize state vector - initially infected individuals will be updated from 0 -> 1
  
  // loop through individuals and determine if they are infected
  for (int i = 0; i < n; ++i){
    int age_i = ages[i];
    double prob_inf_i = prop_inf[age_i - 1];
    double prob_rec_i = prop_rec[age_i - 1];
    //double prob_inf_i2 = prob_inf[age_i - 1] * 8;
    NumericVector rdm = runif(2);
    
    if(rdm[0] < prob_inf_i){ state[i] = 3;} // infected, symptomatic
    if(rdm[1] < prob_rec_i){ state[i] = 7;} // recovered
    
  }
  
  return(state);
}
  

/*** R
n = 10000
age_vec <- sample(1:9, n, TRUE)
prop_inf_vec <- c(rep(0.01/9, 9))
prop_rec_vec <-  c(rep(0.1/9, 9))

initial_infections <- determine_intitial_states(ages = age_vec, 
                                                prop_inf = prop_inf_vec,
                                                prop_rec = prop_rec_vec)

sum(initial_infections[initial_infections == 3])
sum(initial_infections[initial_infections == 7])

*/
