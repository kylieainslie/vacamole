#include <Rcpp.h>
using namespace Rcpp;

// Determine total number of community contacts
// @param contact_matrix Contact matrix

// [[Rcpp::export]]
List get_total_community_contacts(NumericMatrix contact_matrix){

  NumericVector total_comm_contacts = contact_matrix.nrow();
  NumericVector contact_decimal = contact_matrix.nrow();
  NumericMatrix contact_probs(contact_matrix.nrow());

  for (int c = 0; c < contact_matrix.nrow(); ++c){
    total_comm_contacts[c] = sum(contact_matrix(c,_));
    contact_decimal[c] = total_comm_contacts[c] - floor(total_comm_contacts[c]);
    contact_probs(c,_) = contact_matrix(c,_)/total_comm_contacts[c];
  }
  
  List rtn;
  rtn["total_comm_contacts"] = total_comm_contacts;
  rtn["contact_decimal"] = contact_decimal;
  rtn["contact_probabilities"] = contact_probs;
  return rtn;
}
// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
contact_matrices_all <- read.delim("../inst/extdata/data/S2_contact_matrices_withPico3.tsv")
contact_matrix_april2020 <- contact_matrices_all %>%
  filter(survey == "April 2020") %>%
  filter(contact_type == "community") %>%
  select(-survey, -contact_type) %>%
  pivot_wider(., names_from = cont_age, values_from = m_est)

contact_matrix_input <- as.matrix(contact_matrix_april2020[,-1])
rownames(contact_matrix_input) <- contact_matrix_april2020$part_age

x <- get_total_community_contacts(contact_matrix_input)
*/
