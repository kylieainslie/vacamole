#include <Rcpp.h>
#include "helpers.h"

using namespace Rcpp;

// Determine total number of community contacts
// @param contact_matrix Contact matrix
// @param ids Vector of individual ids
// @param hh_ids Vector of household ids for each individual
// @param ages Vector of ages for each individual

// [[Rcpp::export]]
List get_contacts_2(List contact_info,
                             NumericVector ids_to_get_contacts_for,
                             NumericVector ids, 
                             NumericVector hh_ids, 
                             NumericVector ages,
                             NumericVector postcodes, 
                             double prop_contacts_within_postcode,
                             NumericVector household_isolation){
  
  // define function arguments
  int n = ids_to_get_contacts_for.length();
  IntegerVector age_groups = as<IntegerVector>(unique(ages));
  
  // determine number of contacts based on contact matrix for each age group
  NumericVector total_comm_contacts = contact_info["total_comm_contacts"];
  NumericVector contact_decimal = contact_info["contact_decimal"];
  NumericMatrix contact_probs = contact_info["contact_probabilities"];
  int mean_contacts = ceil(mean(total_comm_contacts));
  
  // make matrix of size n * mean_contacts x 2 to store contact pairs
  // try later with std::vector<int> and push_back()
  //NumericMatrix contact_pairs(n*mean_contacts, 2); // column 1 = node, column 2 = edge
  NumericVector nodes(n*mean_contacts);
  NumericVector edges(n*mean_contacts);
  int ind = 0; // to keep track of where to put next contact pair
  
  // loop over individuals
  for (int i = 0; i < n; ++i){
    //Rcpp::Rcout << "Person ID: "<< i+1 << std::endl;
    int id_i = ids_to_get_contacts_for[i]; // id
    
    if(household_isolation[id_i - 1] == 1){continue;}
    // if a person is in household isolation, they don't make any community contacts
    
    // get some info about person i
    int age_i = ages[id_i - 1];        // age
    int pc_i = postcodes[id_i - 1];    // postcode
    int hh_i = hh_ids[id_i - 1];       // household id
    NumericVector contact_probs_i = contact_probs(_,age_i - 1); //contact probabilities by age group
    int comm_contacts_i = total_comm_contacts[age_i-1]; // num contacts
    //Rcpp::Rcout << "Contact probs: "<< contact_probs_i << std::endl;
    NumericVector rdm = runif(1);
    if(rdm[0] < contact_decimal[age_i-1]){ comm_contacts_i += 1; }
    // determine proportion of contacts in postcode and proportion of contacts outside postcode
    int n_in_pc_contacts = floor(prop_contacts_within_postcode * comm_contacts_i);
    int n_out_pc_contacts = comm_contacts_i - n_in_pc_contacts;
    //Rcpp::Rcout << "n in pc contacts: "<< n_in_pc_contacts << std::endl;
    //Rcpp::Rcout << "n out pc contacts: "<< n_out_pc_contacts << std::endl;
    
    // remove contacts in household and isolation
    IntegerVector not_in_hh_contacts = which_cpp(hh_ids, hh_i, "notin");
    IntegerVector not_hh_iso_contacts = which_cpp(household_isolation, 0, "in");
    IntegerVector not_hh_contacts = intersect(not_hh_iso_contacts, not_in_hh_contacts);
    
    // determine if there are already contacts filled in from other individuals
    IntegerVector prior_contacts_i = which_cpp(edges, id_i, "in");
    int n_prior_contacts_i = prior_contacts_i.length();
    //Rcpp::Rcout << "n prior contacts: "<< n_prior_contacts_i << std::endl;
    
    // declare empty vectors
    IntegerVector not_prior_i_contacts;
    NumericVector not_prior_i_contacts_ages;
    NumericVector not_prior_i_contacts_pcs;
    
    if(n_prior_contacts_i > 0){
      
      if (n_prior_contacts_i > comm_contacts_i){ continue; } 
      NumericVector ids_of_prior_i_contacts = subset_vec(prior_contacts_i, nodes);

      //Rcpp::Rcout << "Prior contacts for person i: "<< ids_of_prior_i_contacts << std::endl;
      
      for(int c = 0; c < n_prior_contacts_i; ++c){
        int c_index = ids_of_prior_i_contacts[c] - 1;
        int pc_c = postcodes[c_index];
        if( pc_c == pc_i){ n_in_pc_contacts -= 1;
        } else { n_out_pc_contacts -= 1; }
      }
      
      // get possible contacts for person i 
      IntegerVector pos_of_prior_i_contacts = as<IntegerVector>(ids_of_prior_i_contacts) - 1;
      IntegerVector not_prior_i_contacts = setdiff(not_hh_contacts, pos_of_prior_i_contacts);
      NumericVector not_prior_i_contacts_ages = subset_vec(not_prior_i_contacts, ages);
      NumericVector not_prior_i_contacts_pcs = subset_vec(not_prior_i_contacts, postcodes);
      
      // sample age groups to be contacted
      if ( n_in_pc_contacts > 0 ){
        IntegerVector possible_in_pc_contacts = which_cpp(not_prior_i_contacts_pcs, pc_i, "in");
        NumericVector possible_in_pc_ages = subset_vec(possible_in_pc_contacts, not_prior_i_contacts_ages);
        NumericVector possible_in_pc_age_prob_vec = make_age_prob_vec(as<IntegerVector>(possible_in_pc_ages), as<NumericVector>(contact_probs_i));
        if (n_in_pc_contacts > possible_in_pc_ages.length()) {n_in_pc_contacts = possible_in_pc_ages.length();}
        NumericVector ages_to_contact_in_pc = sample(possible_in_pc_ages, n_in_pc_contacts, false, possible_in_pc_age_prob_vec);
        //Rcpp::Rcout << "Possible in pc ages: "<< possible_in_pc_ages << std::endl;
        //Rcpp::Rcout << "Ages to contact in pc: "<< ages_to_contact_in_pc << std::endl;
        //Rcpp::Rcout << "Possible in pc contacts: "<< possible_in_pc_contacts << std::endl;
        IntegerVector contacts_in_pc = match_and_find(as<IntegerVector>(possible_in_pc_ages), 
                                                      as<IntegerVector>(ages_to_contact_in_pc), 
                                                      possible_in_pc_contacts);
        
        // add contacts for person i to contact_pairs matrix
        for (int j = 0; j < contacts_in_pc.length(); j++){
          nodes[ind] = id_i;
          edges[ind] = contacts_in_pc[j];
          ind += 1;
        }
      }
      
      if ( n_out_pc_contacts > 0 ){
        IntegerVector possible_out_pc_contacts = which_cpp(not_prior_i_contacts_pcs, pc_i, "notin");
        NumericVector possible_out_pc_ages = subset_vec(possible_out_pc_contacts, not_prior_i_contacts_ages);
        NumericVector possible_out_pc_age_prob_vec = make_age_prob_vec(as<IntegerVector>(possible_out_pc_ages), contact_probs_i);
        if (n_out_pc_contacts > possible_out_pc_ages.length()) {n_out_pc_contacts = possible_out_pc_ages.length();}
        NumericVector ages_to_contact_out_pc = sample(possible_out_pc_ages, n_out_pc_contacts, false, possible_out_pc_age_prob_vec);
        //Rcpp::Rcout << "Ages to contact out pc: "<< ages_to_contact_out_pc << std::endl;
        IntegerVector contacts_out_pc = match_and_find(as<IntegerVector>(possible_out_pc_ages), 
                                                       as<IntegerVector>(ages_to_contact_out_pc), 
                                                       possible_out_pc_contacts);
        
        // add contacts for person i to contact_pairs matrix
        for (int k = 0; k < contacts_out_pc.length(); k++){
          nodes[ind] = id_i;
          edges[ind] = contacts_out_pc[k];
          ind += 1;
        }
      }
    } else {
      // get possible contacts for person i
      IntegerVector not_prior_i_contacts = not_hh_contacts;
      //Rcpp::Rcout << "Not prior i contacts length: "<< not_prior_i_contacts.length() << std::endl;
      //Rcpp::Rcout << "Ages length: "<< ages.length() << std::endl;
      NumericVector not_prior_i_contacts_ages = subset_vec(not_prior_i_contacts, ages);
      NumericVector not_prior_i_contacts_pcs = subset_vec(not_prior_i_contacts, postcodes);
      
      // in pc contacts
      IntegerVector possible_in_pc_contacts = which_cpp(not_prior_i_contacts_pcs, pc_i, "in");
      NumericVector possible_in_pc_ages = subset_vec(possible_in_pc_contacts, not_prior_i_contacts_ages);
      NumericVector possible_in_pc_age_prob_vec = make_age_prob_vec(as<IntegerVector>(possible_in_pc_ages), as<NumericVector>(contact_probs_i));
      if (n_in_pc_contacts > possible_in_pc_ages.length()) {n_in_pc_contacts = possible_in_pc_ages.length();}
      NumericVector ages_to_contact_in_pc = sample(possible_in_pc_ages, n_in_pc_contacts, false, possible_in_pc_age_prob_vec);
      //Rcpp::Rcout << "Possible in pc ages: "<< possible_in_pc_ages << std::endl;
      //Rcpp::Rcout << "Ages to contact in pc: "<< ages_to_contact_in_pc << std::endl;
      //Rcpp::Rcout << "Possible in pc contacts: "<< possible_in_pc_contacts << std::endl;
      IntegerVector contacts_in_pc = match_and_find(as<IntegerVector>(possible_in_pc_ages), 
                                                    as<IntegerVector>(ages_to_contact_in_pc), 
                                                    possible_in_pc_contacts);
      
      // add contacts for person i to contact_pairs matrix
      for (int j = 0; j < contacts_in_pc.length(); j++){
        nodes[ind] = id_i;
        edges[ind] = contacts_in_pc[j];
        ind += 1;
      }
      
      // out pc contacts
      IntegerVector possible_out_pc_contacts = which_cpp(not_prior_i_contacts_pcs, pc_i, "notin");
      NumericVector possible_out_pc_ages = subset_vec(possible_out_pc_contacts, not_prior_i_contacts_ages);
      NumericVector possible_out_pc_age_prob_vec = make_age_prob_vec(as<IntegerVector>(possible_out_pc_ages), contact_probs_i);
      if (n_out_pc_contacts > possible_out_pc_ages.length()) {n_out_pc_contacts = possible_out_pc_ages.length();}
      NumericVector ages_to_contact_out_pc = sample(possible_out_pc_ages, n_out_pc_contacts, false, possible_out_pc_age_prob_vec);
      //Rcpp::Rcout << "Possible out pc ages: "<< possible_out_pc_ages << std::endl;
      //Rcpp::Rcout << "Ages to contact out pc: "<< ages_to_contact_out_pc << std::endl;
      //Rcpp::Rcout << "Possible out pc contacts: "<< possible_out_pc_contacts << std::endl;
      IntegerVector contacts_out_pc = match_and_find(as<IntegerVector>(possible_out_pc_ages), 
                                                     as<IntegerVector>(ages_to_contact_out_pc), 
                                                     possible_out_pc_contacts);
      
      // add contacts for person i to contact_pairs matrix
      for (int k = 0; k < contacts_out_pc.length(); k++){
        nodes[ind] = id_i;
        edges[ind] = contacts_out_pc[k];
        ind += 1;
      }
    }
    
  } // end loop of i
     
  // output
  List rtn = List::create(Named("nodes") = nodes,
                          Named("edges") = edges);
  return rtn;
}

/*** R
library(dplyr)
library(tidyr)

contact_matrices_all <- read.delim("../inst/extdata/data/S2_contact_matrices_withPico3.tsv")
contact_matrix_april2020 <- contact_matrices_all %>%
  filter(survey == "April 2020") %>%
  filter(contact_type == "community") %>%
  select(-survey, -contact_type) %>%
  pivot_wider(., names_from = cont_age, values_from = m_est)


contact_matrix_input <- as.matrix(contact_matrix_april2020[-1,-1])
rownames(contact_matrix_input) <- contact_matrix_april2020$part_age[-1]

contact_info <- get_total_community_contacts(contact_matrix_input)

n <- 20
id <- 1:n
get_contact_ids = sample(id, 5)
hh_id <- sort(sample(1:6, n, replace = TRUE))
pc_id <- c(rep(1001, n/2), rep(1002, n/2))
ages <- sample(1:9, n, replace = TRUE)
hh_isolation_ind <- sample(c(0,1), n, prob = c(0.95, 0.05), replace = TRUE)

contacts <- get_contacts_2(contact_info = contact_info,
                         ids_to_get_contacts_for = get_contact_ids,  
                         ids = id,
                         hh_ids = hh_id, 
                         ages = ages,
                         postcodes = pc_id, 
                         prop_contacts_within_postcode = .85,
                         household_isolation = hh_isolation_ind
)
 
*/
