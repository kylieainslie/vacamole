#include <Rcpp.h>
#include "get_total_community_contacts.h"
#include "helpers.h"

using namespace Rcpp;

// Determine total number of community contacts
// @param contact_matrix Contact matrix
// @param ids Vector of individual ids
// @param hh_ids Vector of household ids for each individual
// @param ages Vector of ages for each individual

// [[Rcpp::export]]
NumericMatrix get_contacts(NumericMatrix contact_mat, 
                           NumericVector ids, 
                           NumericVector hh_ids, 
                           NumericVector ages,
                           NumericVector postcodes, 
                           double prop_contacts_within_postcode,
                           NumericVector household_isolation){
  
  // define function arguments
  NumericVector pop_ids = ids;
  int n = pop_ids.length();
  NumericVector hh_id_vec = hh_ids;
  NumericVector age_vec = ages;
  NumericVector age_groups = unique(ages);
  NumericVector postcode_id_vec = postcodes;
  NumericMatrix contacts(n);
  NumericVector hh_isolation = household_isolation;
  
  // determine number of contacts based on contact matrix for each age group
  NumericMatrix contact_matrix = contact_mat;
  List contact_totals = get_total_community_contacts(contact_matrix);
  NumericVector total_comm_contacts = contact_totals["total_comm_contacts"];
  NumericVector contact_decimal = contact_totals["contact_decimal"];
  NumericMatrix contact_probs = contact_totals["contact_probabilities"];
  
  // loop over individuals
  for (int i = 0; i < n; ++i){
    //Rcpp::Rcout << "Person ID: "<< pop_ids[i] << std::endl;
    
    int age_i = age_vec[i]; // individual i's age
    int pc_i = postcode_id_vec[i]; // individual i's postcode
    //NumericVector pc_i_members = in_group(postcode_id_vec, pop_ids, pc_i);
    
    // all household members
    int hh_i = hh_id_vec[i];
    //Rcpp::Rcout << "Household ID for person i: "<< hh_i << std::endl;
    NumericVector hh_i_members = in_group(hh_id_vec, pop_ids, hh_i);
    //IntegerVector hh_i_members_minus_i = which_cpp(hh_i_members, i+1, "notin");
    //Rcpp::Rcout << "IDs of members of i's household: " << hh_i_members << std::endl;
    
    // if a person is in isolation, then they only make contacts with their household
    // change value of contacts 0 -> 1 for household members
    if (household_isolation[i] == 1){
      for (int h = 0; h < hh_i_members.length(); ++h){
        if(hh_i_members[h] == pop_ids[i]){continue;
        } else {contacts(i,hh_i_members[h]) = 1;}
      }
    } else {
    // fixed number of community members - use kernel for # near and far neighbors
    int tot_comm_contacts_i = total_comm_contacts[age_i-1];
    NumericVector rdm = runif(1);
    if(rdm[0] < contact_decimal[age_i-1]){ tot_comm_contacts_i += 1; }
    //Rcpp::Rcout << "Total community contacts for person i: "<< tot_comm_contacts_i << std::endl;
    
    // determine proportion of contacts in postcode and proportion of contacts outside postcode
    int num_in_pc_contacts = floor(prop_contacts_within_postcode * tot_comm_contacts_i);
    int num_out_pc_contacts = tot_comm_contacts_i - num_in_pc_contacts;
    //Rcpp::Rcout << "Community contacts in postcode for person i: "<< num_in_pc_contacts << std::endl;
    //Rcpp::Rcout << "Community contacts outside postcode for person i: "<< num_out_pc_contacts << std::endl;
    
    // determine if there are already contacts filled in from other individuals
    int num_prior_contacts_i = sum(contacts(i,_));
    //Rcpp::Rcout << "Prior community contacts for person i: "<<num_prior_contacts_i << std::endl;
    
    
    if(num_prior_contacts_i > 0){
      IntegerVector prior_contacts_i = which_cpp(contacts(i,_),1, "in");
      //Rcpp::Rcout << "Prior contacts for person i: "<< prior_contacts_i << std::endl;
      
      for(int c = 0; c < num_prior_contacts_i; ++c){
        int c_index = prior_contacts_i[c];
        if((postcode_id_vec[c_index] == postcode_id_vec[i]) & 
           (hh_id_vec[c] != hh_i)){
          num_in_pc_contacts -= 1;
          //Rcpp::Rcout << "New number in postcode contacts for person i: "<< num_in_pc_contacts << std::endl;
          
        } else if (postcode_id_vec[c_index] != postcode_id_vec[i]){
          num_out_pc_contacts -= 1;
          //Rcpp::Rcout << "New number outside postcode contacts for person i: "<< num_out_pc_contacts << std::endl;
        }
      }
    }
      
    // change value of contacts 0 -> 1 for household members
    for (int h = 0; h < hh_i_members.length(); ++h){
      if(hh_i_members[h] == pop_ids[i]){continue;
      } else {contacts(i,hh_i_members[h]) = 1;}
    }
    
    if ( (num_in_pc_contacts < 1) & (num_out_pc_contacts < 1) ){continue;}
      
      // get possible contacts for person i (excluding prior contacts) 
      // exclude people in isolation
      IntegerVector not_hh_iso_contacts = which_cpp(household_isolation, 0, "in");
      IntegerVector not_prior_i_contacts = which_cpp(contacts(i,_), 0, "in");
      IntegerVector possible_i_contacts = intersect(not_hh_iso_contacts, not_prior_i_contacts);
      // need to exclude people in isolation from this!
      
      // loop through possible contacts to determine if a contact is made
      int counter = 0;
      int in_pc_counter = 0;
      int out_pc_counter = 0;
      
      while (counter < (num_in_pc_contacts + num_out_pc_contacts)){
        for (int pc = 0; pc < possible_i_contacts.length(); ++pc){ // find contacts until total number of contacts for that day is reached
          NumericVector rdm2 = runif(2);
          // inside postcode contacts
          if (in_pc_counter < num_in_pc_contacts){
            if((postcode_id_vec[pc] == pc_i) & 
               (rdm2[0] < contact_probs(age_i-1,age_vec[pc]-1)) //&
               //(sum(contacts(possible_i_contacts[pc],_)) < total_comm_contacts[age_vec[pc]-1])
                 ){
              contacts(i,possible_i_contacts[pc]) = 1; // make contact for person i
              contacts(possible_i_contacts[pc],i) = 1; // reciprocal contact
              in_pc_counter += 1;
              counter += 1;
              //Rcpp::Rcout << "While loop counter: "<< counter << std::endl;
              //Rcpp::Rcout << "In postcode counter: "<< in_pc_counter << std::endl;
            } else {continue;}
          }
          // outside postcode contacts
          else if (out_pc_counter < num_out_pc_contacts){
            if((postcode_id_vec[pc] != pc_i) & 
               (rdm2[1] < contact_probs(age_i-1,age_vec[pc]-1)) //&
               //(sum(contacts(possible_i_contacts[pc],_)) < total_comm_contacts[age_vec[pc]-1])
               ){
              contacts(i,possible_i_contacts[pc]) = 1;
              contacts(possible_i_contacts[pc],i) = 1;
              out_pc_counter += 1;
              counter += 1;
              //Rcpp::Rcout << "While loop counter: "<< counter << std::endl;
              //Rcpp::Rcout << "Out postcode counter: "<< out_pc_counter << std::endl;
            } else {continue;}
          } else {continue;}
        } // end for loop
      } // end while loop
    } // end ifelse statement for household isolation
  } // end loop of i
     
  
  // output
  return contacts;
}

/*** R
pop <- readRDS("/s-schijf/ainsliek/vac-a-mole/inst/extdata/data/test_pop_100.rds")
contact_matrix_wo_hh <- readRDS("/s-schijf/ainsliek/vac-a-mole/inst/extdata/data/Contactmatrix_withoutHH_2020-12-04.rds")
hh_isolation_ind <- sample(c(0,1), dim(pop)[1], prob = c(0.95, 0.05), replace = TRUE)

contacts <- get_contacts(contact_mat = contact_matrix_wo_hh, 
                         ids = pop[,"ID"],
                         hh_ids = pop[,"Household_ID"], 
                         ages = pop[,"Age_Group"],
                         postcodes = pop[,"Postcode"], 
                         prop_contacts_within_postcode = .85,
                         household_isolation = hh_isolation_ind
)
 
*/
