#include <Rcpp.h>
using namespace Rcpp;

#ifndef GET_CONTACTS_H
#define GET_CONTACTS_H

NumericMatrix get_contacts(NumericMatrix contact_mat, 
                           NumericVector ids, 
                           NumericVector hh_ids, 
                           NumericVector ages,
                           NumericVector postcodes, 
                           double prop_contacts_within_postcode,
                           NumericVector household_isolation);
  

#endif