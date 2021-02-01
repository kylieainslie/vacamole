#include <Rcpp.h>
using namespace Rcpp;

#ifndef GET_CONTACTS_2_H
#define GET_CONTACTS_2_H

List get_contacts_2(List contact_info,
                             NumericVector ids_to_get_contacts_for,
                             NumericVector ids, 
                             NumericVector hh_ids, 
                             NumericVector ages,
                             NumericVector postcodes, 
                             double prop_contacts_within_postcode,
                             NumericVector household_isolation);
  

#endif