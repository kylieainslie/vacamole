#include <Rcpp.h>
using namespace Rcpp;

#ifndef MAKE_CONTACTS_H
#define MAKE_CONTACTS_H

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
                   int time_step);
  
#endif