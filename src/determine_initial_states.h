#include <Rcpp.h>
using namespace Rcpp;

#ifndef DETERMINE_INITIAL_STATES_H
#define DETERMINE_INITIAL_STATES_H

NumericVector determine_intitial_states(NumericVector ages, 
                                        NumericVector prop_inf,
                                        NumericVector prop_rec);
  
#endif