#include <Rcpp.h>
using namespace Rcpp;

#ifndef CREATE_POP_W_HOUSEHOLDS_H
#define CREATE_POP_W_HOUSEHOLDS_H

NumericMatrix create_pop_w_households(const int pop_size,
                                      const NumericVector age_dist,
                                      const DataFrame hh_info,
                                      const int n_postcodes,
                                      const NumericVector perc_high_risk,
                                      const double perc_vac_hesitant);

#endif