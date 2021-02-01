#include <Rcpp.h>
using namespace Rcpp;

#ifndef HELPERS_H
#define HELPERS_H

NumericVector in_group(NumericVector x, NumericVector y, int id);

NumericVector notin_group(NumericVector x, NumericVector y, int id);

NumericVector subset_vec(IntegerVector x, NumericVector y);

int my_sum(NumericVector x);

IntegerVector which_cpp( NumericVector x, int val, String type);

NumericVector rep_cpp(NumericVector x, NumericVector y);

IntegerVector match_and_find(IntegerVector x, IntegerVector y, IntegerVector z);

NumericVector make_age_prob_vec(IntegerVector ages, NumericVector probabilities);

#endif