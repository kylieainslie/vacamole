#include <Rcpp.h>
using namespace Rcpp;

#ifndef VACCINATE_H
#define VACCINATE_H

List vaccinate(NumericVector num_vaccines,
               NumericVector ages,
               NumericVector high_risk_status,
               NumericVector vaccine_hesitancy,
               NumericVector vac_status,
               NumericVector vac_type,
               int time_step,
               NumericVector time_dose1,
               NumericVector time_dose2,
               NumericVector protection_status,
               int time_btw_doses,
               int time_to_protection,
               NumericMatrix age_group_vac_order);
#endif