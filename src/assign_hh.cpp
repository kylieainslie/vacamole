#include <Rcpp.h>
using namespace Rcpp;

// Assign people to households function
//
// @param pop_size Number of individuals in the population
// @param hh_info_dat dataframe with the following columns:
//        1. household size (number of people in that type of household)
//        2. household composition (single adult, couple w/o children, couple w/ 1 child, etc.)
//        3. proportion of households of that size and type

// [[Rcpp::export]]
DataFrame assign_hh(const int pop_size,
                    const DataFrame hh_info_mat
){
  
  int n_hh = ceil(pop_size / Rcpp::sum(as<NumericVector>(hh_info_mat[0]) * as<NumericVector>(hh_info_mat[2])));
  // Rcpp::Rcout << n_hh << std::endl;
  // Rcpp::Rcout << as<CharacterVector>(hh_info_mat[1]) << std::endl;
  
  CharacterVector pop_hh_vec = Rcpp::as<CharacterVector>(Rcpp::sample(as<CharacterVector>(hh_info_mat[1]), n_hh, true, as<NumericVector>(hh_info_mat[2])));
  // Rcpp::Rcout << pop_hh_vec << std::endl;
  
  // household IDs
  IntegerVector hh_ids = seq(1, n_hh);
  
  IntegerVector hh_size (n_hh);
  for (int i = 0; i < n_hh; ++i){
    // assign household size to each household
    if (pop_hh_vec[i] == "Single"){ hh_size[i] = 1;
    } else if ((pop_hh_vec[i] == "Couple") | (pop_hh_vec[i] == "Single, 1 Child")) { hh_size[i] = 2;
    } else if ((pop_hh_vec[i] == "Couple, 1 Child") | (pop_hh_vec[i] == "Single, 2 Children")) { hh_size[i] = 3;
    } else if ((pop_hh_vec[i] == "Couple, 2 Children") | (pop_hh_vec[i] == "Single, 3 Children")) { hh_size[i] = 4;
    } else if (pop_hh_vec[i] == "Couple, 3 Children") { hh_size[i] = 5;
    } else if (pop_hh_vec[i] == "Multi-person (other)") { hh_size[i] = 10;
    }
  }
  
  // create data frame for output
  DataFrame df_out = DataFrame::create( Named("hh_id") = hh_ids,
                                        Named("hh_type") = pop_hh_vec,
                                        Named("hh_size") = hh_size
  );
  
  return(df_out);
}

/*** R
n <- 100
hh_info_dat <- data.frame(hh_size = c(1,2,2,3,3,4,4,5,10), 
                          hh_type = c("Single", 
                                      "Couple","Single, 1 Child",
                                      "Couple, 1 Child", "Single, 2 Children",
                                      "Couple, 2 Children", "Single, 3 Children",
                                      "Couple, 3 Children","Multi-person (other)"),
                          hh_prob = c(0.385, 0.283, 0.045, 0.095, 0.022, 0.114, 0.007, 0.045, 0.005))
assign_hh(pop_size = n, hh_info_mat = hh_info_dat)
*/
