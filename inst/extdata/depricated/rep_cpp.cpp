#include <Rcpp.h>
using namespace Rcpp;

// Rep function that imitates rep with times argument in R
// @param x vector of inputs you want repeated
// @param y vector of time you want each element in x repeated
// [[Rcpp::export]]
NumericVector rep_cpp(NumericVector x, NumericVector y) {
  int n = y.size();
  NumericVector myvector(sum(y));
  int ind = 0;
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < y(i); ++j) {
      myvector(ind) = x[i];
      ind = ind + 1;
    }
  }
  return myvector;
}


/*** R
x <- 1:5
y<- c(1,2,3,2,1)
rep_cpp(x,y)
*/
