#include <Rcpp.h>
using namespace Rcpp;

// Helper function for subsetting members of same group
// [[Rcpp::export]]
NumericVector in_group(NumericVector x, NumericVector y, int id) {
  return y[x == id];
}

// [[Rcpp::export]]
NumericVector notin_group(NumericVector x, NumericVector y, int id) {
  return y[x != id];
}

// [[Rcpp::export]]
NumericVector subset_vec(IntegerVector positions, NumericVector y) {
    return y[positions];
}

// [[Rcpp::export]]
int my_sum(NumericVector x) {
  NumericVector x2 = x[!is_na(x)];
  return sum(x2);
}

// [[Rcpp::export]]
IntegerVector which_cpp( NumericVector x, int val, String type) {
  
  int nx = x.size();
  std::vector<int> y;
  y.reserve(nx);
  
  for(int i = 0; i < nx; i++) {
    if(type == "in"){
      if (x[i] == val) y.push_back(i);
    } else if (type == "notin"){
      if (x[i] != val) y.push_back(i); 
    }
  }
  
  return wrap(y);
}

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

// [[Rcpp::export]]
IntegerVector match_and_find(IntegerVector x, IntegerVector y, IntegerVector z){
  
  IntegerVector my_table_counts = table(y);
  IntegerVector my_table_names = unique(y);
  IntegerVector match_vec = match(x, my_table_names);
  IntegerVector rtn(y.length());
  int ind = 0;
  for (int i = 0; i < my_table_counts.length(); ++i){
    IntegerVector tmp = which_cpp(as<NumericVector>(match_vec), i+1, "in");
    //Rcpp::Rcout << "tmp: "<< tmp << std::endl;
    if (tmp.length() == 1){
      rtn[ind] = z[tmp[0]];
      ind += 1;
    } else {
      if(my_table_counts[i] > tmp.length()){ my_table_counts[i] = tmp.length();}
      IntegerVector tmp_sample = as<IntegerVector>(sample(tmp, my_table_counts[i]));
      for (int j = 0; j < tmp_sample.length(); j++){
        rtn[ind] = z[tmp_sample[j]];
        ind += 1;
      }
    }
  }

  return rtn;
}

// [[Rcpp::export]]
NumericVector make_age_prob_vec(IntegerVector ages, NumericVector probabilities){
  
  NumericVector age_prob_vec(ages.length());
  
  for (int i = 0; i < ages.length(); ++i){
    age_prob_vec[i] = probabilities[ages[i]-1];
  }
  return age_prob_vec;
}


/*** R
x <- 1:10
y <- 1:10
id = 3

x2 <- 2:5
x3 <- c(1,2,NA)
x4 <- c(1,0,0,0,1,0)

in_group(x, y, id)
notin_group(x, y, id)
subset_vec(x2, y)
my_sum(x3)
which_cpp(x4, 1, "in")
which_cpp(x4, 1, "notin")
rep_cpp(x,y)

# match and find
ages <- sample(1:9, 20, replace = TRUE)
ages_to_contact <- 2:5
ages_to_contact2 <- c(1,2,2,4)
ids <- 1:20

m1 <- match_and_find(ages, ages_to_contact, ids)
m2 <- match_and_find(ages, ages_to_contact2, ids)

# make age prob vec
age_dist <- c(0.103, 0.116, 0.128, 0.122, 0.131, 0.145, 0.121, 0.088, 0.046)
age_prob_vec <- make_age_prob_vec(ages, age_dist)

*/
