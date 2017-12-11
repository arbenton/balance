#include <Rcpp.h>
#include "balance.hpp"

// [[Rcpp::export]]
double pr_cost(double quantity, double inventory,
                     Rcpp::NumericVector demand,
                     Rcpp::NumericVector holding,
                     Rcpp::NumericVector penalty) {
  return balance::pr_cost(quantity, inventory, demand, holding, penalty);
}


// [[Rcpp::export]]
double pr_root(double quantity, double inventory,
                     Rcpp::NumericVector demand,
                     Rcpp::NumericVector holding,
                     Rcpp::NumericVector penalty) {
  return balance::pr_root(quantity, inventory, demand, holding, penalty);
}

// [[Rcpp::export]]
double ls_cost(double quantity, double inventory,
                     Rcpp::NumericVector demand,
                     Rcpp::NumericVector holding,
                     Rcpp::NumericVector penalty,
                     Rcpp::NumericVector ordering) {
  return balance::ls_cost(quantity, inventory, demand, holding, penalty, ordering);
}
