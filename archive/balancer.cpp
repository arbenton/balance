#include <Rcpp.h>
#include "balance.hpp"
/*
// [[Rcpp::export]]
double periodic_review_cost(double quantity, double inventory,
                     Rcpp::NumericVector demand,
                     Rcpp::NumericVector holding,
                     Rcpp::NumericVector penalty) {
  return balance::periodic_review_cost(quantity, inventory, demand, holding, penalty);
}
*/
// [[Rcpp::export]]
Rcpp::NumericVector pr_cost(Rcpp::NumericVector quantity, double inventory,
                     Rcpp::NumericVector demand,
                     Rcpp::NumericVector holding,
                     Rcpp::NumericVector penalty) {
  return balance::pr_cost(quantity, inventory, demand, holding, penalty);
}

// [[Rcpp::export]]
double ls_cost(double quantity, double inventory,
                     Rcpp::NumericVector demand,
                     Rcpp::NumericVector holding,
                     Rcpp::NumericVector penalty,
                     Rcpp::NumericVector ordering) {
  return balance::ls_cost(quantity, inventory, demand, holding, penalty, ordering);
}
