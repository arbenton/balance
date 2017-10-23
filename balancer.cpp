#include <Rcpp.h>
#include "balance.hpp"

// [[Rcpp::export]]
double balancer_cost(double quantity, double inventory,
                     Rcpp::NumericVector demand,
                     Rcpp::NumericVector holding,
                     Rcpp::NumericVector penalty) {
  return balance::pr_cost<Rcpp::NumericVector>(quantity, inventory, demand, holding, penalty);
}
