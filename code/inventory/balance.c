#include <stdio.h>
#include <stdlib.h>
#include <math.h>

inline double pos(double x) {
  return x > 0 ? x : 0;
}

double bal_hcost(double q, double x, double h, double a, double d[], int T) {
  double cum_cost = 0;
  double cum_demand = 0;
  double n_held;
  for (int t = 0; t < T; t++) {
    cum_demand += d[t];
    n_held = pos(q - pos(cum_demand - x));
    if (n_held <= 0) {
      break;
    }
    cum_cost += pow(a, t)*h*n_held;
  }
  return cum_cost;
}

double bal_pcost(double q, double x, double p, double d) {
  double n_back = pos(d - x - q);
  return p*n_back;
}

double bal_base_stock(double q, double x, double h, double p, double a, double d[], double T) {

}
