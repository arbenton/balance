#include <stdio.h>
#include <math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_roots.h>

#define WSPACE_SIZE 10000
#define TOLERANCE 1e-6

typedef struct {
  int horizon;
  double order_quant;
  double starting_inv;
  double* holding_cost;
  double* penalty_cost;
  double* demand_mean;
  double* demand_cov;
} sicp_problem;

double scip_h_holding(double d, void* params) {
  sicp_period* p = (sicp_period*) params; // passed as void* for GSL compliance
  double f =  p->holding_cost*(p->order_quant + p->starting_inv - d)*
              gsl_ran_gaussian_pdf(d - cum_mean, pow(cum_var, .5));
  return f;
}

double scip_holding(double q, void* params) {
  sicp_period* p = (sicp_period*) params; // passed as void* for GSL compliance
  gsl_integration_workspace* w = gsl_integration_workspace_alloc(WSPACE_SIZE);
  p->order_quant = q;
  double delta, result, error;
  double cum_mean = 0;
  double cum_var = 0;
  gsl_function F;
  F.function = &scip_h_holding;
  for (int i = 0; i < p->horizon; ++i) {
      cum_mean += p->demand_mean[i];
      cum_var += pow(p->demand_std[i], 2);
      double p[] = {cum_mean, pow(cum_var, .5)};
      F.params = p;
      gsl_integration_qags(&F, p->starting_inv, p->starting_inv + p->order_quant,
                        0, p->tolerance, 1000, w, &result, &error);
  gsl_integration_workspace_free (w);
  return result;
}

double scip_h_penalty(double d, void* params) {
  sicp_period* p = (sicp_period*) params; // passed as void* for GSL compliance
  double f = p->penalty_cost*(d - p->order_quant - p->starting_inv)*
             gsl_ran_gaussian_pdf(d - p->demand_mean[0], p->demand_std[0]);
  return f;
}

double scip_penalty(double q, void* params) {
  sicp_period* p = (sicp_period*) params; // passed as void* for GSL compliance
  gsl_integration_workspace* w = gsl_integration_workspace_alloc(WSPACE_SIZE);
  p->order_quant = q;
  gsl_function F;
  F.function = &scip_h_penalty;
  F.params = p;
  double result, error;
  gsl_integration_qagiu(&F, p->starting_inv + p->order_quant,
                        0, p->tolerance, 1000, w, &result, &error);
  gsl_integration_workspace_free (w);
  return result;
}

double scip_h_loss(double q, void* params) {
  return scip_holding(q, params) - scip_penalty(q, params);
}

double scip_balance(double q_lo, double q_hi, sicp_period* p) {
  int status;
  int iter = 0;
  double root;
  gsl_function F;
  F.function = &scip_h_loss;
  F.params = p;
  const gsl_root_fsolver_type* T = gsl_root_fsolver_brent;
  gsl_root_fsolver* s = gsl_root_fsolver_alloc(T);
  gsl_root_fsolver_set(s, &F, q_lo, q_hi);
  do {
    iter++;
    status = gsl_root_fsolver_iterate(s);
    root = gsl_root_fsolver_root(s);
    printf("%f\n", root);
    q_lo = gsl_root_fsolver_x_lower(s);
    q_hi = gsl_root_fsolver_x_upper(s);
    status = gsl_root_test_interval(q_lo, q_hi, 0, 0.001);
  } while (status == GSL_CONTINUE && iter < 100);
  gsl_root_fsolver_free(s);
  return root;
}

int main (void) {
  double demand_mean[100], demand_std[100];
  for (int i = 0; i < 100; ++i) {
    demand_mean[i] = 200 + i*10;
    demand_std[i] = 15;
  }

  sicp_period params = {100, 100, 0, 1, 5, demand_mean, demand_std, 1e-6};

  double quantity = scip_balance(0, 500, (void*)&params);
  printf("%6.4f \n", quantity);
  return 0;
}
