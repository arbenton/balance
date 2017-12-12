#include <stdio.h>
#include <math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_roots.h>

#define WSPACE_SIZE 10000
#define TOLERANCE 1e-3
#define MAX_ITER 100

typedef struct {
  double q; // order quantity
  double x; // starting inventory
  double h; // holding cost
  double p; // penalty cost 
  double* d_params; // period demand params
  double* c_params; // cumulative demand params
} sicp_period;

sicp_period* new_sicp_period(double q, double x, double h, double p, double *d_params, double *c_params) {

typedef struct {
  int n;
  double q;
  double x;
  sicp_period* periods;
} sicp_problem;

double sicp_h_holding(double d, void* params) {
  sicp_period* per = (sicp_period*) params;
  // TODO distribution specific parameter handling
  double mean = per->c_params[0];
  double std = per->c_params[1];
  printf("params: %6.4f, %6.4f\n", mean, std);
  double f = per->h*(per->q + per->x - d)*gsl_ran_gaussian_pdf(d - mean, std);
  
  return f;
}

double sicp_holding(double q, void* params) {
  sicp_problem* prob = (sicp_problem*) params; 
  prob->q = q;
  
  for (int i = 0; i < prob->n; ++i) {
    prob->periods[i].q = prob->q;
    prob->periods[i].x = prob->x;
  }
  
  gsl_integration_workspace* w = gsl_integration_workspace_alloc(WSPACE_SIZE);
  gsl_function F;
  F.function = &sicp_h_holding;
  
  double result, error;
  double total_holding = 0;
  
  for (int i = 0; i < prob->n; ++i) {
      F.params = &prob->periods[i];
      gsl_integration_qags(&F, prob->x, prob->x + prob->q, 
                           0, TOLERANCE, WSPACE_SIZE, w, 
                           &result, &error);
      printf("result : %6.4f\n", result);
      if (result < TOLERANCE) {
        break;
      }
      total_holding += result;
  }

  gsl_integration_workspace_free (w);

  return total_holding;
}

/*
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
    q_lo = gsl_root_fsolver_x_lower(s);
    q_hi = gsl_root_fsolver_x_upper(s);
    status = gsl_root_test_interval(q_lo, q_hi, 0, 0.001);
  } while (status == GSL_CONTINUE && iter < 100);
  gsl_root_fsolver_free(s);
  return root;
}
*/
int main (void) {
  
  int n = 100;
  sicp_period* periods = malloc(sizeof(sicp_period)*n);
  for (int i = 0; i < n; ++i) {
    double d_params[] = {100, 5};
    double c_params[] = {(i+1)*100, pow(25*(i+1), .5)};
    sicp_period per = {
      .q = -1, 
      .x = -1, 
      .h =  1, 
      .p =  5, 
      .d_params = d_params, 
      .c_params = c_params
    };
    periods[i] = per;
  }

  for (int i = 0; i < n; ++i) {
    printf("params: %6.4f, %6.4f\n", periods[i].d_params[0], periods[i].d_params[1]);
  }

  sicp_problem problem = {
    .n = n,
    .q = 100,
    .x =   0,
    .periods = periods
  };

  double holding = sicp_holding(100, &problem);  

  printf("holding: %6.4f\n", holding);
 
  return 0;
}
