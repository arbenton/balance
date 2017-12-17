#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>

#include <gsl/gsl_integration.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_roots.h>

#define WSPACE_SIZE 10000
#define TOLERANCE 1e-9
#define MAX_ITER 500

typedef struct {
  double q; // order quantity
  double x; // starting inventory
  double u; // capacity constraint
  double h; // holding cost
  double p; // penalty cost
  double* d_params; // period demand params
  double* c_params; // cumulative demand params
} sicp_period;

sicp_period* new_sicp_period(double u, double h, double p, double *d_params, double *c_params) {
  sicp_period* per = malloc(sizeof(sicp_period));
  per->q = -1;
  per->x = -1;
  per->u = u;
  per->h = h;
  per->p = p;
  per->d_params = malloc(sizeof(double*)*2);
  per->c_params = malloc(sizeof(double*)*2);
  memcpy(per->d_params, d_params, sizeof(double)*2);
  memcpy(per->c_params, c_params, sizeof(double)*2);
  return per;
}

int free_sicp_period(sicp_period* per) {
  free(per->d_params);
  free(per->c_params);
  free(per);
  return 1;
}

typedef struct {
  int n;
  double q;
  double x;
  sicp_period* periods;
  gsl_integration_workspace* w;
} sicp_problem;

sicp_problem* new_sicp_problem(double n, double x, sicp_period* periods) {
    sicp_problem* prob = malloc(sizeof(sicp_problem));
    prob->n = n;
    prob->q = -1;
    prob->x = x;
    prob->periods = periods;
    prob->w = gsl_integration_workspace_alloc(WSPACE_SIZE);
    return prob;
}

double sicp_h_holding(double d, void* params) {
  /* helper for gsl quadrature functions */
  sicp_period* per = (sicp_period*) params;
  // TODO distribution specific parameter handling
  double mean = per->c_params[0];
  double std = per->c_params[1];
  return per->h*(per->q + per->x - d)*gsl_ran_gaussian_pdf(d - mean, std);
}

double sicp_h_holding_der(double d, void* params) {
  /* helper for gsl quadrature functions */
  sicp_period* per = (sicp_period*) params;
  // TODO distribution specific parameter handling
  double mean = per->c_params[0];
  double std = per->c_params[1];
  return per->h*(gsl_cdf_gaussian_P(per->q + per->x - mean, std)
    - gsl_cdf_gaussian_P(per->x - mean, std));
}

double sicp_myopic_holding(double q, void* params) {
  sicp_problem* prob = (sicp_problem*) params;
  prob->q = q;
  prob->periods[0].q = prob->q;
  prob->periods[0].x = prob->x;

  gsl_function F;
  F.function = &sicp_h_holding;

  double result, error;
  F.params = &prob->periods[0];
  gsl_integration_qags(&F, prob->x, prob->x + prob->q,
                       TOLERANCE, TOLERANCE, WSPACE_SIZE, prob->w,
                       &result, &error);

  return result;
}

double sicp_horizon_holding(double q, void* params) {
  sicp_problem* prob = (sicp_problem*) params;
  prob->q = q;
  for (int i = 0; i < prob->n; ++i) {
    prob->periods[i].q = prob->q;
    prob->periods[i].x = prob->x;
  }

  gsl_function F;
  F.function = &sicp_h_holding;

  double result, error;
  double total_holding = 0;

  for (int i = 0; i < prob->n; ++i) {
      F.params = &prob->periods[i];
      gsl_integration_qags(&F, prob->x, prob->x + prob->q,
                           TOLERANCE, TOLERANCE, WSPACE_SIZE, prob->w,
                           &result, &error);
      total_holding += result;
      if (result < TOLERANCE) {
        // early stopping condition is not described in formal algorithm
        // this tolerance should be very low.
        break;
      }
  }

  return total_holding;
}


double sicp_horizon_holding_der(double q, void* params) {
  sicp_problem* prob = (sicp_problem*) params;
  prob->q = q;
  for (int i = 0; i < prob->n; ++i) {
    prob->periods[i].q = prob->q;
    prob->periods[i].x = prob->x;
  }

  double result;
  double total_holding_der = 0;

  for (int i = 0; i < prob->n; ++i) {
      double mean = prob->periods[i].c_params[0];
      double std = prob->periods[i].c_params[1];
      result = prob->periods[i].h*(
          gsl_cdf_gaussian_P(prob->periods[i].q + prob->periods[i].x - mean, std)
        - gsl_cdf_gaussian_P(prob->periods[i].x - mean, std));
      total_holding_der += result;
      if (result < TOLERANCE) {
        // early stopping condition is not described in formal algorithm
        // this tolerance should be very low.
        break;
      }
  }

  return total_holding_der;
}

double sicp_h_penalty(double d, void* params) {
  /* helper for gsl quadrature functions */
  sicp_period* per = (sicp_period*) params;
  // TODO distribution specific parameter handling
  double mean = per->c_params[0];
  double std = per->c_params[1];
  return per->p*(d - per->q - per->x)*gsl_ran_gaussian_pdf(d - mean, std);
}

double sicp_myopic_penalty(double q, void* params) {
  sicp_problem* prob = (sicp_problem*) params;
  prob->q = q;
  for (int i = 0; i < prob->n; ++i) {
    prob->periods[i].q = prob->q;
    prob->periods[i].x = prob->x;
  }

  gsl_function F;
  F.function = &sicp_h_penalty;
  F.params = &prob->periods[0];
  double result, error;
  gsl_integration_qagiu(&F, prob->x + prob->q,
                        0, TOLERANCE, WSPACE_SIZE, prob->w,
                        &result, &error);

  return result;
}

double sicp_myopic_penalty_der(double d, void* params) {
  /* helper for gsl quadrature functions */
  sicp_problem* prob = (sicp_problem*) params;
  // TODO distribution specific parameter handling
  double mean = prob->periods[0].c_params[0];
  double std = prob->periods[0].c_params[1];
  return -prob->periods[0].p*(1 - gsl_cdf_gaussian_P(prob->periods[0].x + prob->periods[0].q - mean, std));
}

double sicp_h_horizon_root(double q, void* params) {
  return sicp_horizon_holding(q, params) - sicp_myopic_penalty(q, params);
}

double sicp_h_horizon_root_der(double q, void* params) {
  return sicp_horizon_holding_der(q, params) - sicp_myopic_penalty_der(q, params);
}

void sicp_h_horizon_root_fdf(double q, void* params, double* y, double* dy) {
  *y = sicp_h_horizon_root(q, params);
  *dy = sicp_h_horizon_root_der(q, params);
}

double sicp_h_horizon_max(double q, void* params) {
  double holding = sicp_horizon_holding(q, params);
  double penalty = sicp_myopic_penalty(q, params);
  return holding > penalty ? holding : penalty;
}


double sicp_h_myopic_root(double q, void* params) {
  return sicp_myopic_holding(q, params) - sicp_myopic_penalty(q, params);
}

double sicp_h_myopic_sum(double q, void* params) {
  double holding = sicp_myopic_holding(q, params);
  double penalty = sicp_myopic_penalty(q, params);
  return holding + penalty;
}

double sicp_h_horizon_sum(double q, void* params) {
  double holding = sicp_horizon_holding(q, params);
  double penalty = sicp_myopic_penalty(q, params);
  return holding + penalty;
}

double sicp_h_horizon_sum_der(double q, void* params) {
  double holding = sicp_horizon_holding_der(q, params);
  double penalty = sicp_myopic_penalty_der(q, params);
  return holding + penalty;
}


double sicp_myopic(double q_lo, double q_hi, sicp_problem* prob) {
  int status;
  int iter = 0;
  double root;
  gsl_function F;
  F.function = &sicp_h_myopic_root;
  F.params = prob;
  const gsl_root_fsolver_type* T = gsl_root_fsolver_brent;
  gsl_root_fsolver* s = gsl_root_fsolver_alloc(T);
  gsl_root_fsolver_set(s, &F, q_lo, q_hi);
  do {
    iter++;
    status = gsl_root_fsolver_iterate(s);
    root = gsl_root_fsolver_root(s);
    q_lo = gsl_root_fsolver_x_lower(s);
    q_hi = gsl_root_fsolver_x_upper(s);
    status = gsl_root_test_interval(q_lo, q_hi, 0, TOLERANCE);
  } while (status == GSL_CONTINUE && iter < MAX_ITER);
  gsl_root_fsolver_free(s);

  return root;
}

double sicp_balance(double root, sicp_problem* prob) {
  int status;
  int iter = 0;
  gsl_function_fdf FDF;
  FDF.f = &sicp_h_horizon_root;
  FDF.df = &sicp_h_horizon_root_der;
  FDF.fdf = &sicp_h_horizon_root_fdf;
  FDF.params = prob;
  const gsl_root_fdfsolver_type *T = gsl_root_fdfsolver_newton;
  gsl_root_fdfsolver *s = gsl_root_fdfsolver_alloc (T);
  gsl_root_fdfsolver_set(s, &FDF, root);
  double x0, x;
  x = root;
  do {
      iter++;
      status = gsl_root_fdfsolver_iterate (s);
      x0 = x;
      x = gsl_root_fdfsolver_root (s);
      status = gsl_root_test_delta (x, x0, 0, TOLERANCE);
      if (status == GSL_SUCCESS)
        break;
    }
  while (status == GSL_CONTINUE && iter < MAX_ITER);

  gsl_root_fdfsolver_free (s);

  return x;
}

void multivariate2cumulative(double* mean, double** cov, double* c_mean, double* c_var, int n) {
  c_mean[0] = mean[0];
  c_var[0] = cov[0][0];
  for (int i = 1; i < n; ++i) {
    c_mean[i] = c_mean[i - 1] + mean[i];
    c_var[i] = c_var[i - 1] + cov[i][i];
    for (int j = 0; j < i; ++j)
      c_var[i] += cov[i][j] + cov[j][i];
  }
}

int main (void) {
  int n = 100;
  sicp_period* periods = malloc(sizeof(sicp_period)*n);
  for (int i = 0; i < n; ++i) {
    double d_params[] = {1000, 300};
    double c_params[] = {(i+1)*1000, pow(300*300*(i+1), .5)};
    periods[i] = *new_sicp_period(-1, 1, 1.1, d_params, c_params);
  }

  sicp_problem* problem = new_sicp_problem(n, 0, periods);

  printf("q,h,dh,l,dl,f,df,my,s,m\n");
  for (int q = 0; q < 10000; q++) {
    printf("%4d,%9.6f,%9.6f,%9.6f,%9.6f,%9.6f,%9.6f,%9.6f,%9.6f,%9.6f\n", q,
    sicp_horizon_holding(q, problem),
    sicp_horizon_holding_der(q, problem),
    sicp_myopic_penalty(q, problem),
    sicp_myopic_penalty_der(q, problem),
    sicp_h_horizon_root(q, problem),
    sicp_h_horizon_root_der(q, problem),
    sicp_h_myopic_sum(q,problem),
    sicp_h_horizon_sum(q,problem),
    sicp_h_horizon_max(q, problem));
  }


//  printf("last solution found:  %9.6f\n", opt_q);
//  printf("average milliseconds to compute: %9.6f\n", time_spent*1000/1000);

  return 1;
}
