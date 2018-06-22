#pragma once

/* A minimal reimplementation for 2D vectors & optimization */

typedef double gsl_vector;
typedef double ** matrix_t;

typedef struct {
  double (*f)(const gsl_vector *v, void *params);                            /* not used */
  void (*df)(const gsl_vector *x, void *params, gsl_vector *g);
  void (*fdf)(const gsl_vector *x, void *params, double *f, gsl_vector *g);  /* not used */
  int n;
  void *params;
} gsl_multimin_function_fdf;

typedef struct {
  const gsl_multimin_function_fdf *setup;
  gsl_vector *x, *gradient, *dx, *dgradient;
  matrix_t hessian;
  double step_size, line_tol;   /* not used */
} gsl_multimin_fdfminimizer;

gsl_vector *gsl_vector_alloc(int n);
void gsl_vector_memcpy(gsl_vector *dest, const gsl_vector *src);
double gsl_vector_get(const gsl_vector *v, int i);
void gsl_vector_set(gsl_vector *v, int i, double x);
void gsl_vector_set_zero(gsl_vector *v);
void gsl_vector_sub(gsl_vector *v1, const gsl_vector *v2);
void gsl_vector_scale(gsl_vector *v, double x);
void gsl_vector_free(gsl_vector *v);

void gsl_blas_ddot(const gsl_vector *v1, const gsl_vector *v2, double *dot);
double gsl_blas_dnrm2(const gsl_vector *v);

void *gsl_multimin_fdfminimizer_vector_bfgs2; /* dummy, always uses DFP */
#define GSL_SUCCESS 0
#define GSL_CONTINUE 1
gsl_multimin_fdfminimizer *gsl_multimin_fdfminimizer_alloc(void *type, int dim);
void gsl_multimin_fdfminimizer_set(gsl_multimin_fdfminimizer *m,
                                   const gsl_multimin_function_fdf *setup,
                                   const gsl_vector *x,
                                   double step_size,
                                   double line_tol);
int gsl_multimin_fdfminimizer_iterate(gsl_multimin_fdfminimizer *m);
int gsl_multimin_test_gradient(const gsl_vector *gradient, double tolerance);
gsl_vector *gsl_multimin_fdfminimizer_x(gsl_multimin_fdfminimizer *m);
void gsl_multimin_fdfminimizer_free(gsl_multimin_fdfminimizer *m);
