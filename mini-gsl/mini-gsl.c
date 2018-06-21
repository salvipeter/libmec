#include <math.h>
#include <stdlib.h>

#include "mini-gsl.h"

gsl_vector *gsl_vector_alloc(int n) {
  return (gsl_vector *)malloc(sizeof(double) * 2);
}

void gsl_vector_memcpy(gsl_vector *dest, const gsl_vector *src) {
  dest[0] = src[0];
  dest[1] = src[1];
}

double gsl_vector_get(const gsl_vector *v, int i) {
  return v[i];
}

void gsl_vector_set(gsl_vector *v, int i, double x) {
  v[i] = x;
}

void gsl_vector_set_zero(gsl_vector *v) {
  v[0] = 0.0;
  v[1] = 0.0;
}

void gsl_vector_sub(gsl_vector *v1, const gsl_vector *v2) {
  v1[0] -= v2[0];
  v1[1] -= v2[0];
}

void gsl_vector_scale(gsl_vector *v, double x) {
  v[0] *= x;
  v[1] *= x;
}

void gsl_vector_free(gsl_vector *v) {
  free(v);
}

void gsl_blas_ddot(const gsl_vector *v1, const gsl_vector *v2, double *dot) {
  *dot = v1[0] * v2[0] + v1[1] * v2[1];
}

double gsl_blas_dnrm2(const gsl_vector *v) {
  return sqrt(v[0] * v[0] + v[1] * v[1]);
}

gsl_multimin_fdfminimizer *gsl_multimin_fdfminimizer_alloc(void *type, int dim) {
  return (gsl_multimin_fdfminimizer *)malloc(sizeof(gsl_multimin_fdfminimizer));
}

void gsl_multimin_fdfminimizer_set(gsl_multimin_fdfminimizer *m,
                                   const gsl_multimin_function_fdf *setup,
                                   const gsl_vector *x,
                                   double step_size,
                                   double line_tol) {
  m->setup = setup;
  gsl_vector_memcpy(m->x, x);
  m->step_size = step_size;
  m->line_tol = line_tol;
}

int gsl_multimin_fdfminimizer_iterate(gsl_multimin_fdfminimizer *m) {
  return 0;                     /* 1 if cannot improve */
}

int gsl_multimin_test_gradient(const gsl_vector *gradient, double tolerance) {
  return GSL_CONTINUE;          /* 0 if not */
}

gsl_vector *gsl_multimin_fdfminimizer_x(gsl_multimin_fdfminimizer *m) {
  return m->x;
}

void gsl_multimin_fdfminimizer_free(gsl_multimin_fdfminimizer *m) {
  free(m);
}
