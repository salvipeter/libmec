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

static
matrix_t matrix_alloc() {
  matrix_t m = (matrix_t)malloc(sizeof(double *) * 2);
  m[0] = (double *)malloc(sizeof(double) * 2);
  m[1] = (double *)malloc(sizeof(double) * 2);
  return m;
}

static
void matrix_free(matrix_t m) {
  free(m[0]);
  free(m[1]);
  free(m);
}

static
void matrix_set_identity(matrix_t m) {
  m[0][0] = 1.0; m[0][1] = 0.0;
  m[1][0] = 0.0; m[1][1] = 1.0;
}

static
void outer_product(const gsl_vector *u, const gsl_vector *v, matrix_t m) {
  m[0][0] = u[0] * v[0]; m[0][1] = u[0] * v[1];
  m[1][0] = u[1] * v[0]; m[1][1] = u[1] * v[1];
}

static
void matrix_mul_vector(const matrix_t m, const gsl_vector *v, gsl_vector *result) {
  result[0] = m[0][0] * v[0] + m[0][1] * v[1];
  result[1] = m[1][0] * v[0] + m[1][1] * v[1];
}

static
void matrix_add(matrix_t m1, const matrix_t m2) {
  m1[0][0] += m2[0][0]; m1[0][1] += m2[0][1];
  m1[1][0] += m2[1][0]; m1[1][1] += m2[1][1];
}

static
void matrix_div(matrix_t m, double x) {
  m[0][0] /= x; m[0][1] /= x;
  m[1][0] /= x; m[1][1] /= x;
}

gsl_multimin_fdfminimizer *gsl_multimin_fdfminimizer_alloc(void *type, int dim) {
  gsl_multimin_fdfminimizer *m =
    (gsl_multimin_fdfminimizer *)malloc(sizeof(gsl_multimin_fdfminimizer));
  m->x = gsl_vector_alloc(2);
  m->gradient = gsl_vector_alloc(2);
  m->dx = gsl_vector_alloc(2);
  m->dgradient = gsl_vector_alloc(2);
  m->hessian = matrix_alloc();
  return m;
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
  matrix_set_identity(m->hessian);
  m->setup->df(m->x, m->setup->params, m->gradient);
}

static
void dfp(const gsl_vector *dx, const gsl_vector *dg, matrix_t h) {
  double denom;
  matrix_t m = matrix_alloc();
  gsl_blas_ddot(dx, dg, &denom);
  outer_product(dx, dx, m);
  matrix_div(m, denom);
  gsl_vector *v = gsl_vector_alloc(2);
  matrix_mul_vector(h, dg, v);
  matrix_add(h, m);
  outer_product(v, v, m);
  gsl_blas_ddot(dg, v, &denom);
  gsl_vector_free(v);
  matrix_div(m, -denom);
  matrix_add(h, m);
  matrix_free(m);
}

int gsl_multimin_fdfminimizer_iterate(gsl_multimin_fdfminimizer *m) {
  matrix_mul_vector(m->hessian, m->gradient, m->dx); /* dx = dx_k - dx_k+1 (!) */
  gsl_vector_sub(m->x, m->dx);
  gsl_vector_memcpy(m->dgradient, m->gradient);
  m->setup->df(m->x, m->setup->params, m->gradient);
  gsl_vector_sub(m->dgradient, m->gradient);         /* dg = dg_k - dg_k+1 (!) */
  dfp(m->dx, m->dgradient, m->hessian);
  return 0;
}

int gsl_multimin_test_gradient(const gsl_vector *gradient, double tolerance) {
  if (gsl_blas_dnrm2(gradient) > tolerance)
    return GSL_CONTINUE;
  return GSL_SUCCESS;
}

gsl_vector *gsl_multimin_fdfminimizer_x(gsl_multimin_fdfminimizer *m) {
  return m->x;
}

void gsl_multimin_fdfminimizer_free(gsl_multimin_fdfminimizer *m) {
  matrix_free(m->hessian);
  gsl_vector_free(m->dgradient);
  gsl_vector_free(m->dx);
  gsl_vector_free(m->gradient);
  gsl_vector_free(m->x);
  free(m);
}
