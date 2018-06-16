#include <math.h>
#include <string.h>

#include <gsl/gsl_blas.h>
#include <gsl/gsl_multimin.h>

#define EPSILON 1e-8

typedef struct {
  int n;
  gsl_vector **points, *p, **deviations, *x;
  double *edge, *prior;
} mec_t;

mec_t *mec_init(int n, const double *points) {
  mec_t *mec = (mec_t *)malloc(sizeof(mec_t));
  mec->n = n;
  mec->points = (gsl_vector **)malloc(n * sizeof(gsl_vector *));
  for (int i = 0; i < n; ++i) {
    mec->points[i] = gsl_vector_alloc(2);
    memcpy(gsl_vector_ptr(mec->points[i], 0), points + 2 * i, 2 * sizeof(double));
  }
  mec->p = gsl_vector_alloc(2);
  mec->deviations = (gsl_vector **)malloc(n * sizeof(gsl_vector *));
  for (int i = 0; i < n; ++i)
    mec->deviations[i] = gsl_vector_alloc(2);
  mec->x = gsl_vector_alloc(2);
  mec->edge = (double *)malloc(n * sizeof(double));
  mec->prior = (double *)malloc(n * sizeof(double));
  return mec;
}

static
double distance(const gsl_vector *a, const gsl_vector *b) {
  gsl_vector *tmp = gsl_vector_alloc(2);
  gsl_vector_memcpy(tmp, a);
  gsl_vector_sub(tmp, b);
  double result = gsl_blas_dnrm2(tmp);
  gsl_vector_free(tmp);
  return result;
}

static
double f(const gsl_vector *x, void *params) {
  mec_t *mec = (mec_t *)params;
  int n = mec->n;
  double result = 0.0;
  for (int i = 0; i < n; ++i) {
    double dot;
    gsl_blas_ddot(x, mec->deviations[i], &dot);
    result += mec->prior[i] * exp(-dot);
  }
  return log(result);
}

static
void df(const gsl_vector *x, void *params, gsl_vector *g) {
  mec_t *mec = (mec_t *)params;
  int n = mec->n;
  double sum = 0.0;
  gsl_vector_set_zero(g);
  for (int i = 0; i < n; ++i) {
    double dot;
    gsl_blas_ddot(x, mec->deviations[i], &dot);
    double zi = mec->prior[i] * exp(-dot);
    sum += zi;
    for (int d = 0; d < 2; ++d)
      gsl_vector_set(g, d, gsl_vector_get(g, d) - zi * gsl_vector_get(mec->deviations[i], d));
  }
  gsl_vector_scale(g, 1.0 / sum);
}

static
void fdf(const gsl_vector *x, void *params, double *f, gsl_vector *g) {
  mec_t *mec = (mec_t *)params;
  int n = mec->n;
  *f = 0.0;
  gsl_vector_set_zero(g);
  for (int i = 0; i < n; ++i) {
    double dot;
    gsl_blas_ddot(x, mec->deviations[i], &dot);
    double zi = mec->prior[i] * exp(-dot);
    *f += zi;
    for (int d = 0; d < 2; ++d)
      gsl_vector_set(g, d, gsl_vector_get(g, d) - zi * gsl_vector_get(mec->deviations[i], d));
  }
  gsl_vector_scale(g, 1.0 / *f);
  *f = log(*f);
}

void mec_eval(mec_t *mec, const double *point, double *coordinates) {
  int n = mec->n;
  gsl_vector **points = mec->points;
  memcpy(gsl_vector_ptr(mec->p, 0), point, 2 * sizeof(double));

  /* Compute prior */
  int small = 0;
  for (int i = 0; i < n; ++i) {
    int im = (i + n - 1) % n;
    mec->edge[i] =
      distance(mec->p, points[i]) + distance(mec->p, points[im]) - distance(points[i], points[im]);
  }
  double sum = 0.0;
  for (int i = 0; i < n; ++i) {
    int ip = (i + 1) % n;
    mec->prior[i] = 1.0;
    for (int j = 0; j < n; ++j)
      if (j != i && j != ip)
        mec->prior[i] *= mec->edge[j];
    sum += mec->prior[i];
  }
  if (sum >= EPSILON)
    for (int i = 0; i < n; ++i) {
      mec->prior[i] /= sum;
      if (mec->prior[i] < EPSILON)
        ++small;
    }
  else {
    for (int i = 0; i < n; ++i) {
      int ip = (i + 1) % n;
      if (mec->edge[i] < EPSILON && mec->edge[ip] < EPSILON)
        coordinates[i] = 1.0;
      else
        coordinates[i] = 0.0;
    }
    return;
  }

  /* Special case: point on edge */
  if (small == n - 2) {
    for (int i = 0; i < n; ++i)
      if (mec->prior[i] >= EPSILON) {
        int ip = (i + 1) % n;
        if (mec->prior[ip] >= EPSILON)
          coordinates[i] = distance(mec->p, points[ip]) / distance(points[i], points[ip]);
        else {
          int im = (i + n - 1) % n;
          coordinates[i] = distance(mec->p, points[im]) / distance(points[i], points[im]);
        }
      } else
        coordinates[i] = 0.0;
    return;
  }

  /* Setup deviations */
  for (int i = 0; i < n; ++i) {
    gsl_vector_memcpy(mec->deviations[i], points[i]);
    gsl_vector_sub(mec->deviations[i], mec->p);
  }

  /* Setup the minimizer */
  gsl_vector_set_zero(mec->x);
  gsl_multimin_fdfminimizer *minimizer =
    gsl_multimin_fdfminimizer_alloc(gsl_multimin_fdfminimizer_vector_bfgs2, 2);
  gsl_multimin_function_fdf setup;
  setup.f = f;
  setup.df = df;
  setup.fdf = fdf;
  setup.n = 2;
  setup.params = (void *)mec;
  df(mec->x, mec, mec->p);                                     /* p is not used any more */
  double step_size = gsl_blas_dnrm2(mec->p);                   /* as good a guess as any */
  double line_tol = 0.1;                                       /* line search tolerance */
  gsl_multimin_fdfminimizer_set(minimizer, &setup, mec->x, step_size, line_tol);

  /* Iterate */
  int status;
  do {
    status = gsl_multimin_fdfminimizer_iterate(minimizer);
    if (status)                 /* cannot improve */
      break;
    status = gsl_multimin_test_gradient(minimizer->gradient, EPSILON);
  } while (status == GSL_CONTINUE);

  /* Compute the coordinates */
  gsl_vector *lambda = gsl_multimin_fdfminimizer_x(minimizer);
  sum = 0.0;
  for (int i = 0; i < n; ++i) {
    double dot;
    gsl_blas_ddot(lambda, mec->deviations[i], &dot);
    coordinates[i] = mec->prior[i] * exp(-dot);
    sum += coordinates[i];
  }
  for (int i = 0; i < n; ++i)
    coordinates[i] /= sum;

  /* Free the minimizer */
  gsl_multimin_fdfminimizer_free(minimizer);
}

void mec_free(mec_t *mec) {
  free(mec->prior);
  free(mec->edge);
  gsl_vector_free(mec->x);
  for (int i = 0; i < mec->n; ++i)
    gsl_vector_free(mec->deviations[i]);
  free(mec->deviations);
  gsl_vector_free(mec->p);
  for (int i = 0; i < mec->n; ++i)
    gsl_vector_free(mec->points[i]);
  free(mec->points);
  free(mec);
}
