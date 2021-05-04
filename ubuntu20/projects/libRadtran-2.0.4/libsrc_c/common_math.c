#include <stdio.h>
#include <float.h>
#include <math.h>

#include "common_math.h"
#include "errors.h"

/***********************************************************************************/
/* Functions: v_mult, v_sub, vn_mult                                      @62_30i@ */
/* Description:                                                                    */
/*   short cuts for vector*vector -> number, vector - vector -> vector,            */
/*                  vector*number -> vector                                        */
/*                                                                                 */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author: Robert Buras                                                            */
/*                                                                        @i62_30@ */
/***********************************************************************************/
int v_cross_product_norm (const double a[3], const double b[3], double c[3]) {
  double n = 0.0;
  int    i = 0;

  c[0] = a[1] * b[2] - a[2] * b[1];
  c[1] = a[2] * b[0] - a[0] * b[2];
  c[2] = a[0] * b[1] - a[1] * b[0];

  /* Normalization*/
  for (i = 0; i < 3; i++)
    n += c[i] * c[i];

  /* a and b are parallel */
  if (n == 0.0)
    return 1;

  n = sqrt (n);

  for (i = 0; i < 3; i++)
    c[i] /= n;

  return 0;
}

/* const problem in C for a, see header file for details! */
void mat_v_mult (size_t N, double** a, const double* b, double* c) {
  int i = 0, j = 0;

  for (i = 0; i < N; i++) {
    c[i] = 0.0;
    for (j = 0; j < N; j++)
      c[i] += a[i][j] * b[j];
  }
}

void v_mult (const double a[3], const double b[3], double* c) {
  *c = vec3_dot_product (a, b);
}

void v_sub (const double a[3], const double b[3], double c[3]) {
  c[0] = a[0] - b[0];
  c[1] = a[1] - b[1];
  c[2] = a[2] - b[2];
}

void v_add (const double a[3], const double b[3], double c[3]) {
  c[0] = a[0] + b[0];
  c[1] = a[1] + b[1];
  c[2] = a[2] + b[2];
}

void vn_mult (const double a[3], const double* b, double c[3]) {
  c[0] = a[0] * *b;
  c[1] = a[1] * *b;
  c[2] = a[2] * *b;
}

// Cross product, right hand rule, a(thumb), b(pointing finger)
void v_cross_product (const double a[3], const double b[3], double c[3]) {
  c[0] = a[1] * b[2] - a[2] * b[1];
  c[1] = a[2] * b[0] - a[0] * b[2];
  c[2] = a[0] * b[1] - a[1] * b[0];
}

int v_mult_mu (const double a[3], const double b[3], double* c) {
  *c = vec3_dot_product (a, b);

  if (*c > 1.0) {
    if (*c > 1.0 + DBL_EPSILON)
      *c = NOT_A_NUMBER;
    *c = 1.0;
  }

  if (*c < -1.0) {
    if (*c < -1.0 - DBL_EPSILON)
      *c = NOT_A_NUMBER;
    *c = -1.0;
  }

  if (isnan (*c)) {
    fprintf (stderr, "Error! v_mult_mu([%f %f %f] [%f %f %f]) has unphysical result %e \n", a[0], a[1], a[2], b[0], b[1], b[2], *c);
    CHKERR (-1);
  }
  return 0;
}

void v_neg (double* a) {
  a[0] = -a[0];
  a[1] = -a[1];
  a[2] = -a[2];
}

double cosd (double angle) {
  return cos (angle * M_PI / 180.0);
}

double sind (double angle) {
  return sin (angle * M_PI / 180.0);
}

double tand (double angle) {
  // tan is extremely slow in current gcc version, CE 180824
  // return tan ( angle * PI / 180.0 );
  return sin (angle * M_PI / 180.0) / cos (angle * M_PI / 180.0);
}

double acosd (double rad) {
  return acos (rad) * 180.0 / M_PI;
}

double asind (double rad) {
  return asin (rad) * 180.0 / M_PI;
}

double atand (double rad) {
  return atan (rad) * 180.0 / M_PI;
}

/* return -1 or +1 according to the sign of floating point variable x, 0 if the value is zero*/
int sign (double x) {
  if (x > 0)
    return 1;
  if (x < 0)
    return -1;
  return 0;
}

/* applies out[i] = fabs(inp[i]) on all elements of inp[N] */
void vec_fabs (const size_t N, const double* inp, double* out) {
  for (size_t i = 0; i < N; ++i)
    out[i] = fabs (inp[i]);
}

/* Dot product of two vecs a, b with size N */
int vec_dot_product (const size_t N, const double* a, const double* b, double* c) {
  *c = 0;
  for (size_t i = 0; i < N; ++i)
    *c += a[i] * b[i];
  return 0;
}

double vec3_dot_product (const double a[3], const double b[3]) {
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

/* Compute 2-Norm of vector */
int vec_norm2 (const size_t N, const double* v, double* norm) {
  int ierr = vec_dot_product (N, v, v, norm);
  CHKERR (ierr);
  CHKERR (*norm < 0);
  *norm = sqrt (*norm);
  return 0;
}
double vec3_norm2 (const double v[3]) {
  return sqrt (vec3_dot_product (v, v));
}

/* Computes y = alpha y */
int vec3_scale (const double alpha, double y[3]) {
  if (alpha == 1)
    return 0;
  y[0] = alpha * y[0];
  y[1] = alpha * y[1];
  y[2] = alpha * y[2];
  return 0;
}
int vec_scale (size_t N, const double alpha, double* y) {
  if (alpha == 1)
    return 0;
  for (size_t i = 0; i < N; ++i)
    y[i] = alpha * y[i];
  return 0;
}

/* Computes y = alpha x + y */
int vec_AXPY (size_t N, const double alpha, const double* x, double* y) {
  if (alpha == 0)
    return 0;
  for (size_t i = 0; i < N; ++i)
    y[i] = alpha * x[i] + y[i];
  return 0;
}

/* Computes y = alpha x + beta y */
int vec_AXPBY (size_t N, const double alpha, const double* x, const double beta, double* y) {
  if (beta == 1)
    return vec_AXPY (N, alpha, x, y);
  if (alpha == 0)
    return vec_scale (N, beta, y);
  for (size_t i = 0; i < N; ++i)
    y[i] = alpha * x[i] + beta * y[i];
  return 0;
}

/* Computes z = alpha x + beta y + gamma z
     The implementation is optimized for alpha of 1 and gamma of 1 or 0
*/
int vec_AXPBYPCZ (size_t        N,
                  const double  alpha,
                  const double* x,
                  const double  beta,
                  const double* y,
                  const double  gamma,
                  double*       z) {
  if (alpha == 1) {
    if (gamma == 0) {
      for (size_t i = 0; i < N; ++i)
        z[i] = x[i] + beta * y[i];
    } else if (gamma == 1) {
      for (size_t i = 0; i < N; ++i)
        z[i] = x[i] + beta * y[i] + z[i];
    } else {
      for (size_t i = 0; i < N; ++i)
        z[i] = x[i] + beta * y[i] + gamma * z[i];
    }
  } else {
    if (gamma == 0) {
      for (size_t i = 0; i < N; ++i)
        z[i] = alpha * x[i] + beta * y[i];
    } else if (gamma == 1) {
      for (size_t i = 0; i < N; ++i)
        z[i] = alpha * x[i] + beta * y[i] + z[i];
    } else {
      for (size_t i = 0; i < N; ++i)
        z[i] = alpha * x[i] + beta * y[i] + gamma * z[i];
    }
  }
  return 0;
}

/* Compares 2 floats if they are equal within a specified range, if precision==0, use the machine precision */
int approx (const double a, const double b, const double precision) {
  double epsilon = precision == 0 ? DBL_EPSILON : precision;
  return fabs (a - b) <= ((fabs (a) < fabs (b) ? fabs (b) : fabs (a)) * epsilon);
}

/* Searches for index of max value in array
    - use ints to allow for negative search indices
    - returns error 1 if no maximum is found or -1 if ambiguous results
*/
int max_loc (const double* arr, int start, int end, int* ierr) {
  *ierr      = 1;
  int maxloc = start;
  if (arr[start] >= arr[end])
    *ierr = 0;
  for (int i = start + 1; i <= end; ++i) {
    if (arr[i] >= arr[maxloc]) {
      *ierr = 0;
      if (arr[i] == arr[maxloc])
        *ierr = -1;
      maxloc = i;
    }
  }
  return maxloc;
}

/* Inplace reverse array entries */
int reverse_1d_f4 (size_t N, float* arr) {
  // dont do stack allocations if temp array is bigger than approx 1MB
  if ((N / 2) * sizeof (float) > 1e6)
    return 1;

  float  tmp[N / 2];
  size_t i;
  for (i = 0; i < N / 2; i++)
    tmp[i] = arr[i];

  for (i = 0; i < N / 2; i++)
    arr[i] = arr[N - i - 1];

  for (i = 0; i < N / 2; i++)
    arr[N - i - 1] = tmp[i];

  return 0;
}
int reverse_1d_f8 (size_t N, double* arr) {
  // dont do stack allocations if temp array is bigger than approx 1MB
  if ((N / 2) * sizeof (double) > 1e6)
    return 1;

  double tmp[N / 2];
  size_t i;
  for (i = 0; i < N / 2; i++)
    tmp[i] = arr[i];

  for (i = 0; i < N / 2; i++)
    arr[i] = arr[N - i - 1];

  for (i = 0; i < N / 2; i++)
    arr[N - i - 1] = tmp[i];

  return 0;
}

// count all true/false values in arr[N]
size_t vec_count_true (const size_t N, const int* arr) {
  size_t cnt = 0;
  for (size_t i = 0; i < N; ++i)
    if (arr[i])
      ++cnt;
  return cnt;
}
size_t vec_count_false (const size_t N, const int* arr) {
  size_t cnt = 0;
  for (size_t i = 0; i < N; ++i)
    if (!arr[i])
      ++cnt;
  return cnt;
}

// Determine distance/edge length between two points with N dimensions
double vec_distance_between_two_points (const size_t N, const double* p1, const double* p2) {
  double distance = 0;
  for (size_t i = 0; i < N; ++i)
    distance += (p2[i] - p1[i]) * (p2[i] - p1[i]);
  return sqrt (distance);
}

// Use Herons Formula to determine the area of a triangle given the 3 edge lengths
double triangle_area_by_edgelengths (const double e1, const double e2, const double e3) {
  const double p = (e1 + e2 + e3) / 2;
  return sqrt (p * (p - e1) * (p - e2) * (p - e3));
}

double triangle_area_by_vertices (const double p1[3], const double p2[3], const double p3[3]) {
  const double e1 = vec_distance_between_two_points (3, p2, p3);
  const double e2 = vec_distance_between_two_points (3, p3, p1);
  const double e3 = vec_distance_between_two_points (3, p1, p2);
  return triangle_area_by_edgelengths (e1, e2, e3);
}
