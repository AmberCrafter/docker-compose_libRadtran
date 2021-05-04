/************************************************************************
 * $Id$
 ************************************************************************/

/*--------------------------------------------------------------------
 * This subroutine was implemented for the two-stream model of
 * Robert Buras (rodents)
 * ulrike 16.06.2010
 *--------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>

#include <math.h>
#define SWAP(a, b)                                                                                                                 \
  {                                                                                                                                \
    dum = (a);                                                                                                                     \
    (a) = (b);                                                                                                                     \
    (b) = dum;                                                                                                                     \
  }
#define TINY 1.0e-20

#include "bandec.h"

/***************************************************************************/
/* Given an n x n band diagonal matrix A with m1 subdiagonal rows and      */
/* m2 superdiagonal rows compactly stored in the array		           */
/* a[0..n-1][0..4] as described in the comment for routine banmul,         */
/* this routine constructs an LU decomposition of a rowwise permuation     */
/* of A. The upper triangular matrix replaces a, while the lower           */
/* triangular matrix is returned in al[0..n-1][0..1]. indx[0..n-1] is      */
/* an output vector which records the row permutation effected by the      */
/* partial pivoting; d is output as +-1 depending on whether the           */
/* number of row interchanges was even or odd, respectively. This          */
/* routine is used in combination with banbks to solve band-diagonal       */
/* sets of equations.                                                      */
/*                                                                         */
/* Indices changed from Fortran convention [1..n] to C convention [0..n-1] */
/***************************************************************************/

void bandec5d (double** a, long n, double** al, long indx[]) {
  long   i, j, k, l;
  double dum;

  l = 2;
  for (i = 0; i < 2; i++) { /* Rearrange the storage a bit. */
    for (j = 2 - i; j < 5; j++)
      a[i][j - l] = a[i][j];
    l--;
    for (j = 4 - l; j < 5; j++)
      a[i][j] = 0.0;
  }

  l = 2;
  for (k = 0; k < n; k++) { /* For each row ... */
    dum = a[k][0];
    i   = k;
    if (l < n)
      l++;
    for (j = k + 1; j < l; j++) { /* Find the pivot element. */
      if (fabs (a[j][0]) > fabs (dum)) {
        dum = a[j][0];
        i   = j;
      }
    }
    indx[k] = i;
    if (dum == 0.0)
      a[k][0] = TINY;
    /* Matrix is algorithmically singular, but proceed anyway
       with TINY pivot (desirable in some applications). */
    if (i != k) { /* Interchange rows. */
      for (j = 0; j < 5; j++)
        SWAP (a[k][j], a[i][j]);
    }
    for (i = k + 1; i < l; i++) { /* Do the elimination. */
      dum              = a[i][0] / a[k][0];
      al[k][i - k - 1] = dum;
      for (j = 1; j < 5; j++)
        a[i][j - 1] = a[i][j] - dum * a[k][j];
      a[i][4] = 0.0;
    }
  }
}

/***************************************************************************/
/* Given the arrays a, al, and indx as returned from bandec, and given     */
/* a right-hand side vector b[0..n-1], solves the band diagonal linear     */
/* equations Ax=b. The solution vector x overwrites b[0..n-1]. The         */
/* other input arrays are not modified, and can be left in place for       */
/* successive calls with different right-hand sides.                       */
/* Indices changed from Fortran convention [1..n] to C convention [0..n-1] */
/***************************************************************************/

void banbks5d (double** a, long n, double** al, long indx[], double b[]) {
  long   i = 0, k = 0, l = 0;
  double dum;

  l = 2;
  for (k = 0; k < n; k++) { /* Forward substitution, unscrambling the
			  permuted rows as we go. */
    if (indx[k] != k)
      SWAP (b[k], b[i]);
    if (l < n)
      l++;
    for (i = k + 1; i < l; i++)
      b[i] -= al[k][i - k - 1] * b[k];
  }

  l = 1;
  for (i = n - 1; i >= 0; i--) { /* Backsubstitution. */
    for (k = 1; k < l; k++)
      b[i] -= a[i][k] * b[k + i];
    b[i] /= a[i][0];
    if (l < 5)
      l++;
  }
}

/***********************************************************************/
/* Given an n x n band diagonal matrix A with m1 subdiagonal rows and  */
/* m2 superdiagonal rows compactly stored in the array		       */
/* a[1..n][1..m1+m2+1] as described in the comment for routine banmul, */
/* this routine constructs an LU decomposition of a rowwise permuation */
/* of A. The upper triangular matrix replaces a, while the lower       */
/* triangular matrix is returned in al[1..n][1..m1]. indx[1..n] is     */
/* an output vector which records the row permutation effected by the  */
/* partial pivoting; d is output as +-1 depending on whether the       */
/* number of row interchanges was even or odd, respectively. This      */
/* routine is used in combination with banbks to solve band-diagonal   */
/* sets of equations.                                                  */
/***********************************************************************/

void bandec (float** a, unsigned long n, int m1, int m2, float** al, unsigned long indx[], float* d) {
  unsigned long i, j, k, l;
  int           mm;
  float         dum;

  mm = m1 + m2 + 1;
  l  = m1;
  for (i = 1; i <= m1; i++) {
    for (j = m1 + 2 - i; j <= mm; j++)
      a[i][j - l] = a[i][j];
    l--;
    for (j = mm - l; j <= mm; j++)
      a[i][j] = 0.0;
  }
  *d = 1.0;
  l  = m1;
  for (k = 1; k <= n; k++) {
    dum = a[k][1];
    i   = k;
    if (l < n)
      l++;
    for (j = k + 1; j <= l; j++) {
      if (fabs (a[j][1]) > fabs (dum)) {
        dum = a[j][1];
        i   = j;
      }
    }
    indx[k] = i;
    if (dum == 0.0)
      a[k][1] = TINY;
    if (i != k) {
      *d = -(*d);
      for (j = 1; j <= mm; j++)
        SWAP (a[k][j], a[i][j])
    }
    for (i = k + 1; i <= l; i++) {
      dum          = a[i][1] / a[k][1];
      al[k][i - k] = dum;
      for (j = 2; j <= mm; j++)
        a[i][j - 1] = a[i][j] - dum * a[k][j];
      a[i][mm] = 0.0;
    }
  }
}

/***********************************************************************/
/* Given the arrays a, al, and indx as returned from bandec, and given */
/* a right-hand side vector b[1..n], solves the band diagonal linear   */
/* equations Ax=b. The solution vector x overwrites b[1..n]. The       */
/* other input arrays are not modified, and can be left in place for   */
/* successive calls with different right-hand sides.                   */
/***********************************************************************/

void banbks (float** a, unsigned long n, int m1, int m2, float** al, unsigned long indx[], float b[]) {
  unsigned long i, k, l;
  int           mm;
  float         dum;

  mm = m1 + m2 + 1;
  l  = m1;
  for (k = 1; k <= n; k++) {
    i = indx[k];
    if (i != k)
      SWAP (b[k], b[i])
    if (l < n)
      l++;
    for (i = k + 1; i <= l; i++)
      b[i] -= al[k][i - k] * b[k];
  }
  l = 1;
  for (i = n; i >= 1; i--) {
    dum = b[i];
    for (k = 2; k <= l; k++)
      dum -= a[i][k] * b[k + i - 1];
    b[i] = dum / a[i][1];
    if (l < mm)
      l++;
  }
}

void matrix_to_band_form (double** M, long n, int m1, int m2, double** Mb) {
  long k, l;
  for (k = 0; k < n; k++) {
    for (l = -m1; l <= m2; l++) {
      if (k + l >= 0 && k + l < n) {
        Mb[k][m1 + l] = M[k][k + l];
      }
    }
  }
}

void bandec_M (double** au, long n, int m1, int m2, double** al, long indx[]) {
  long   i, j, k, l;
  int    mm;
  double dum;

  mm = m1 + m2 + 1;
  l  = m1;
  for (i = 0; i < m1; i++) {
    for (j = m1 - i; j < mm; j++)
      au[i][j - l] = au[i][j];
    l--;
    for (j = m1 + m2 - l; j < mm; j++)
      au[i][j] = 0.0;
  }

  l = m1;
  for (k = 0; k < n; k++) {
    dum = au[k][0];
    i   = k;
    if (l < n)
      l++;
    for (j = k + 1; j < l; j++) {
      if (fabs (au[j][0]) > fabs (dum)) {
        dum = au[j][0];
        i   = j;
      }
    }
    indx[k] = i;
    if (dum == 0.0)
      au[k][0] = TINY;
    if (i != k) {
      for (j = 0; j < mm; j++)
        SWAP (au[k][j], au[i][j]);
    }
    for (i = k + 1; i < l; i++) {
      dum              = au[i][0] / au[k][0];
      al[k][i - k - 1] = dum;
      for (j = 1; j < mm; j++)
        au[i][j - 1] = au[i][j] - dum * au[k][j];
      au[i][m1 + m2] = 0.0;
    }
  }
}

void bandec11d (double** a, long n, double** al, long indx[]) {
  long   i, j, k, l;
  double dum;

  l = 5;
  for (i = 0; i < 5; i++) { /* Rearrange the storage a bit. */
    for (j = 5 - i; j < 11; j++)
      a[i][j - l] = a[i][j];
    l--;
    for (j = 10 - l; j < 11; j++)
      a[i][j] = 0.0;
  }

  l = 5;
  for (k = 0; k < n; k++) { /* For each row ... */
    dum = a[k][0];
    i   = k;
    if (l < n)
      l++;
    for (j = k + 1; j < l; j++) { /* Find the pivot element. */
      if (fabs (a[j][0]) > fabs (dum)) {
        dum = a[j][0];
        i   = j;
      }
    }
    indx[k] = i;
    if (dum == 0.0)
      a[k][0] = TINY;
    /* Matrix is algorithmically singular, but proceed anyway
       with TINY pivot (desirable in some applications). */
    if (i != k) { /* Interchange rows. */
      for (j = 0; j < 11; j++)
        SWAP (a[k][j], a[i][j]);
    }
    for (i = k + 1; i < l; i++) { /* Do the elimination. */
      dum              = a[i][0] / a[k][0];
      al[k][i - k - 1] = dum;
      for (j = 1; j < 11; j++)
        a[i][j - 1] = a[i][j] - dum * a[k][j];
      a[i][10] = 0.0;
    }
  }
}
void banbks_M (double** a, long n, int m1, int m2, double** al, long indx[], double b[]) {
  long  i = 0, k = 0, l = 0;
  int   mm;
  float dum;
  mm = m1 + m2 + 1;

  l = m1;
  for (k = 0; k < n; k++) {
    i = indx[k];
    if (i != k)
      SWAP (b[k], b[i]);
    if (l < n)
      l++;
    for (i = k + 1; i < l; i++)
      b[i] -= al[k][i - k - 1] * b[k];
  }
  l = 1;
  for (i = n - 1; i >= 0; i--) {
    dum = b[i];
    for (k = 1; k < l; k++)
      dum -= a[i][k] * b[k + i];
    b[i] = dum / a[i][0];
    if (l < mm)
      l++;
  }
}

void banbks11d (double** a, long n, double** al, long indx[], double b[]) {
  long   i = 0, k = 0, l = 0;
  double dum;

  l = 5;
  for (k = 0; k < n; k++) {  //Forward substitution, unscrambling the permuted rows as we go.
    if (indx[k] != k)
      SWAP (b[k], b[i]);
    if (l < n)
      l++;
    for (i = k + 1; i < l; i++)
      b[i] -= al[k][i - k - 1] * b[k];
  }

  l = 1;
  for (i = n - 1; i >= 0; i--) {  // Backsubstitution.
    for (k = 1; k < l; k++)
      b[i] -= a[i][k] * b[k + i];
    b[i] /= a[i][0];
    if (l < 11)
      l++;
  }
}

#undef SWAP
#undef TINY
