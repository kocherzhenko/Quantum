// Matrix and Vector Initialization Routines
// Adapted from "Numerical Recipies in C"

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stddef.h>

#define NR_END 1
#define FREE_ARG char*

typedef struct { double amp; double ph; double re; double im; } complex_t;

// Allocate a double matrix with subscript range m[nr1..nrh][nc1..nch]
double **dmatrix(long nr1, long nrh, long nc1, long nch)
{
  long i, nrow=nrh-nr1+1, ncol = nch-nc1+1;
  double **m;

  m = (double **) malloc((size_t)((nrow+NR_END)*sizeof(double)));
  if (!m) { printf("Allocation failed in dmatrix\n"); }
  m += NR_END; m -= nr1;
  m[nr1] = (double *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(double)));
  if (!m[nr1]) { printf("Allocation failed in dmatrix\n"); }
  m[nr1] += NR_END; m[nr1] -= nc1;
  for (i=nr1+1; i<=nrh; i++) m[i]=m[i-1]+ncol;

  return m;
}

// Free a double matrix with subscript range m[nr1..nrh][nc1..nch]
void free_dmatrix(double **m, long nr1, long nrh, long nc1, long nch)
{
  free((FREE_ARG) (m[nr1]+nc1-NR_END));
  free((FREE_ARG) (m+nr1-NR_END));
}

// Allocate a double vector with subscript range v[n1..nh]
double *dvector(long n1, long nh)
{
  double *v;

  v = (double *) malloc((size_t) ((nh-n1+1+NR_END)*sizeof(double)));
  if (!v) printf("Allocation failure in dvector\n");
  return v-n1+NR_END;
}

// Free a double vector with subscript range v[n1..nh]
void free_dvector(double *v, long n1, long nh)
{
  free((FREE_ARG) (v+n1-NR_END));
}

// Allocate an integer vector with subscript range v[n1..nh]
int *ivector(long n1, long nh)
{
  int *v;

  v = (int *) malloc((size_t) ((nh-n1+1+NR_END)*sizeof(int)));
  if (!v) printf("Allocation failure in ivector\n");
  return v-n1+NR_END;
}

// Free an integer vector with subscript range v[n1..nh]
void free_ivector(int *v, long n1, long nh)
{
  free((FREE_ARG) (v+n1-NR_END));
}

// Allocate a complex_t matrix with subscript range m[nr1..nrh][nc1..nch]
complex_t **cmatrix(long nr1, long nrh, long nc1, long nch)
{
  long i, nrow=nrh-nr1+1, ncol = nch-nc1+1;
  complex_t **m;

  m = (complex_t **) malloc((nrow+NR_END)*sizeof(complex_t));
  if (!m) { printf("Allocation failed in cmatrix\n"); }
  m += NR_END; m -= nr1;
  m[nr1] = (complex_t *) malloc((nrow*ncol+NR_END)*sizeof(complex_t));
  if (!m[nr1]) { printf("Allocation failed in cmatrix\n"); }
  m[nr1] += NR_END; m[nr1] -= nc1;
  for (i=nr1+1; i<=nrh; i++) m[i]=m[i-1]+ncol;

  return m;
}

// Free a complex_t matrix with subscript range m[nr1..nrh][nc1..nch]
void free_cmatrix(complex_t **m, long nr1, long nrh, long nc1, long nch)
{
  free((FREE_ARG) (m[nr1]+nc1-NR_END));
  free((FREE_ARG) (m+nr1-NR_END));
}

// Allocate a complex_t vector with subscript range v[n1..nh]
complex_t *cvector(long n1, long nh)
{
  complex_t *v;

  v = (complex_t *) malloc((size_t) ((nh-n1+1+NR_END)*sizeof(complex_t)));
  if (!v) printf("Allocation failure in cvector\n");
  return v-n1+NR_END;
}

// Free a complex_t vector with subscript range v[n1..nh]
void free_cvector(complex_t *v, long n1, long nh)
{
  free((FREE_ARG) (v+n1-NR_END));
}
