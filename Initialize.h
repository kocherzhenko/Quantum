// Matrix and Vector Initialization Routines
// Adapted from "Numerical Recipies in C"

typedef struct { double amp; double ph; double re; double im; } complex_t;

// Allocate a double matrix with subscript range m[nr1..nrh][nc1..nch]
double **dmatrix(long nr1, long nrh, long nc1, long nch);

// Free a double matrix with subscript range m[nr1..nrh][nc1..nch]
void free_dmatrix(double **m, long nr1, long nrh, long nc1, long nch);

// Allocate a double vector with subscript range v[n1..nh]
double *dvector(long n1, long nh);

// Free a double vector with subscript range v[n1..nh]
void free_dvector(double *v, long n1, long nh);

// Allocate an integer vector with subscript range v[n1..nh]
int *ivector(long n1, long nh);

// Free an integer vector with subscript range v[n1..nh]
void free_ivector(int *v, long n1, long nh);

// Allocate a complex_t matrix with subscript range m[nr1..nrh][nc1..nch]
complex_t **cmatrix(long nr1, long nrh, long nc1, long nch);

// Free a complex_t matrix with subscript range m[nr1..nrh][nc1..nch]
void free_cmatrix(complex_t **m, long nr1, long nrh, long nc1, long nch);

// Allocate a complex_t vector with subscript range v[n1..nh]
complex_t *cvector(long n1, long nh);

// Free a complex_t vector with subscript range v[n1..nh]
void free_cvector(complex_t *v, long n1, long nh);
