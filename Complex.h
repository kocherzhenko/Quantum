/* --------------------------------------
        Complex number operations
---------------------------------------*/
#define PI 3.14159265358

// Convert real and imaginary parts of complex number to amplitude and phase 
complex_t ReIm2AmpPh(complex_t c);

// Find product of 2 complex numbers
complex_t prod(complex_t c1, complex_t c2);

// Add 2 complex numbers
complex_t add(complex_t c1, complex_t c2);

// Multiply a complex number by a real number
complex_t prod(complex_t c, double d);

// Add a complex number and a real number
complex_t add(complex_t c, double d);

// Find the inverse of a complex number
complex_t inv(complex_t c);

// Find negative of a complex number
complex_t neg(complex_t c);

// Find the complex conjugate
complex_t conj(complex_t c);

// Find the square root of a complex number
complex_t csqrt(complex_t c);

// Find the exponential of a complex number
complex_t cexp(complex_t c);
