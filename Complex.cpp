#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "Initialize.h"

#define PI 3.14159265358
#define EXP 2.7182818285

/* --------------------------------------
        Complex number operations
---------------------------------------*/

// Convert real and imaginary parts of complex number to amplitude and phase 
complex_t ReIm2AmpPh(complex_t c)
{ 
  c.amp = sqrt(c.re*c.re+c.im*c.im);
  if (c.re == 0 && c.im == 0)
  {
    c.ph = 0.;
  }
  else if (c.re == 0)
  {
    c.ph = PI/2.;
  }
  else
  {
    c.ph = atan(c.im/c.re);
    if (c.re < 0.0)
    {
      if (c.im >= 0) { c.ph += PI; }
      else { c.ph -= PI; }
    }
  } 
  return c;
}

// Find product of 2 complex numbers
complex_t prod(complex_t c1, complex_t c2)
{
  complex_t prod;
  prod.re = c1.re*c2.re - c1.im*c2.im;
  prod.im = c1.re*c2.im + c1.im*c2.re;
  prod = ReIm2AmpPh(prod);
  return prod;
}

// Add 2 complex numbers
complex_t add(complex_t c1, complex_t c2)
{
  complex_t add;
  add.re = c1.re + c2.re;
  add.im = c1.im + c2.im;
  add = ReIm2AmpPh(add);
  return add;
}

// Multiply a complex number by a real number
complex_t prod(complex_t c, double d)
{
  complex_t prod;
  prod.re = c.re*d;
  prod.im = c.im*d;
  prod = ReIm2AmpPh(prod);
  return prod;
}

// Add a complex number and a real number
complex_t add(complex_t c, double d)
{
  complex_t add;
  add.re = c.re + d;
  add.im = c.im;
  add = ReIm2AmpPh(add);
  return add;
}

// Find the inverse of a complex number
complex_t inv(complex_t c)
{
  complex_t i;
  i.re = c.re / (c.re*c.re + c.im*c.im);
  i.im = -c.im / (c.re*c.re + c.im*c.im);
  i = ReIm2AmpPh(i);
  return i;
}

// Find negative of a complex number
complex_t neg(complex_t c)
{
  complex_t n;
  n.re = -c.re;
  n.im = -c.im;
  n = ReIm2AmpPh(n);
  return n;
}

// Find the complex conjugate
complex_t conj(complex_t c)
{
  complex_t c1;
  c1.re = c.re;
  c1.im = -c.im;
  c1 = ReIm2AmpPh(c1);
  return c1;
}

// Find the square root of a complex number
complex_t csqrt(complex_t c)
{
  double r, phi;
  complex_t cs; 
  cs.amp = sqrt(c.amp);
  cs.ph = c.ph/2.0;
  cs.re = cs.amp*cos(cs.ph);
  cs.im = cs.amp*sin(cs.ph);
  return cs;
}

// Find the exponential of a complex number
complex_t cexp(complex_t c)
{
  complex_t ce;
  ce.re = pow(EXP, c.re) * cos(c.im);
  ce.im = pow(EXP, c.re) * sin(c.im);
  ce = ReIm2AmpPh(ce);
  return ce;
}
