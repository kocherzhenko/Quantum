#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include "Initialize.h"
#include "Complex.h"

#define PI 3.14159265358

/*---------------------------------
         Set up gates
---------------------------------*/
//
// 1-qubit gates 
//

// Identity
complex_t** Identity(complex_t **I)
{
  int i, j;

  I[0][0].re = I[1][1].re = 1.;
  I[0][1].re = I[1][0].re = 0.;
  for (i=0; i<=1; i++)
  {
    for (j=0; j<=1; j++)
    {
      I[i][j].im = 0.;
      I[i][j] = ReIm2AmpPh(I[i][j]);
    }
  }
  return I;
}

// Hadamard
complex_t** Hadamard(complex_t **Had)
{
  int i, j;

  Had[0][0].re = Had[0][1].re = Had[1][0].re = 1./sqrt(2.);
  Had[1][1].re = -1./sqrt(2.);
  for (i=0; i<=1; i++)
  {
    for (j=0; j<=1; j++)
    {
      Had[i][j].im = 0.;
      Had[i][j] = ReIm2AmpPh(Had[i][j]);
    }
  }
  return Had;
}

// PauliX
complex_t** PauliX(complex_t **PauX)
{
  int i, j;

  PauX[0][0].re = PauX[1][1].re = 0.;
  PauX[0][1].re = PauX[1][0].re = 1.;
  for (i=0; i<=1; i++)
  {
    for (j=0; j<=1; j++)
    {
      PauX[i][j].im = 0.;
      PauX[i][j] = ReIm2AmpPh(PauX[i][j]);
    }
  }
  return PauX;
}

// PauliY
complex_t** PauliY(complex_t **PauY)
{
  int i, j;

  PauY[0][0].im = PauY[1][1].im = 0.;
  PauY[0][1].im = -1.;
  PauY[1][0].im = 1.;
  for (i=0; i<=1; i++)
  {
    for (j=0; j<=1; j++)
    {
      PauY[i][j].re = 0.;
      PauY[i][j] = ReIm2AmpPh(PauY[i][j]);
    }
  }
  return PauY;
}

// PauliZ
complex_t** PauliZ(complex_t **PauZ)
{
  int i, j;

  PauZ[0][0].re = 1.;
  PauZ[1][1].re = -1.;
  PauZ[0][1].re = PauZ[1][0].re = 0.;
  for (i=0; i<=1; i++)
  {
    for (j=0; j<=1; j++)
    {
      PauZ[i][j].im = 0.;
      PauZ[i][j] = ReIm2AmpPh(PauZ[i][j]);
    }
  }
  return PauZ;
}

// RotX
complex_t** RX(complex_t **RotX, double theta)
{
  int i, j;

  theta *= (PI/180.);
  RotX[0][0].re = RotX[1][1].re = cos(theta/2.);
  RotX[0][1].im = RotX[1][0].im = -sin(theta/2.);
  RotX[0][0].im = RotX[1][1].im = 0.;
  RotX[0][1].re = RotX[1][0].re = 0.;
  for (i=0; i<=1; i++)
  {
    for (j=0; j<=1; j++)
    {
      RotX[i][j] = ReIm2AmpPh(RotX[i][j]);
    }
  }
  return RotX;
}

// RotY
complex_t** RY(complex_t **RotY, double theta)
{
  int i, j;

  theta *= (PI/180.);
  RotY[0][0].re = RotY[1][1].re = cos(theta/2.);
  RotY[0][1].re = -sin(theta/2.);
  RotY[1][0].re = sin(theta/2.);
  for (i=0; i<=1; i++)
  {
    for (j=0; j<=1; j++)
    { 
      RotY[i][j].im = 0.;
      RotY[i][j] = ReIm2AmpPh(RotY[i][j]);
    }
  }
  return RotY;
}

// RotZ
complex_t** RZ(complex_t **RotZ, double theta)
{
  int i, j;

  theta *= (PI/180.);
  RotZ[0][0].re = RotZ[1][1].re = cos(theta/2.);
  RotZ[0][0].im = -sin(theta/2.);
  RotZ[1][1].im = sin(theta/2.);
  RotZ[0][1].re = RotZ[0][1].im = RotZ[1][0].re = RotZ[1][0].im = 0.;
  for (i=0; i<=1; i++)
  {
    for (j=0; j<=1; j++)
    { 
      RotZ[i][j] = ReIm2AmpPh(RotZ[i][j]);
    }
  }
  return RotZ;
}

// T gate if dagger = false, T^dagger gate if dagger = true
complex_t** TGate(complex_t **T, bool dagger)
{
  int i, j;

  T[0][0].re = 1.;
  T[0][1].re = T[1][0].re = 0.;
  T[0][0].im = T[0][1].im = T[1][0].im = 0.;
  T[1][1].re = cos(PI/4.);
  if (dagger == false)
  {
    T[1][1].im = sin(PI/4.);
  }
  else
  {
    T[1][1].im = -sin(PI/4.);
  }
  for (i=0; i<=1; i++)
  {
    for (j=0; j<=1; j++)
    { 
      T[i][j] = ReIm2AmpPh(T[i][j]);
    }
  }
  return T;
}

// S gate if dagger = false, S^dagger gate if dagger = true
complex_t** SGate(complex_t **S, bool dagger)
{
  int i, j;
  S[0][0].re = 1.;
  if (dagger == false)
  {
    S[1][1].im = 1.;
  }
  else
  {
    S[1][1].im = -1.;
  }
  S[0][1].re = S[1][0].re = S[1][1].re = 0.;
  S[0][0].im = S[0][1].im = S[1][0].im = 0.;
  for (i=0; i<=1; i++)
  {
    for (j=0; j<=1; j++)
    { 
      S[i][j] = ReIm2AmpPh(S[i][j]);
    }
  }
  return S;
}

// Projector 0
complex_t** Proj0(complex_t **P0)
{
  int i, j;

  for (i=0; i<=1; i++)
  {
    for (j=0; j<=1; j++)
    {
      P0[i][j].re = 0.;
      P0[i][j].im = 0.;
      P0[0][0].re = 1.;
      P0[i][j] = ReIm2AmpPh(P0[i][j]);
    }
  }
  return P0;
}

// Projector 1
complex_t** Proj1(complex_t **P1)
{
  int i, j;

  for (i=0; i<=1; i++)
  {
    for (j=0; j<=1; j++)
    {
      P1[i][j].re = 0.;
      P1[i][j].im = 0.;
      P1[1][1].re = 1.;
      P1[i][j] = ReIm2AmpPh(P1[i][j]);
    }
  }
  return P1;
}
