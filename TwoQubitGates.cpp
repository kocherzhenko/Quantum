#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "Initialize.h"
#include "Complex.h"
#include "Matrix.h"
#include "OneQubitGates.h"

//
// Controlled gates 
//

// CNOT
complex_t** CNOT(complex_t **cNOT, int qubit, int cqubit, int nqubits)
{
  int i, j;
  complex_t **temp0, **temp1, **temp2, **temp;

  if (qubit != cqubit)
  {
    temp0 = cmatrix(0,1, 0,1);
    temp = cmatrix(0,1, 0,1);
    if (cqubit != 0)
    {
      temp0 = Identity(temp0);
    }
    else
    {
      temp0 = Proj0(temp0);
    }
    for (i=1; i<nqubits; i++)
    {
      temp1 = cmatrix(0,pow(2,i+1)-1, 0,pow(2,i+1)-1);
      if (i != cqubit)
      {
        temp = Identity(temp);
      }
      else
      {
        temp = Proj0(temp);
      }
      temp1 = direct(temp1, temp0, temp, pow(2,i), 2);
      free_cmatrix(temp0, 0,pow(2,i)-1, 0,pow(2,i)-1);
      temp0 = temp1;
    }

    temp0 = cmatrix(0,1, 0,1);
    if (cqubit != 0 && qubit != 0)
    {
      temp0 = Identity(temp0);
    }
    else if (qubit != 0)
    {
      temp0 = Proj1(temp0);
    }
    else
    {
      temp0 = PauliX(temp0);
    }
    for (i=1; i<nqubits; i++)
    {
      temp2 = cmatrix(0,pow(2,i+1)-1, 0,pow(2,i+1)-1);
      if (i != cqubit && i != qubit)
      {
        temp = Identity(temp);
      }
      else if (i != qubit)
      {
        temp = Proj1(temp);
      }
      else
      {
        temp = PauliX(temp);
      }
      temp2 = direct(temp2, temp0, temp, pow(2,i), 2);
      free_cmatrix(temp0, 0,pow(2,i)-1, 0,pow(2,i)-1);
      temp0 = temp2;
    }
    free_cmatrix(temp, 0,1, 0,1);
    
    for (i=0; i<pow(2,nqubits); i++)
    {
      for (j=0; j<pow(2,nqubits); j++)
      {
        cNOT[i][j] = add(temp1[i][j],temp2[i][j]);
      }
    }

    free_cmatrix(temp1, 0,pow(2,nqubits)-1, 0,pow(2,nqubits)-1);
    free_cmatrix(temp2, 0,pow(2,nqubits)-1, 0,pow(2,nqubits)-1);
  }
  else
  {
    printf("Error: qubit and cqubit cannot be the same");
    cNOT = NULL;
  }

  return cNOT;
}

// CPHASE
complex_t** CPHASE(complex_t **cPhase, int qubit, int cqubit, int nqubits)
{
  int i, j;
  complex_t **temp0, **temp1, **temp2, **temp;

  if (qubit != cqubit)
  {
    temp1 = cmatrix(0,pow(2,nqubits)-1, 0,pow(2,nqubits)-1);
    for (i=0; i<pow(2,nqubits); i++)
    {
      for (j=0; j<pow(2,nqubits); j++)
      {
        if (i != j) { temp1[i][j].re = 0.; }
        else { temp1[i][j].re = 1.; }
        temp1[i][j].im = 0.;
        temp1[i][j] = ReIm2AmpPh(temp1[i][j]);
      }
    }

    temp0 = cmatrix(0,1, 0,1);
    temp = cmatrix(0,1, 0,1);
    if (cqubit != 0 && qubit != 0)
    {
      temp0 = Identity(temp0);
    }
    else
    {
      temp0 = Proj1(temp0);
    }
    for (i=1; i<nqubits; i++)
    {
      temp2 = cmatrix(0,pow(2,i+1), 0,pow(2,i+1));
      if (i != cqubit && i != qubit)
      {
        temp = Identity(temp);
      }
      else
      {
        temp = Proj1(temp);
      }
      temp2 = direct(temp2, temp0, temp, pow(2,i), 2);
      free_cmatrix(temp0, 0,pow(2,i), 0,pow(2,i));
      temp0 = temp2;
    }
    free_cmatrix(temp, 0,1, 0,1);

    for (i=0; i<pow(2,nqubits); i++)
    {
      for (j=0; j<pow(2,nqubits); j++)
      {
        cPhase[i][j] = add(temp1[i][j],prod(temp2[i][j],-2.));
      }
    }

    free_cmatrix(temp1, 0,pow(2,nqubits), 0,pow(2,nqubits));
    free_cmatrix(temp2, 0,pow(2,nqubits), 0,pow(2,nqubits));
  }
  else
  {
    printf("Error: qubit and cqubit cannot be the same");
    cPhase = NULL;
  }

  return cPhase;
}

// SWAP
complex_t** SWAP(complex_t **swap, int qubit1, int qubit2, int nqubits)
{
  int i, j;
  complex_t **temp1, **temp2, **temp3, **temp;
  
  temp1 = cmatrix(0,pow(2,nqubits)-1, 0,pow(2,nqubits)-1);
  temp1 = CNOT(temp1, qubit1, qubit2, nqubits);
  temp2 = cmatrix(0,pow(2,nqubits)-1, 0,pow(2,nqubits)-1);
  temp2 = CNOT(temp2, qubit2, qubit1, nqubits);
  temp3 = cmatrix(0,pow(2,nqubits)-1, 0,pow(2,nqubits)-1);
  temp3 = CNOT(temp3, qubit1, qubit2, nqubits);
  temp = cmatrix(0,pow(2,nqubits)-1, 0,pow(2,nqubits)-1);
  temp = matrprod(temp, temp1, temp2, pow(2,nqubits));
  free_cmatrix(temp1, 0,pow(2,nqubits)-1, 0,pow(2,nqubits)-1);
  free_cmatrix(temp2, 0,pow(2,nqubits)-1, 0,pow(2,nqubits)-1);
  swap = matrprod(swap, temp, temp3, pow(2,nqubits));
  free_cmatrix(temp3, 0,pow(2,nqubits)-1, 0,pow(2,nqubits)-1);
  free_cmatrix(temp, 0,pow(2,nqubits)-1, 0,pow(2,nqubits)-1);

  return swap;
}
