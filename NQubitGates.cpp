#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <stdbool.h>
#include "Initialize.h"
#include "Complex.h"
#include "Matrix.h"
#include "OneQubitGates.h"
#include "TwoQubitGates.h"

typedef struct { char name; int qubit; int cqubit; int cqubit2; double theta; } gate;

/*-----------------------------------------------
       Construct Multi-Qubit Gates
-----------------------------------------------*/

complex_t** getGate(gate g, int i, complex_t **temp)
{
  temp = cmatrix(0,1, 0,1);
  if (g.qubit != i)
  {
    temp = Identity(temp);
  } 
  else
  {
    switch (g.name)
    {
      case 'I':
        temp = Identity(temp);
        break;
      case 'H':
        temp = Hadamard(temp);
        break;
      case 'X':
        temp = PauliX(temp);
        break;
      case 'Y':
        temp = PauliY(temp);
        break;
      case 'Z':
        temp = PauliZ(temp);
        break;
      case 'x':
        temp = RX(temp, g.theta);
        break;
      case 'y':
        temp = RY(temp, g.theta);
        break;
      case 'z':
        temp = RZ(temp, g.theta);
        break;
      case 't':
        temp = TGate(temp, false);
        break;
      case 'd':
        temp = TGate(temp, true);
        break;
      case 's':
        temp = SGate(temp, false);
        break;
      case 'a':
        temp = SGate(temp, true);
        break;
      default:
        printf("Invalid gate\n");
    }  
  }

  return temp;
}

// Arguments:
//    cmatr - resultant matrix for direct product of all 1 and 2 cubit operations;
//    qubit1 - number of operator 
complex_t** getNqubitGate(complex_t **cmatr, gate g, int qubitNum)
{
  int i, j;
  int k, l;
  complex_t **temp, **temp0, **temp1;
  
  temp = cmatrix(0,1, 0,1);

  if (g.name == 'I' || g.name == 'H' || g.name == 'X' || g.name == 'Y' || g.name == 'Z' || g.name == 'x' || g.name == 'y' || g.name == 'z' || g.name == 't' || g.name == 'd' || g.name == 's' || g.name == 'a')
  {
//    printf(" %c %6.3f %d %d \n", g.name, g.theta, g.qubit, g.cqubit);
    temp0 = cmatrix(0,1, 0,1);
    temp0 = getGate(g, 0, temp0);
/*
    for (k=0; k<2; k++)
    {
      for (l=0; l<2; l++)
      {
        printf("%6.3f+i%6.3f ", temp0[k][l].re, temp0[k][l].im);
      }
    } 
*/    
    for (i=1; i<qubitNum; i++)
    {
      temp = getGate(g, i, temp);
/*
      for (k=0; k<2; k++)
      {
        for (l=0; l<2; l++)
        {
          printf("%6.3f+i%6.3f ", temp[k][l].re, temp[k][l].im);
        }
      } 
*/
      temp1 = cmatrix(0,pow(2,i+1)-1, 0,pow(2,i+1)-1);
      temp1 = direct(temp1, temp0, temp, pow(2,i), 2);
      free_cmatrix(temp0, 0,pow(2,i)-1, 0,pow(2,i)-1);
      temp0 = temp1;
    }
    for (k=0; k<pow(2,qubitNum); k++)
    {
      for (l=0; l<pow(2,qubitNum); l++)
      {
        cmatr[k][l] = temp0[k][l];
      }
    }

    free_cmatrix(temp1, 0,pow(2,qubitNum)-1, 0,pow(2,qubitNum)-1);
  }
  else if (g.name == 'C' || g.name == 'c' || g.name == 'S')
  {
//    printf(" %c %6.3f %d %d \n", g.name, g.theta, g.qubit, g.cqubit);
    if (g.name == 'C')
    {
      cmatr = CNOT(cmatr, g.qubit, g.cqubit, qubitNum);
    }
    else if (g.name == 'c')
    {
      cmatr = CPHASE(cmatr, g.qubit, g.cqubit, qubitNum);
    }
    else if (g.name == 'S')
    {
      cmatr = SWAP(cmatr, g.qubit, g.cqubit, qubitNum);
    }
  }
  else
  {
    printf("Error: undefined gate.\n");
  }

  free_cmatrix(temp, 0,1, 0,1);
  return cmatr;
}

//
// 3-qubit gates
//

// Toffoli gate
complex_t** Toffoli(complex_t **cmatr, gate g, int nqubits)
{
  gate g1;
  complex_t **temp1, **temp2, **temp, **hold;
  int i, j, size;

  size = pow(2,nqubits);
  temp1 = cmatrix(0,size-1, 0,size-1);
  temp2 = cmatrix(0,size-1, 0,size-1);
  temp = cmatrix(0,size-1, 0,size-1);

// apply Hadamard to g.qubit
  g1.name = 'H';
  g1.theta = 0.;
  g1.qubit = g.qubit;
  g1.cqubit = nqubits;
  g1.cqubit2 = nqubits+1;
  temp1 = getNqubitGate(temp1, g1, nqubits);

// apply CNOT to g.qubit with g.cqubit as control
  g1.name = 'C';
  g1.theta = 0.;
  g1.qubit = g.qubit;
  g1.cqubit = g.cqubit;
  g1.cqubit2 = nqubits+1;
  temp2 = getNqubitGate(temp2, g1, nqubits);

// take product of first 2 gates, point temp1 to result, point temp to former temp1 to store next product
  temp = matrprod(temp, temp1, temp2, size);
  hold = temp1;
  temp1 = temp;
  temp = hold;

// apply T^dagger to g.qubit
  g1.name = 'd';
  g1.theta = 0.;
  g1.qubit = g.qubit;
  g1.cqubit = nqubits;
  g1.cqubit2 = nqubits+1;
  temp2 = getNqubitGate(temp2, g1, nqubits);

  temp = matrprod(temp, temp1, temp2, size);
  hold = temp1;
  temp1 = temp;
  temp = hold;

// apply CNOT to g.qubit with g.cqubit2 as control
  g1.name = 'C';
  g1.theta = 0.;
  g1.qubit = g.qubit;
  g1.cqubit = g.cqubit2;
  g1.cqubit2 = nqubits;
  temp2 = getNqubitGate(temp2, g1, nqubits);
   
  temp = matrprod(temp, temp1, temp2, size);
  hold = temp1;
  temp1 = temp;
  temp = hold;

// apply T to g.qubit
  g1.name = 't';
  g1.theta = 0.;
  g1.qubit = g.qubit;
  g1.cqubit = nqubits;
  g1.cqubit2 = nqubits+1;
  temp2 = getNqubitGate(temp2, g1, nqubits);

  temp = matrprod(temp, temp1, temp2, size);
  hold = temp1;
  temp1 = temp;
  temp = hold;

// apply CNOT to g.qubit with g.cqubit as control
  g1.name = 'C';
  g1.theta = 0.;
  g1.qubit = g.qubit;
  g1.cqubit = g.cqubit;
  g1.cqubit2 = nqubits+1;
  temp2 = getNqubitGate(temp2, g1, nqubits);
  
  temp = matrprod(temp, temp1, temp2, size);
  hold = temp1;
  temp1 = temp;
  temp = hold;

// apply T^dagger to g.qubit
  g1.name = 'd';
  g1.theta = 0.;
  g1.qubit = g.qubit;
  g1.cqubit = nqubits;
  g1.cqubit2 = nqubits+1;
  temp2 = getNqubitGate(temp2, g1, nqubits);

  temp = matrprod(temp, temp1, temp2, size);
  hold = temp1;
  temp1 = temp;
  temp = hold;

// apply CNOT to g.qubit with g.cqubit2 as control
  g1.name = 'C';
  g1.theta = 0.;
  g1.qubit = g.qubit;
  g1.cqubit = g.cqubit2;
  g1.cqubit2 = nqubits;
  temp2 = getNqubitGate(temp2, g1, nqubits);
   
  temp = matrprod(temp, temp1, temp2, size);
  hold = temp1;
  temp1 = temp;
  temp = hold;

// apply T to g.cqubit
  g1.name = 't';
  g1.theta = 0.;
  g1.qubit = g.cqubit;
  g1.cqubit = nqubits;
  g1.cqubit2 = nqubits+1;
  temp2 = getNqubitGate(temp2, g1, nqubits);

  temp = matrprod(temp, temp1, temp2, size);
  hold = temp1;
  temp1 = temp;
  temp = hold;

// apply T to g.cubit
  g1.name = 't';
  g1.theta = 0.;
  g1.qubit = g.qubit;
  g1.cqubit = nqubits;
  g1.cqubit2 = nqubits+1;
  temp2 = getNqubitGate(temp2, g1, nqubits);

  temp = matrprod(temp, temp1, temp2, size);
  hold = temp1;
  temp1 = temp;
  temp = hold;

// apply CNOT to g.cqubit with g.cqubit2 as control
  g1.name = 'C';
  g1.theta = 0.;
  g1.qubit = g.cqubit;
  g1.cqubit = g.cqubit2;
  g1.cqubit2 = nqubits;
  temp2 = getNqubitGate(temp2, g1, nqubits);
   
  temp = matrprod(temp, temp1, temp2, size);
  hold = temp1;
  temp1 = temp;
  temp = hold;

// apply Hadamard to g.cubit
  g1.name = 'H';
  g1.theta = 0.;
  g1.qubit = g.qubit;
  g1.cqubit = nqubits;
  g1.cqubit2 = nqubits+1;
  temp2 = getNqubitGate(temp2, g1, nqubits);
   
  temp = matrprod(temp, temp1, temp2, size);
  hold = temp1;
  temp1 = temp;
  temp = hold;

// apply T to g.cqubit2
  g1.name = 't';
  g1.theta = 0.;
  g1.qubit = g.cqubit2;
  g1.cqubit = nqubits;
  g1.cqubit2 = nqubits+1;
  temp2 = getNqubitGate(temp2, g1, nqubits);

  temp = matrprod(temp, temp1, temp2, size);
  hold = temp1;
  temp1 = temp;
  temp = hold;

// apply T^dagger to g.cqubit
  g1.name = 'd';
  g1.theta = 0.;
  g1.qubit = g.cqubit;
  g1.cqubit = nqubits;
  g1.cqubit2 = nqubits+1;
  temp2 = getNqubitGate(temp2, g1, nqubits);

  temp = matrprod(temp, temp1, temp2, size);
  hold = temp1;
  temp1 = temp;
  temp = hold;

// apply CNOT to g.cqubit with g.cqubit2 as control
  g1.name = 'C';
  g1.theta = 0.;
  g1.qubit = g.cqubit;
  g1.cqubit = g.cqubit2;
  g1.cqubit2 = nqubits;
  temp2 = getNqubitGate(temp2, g1, nqubits);
   
  temp = matrprod(temp, temp1, temp2, size);

  for (i=0; i<size; i++)
  {
    for (j=0; j<size; j++)
    {
      cmatr[i][j] = temp[i][j];
    }
  }

  free_cmatrix(temp1, 0,size-1, 0,size-1);
  free_cmatrix(temp2, 0,size-1, 0,size-1);
  free_cmatrix(temp, 0,size-1, 0,size-1);

  return cmatr;
}

// Fredkin (CSWAP) gate
complex_t** Fredkin(complex_t **cmatr, gate g, int nqubits)
{
  gate g1;
  complex_t **temp1, **temp2, **temp, **hold;
  int i, j, size;

  size = pow(2,nqubits);
  temp1 = cmatrix(0,size-1, 0,size-1);
  temp2 = cmatrix(0,size-1, 0,size-1);
  temp = cmatrix(0,size-1, 0,size-1);
  
  g1.name = 'T';
  g1.theta = 0.;
  g1.qubit = g.qubit;
  g1.cqubit = g.cqubit;
  g1.cqubit2 = g.cqubit2;
  temp1 = Toffoli(temp1, g1, nqubits);

  g1.name = 'T';
  g1.theta = 0.;
  g1.qubit = g.cqubit;
  g1.cqubit = g.qubit;
  g1.cqubit2 = g.cqubit2;
  temp2 = Toffoli(temp2, g1, nqubits);
  
  temp = matrprod(temp, temp1, temp2, size);
  hold = temp1;
  temp1 = temp;
  temp = hold;
  
  g1.name = 'T';
  g1.theta = 0.;
  g1.qubit = g.qubit;
  g1.cqubit = g.cqubit;
  g1.cqubit2 = g.cqubit2;
  temp2 = Toffoli(temp2, g1, nqubits);
  
  temp = matrprod(temp, temp1, temp2, size);

  for (i=0; i<size; i++)
  {
    for (j=0; j<size; j++)
    {
      cmatr[i][j] = temp[i][j];
    }
  }
  
  free_cmatrix(temp1, 0,size-1, 0,size-1);
  free_cmatrix(temp2, 0,size-1, 0,size-1);
  free_cmatrix(temp, 0,size-1, 0,size-1);

  return cmatr;
}

// Measure qubit qnum in the z-basis
complex_t* measureZ(complex_t *pStates, int qnum, int nqubits)
{
  int i, mask, nstates, flag;
  double prob0 = 0., prob;

  nstates = (int) pow(2,nqubits);

  if (qnum >= 0 && qnum < nqubits)
  {
    mask = (int) pow(2,nqubits-qnum-1);
    for (i=0; i<nstates; i++)
    {
      if ( (i & mask) == 0 )
      {
        prob0 += pow(pStates[i].amp, 2);
      }
    } 

    prob = ( (double) rand() ) / ( (double) RAND_MAX ); 
    printf("Prob0 = %10.8f, random = %10.8f\n", prob0, prob);
    if (prob < prob0) flag = 0;
    else flag = mask;

// Select N-qubit states with measured value for the qnum-th qubit and renormalize probabilities
    for (i=0; i<nstates; i++)
    {
      if ( (i & mask) == flag )
      {
        if (flag == 0) pStates[i].amp /= sqrt(prob0);
        else pStates[i].amp /= sqrt(1.-prob0);
        pStates[i].re = pStates[i].amp*cos(pStates[i].ph);
        pStates[i].im = pStates[i].amp*sin(pStates[i].ph);
      }
      else
      {
        pStates[i].amp = pStates[i].ph = pStates[i].re = pStates[i].im = 0.;
      }
    }

  }
  return pStates;
}
