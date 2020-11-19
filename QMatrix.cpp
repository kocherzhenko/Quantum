#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <stdbool.h>
#include <string.h>
#include "Initialize.h"
#include "Complex.h"
#include "Matrix.h"
#include "OneQubitGates.h"
#include "TwoQubitGates.h"
#include "NQubitGates.h"
#include "StateSetup.h"
#include "GateSetup.h"


// Print binary representation of a number
int printBin(int i, int nqubits)
{
  int j, mask;
    
  for (j=0; j<nqubits; j++)
  {
    mask = (int) pow(2,j); // set up the j-th qubit in state "1", and the rest of qubits in state "0" 
    if ( (i & mask) == 0 ) // if the j-th qubit in the i-th total system state is "0", print "0"
    {
      printf("0");
    }
    else // if the j-th qubit in the i-th total system state is "1", print "1"
    {
      printf("1");
    }
  }

  return 0;
}


int main(int argc, char **argv)
{
  FILE *fp;
  char *filename; //= "InitialQubitStates.in";
  char *filename1; //= "GateSequence.in";

  int i, j, k;
  int nqubits;
  int nstates;
  int ngates;
  int niter;
  double sum;
  double **pQubits;
  complex_t *pStates;

  complex_t **I, **Had;
  complex_t **PauX, **PauY, **PauZ;
  complex_t **RotX, **RotY, **RotZ;
  double theta;
  complex_t **CNOT, **CPhase; 
  complex_t **cmatr;
  gate g, *allgates;
  
  time_t t; 
  srand((unsigned) time(&t)); 

  if (argc != 4)
  {
    printf("Error: The program requires two input files: initial qubit states and gate sequence, as well as the number of iterations to be performed..\n");
    return 1;
  }
  filename = argv[1];
  filename1 = argv[2];
  niter = atoi(argv[3]);
//  printf(" %s %s %s %s \n", argv[0], filename, );

// Count qubits
  nqubits = countQubits(filename);
  nstates = (int) pow(2., (double) nqubits);
  printf("Number of qubits is %d, number of states is %d.\n", nqubits, nstates);

// Get initial single-qubit states
  pQubits = dmatrix(0,nqubits-1, 0,3);
  pQubits = getInitialQubitStates(nqubits, pQubits, filename);
  if (pQubits == NULL)
  {
    return 1;
  }
  /*
  else
  {
    for (i = 0; i < nqubits; i++)
    {
      printf(" %lf  %lf %lf %lf \n", pQubits[i][0], pQubits[i][1], pQubits[i][2], pQubits[i][3]);
    }
  }
  */

  for (k=0; k<niter; k++)
  {
    printf("\nIteration %d.\n\n", k+1);
    // Construct n-qubit states
    pStates = cvector(0,nstates-1);
    pStates = getStates(pStates, pQubits, nqubits);
    printf("Initial state:\n");
    for (i=0; i<nstates; i++)
    {
      sum += pow(pStates[i].amp, 2); 
      if (pow(pStates[i].amp,2) > CUTOFF)
      {
        printBin(i, nqubits);
        printf(" %5d  %6.3f  %6.3f %6.3f \n", i, pStates[i].re, pStates[i].im, pow(pStates[i].amp, 2));
      }
    } 
    printf("\n");
 
    cmatr = cmatrix(0,nstates-1, 0,nstates-1);
 
    // Allocate a vector of gate structures with subscript range [0..ngates]
    ngates = countGates(filename1);
    
    allgates = (gate*) malloc( ngates * sizeof(gate) );
    if (!allgates)
    {
      printf("Allocation failure for gate array.\n");
    }
    
    allgates = getAllGates(allgates, filename1, nqubits);
    /*
    for (i=0; i<ngates; i++)
    {
      printf(" %5c %5d %5d %5d %6.3f \n", allgates[i].name, allgates[i].qubit, allgates[i].cqubit, allgates[i].cqubit2, allgates[i].theta);
    }
    */
    
    // Apply the sequence of gates to the N-qubit state  
    for (i=0; i<ngates; i++)
    {
      if (allgates[i].name != 'T' && allgates[i].name != 'F' && allgates[i].name != 'M')
      {
        cmatr = getNqubitGate(cmatr, allgates[i], nqubits);
        pStates = applyGate(cmatr, pStates, nstates); 
      }
      else if (allgates[i].name == 'T')
      {
        cmatr = Toffoli(cmatr, allgates[i], nqubits);
        pStates = applyGate(cmatr, pStates, nstates); 
      }
      else if (allgates[i].name == 'F')
      {
        cmatr = Fredkin(cmatr, allgates[i], nqubits);
        pStates = applyGate(cmatr, pStates, nstates); 
      }
      else // measure qubit
      {
        pStates = measureZ(pStates, allgates[i].qubit, nqubits);
      }
 
      printf("Gate %d: %c \n", i+1, allgates[i].name);
      
    /*
      for (j=0; j<nstates; j++)
      {
        for (k=0; k<nstates; k++)
        {
          printf(" %6.3f ", cmatr[j][k].re);
        }
        printf("\n");
      }
    */
 
      sum = 0.;
      printf("State after gate %d: \n", i+1);
      for (j=0; j<nstates; j++)
      {
        sum += pow(pStates[j].amp, 2); 
        if (pow(pStates[j].amp,2) > CUTOFF)
        {
          printBin(j, nqubits);
          printf(" %5d  %6.3f  %6.3f %6.3f \n", j, pStates[j].re, pStates[j].im, pow(pStates[j].amp, 2));
        }
      } 
      printf("\n");
      if (abs(1.-sum) > CUTOFF)
      {
        printf("Error, normalization = %6.3f != 1.\n", sum);
        return 1;
      }
    } 
    /*  
    for (i=0; i<nstates; i++)
    {
      printf(" %5d    %6.3f  %6.3f  %6.3f  %6.3f \n", i, pStates[i].amp, pStates[i].ph, pStates[i].re, pStates[i].im);
    } 
    */ 
  }

  free(allgates);
  free_cmatrix(cmatr, 0,nstates-1, 0,nstates-1);
  free_cvector(pStates, 0,nstates-1); 

  return 0;
}
