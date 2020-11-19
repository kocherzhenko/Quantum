#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "Initialize.h"

#define PI 3.14159265358
#define CUTOFF 1.e-4

/* ------------------------------------------------
         Qubit & state setup
-------------------------------------------------*/

// Count the number of qubits in the input
// Arguments:
//   *filename - name of file with initial qubit states (2 coefficients per qubit)
int countQubits(char *filename)
{
  FILE *fp;
  int nqubits = 0;
  double state;

  fp = fopen(filename, "r");
  if (fp)
  {
    while (fscanf(fp, "%lf", &state) != EOF)
    {
      nqubits++;
    }
    if (nqubits % 3 == 0) // There must be 3 coeffs for the state of each qubit
    {
       nqubits /= 3; 
       printf("%d qubits.\n", nqubits);
    }
    else
    {
       printf("Input error.\n");
       return 1;
    }
  }
  else
  {
    nqubits = 0;
    printf("Error opening file.\n");
  }
  fclose(fp);
  
  return nqubits;
}

// Read in initial qubit states from file
// Arguments:
//   nqubits - number of qubits
//   **pQubits - matrix with qubit state coefficients
//               (row - # of qubit,
//                1st column - amplitude of "0" state, 2nd column - phase of "0" state,
//                3rd column - amplitude of "1" state, 4th column - phase of "1" state)
//   *filename - name of file with initial qubit states (2 coefficients per qubit)
double** getInitialQubitStates(int nqubits, double **pQubits, char *filename)
{
  FILE *fp;
  int i, j;
  float sum;
  double temp;

  fp = fopen(filename, "r");
  for (i = 0; i < 3*nqubits; i++)
  {
    fscanf(fp, "%lf", &temp); 
    if (i % 3 == 0)
    {
      pQubits[i/3][i%3] = temp;
      pQubits[i/3][i%3 + 1] = 0.;
    }
    else
    {
      pQubits[i/3][i%3 + 1] = temp;
    }
  }
  fclose (fp);

  for (i = 0; i < nqubits; i++)
  {
    // Check that the single-qubit states are normalized
    sum = pow(pQubits[i][0], 2) + pow(pQubits[i][2], 2);
    if ( abs(sum - 1.0) > CUTOFF )
    {
       printf("Error: qubit %d state not normalized.\n", i);
       free_dmatrix(pQubits, 0,nqubits-1, 0,1);
       return NULL;
    }
    // Check that the phase is in the (0, 360) degrees range, if not, adjust
    while ( pQubits[i][1] < 0. )
    {
      pQubits[i][1] += 360.;
    }
    while ( pQubits[i][3] < 0. )
    {
      pQubits[i][3] += 360.;
    }
    while ( pQubits[i][1] >= 360. )
    {
      pQubits[i][1] -= 360.;
    }
    while ( pQubits[i][3] >= 360. )
    {
      pQubits[i][3] -= 360.;
    }
    // Convert phases to radians
    pQubits[i][1] *= (PI/180.);
    pQubits[i][3] *= (PI/180.);
    if (pQubits[i][0] < CUTOFF)
    {
      pQubits[i][1] = 0.;
    }
    if (pQubits[i][2] < CUTOFF)
    {
      pQubits[i][3] = 0.;
    } 
  }

  return pQubits;
}

// Set up states of the entire N-qubit system
complex_t *getStates(complex_t *pStates, double **pQubits, int nqubits)
{
  int nstates;
  int i, j;
  int mask;
  
  nstates = (int) pow(2,nqubits);
  for (i=0; i<nstates; i++)
  {
    pStates[i].amp = 1.;
    pStates[i].ph = 0.;
    for (j=0; j<nqubits; j++)
    {
      mask = (int) pow(2,j); // set up the j-th qubit in state "1", and the rest of qubits in state "0" 
      if ( (i & mask) == 0 ) // if the j-th qubit in the i-th total system state is "0", use the amplitude and phase for its "0" state
      {
        pStates[i].amp *= pQubits[j][0];
        pStates[i].ph += pQubits[j][1];
      }
      else // if the j-th qubit in the i-th total system state is "1", use the amplitude and phase for its "1" state
      {
        pStates[i].amp *= pQubits[j][2];
        pStates[i].ph += pQubits[j][3];
      }
    }
    while (pStates[i].ph > 2.*PI)
    {
      pStates[i].ph -= 2.*PI;
    }
    if (pStates[i].amp < CUTOFF)
    {
      pStates[i].ph = 0.;
    }
    pStates[i].re = pStates[i].amp * cos(pStates[i].ph);
    pStates[i].im = pStates[i].amp * sin(pStates[i].ph);
  }

  return pStates;
}
