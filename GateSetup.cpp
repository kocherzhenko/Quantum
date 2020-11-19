#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "Initialize.h"
#include "Complex.h"
#include "NQubitGates.h"
#include "StateSetup.h"

complex_t *applyGate(complex_t **cmatr, complex_t *pStates, int nstates)
{
  int i, j;
  complex_t *temp;

  temp = cvector(0,nstates-1);
  for (i=0; i<nstates; i++)
  {
    temp[i].re = temp[i].im = 0.;
  }

  for (i=0; i<nstates; i++)
  {
    for (j=0; j<nstates; j++)
    {
      temp[i] = add(temp[i], prod(cmatr[i][j], pStates[j]));
    }    
  }

  for (i=0; i<nstates; i++)
  {
    pStates[i] = temp[i];
  }
  free_cvector(temp, 0,nstates-1);
 
  return pStates;
}

// Count the number of gates
int countGates(char *filename)
{
  FILE *fp1;
  char c;
  int gatenum = 0;

  fp1 = fopen(filename, "r");
  if (fp1)
  {
    c = getc(fp1);
    while (c != EOF)
    {
// Read in one line in the file (1 gate)
      if (c != '#') gatenum++;
      while (c != '\n')
      {
        c = getc(fp1);
      }
      c = getc(fp1);
    }
  }
  else
  {
    printf("Error opening file %s.\n", filename);
  }
  fclose(fp1);

  return gatenum;
}


// Read in the sequence of gates
gate* getAllGates(gate* allgates, char *filename, int nqubits)
{
  FILE *fp1;
  char onegate[128], element[32];
  char c;
  int i, j, elnum;
  gate g;
  int qubit, cqubit, cqubit2;
  double theta = -100., dnum;
  int gatenum;

  fp1 = fopen(filename, "r");
  if (fp1)
  {
    c = getc(fp1);
    gatenum = 0;
    while (c != EOF)
    {
// Detect comments
      if (c == '#')
      {
        while (c != '\n') c = getc(fp1);
      }
// Read in one line in the file (1 gate)
      i = 0;
      while (c != '\n')
      {
        onegate[i] = c;
        i++;
        c = getc(fp1);
      }
      onegate[i] = '\0';
//      printf("%s\n", onegate);
      c = getc(fp1);
// Ignore comments
      if (onegate[0] == '\0') continue;

// Break down the parameters for this gate
      i = 0;
      elnum = 0;
      g.name = 'U';
      g.qubit = -1;
      g.cqubit = -1;
      g.cqubit2 = -1;
      g.theta = -100;
      while (onegate[i] != '\0')
      {
        while (onegate[i] == ' ') // Treat multiple spaces as one
        {
          i++;
        }
        if (onegate[i] == 0) // Ignore spaces at the end of a line
        {
          continue;
        }
        j = 0;
        while ((onegate[i] != ' ') && (onegate[i] != '\0')) // Split each line into elements of gate
        {
          element[j] = onegate[i];
          i++;
          j++;
        }
        element[j] = '\0';
// Detect gate type
        if (elnum == 0)
        {
          if ( strncmp(element, "I", 5) == 0 )
          {
            g.name = 'I';
          }
          else if ( strncmp(element, "M", 5) == 0 )
          {
            g.name = 'M';
          }
          else if ( strncmp(element, "HAD", 5) == 0 )
          {
            g.name = 'H';
          }
          else if ( strncmp(element, "PX", 5) == 0 )
          {
            g.name = 'X';
          }
          else if ( strncmp(element, "PY", 5) == 0 )
          {
            g.name = 'Y';
          }
          else if ( strncmp(element, "PZ", 5) == 0 )
          {
            g.name = 'Z';
          }
          else if ( strncmp(element, "RX", 5) == 0 )
          {
            g.name = 'x';
          }
          else if ( strncmp(element, "RY", 5) == 0 )
          {
            g.name = 'y';
          }
          else if ( strncmp(element, "RZ", 5) == 0 )
          {
            g.name = 'z';
          }
          else if ( strncmp(element, "T", 5) == 0 )
          {
            g.name = 't';
          }
          else if ( strncmp(element, "S", 5) == 0 )
          {
            g.name = 's';
          }
          else if ( strncmp(element, "P0", 5) == 0 )
          {
            g.name = '0';
          }
          else if ( strncmp(element, "P1", 5) == 0 )
          {
            g.name = '1';
          }
          else if ( strncmp(element, "CX", 5) == 0 )
          {
            g.name = 'C';
          }
          else if ( strncmp(element, "CZ", 5) == 0 )
          {
            g.name = 'c';
          }
          else if ( strncmp(element, "SW", 5) == 0 )
          {
            g.name = 'S';
          }
          else if ( strncmp(element, "TF", 5) == 0 )
          {
            g.name = 'T';
          }
          else if ( strncmp(element, "FR", 5) == 0 )
          {
            g.name = 'F';
          }
          else
          {
            printf("Invalid gate %d: %s.\n", gatenum, element);
            return NULL;
          } 
        }
// Read in the qubit to which the gate is applied
        else if (elnum == 1)
        {
          qubit = atoi(element);
          dnum = strtod(element, NULL);
          if (qubit == 0 || abs(qubit - dnum) > CUTOFF || qubit < 1 || qubit > nqubits) // Check that qubit is an integer between 1 and nqubits
          {
            printf("Invalid qubit number for gate %d: %s.\n", gatenum, element);
            return NULL;
          }
          g.qubit = qubit - 1;
        }
// Read in the control qubit for CNOT, CPHASE, and Toffoli, the second qubit for SWAP, theta for RX, RY, and RZ, or dagger/no dagger for T and S
        else if (elnum == 2)
        {
          if (g.name == 'I' || g.name == 'M' || g.name == 'N' || g.name == 'H' || g.name == 'X' || g.name == 'Y' || g.name == 'Z' || g.name == '0' || g.name == '1')
          {
            printf("Unexpected parameter for gate %d: %s.\n", gatenum, element);
            return NULL;
          }
          else if (g.name == 't' || g.name == 's')
          {
            if (strncmp(element, "false", 8) == 0)
            {
              continue;
            }
            else if (strncmp(element, "true", 8) == 0)
            {
              if (g.name == 't')
              {
                g.name = 'd';
              }
              else
              {
                g.name = 'a';
              }
            }
            else
            {
              printf("Not a boolean value for gate %d: %s.\n", gatenum, element);
              return NULL;
            }
          }
          else if (g.name == 'x' || g.name == 'y' || g.name == 'z')
          {
            theta = strtod(element, NULL);
            if ( (theta > 360.) || (theta < 0.) )
            {
              printf("Invalid angle for gate %d: %s.\n", gatenum, element);
              return NULL;
            }
            g.theta = theta*PI/180.;
            theta = -100.;
          }
          else if (g.name == 'C' || g.name == 'c' || g.name == 'S' || g.name == 'T' || g.name == 'F')
          {
            qubit = atoi(element);
            dnum = strtod(element, NULL);
            // Check that cqubit is an integer between 1 and nqubits and is not the same as qubit
            if (qubit == 0 || abs(qubit - dnum) > CUTOFF || qubit < 1 || qubit > nqubits || qubit-1 == g.qubit)
            {
              printf("Invalid qubit number for gate %d: %s.\n", gatenum, element);
              return NULL;
            }
            g.cqubit = qubit - 1;
          }
        }
// Read in the control control qubit for Toffoli and Fredkin gates
        else if (elnum == 3)
        {
          if (g.name == 'T' || g.name == 'F')
          {
            qubit = atoi(element);
            dnum = strtod(element, NULL);
            // Check that cqubit2 is an integer between 1 and nqubits and is not the same as qubit or cqubit
            if (qubit == 0 || abs(qubit - dnum) > CUTOFF || qubit < 1 || qubit > nqubits || qubit-1 == g.qubit || qubit-1 == g.cqubit)
            {
              printf("Invalid qubit number for gate %d: %s.\n", gatenum, element);
              return NULL;
            }
            g.cqubit2 = qubit - 1;
          }
          else
          {
            printf("Unexpected parameter for gate %d: %s.\n", gatenum, element);
            return NULL;
          }
        }
        else
        {
          printf("Unexpected parameter for gate %d: %s.\n", gatenum, element);
          return NULL;
        }
        if (onegate[i] == '\0')
        {
          break;
        }
        elnum++;
      }

// Check that all required parameters for the gate have been set
      if (g.qubit == -1)
      {
        printf("Error: not enough parameters for gate %d.\n", gatenum);
        return NULL;
      }
      if ( (g.name == 'C' || g.name == 'c' || g.name == 'S' || g.name == 'T') && (g.cqubit == -1) )
      {
        printf("Error: not enough parameters for gate %d.\n", gatenum);
        return NULL;
      }
      if ( (g.name == 'x' || g.name == 'y' || g.name == 'z') && (g.theta < 0) )
      {
        printf("Error: not enough parameters for gate %d.\n", gatenum);
        return NULL;
      }
      if ( (g.name == 'T' || g.name == 'F') && (g.cqubit2 == -1) )
      {
        printf("Error: not enough parameters for gate %d.\n", gatenum);
        return NULL;
      }
//      printf(" %5c %5d %5d %5d %6.3f \n", g.name, g.qubit, g.cqubit, g.cqubit2, g.theta);
      allgates[gatenum] = g;
      gatenum++;
    }
  }
  else
  {
    printf("Couldn't open file %s.\n", filename);
  }

  fclose(fp1);
  

  return allgates;
}
