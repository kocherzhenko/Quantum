#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "Initialize.h"
#include "Complex.h"

/* ------------------------------------------------
         Matrix operations
-------------------------------------------------*/

// Matrix product
complex_t** matrprod(complex_t** product, complex_t** matr1, complex_t** matr2, int size)
{
  int i,j,k;

  for (i=0; i<size; i++)
  {
    for (j=0; j<size; j++)
    {
      product[i][j].re = 0.;
      product[i][j].im = 0.;
    }
  }
  
  for (i=0; i<size; i++)
  {
    for (j=0; j<size; j++)
    {
      for (k=0; k<size; k++)
      {
        product[i][j] = add( product[i][j], prod(matr1[i][k],matr2[k][j]) );
      }
    }
  }

  return product;
}

// Direct product
// Arguments:
//    product - diect product of square matrices matr1 and matr2 (of sizes size1 and size 2, respectively)
complex_t** direct(complex_t** product, complex_t** matr1, complex_t** matr2, int size1, int size2)
{
  int i1, j1, i2, j2;

  for (i1=0; i1<size1; i1++)
  {
    for (j1=0; j1<size1; j1++)
    {
       for (i2=0; i2<size2; i2++)
       {
         for (j2=0; j2<size2; j2++)
         {
           product[i1*size2+i2][j1*size2+j2] = prod(matr1[i1][j1],matr2[i2][j2]);
         }
       }
    }
  }
  
  return product;
}
