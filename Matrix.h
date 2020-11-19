/* ------------------------------------------------
         Matrix operations
-------------------------------------------------*/

// Matrix product
complex_t** matrprod(complex_t** product, complex_t** matr1, complex_t** matr2, int size);

// Direct product
// Arguments:
//    product - diect product of square matrices matr1 and matr2 (of sizes size1 and size 2, respectively)
complex_t** direct(complex_t** product, complex_t** matr1, complex_t** matr2, int size1, int size2);
