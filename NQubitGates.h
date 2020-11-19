typedef struct { char name; int qubit; int cqubit; int cqubit2; double theta; } gate;

/*-----------------------------------------------
       Construct Multi-Qubit Gates
-----------------------------------------------*/

complex_t** getGate(gate g, int i, complex_t **temp);

// Arguments:
//    cmatr - resultant matrix for direct product of all 1 and 2 cubit operations;
//    qubit1 - number of operator 
complex_t** getNqubitGate(complex_t **cmatr, gate g, int qubitNum);

//
// 3-qubit gates
//

// Toffoli gate
complex_t** Toffoli(complex_t **cmatr, gate g, int nqubits);

// Fredkin (CSWAP) gate
complex_t** Fredkin(complex_t **cmatr, gate g, int nqubits);


// Measure qubit qnum in the z-basis
complex_t* measureZ(complex_t *pStates, int qnum, int nqubits);
