//
// Controlled gates 
//

// CNOT
complex_t** CNOT(complex_t **cNOT, int qubit, int cqubit, int nqubits);

// CPHASE
complex_t** CPHASE(complex_t **cPhase, int qubit, int cqubit, int nqubits);

// SWAP
complex_t** SWAP(complex_t **swap, int qubit1, int qubit2, int nqubits);
