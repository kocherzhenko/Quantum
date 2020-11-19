// Set up and apply gates

// Apply an N-qubit gate to an N-qubit state vector
complex_t *applyGate(complex_t **cmatr, complex_t *pStates, int nstates);

// Count the number of gates in the input
int countGates(char *filename);

// Read in the sequence of gates for the input file
gate* getAllGates(gate* allgates, char *filename, int nqubits);
