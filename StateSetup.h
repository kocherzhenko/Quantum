/* ------------------------------------------------
         Qubit & state setup
-------------------------------------------------*/

#define CUTOFF 1.e-4

// Count the number of qubits in the input
// Arguments:
//   *filename - name of file with initial qubit states (2 coefficients per qubit)
int countQubits(char *filename);

// Read in initial qubit states from file
// Arguments:
//   nqubits - number of qubits
//   **pQubits - matrix with qubit state coefficients
//               (row - # of qubit,
//                1st column - amplitude of "0" state, 2nd column - phase of "0" state,
//                3rd column - amplitude of "1" state, 4th column - phase of "1" state)
//   *filename - name of file with initial qubit states (2 coefficients per qubit)
double** getInitialQubitStates(int nqubits, double **pQubits, char *filename);

// Set up states of the entire N-qubit system
complex_t *getStates(complex_t *pStates, double **pQubits, int nqubits);

