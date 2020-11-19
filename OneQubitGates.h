/*---------------------------------
         Set up gates
---------------------------------*/
//
// 1-qubit gates 
//

// Identity
complex_t** Identity(complex_t **I);

// Hadamard
complex_t** Hadamard(complex_t **Had);

// PauliX
complex_t** PauliX(complex_t **PauX);

// PauliY
complex_t** PauliY(complex_t **PauY);

// PauliZ
complex_t** PauliZ(complex_t **PauZ);

// RotX
complex_t** RX(complex_t **RotX, double theta);

// RotY
complex_t** RY(complex_t **RotY, double theta);

// RotZ
complex_t** RZ(complex_t **RotZ, double theta);

// T gate if dagger = false, T^dagger gate if dagger = true
complex_t** TGate(complex_t **T, bool dagger);

// S gate if dagger = false, S^dagger gate if dagger = true
complex_t** SGate(complex_t **S, bool dagger);

// Projector 0
complex_t** Proj0(complex_t **P0); 

// Projector 1
complex_t** Proj1(complex_t **P1);
