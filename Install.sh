g++ -c Initialize.cpp
g++ -c Complex.cpp
g++ -c Matrix.cpp
g++ -c OneQubitGates.cpp
g++ -c TwoQubitGates.cpp
g++ -c NQubitGates.cpp
g++ -c StateSetup.cpp
g++ -c GateSetup.cpp
g++ QMatrix.cpp Initialize.o Complex.o Matrix.o OneQubitGates.o TwoQubitGates.o NQubitGates.o StateSetup.o GateSetup.o -o qgc
rm *.o
