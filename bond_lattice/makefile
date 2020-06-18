all: chain_lattice.cpp
	g++ -std=c++11 -Wall -O3 -funroll-loops -fPIC -O2 -c chain_lattice.cpp -o Backend.o
	g++ -std=c++11 -shared -o libmcsim.so Backend.o
