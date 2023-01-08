#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <fstream>
#include <iostream>

#include "lattice.h"
#include "mpi.h"

int main(int argc, char* argv[]) {
	std::map<std::string, double> parameters = {
		{"size", 600},	  {"iterations", 50}, {"rho", 0.7},	  {"beta", 2.0},
		{"alpha", 0.8},	  {"theta", 0.011},	  {"kappa", 0.1}, {"mu", 0.01},
		{"gamma", 0.001}, {"sigma", 0.001}};

	size_t iterations = parameters["iterations"];

	Lattice lattice(parameters["size"], parameters);

	MPI_Init(&argc, &argv);

	// printf("Main\n");
	// fflush(stdout);

	for (size_t i = 0; i < iterations; i++) {
		std::cout << "Iteration: " + std::to_string(i) + "/" +
						 std::to_string(iterations)
				  << "\r" << std::flush;

		// printf("In main\n");
		// fflush(stdout);

		lattice.diffuse();
		lattice.freeze();
		lattice.attach();
		lattice.melt();
		lattice.addNoise();

		if (i % 3 == 0)
			lattice.saveTo("frames/snowflake" + std::string("_") +
						   std::to_string(i));
	}
	lattice.saveTo("frames/snowflake");

	MPI_Finalize();

	return 0;
}
