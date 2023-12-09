#include <iostream>
#include "solver.h"


int main(int argc, char* argv[])
{
	Solver solver("../meshes/mesh.txt");
	solver.solve();
	solver.saveToFile("../results/result.txt");
	solver.resultsReport();
	return 0;
}
