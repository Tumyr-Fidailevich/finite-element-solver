#include "solver.h"


int main(int argc, char* argv[])
{
	Solver solver("D:/finite-element-solver/meshes/mesh8.txt");
	solver.solve();
	solver.saveToFile("D:/finite-element-solver/results/mesh8_result_1e6.txt");
//	std::cin.get();
	return 0;
}
