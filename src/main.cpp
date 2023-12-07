#include <iostream>
#include "solver.h"


int main(int argc, char* argv[])
{
	Solver solver("../mesh.txt");
	solver.solve();
	solver.saveToFile("../result.txt");
	std::cout << solver.getDisplacements() << std::endl;
	return 0;
}
