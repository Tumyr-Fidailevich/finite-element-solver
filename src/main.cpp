#include "solver.h"


int main(int argc, char* argv[])
{
	PlanarSolver solver("D:/finite-element-solver/meshes/planar_mesh.txt");
	solver.solve();
	solver.saveToFile("D:/finite-element-solver/results/planar_mesh_result.txt", PlanarSolver::ScaleStatus::On);
	std::cin.get();
	return 0;
}