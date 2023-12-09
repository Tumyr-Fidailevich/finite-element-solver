#ifndef FINITE_ELEMENT_SOLVER_SOLVER_H
#define FINITE_ELEMENT_SOLVER_SOLVER_H


#include "pch.h"
#include "finite_element.h"

class Solver
{
public:
	explicit Solver(const std::string &path);

	~Solver();

	void saveToFile(const std::string &path);

	void calculateGlobalStiffnessMatrix();

	void solve();

	void showDisplacements() const;

	void showLocalStiffnessMatrix() const;

	void showRotateMatrix() const;

	void showGlobalStiffnessMatrix() const;

	void resultsReport() const;

private:

	void initializeFromFile(const std::string &path);

	void setBoundaryConditions();

	Eigen::Vector2d loads;
	std::vector<Node> nodes;
	std::vector<BeamElement> elements;
	Eigen::SparseMatrix<double> stiffnessMatrixBeforeBoundaryConditions;
	Eigen::SparseMatrix<double> stiffnessMatrix;
	std::unordered_map<int, Material> materials;
	Eigen::Vector2d displacements;
};

#endif