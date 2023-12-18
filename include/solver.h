#ifndef FINITE_ELEMENT_SOLVER_SOLVER_H
#define FINITE_ELEMENT_SOLVER_SOLVER_H


#include "pch.h"
#include "finite_element.h"

class Solver
{
public:
	explicit Solver(const std::string &path);

	~Solver();

	void saveToFile(const std::string &path) const;

	void calculateGlobalStiffnessMatrix();

	void solve();

	void showDisplacements() const;

	void showLocalStiffnessMatrix() const;

	void showRotateMatrix() const;

	void showGlobalStiffnessMatrix() const;

	void showQuasiDiagonalStiffnessMatrix() const;

	void showCorrespondenceMatrix() const;

	void resultsReport() const;

private:

	void initializeFromFile(const std::string &path);

	void setConstraint();

	Eigen::VectorXd loads;
	std::vector<Node> nodes;
	std::vector<BeamElement> elements;
	Eigen::SparseMatrix<double> stiffnessMatrixBeforeBoundaryConditions;
	Eigen::SparseMatrix<double> quasiDiagonalStiffnessMatrix;
	Eigen::SparseMatrix<double> correspondenceMatrix;
	Eigen::SparseMatrix<double> stiffnessMatrix;
	std::unordered_map<int, Material> materials;
	Eigen::VectorXd displacements;
};

#endif