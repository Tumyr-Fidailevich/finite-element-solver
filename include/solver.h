#ifndef BEAM_ELEMENT_SOLVER_SOLVER_H
#define BEAM_ELEMENT_SOLVER_SOLVER_H


#include "pch.h"
#include "beam_element.h"
#include "planar_element.h"


class BeamSolver
{
public:
	explicit BeamSolver(const std::string &path);

	~BeamSolver();

	void saveToFile(const std::string &path) const;

	void calculateGlobalStiffnessMatrix();

	void solve();

	void showDisplacements() const;

	void showLocalStiffnessMatrix() const;

	void showRotateMatrix() const;

	void showGlobalStiffnessMatrix() const;

	void showQuasiDiagonalStiffnessMatrix() const;

	void showCorrespondenceMatrix() const;

	void showNodalForces() const;

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


class PlanarSolver
{
public:
	enum class ScaleStatus : short
	{
		On,
		Off
	};

	explicit PlanarSolver(const std::string &path);

	~PlanarSolver();

	void saveToFile(const std::string &path, ScaleStatus scale = ScaleStatus::Off) const;

	void calculateGlobalStiffnessMatrix();

	void setCorrespondenceMatrix();

	void solve();

	void showDisplacements(std::ostream &stream = std::cout, ScaleStatus scale = ScaleStatus::Off) const;

	void showElementsGlobalStiffnessMatrix(std::ostream &stream = std::cout, ScaleStatus scale = ScaleStatus::Off) const;

	void showCorrespondenceMatrix(std::ostream &stream = std::cout) const;

	void showSystemStiffnessMatrix(std::ostream &stream = std::cout, ScaleStatus scale = ScaleStatus::Off) const;

	void showSystemStiffnessMatrixBeforeApplyingBoundaryConditions(std::ostream &stream = std::cout, ScaleStatus scale = ScaleStatus::Off) const;

	void showNodalForces(std::ostream &stream = std::cout) const;

	void showDeformations(std::ostream &stream = std::cout, ScaleStatus scale = ScaleStatus::Off) const;

	void showStrains(std::ostream &stream = std::cout, ScaleStatus scale = ScaleStatus::Off) const;

	void showShapeFunctions(std::ostream &stream = std::cout, ScaleStatus scale = ScaleStatus::Off) const;

	void resultsReport(std::ostream &stream = std::cout, ScaleStatus scale = ScaleStatus::Off) const;

private:

	void initializeFromFile(const std::string &path);

	void setConstraint();

	Eigen::VectorXd loads;
	std::vector<Node> nodes;
	std::vector<PlanarElement> elements;
	Eigen::MatrixX<double> stiffnessMatrixBeforeBoundaryConditions;
	Eigen::MatrixX<std::size_t> correspondenceMatrix;
	Eigen::MatrixX<double> stiffnessMatrix;
	std::unordered_map<int, Material> materials;
	Eigen::VectorXd displacements;
};

#endif