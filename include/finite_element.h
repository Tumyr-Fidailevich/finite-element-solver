#ifndef FINITE_ELEMENT_SOLVER_FINITE_ELEMENT_H
#define FINITE_ELEMENT_SOLVER_FINITE_ELEMENT_H

#include "pch.h"

struct Node
{
	double x = 0;
	double y = 0;
	enum ConstraintType
	{
		None = 0,
		X = 1 << 0,
		Y = 1 << 1,
		Theta = 1 << 2
	};
	int constraint = None;
};

struct Material
{
	double youngModulus;
	double inertiaMoment;
	double area;
};

struct BeamElement
{

	void calculateGlobalStiffnessMatrix(std::vector<Node> &nodes, std::unordered_map<int, Material> &materials);

	void calculateRotateMatrix(std::vector<Node> &nodes, std::unordered_map<int, Material> &materials);

	void calculateLocalStiffnessMatrix(std::vector<Node> &nodes, std::unordered_map<int, Material> &materials);

	std::pair<int, int> nodesIds;
	int materialId;
	Eigen::SparseMatrix<double> localStiffnessMatrix{6, 6};
	Eigen::SparseMatrix<double> rotateMatrix{6, 6};
	Eigen::SparseMatrix<double> globalStiffnessMatrix{6, 6};
};


#endif