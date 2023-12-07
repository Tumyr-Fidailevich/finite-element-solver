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

	Eigen::SparseMatrix<double>
	calculateGlobalStiffnessMatrix(std::vector<Node> &nodes, std::unordered_map<int, Material> &materials) const;

	Eigen::SparseMatrix<double>
	calculateRotateMatrix(std::vector<Node> &nodes, std::unordered_map<int, Material> &materials) const;

	Eigen::SparseMatrix<double>
	calculateLocalStiffnessMatrix(std::vector<Node> &nodes, std::unordered_map<int, Material> &materials) const;

	std::pair<int, int> nodesIds;
	int materialId;
};


#endif