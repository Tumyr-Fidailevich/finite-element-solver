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
	static const short DOF = 3;
	static const short ORDER = 2;
	Eigen::SparseMatrix<double> localStiffnessMatrix{DOF * ORDER, DOF * ORDER};
	Eigen::SparseMatrix<double> rotateMatrix{DOF * ORDER, DOF * ORDER};
	Eigen::SparseMatrix<double> globalStiffnessMatrix{DOF * ORDER, DOF * ORDER};
};


#endif