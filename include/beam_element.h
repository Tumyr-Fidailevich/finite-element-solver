#ifndef FINITE_ELEMENT_SOLVER_BEAM_ELEMENT_H
#define FINITE_ELEMENT_SOLVER_BEAM_ELEMENT_H

#include "pch.h"
#include "material.h"
#include "node.h"


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