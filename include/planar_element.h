#ifndef FINITE_ELEMENT_SOLVER_PLANAR_ELEMENT_H
#define FINITE_ELEMENT_SOLVER_PLANAR_ELEMENT_H

#include "pch.h"
#include "material.h"
#include "node.h"


struct PlanarElement
{
	static const short DOF = 2;
	static const short NODES_NUMBER = 4;
	static const short ORDER = 1;
	void calculateStiffnessMatrix(std::vector<Node> &nodes, std::unordered_map<int, Material> &materials);

	enum class NodePosition : short
	{
		LowerLeft,
		LowerRight,
		UpperLeft,
		UpperRight
	};

	int materialId;
	double thickness;
	std::map<NodePosition, std::size_t> posToID;
	Eigen::MatrixX <double> stiffnessMatrix{DOF * NODES_NUMBER, DOF * NODES_NUMBER};
};

#endif
