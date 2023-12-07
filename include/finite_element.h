#ifndef FINITE_ELEMENT_SOLVER_STRUCTS_H
#define FINITE_ELEMENT_SOLVER_STRUCTS_H

#include "pch.h"

struct Element
{

	void calculateStiffnessMatrix(const Eigen::Matrix3f& D, std::vector<Eigen::Triplet<float> >& triplets);

	Eigen::Matrix<float, 3, 6> B;
	int nodesIds[3];
};

struct Constraint
{
	enum Type
	{
		UX = 1 << 0,
		UY = 1 << 1,
		UXY = UX | UY
	};
	int node;
	Type type;
};


#endif