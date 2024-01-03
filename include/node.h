#ifndef FINITE_ELEMENT_SOLVER_NODE_H
#define FINITE_ELEMENT_SOLVER_NODE_H

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

#endif //FINITE_ELEMENT_SOLVER_NODE_H
