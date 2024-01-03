#ifndef FINITE_ELEMENT_SOLVER_MATERIAL_H
#define FINITE_ELEMENT_SOLVER_MATERIAL_H

struct Material
{
	double youngModulus;
	double inertiaMoment;
	double area;
	double poissonRatio;
};

#endif
