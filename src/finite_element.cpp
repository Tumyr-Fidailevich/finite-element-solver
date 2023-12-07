#include "finite_element.h"

Eigen::SparseMatrix<double>
BeamElement::calculateGlobalStiffnessMatrix(std::vector<Node> &nodes, std::unordered_map<int, Material> &materials) const
{
	auto rotateMatrix = calculateRotateMatrix(nodes, materials);
	auto localStiffnessMatrix = calculateLocalStiffnessMatrix(nodes, materials);
	return rotateMatrix.transpose() * localStiffnessMatrix * rotateMatrix;
}

Eigen::SparseMatrix<double>
BeamElement::calculateLocalStiffnessMatrix(std::vector<Node> &nodes, std::unordered_map<int, Material> &materials) const
{
	Eigen::SparseMatrix<double> stiffnessMatrix(6, 6);
	Node &left = nodes[nodesIds.first];
	Node &right = nodes[nodesIds.second];
	Material &material = materials[materialId];
	double length = std::sqrt((right.x - left.x) * (right.x - left.x) + (right.y - left.y) * (right.y - left.y));
	stiffnessMatrix.insert(0, 0) = material.youngModulus * material.area / length;
	stiffnessMatrix.insert(0, 3) = -material.youngModulus * material.area / length;
	stiffnessMatrix.insert(3, 3) = material.youngModulus * material.area / length;
	stiffnessMatrix.insert(3, 0) = -material.youngModulus * material.area / length;

	stiffnessMatrix.insert(1, 1) = 12 * material.youngModulus * material.inertiaMoment / (length * length * length);
	stiffnessMatrix.insert(1, 2) = 6 * material.youngModulus * material.inertiaMoment / (length * length);
	stiffnessMatrix.insert(2, 1) = 6 * material.youngModulus * material.inertiaMoment / (length * length);
	stiffnessMatrix.insert(2, 2) = 4 * material.youngModulus * material.inertiaMoment / length;

	stiffnessMatrix.insert(1, 4) = -12 * material.youngModulus * material.inertiaMoment / (length * length * length);
	stiffnessMatrix.insert(1, 5) = 6 * material.youngModulus * material.inertiaMoment / (length * length);
	stiffnessMatrix.insert(2, 4) = -6 * material.youngModulus * material.inertiaMoment / (length * length);
	stiffnessMatrix.insert(2, 5) = 2 * material.youngModulus * material.inertiaMoment / length;

	stiffnessMatrix.insert(4, 1) = -12 * material.youngModulus * material.inertiaMoment / (length * length * length);
	stiffnessMatrix.insert(4, 2) = -6 * material.youngModulus * material.inertiaMoment / (length * length);
	stiffnessMatrix.insert(5, 1) = 6 * material.youngModulus * material.inertiaMoment / (length * length);
	stiffnessMatrix.insert(5, 2) = 2 * material.youngModulus * material.inertiaMoment / length;

	stiffnessMatrix.insert(4, 4) = 12 * material.youngModulus * material.inertiaMoment / (length * length * length);
	stiffnessMatrix.insert(4, 5) = -6 * material.youngModulus * material.inertiaMoment / (length * length);
	stiffnessMatrix.insert(5, 4) = -6 * material.youngModulus * material.inertiaMoment / (length * length);
	stiffnessMatrix.insert(5, 5) = 4 * material.youngModulus * material.inertiaMoment / length;

	return stiffnessMatrix;
}

Eigen::SparseMatrix<double>
BeamElement::calculateRotateMatrix(std::vector<Node> &nodes, std::unordered_map<int, Material> &materials) const
{
	Eigen::SparseMatrix<double> rotateMatrix(6, 6);
	Node &left = nodes[nodesIds.first];
	Node &right = nodes[nodesIds.second];
	double cos = (right.x - left.x) /
				 (std::sqrt((right.x - left.x) * (right.x - left.x) + (right.y - left.y) * (right.y - left.y)));
	double sin = (right.y - left.y) /
				 (std::sqrt((right.x - left.x) * (right.x - left.x) + (right.y - left.y) * (right.y - left.y)));
	for (int i = 0; i < 3; ++i)
	{
		rotateMatrix.insert(i, 0) = cos;
		rotateMatrix.insert(i, 1) = sin;
		rotateMatrix.insert(i + 3, 0) = -sin;
		rotateMatrix.insert(i + 3, 1) = cos;

		rotateMatrix.insert(i, 2) = 0.0;
		rotateMatrix.insert(i + 3, 2) = 0.0;

		for (int j = 3; j < 6; ++j)
		{
			rotateMatrix.insert(i, j) = 0.0;
			rotateMatrix.insert(i + 3, j) = 0.0;
		}
	}
	return rotateMatrix;
}


