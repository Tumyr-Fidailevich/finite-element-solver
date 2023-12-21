#include "finite_element.h"

void BeamElement::calculateGlobalStiffnessMatrix(std::vector<Node> &nodes, std::unordered_map<int, Material> &materials)
{
	/**
 	* Расчет глобальной матрицы жесткости элемента
	*
	* @param nodes Массив узлов
	* @param materials хэш таблица с материалами
 	*/
	calculateLocalStiffnessMatrix(nodes, materials);
	calculateRotateMatrix(nodes, materials);
	globalStiffnessMatrix = rotateMatrix.transpose() * localStiffnessMatrix * rotateMatrix;
}

void BeamElement::calculateLocalStiffnessMatrix(std::vector<Node> &nodes, std::unordered_map<int, Material> &materials)
{
	/**
 	* Расчет локальной матрицы жесткости элемента
	*
	* @param nodes Массив узлов
	* @param materials хэш таблица с материалами
 	*/
	Node &left = nodes[nodesIds.first];
	Node &right = nodes[nodesIds.second];
	Material &material = materials[materialId];
	double length = std::sqrt((right.x - left.x) * (right.x - left.x) + (right.y - left.y) * (right.y - left.y));
	localStiffnessMatrix.insert(0, 0) = material.youngModulus * material.area / length;
	localStiffnessMatrix.insert(0, 3) = -material.youngModulus * material.area / length;
	localStiffnessMatrix.insert(3, 3) = material.youngModulus * material.area / length;
	localStiffnessMatrix.insert(3, 0) = -material.youngModulus * material.area / length;

	localStiffnessMatrix.insert(1, 1) =
			12 * material.youngModulus * material.inertiaMoment / (length * length * length);
	localStiffnessMatrix.insert(1, 2) = 6 * material.youngModulus * material.inertiaMoment / (length * length);
	localStiffnessMatrix.insert(2, 1) = 6 * material.youngModulus * material.inertiaMoment / (length * length);
	localStiffnessMatrix.insert(2, 2) = 4 * material.youngModulus * material.inertiaMoment / length;

	localStiffnessMatrix.insert(1, 4) =
			-12 * material.youngModulus * material.inertiaMoment / (length * length * length);
	localStiffnessMatrix.insert(1, 5) = 6 * material.youngModulus * material.inertiaMoment / (length * length);
	localStiffnessMatrix.insert(2, 4) = -6 * material.youngModulus * material.inertiaMoment / (length * length);
	localStiffnessMatrix.insert(2, 5) = 2 * material.youngModulus * material.inertiaMoment / length;

	localStiffnessMatrix.insert(4, 1) =
			-12 * material.youngModulus * material.inertiaMoment / (length * length * length);
	localStiffnessMatrix.insert(4, 2) = -6 * material.youngModulus * material.inertiaMoment / (length * length);
	localStiffnessMatrix.insert(5, 1) = 6 * material.youngModulus * material.inertiaMoment / (length * length);
	localStiffnessMatrix.insert(5, 2) = 2 * material.youngModulus * material.inertiaMoment / length;

	localStiffnessMatrix.insert(4, 4) =
			12 * material.youngModulus * material.inertiaMoment / (length * length * length);
	localStiffnessMatrix.insert(4, 5) = -6 * material.youngModulus * material.inertiaMoment / (length * length);
	localStiffnessMatrix.insert(5, 4) = -6 * material.youngModulus * material.inertiaMoment / (length * length);
	localStiffnessMatrix.insert(5, 5) = 4 * material.youngModulus * material.inertiaMoment / length;

}

void BeamElement::calculateRotateMatrix(std::vector<Node> &nodes, std::unordered_map<int, Material> &materials)
{
	/**
 	* Расчет матрицы перехода
	*
	* @param nodes Массив узлов
	* @param materials хэш таблица с материалами
 	*/
	Node &left = nodes[nodesIds.first];
	Node &right = nodes[nodesIds.second];
	double cos = (right.x - left.x) /
				 (std::sqrt((right.x - left.x) * (right.x - left.x) + (right.y - left.y) * (right.y - left.y)));
	double sin = (right.y - left.y) /
				 (std::sqrt((right.x - left.x) * (right.x - left.x) + (right.y - left.y) * (right.y - left.y)));
	rotateMatrix.insert(0, 0) = cos;
	rotateMatrix.insert(1, 1) = cos;
	rotateMatrix.insert(0, 1) = sin;
	rotateMatrix.insert(1, 0) = -sin;
	rotateMatrix.insert(2, 2) = 1;

	rotateMatrix.insert(3, 3) = cos;
	rotateMatrix.insert(4, 4) = cos;
	rotateMatrix.insert(3, 4) = sin;
	rotateMatrix.insert(4, 3) = -sin;
	rotateMatrix.insert(5, 5) = 1;
}


