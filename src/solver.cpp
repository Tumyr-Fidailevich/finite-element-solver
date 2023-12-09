#include "solver.h"

Solver::Solver(const std::string &path)
{
	initializeFromFile(path);
	stiffnessMatrix.resize(3 * nodes.size(), 3 * nodes.size());
}

Solver::~Solver() = default;

void Solver::setBoundaryConditions()
{
	for(int i = 0; i < nodes.size(); ++i)
	{
		if(nodes[i].constraint & Node::ConstraintType::X)
		{
			for(int j = 0; j < stiffnessMatrix.cols(); ++j) stiffnessMatrix.coeffRef(3 * i, j) = 0;
			for(int j = 0; j < stiffnessMatrix.rows(); ++j) stiffnessMatrix.coeffRef(j, 3 * i) = 0;
		}
		if(nodes[i].constraint & Node::ConstraintType::Y)
		{
			for(int j = 0; j < stiffnessMatrix.cols(); ++j) stiffnessMatrix.coeffRef(3 * i + 1, j) = 0;
			for(int j = 0; j < stiffnessMatrix.rows(); ++j) stiffnessMatrix.coeffRef(j, 3 * i + 1) = 0;
		}
		if(nodes[i].constraint & Node::ConstraintType::Theta)
		{
			for(int j = 0; j < stiffnessMatrix.cols(); ++j) stiffnessMatrix.coeffRef(3 * i + 2, j) = 0;
			for(int j = 0; j < stiffnessMatrix.rows(); ++j) stiffnessMatrix.coeffRef(j, 3 * i + 2) = 0;
		}
	}
}

void Solver::calculateGlobalStiffnessMatrix()
{
	auto offset = 0;
	for (auto &element: elements)
	{
		element.calculateGlobalStiffnessMatrix(nodes, materials);
		auto size = element.globalStiffnessMatrix.outerSize();
		for (int k = 0; k < size; ++k)
		{
			for (Eigen::SparseMatrix<double>::InnerIterator it(element.globalStiffnessMatrix, k); it; ++it)
			{
				stiffnessMatrix.insert(offset + it.row(), offset + it.col()) += it.value();
			}
		}
		offset += size / 2;
	}
}

void Solver::solve()
{
	Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
	calculateGlobalStiffnessMatrix();
	setBoundaryConditions();
	solver.compute(stiffnessMatrix);
	displacements = solver.solve(loads);
}

void Solver::initializeFromFile(const std::string &path)
{
	std::fstream file(path);

	int materialsCount;
	file >> materialsCount;
	for (int _ = 0; _ < materialsCount; _++)
	{
		Material material{};
		int id;
		file >> id >> material.youngModulus >> material.inertiaMoment >> material.area;
		materials[id] = material;
	}

	int nodesCount;
	file >> nodesCount;
	for (int i = 0; i < nodesCount; i++)
	{
		Node node;
		file >> node.x >> node.y;
		nodes.emplace_back(node);
	}
	nodes.shrink_to_fit();

	int elementsCount;
	file >> elementsCount;
	for (int i = 0; i < elementsCount; i++)
	{
		BeamElement element;
		file >> element.materialId >> element.nodesIds.first >> element.nodesIds.second;
		elements.push_back(element);
	}
	elements.shrink_to_fit();

	int boundaryConditionsCount;
	file >> boundaryConditionsCount;
	for (int i = 0; i < boundaryConditionsCount; i++)
	{
		int cons;
		file >> cons;
		nodes[i].constraint |= cons == 1 ? Node::ConstraintType::X : 0;
		file >> cons;
		nodes[i].constraint |= cons == 1 ? Node::ConstraintType::Y : 0;
		file >> cons;
		nodes[i].constraint |= cons == 1 ? Node::ConstraintType::Theta : 0;
	}

	int loadsCount;
	file >> loadsCount;
	loads.resize(loadsCount);
	for (int i = 0; i < loadsCount; i++)
	{
		double X, Y, M;
		file >> X >> Y >> M;
		loads[i + 0] = X;
		loads[i + 1] = Y;
		loads[i + 2] = M;
	}
}

void Solver::saveToFile(const std::string &path)
{

}

Eigen::Vector2d Solver::getDisplacements() const
{
	return displacements;
}
