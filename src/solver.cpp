#include "solver.h"

BeamSolver::BeamSolver(const std::string &path)
{
	/**
	 * Инициализирует объект решателя по пути сетки
	 *
	 * @param path Путь к файлу с сеткой
	 */
	initializeFromFile(path);
	stiffnessMatrix.resize(BeamElement::DOF * nodes.size(), BeamElement::DOF * nodes.size());
	quasiDiagonalStiffnessMatrix.resize(elements.size() * BeamElement::DOF * BeamElement::ORDER,
										elements.size() * BeamElement::DOF * BeamElement::ORDER);
	correspondenceMatrix.resize(elements.size() * BeamElement::DOF * BeamElement::ORDER,
								nodes.size() * BeamElement::DOF);
}

BeamSolver::~BeamSolver() = default;

void BeamSolver::setConstraint()
{
	/**
 	* Устанавливает граничные условия
 	*/
	for (int i = 0; i < nodes.size(); ++i)
	{
		if (nodes[i].constraint & Node::ConstraintType::X)
		{
			for (int j = 0; j < stiffnessMatrix.cols(); ++j)
			{
				stiffnessMatrix.coeffRef(3 * i, j) = 0;
				stiffnessMatrix.coeffRef(j, 3 * i) = 0;
			}
			stiffnessMatrix.coeffRef(3 * i, 3 * i) = 1;
		}
		if (nodes[i].constraint & Node::ConstraintType::Y)
		{
			for (int j = 0; j < stiffnessMatrix.cols(); ++j)
			{
				stiffnessMatrix.coeffRef(3 * i + 1, j) = 0;
				stiffnessMatrix.coeffRef(j, 3 * i + 1) = 0;
			}
			stiffnessMatrix.coeffRef(3 * i + 1, 3 * i + 1) = 1;
		}
		if (nodes[i].constraint & Node::ConstraintType::Theta)
		{
			for (int j = 0; j < stiffnessMatrix.cols(); ++j)
			{
				stiffnessMatrix.coeffRef(3 * i + 2, j) = 0;
				stiffnessMatrix.coeffRef(j, 3 * i + 2) = 0;
			}
			stiffnessMatrix.coeffRef(3 * i + 2, 3 * i + 2) = 1;
		}
	}
}

void BeamSolver::calculateGlobalStiffnessMatrix()
{
	/**
 	* Расчет матрицы жесткости системы
 	*/
	auto offset = 0LL;
	for (auto &element: elements)
	{
		element.calculateGlobalStiffnessMatrix(nodes, materials);
		auto size = element.globalStiffnessMatrix.outerSize();
		for (int i = 0; i < size; i++)
		{
			for (int j = 0; j < size; j++)
			{
				quasiDiagonalStiffnessMatrix.coeffRef(offset + i, offset + j) =
						element.globalStiffnessMatrix.coeff(i, j);
			}
		}

		correspondenceMatrix.coeffRef(offset, 3 * element.nodesIds.first) = 1;
		correspondenceMatrix.coeffRef(offset + 1, 3 * element.nodesIds.first + 1) = 1;
		correspondenceMatrix.coeffRef(offset + 2, 3 * element.nodesIds.first + 2) = 1;

		correspondenceMatrix.coeffRef(offset + 3, 3 * element.nodesIds.second) = 1;
		correspondenceMatrix.coeffRef(offset + 4, 3 * element.nodesIds.second + 1) = 1;
		correspondenceMatrix.coeffRef(offset + 5, 3 * element.nodesIds.second + 2) = 1;

		offset += size;
	}
	stiffnessMatrix = correspondenceMatrix.transpose() * quasiDiagonalStiffnessMatrix * correspondenceMatrix;
	stiffnessMatrixBeforeBoundaryConditions = stiffnessMatrix;
}

void BeamSolver::solve()
{
	/**
 	* Запускает решение задачи
 	*/
	Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
	calculateGlobalStiffnessMatrix();
	setConstraint();
	solver.compute(stiffnessMatrix);
	displacements = solver.solve(loads);
}

void BeamSolver::initializeFromFile(const std::string &path)
{
	/**
 	* Инициализирует параметры решателя по исходному файлу
 	*
 	* @param path Путь к файлу с сеткой
 	*/
	std::fstream file(path);
	if (!file.is_open())
	{
		file.close();
		throw std::exception("Cant open mesh file");
	}

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
		nodes.push_back(node);
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
		nodes[i].constraint |= cons == 0 ? Node::ConstraintType::X : 0;
		file >> cons;
		nodes[i].constraint |= cons == 0 ? Node::ConstraintType::Y : 0;
		file >> cons;
		nodes[i].constraint |= cons == 0 ? Node::ConstraintType::Theta : 0;
	}

	int loadsCount;
	file >> loadsCount;
	loads.resize(3 * loadsCount);
	for (int i = 0; i < loadsCount; i++)
	{
		double X, Y, M;
		file >> X >> Y >> M;
		loads[3 * i + 0] = X;
		loads[3 * i + 1] = Y;
		loads[3 * i + 2] = M;
	}
	file.close();
}

void BeamSolver::saveToFile(const std::string &path) const
{
	/**
 	* Сохраняет результаты расчета в выбранный файл
 	*
 	* @param path Путь к файлу с результатами
 	*/
	std::ofstream file(path);
	if (file.is_open())
	{
		file << "X Displacements: " << std::endl;
		for (int i = 0; i < displacements.size(); i += 3) file << displacements[i] << std::endl;
		file << "Y Displacements: " << std::endl;
		for (int i = 1; i < displacements.size(); i += 3) file << displacements[i] << std::endl;
		file << "Angle Displacements: " << std::endl;
		for (int i = 2; i < displacements.size(); i += 3) file << displacements[i] << std::endl;

		for (int i = 0; i < elements.size(); ++i)
		{
			file << i + 1 << " elements local stiffness matrix:" << std::endl;
			file << 1e-6 * Eigen::MatrixXd(elements[i].localStiffnessMatrix) << std::endl;
		}

		for (int i = 0; i < elements.size(); ++i)
		{
			file << i + 1 << " elements rotate matrix:" << std::endl;
			file << Eigen::MatrixXd(elements[i].rotateMatrix) << std::endl;
		}

		for (int i = 0; i < elements.size(); ++i)
		{
			file << i + 1 << " elements global stiffness matrix:" << std::endl;
			file << 1e-6 * Eigen::MatrixXd(elements[i].globalStiffnessMatrix) << std::endl;
		}

		file << "Correspondence matrix: " << std::endl;
		file << Eigen::MatrixXd(correspondenceMatrix) << std::endl;

		file << "Quasi diagonal stiffness matrix: " << std::endl;
		file << 1e-6 * Eigen::MatrixXd(quasiDiagonalStiffnessMatrix) << std::endl;

		file << "Full system matrix before applying boundary conditions: " << std::endl;
		file << 1e-6 * Eigen::MatrixXd(stiffnessMatrixBeforeBoundaryConditions) << std::endl;

		file << "System stiffness matrix: " << std::endl;
		file << 1e-6 * Eigen::MatrixXd(stiffnessMatrix) << std::endl;

		file << "Nodal forces: " << std::endl;
		file << stiffnessMatrixBeforeBoundaryConditions * displacements << std::endl;
	} else
	{
		std::cerr << "Error during opening destination file" << std::endl;
	}
	file.close();
}

void BeamSolver::showDisplacements() const
{
	/**
 	* Отрисовывает перемещения
 	*/
	std::cout << "Displacements: " << std::endl;
	std::cout << displacements << std::endl;
}

void BeamSolver::showLocalStiffnessMatrix() const
{
	/**
 	* Отрисовывает локальную матрицу жесткости каждого элемента
 	*/
	for (int i = 0; i < elements.size(); ++i)
	{
		std::cout << i + 1 << " elements local stiffness matrix:" << std::endl;
		std::cout << Eigen::MatrixXd(elements[i].localStiffnessMatrix) << std::endl;
	}
}

void BeamSolver::showRotateMatrix() const
{
	/**
 	* Отрисовывает матрицу перехода в глобальную систему координат
 	*/
	for (int i = 0; i < elements.size(); ++i)
	{
		std::cout << i + 1 << " elements rotate matrix:" << std::endl;
		std::cout << Eigen::MatrixXd(elements[i].rotateMatrix) << std::endl;
	}
}

void BeamSolver::showGlobalStiffnessMatrix() const
{
	/**
 	* Отрисовывает глобальную матрицу жесткости
 	*/
	for (int i = 0; i < elements.size(); ++i)
	{
		std::cout << i + 1 << " elements global stiffness matrix:" << std::endl;
		std::cout << Eigen::MatrixXd(elements[i].globalStiffnessMatrix) << std::endl;
	}
}

void BeamSolver::showQuasiDiagonalStiffnessMatrix() const
{
	/**
 	* Отрисовывает квазидиагональную матрицу жесткости системы
 	*/
	std::cout << "Quasi diagonal stiffness matrix: " << std::endl;
	std::cout << Eigen::MatrixXd(quasiDiagonalStiffnessMatrix) << std::endl;
}

void BeamSolver::showCorrespondenceMatrix() const
{
	/**
 	* Отрисовывает матрицу соответствия
 	*/
	std::cout << "Correspondence matrix: " << std::endl;
	std::cout << Eigen::MatrixXd(correspondenceMatrix) << std::endl;
}

void BeamSolver::showNodalForces() const
{
	/**
 	* Отрисовывает узловые нагрузки
 	*/
	std::cout << "Nodal forces: " << std::endl;
	std::cout << stiffnessMatrixBeforeBoundaryConditions * displacements << std::endl;
}

void BeamSolver::resultsReport() const
{
	/**
 	* Отрисоывает результаты в консоль
 	*/
	showDisplacements();
	showLocalStiffnessMatrix();
	showRotateMatrix();
	showGlobalStiffnessMatrix();
	showQuasiDiagonalStiffnessMatrix();
	showCorrespondenceMatrix();
	std::cout << "Full system matrix before applying boundary conditions: " << std::endl;
	std::cout << Eigen::MatrixXd(stiffnessMatrixBeforeBoundaryConditions) << std::endl;
	std::cout << "System stiffness matrix: " << std::endl;
	std::cout << Eigen::MatrixXd(stiffnessMatrix);
}
