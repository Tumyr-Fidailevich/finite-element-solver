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

PlanarSolver::PlanarSolver(const std::string &path)
{
	/**
	 * Инициализирует объект решателя по пути сетки
	 *
	 * @param path Путь к файлу с сеткой
	 */
	initializeFromFile(path);
	stiffnessMatrix.resize(PlanarElement::NODES_NUMBER * elements.size(), PlanarElement::NODES_NUMBER * elements.size());
	stiffnessMatrix.setZero();
	correspondenceMatrix.resize(elements.size(), PlanarElement::NODES_NUMBER);
}

PlanarSolver::~PlanarSolver() = default;

void PlanarSolver::saveToFile(const std::string &path, ScaleStatus scale) const
{
	/**
	 * Сохраняет результаты расчетов в файл
	 *
	 * @param path путь к файлу с результатами
	 * @param scale включить/отключить скалирование результатов
	 */
	std::ofstream fileStream(path);

	if (!fileStream.is_open())
	{
		throw std::exception("Cant open mesh file");
	}

	resultsReport(fileStream, scale);
	fileStream.close();
}

void PlanarSolver::calculateGlobalStiffnessMatrix()
{
	/**
	 * Расчет матрицы жесткости системы
	 */
	for (auto &elem: elements)
	{
		elem.calculateStiffnessMatrix(nodes, materials);
	}
	for (int k = 0; k < correspondenceMatrix.rows(); k++)
	{
		auto &elemStiffnessMatrix = elements[k].stiffnessMatrix;
		for (int i = 0; i < correspondenceMatrix.cols(); i++)
		{
			for (int j = 0; j < correspondenceMatrix.cols(); j++)
			{
				auto it1 = 2 * i;
				auto it2 = 2 * i + 1;
				auto jt1 = 2 * j;
				auto jt2 = 2 * j + 1;

				auto ig1 = 2 * correspondenceMatrix(k, i);
				auto ig2 = 2 * correspondenceMatrix(k, i) + 1;
				auto jg1 = 2 * correspondenceMatrix(k, j);
				auto jg2 = 2 * correspondenceMatrix(k, j) + 1;

				stiffnessMatrix(ig1, jg1) = stiffnessMatrix(ig1, jg1) + elemStiffnessMatrix(it1, jt1);
				stiffnessMatrix(ig1, jg2) = stiffnessMatrix(ig1, jg2) + elemStiffnessMatrix(it1, jt2);
				stiffnessMatrix(ig2, jg1) = stiffnessMatrix(ig2, jg1) + elemStiffnessMatrix(it2, jt1);
				stiffnessMatrix(ig2, jg2) = stiffnessMatrix(ig2, jg2) + elemStiffnessMatrix(it2, jt2);
			}
		}
	}
	stiffnessMatrixBeforeBoundaryConditions = stiffnessMatrix;
}

void PlanarSolver::solve()
{
	/**
 	* Запускает решение задачи
 	*/
 	setCorrespondenceMatrix();
	calculateGlobalStiffnessMatrix();
	setConstraint();
	Eigen::ColPivHouseholderQR<Eigen::MatrixX<double>> solver(stiffnessMatrix);
	displacements = solver.solve(loads);
}

void PlanarSolver::initializeFromFile(const std::string &path)
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
		file >> id >> material.youngModulus >> material.inertiaMoment >> material.area >> material.poissonRatio;
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
		PlanarElement element;
		file >> element.materialId  >> element.thickness
			 >> element.posToID[PlanarElement::NodePosition::LowerLeft]
			 >> element.posToID[PlanarElement::NodePosition::LowerRight]
			 >> element.posToID[PlanarElement::NodePosition::UpperRight]
			 >> element.posToID[PlanarElement::NodePosition::UpperLeft];
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
	}

	int loadsCount;
	file >> loadsCount;
	loads.resize(2 * loadsCount);
	for (int i = 0; i < loadsCount; i++)
	{
		file >> loads[2 * i + 0] >> loads[2 * i + 1];
	}
	file.close();
}

void PlanarSolver::setConstraint()
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
				stiffnessMatrix.coeffRef(PlanarElement::DOF * i, j) = 0;
				stiffnessMatrix.coeffRef(j, PlanarElement::DOF * i) = 0;
			}
			stiffnessMatrix.coeffRef(PlanarElement::DOF * i, PlanarElement::DOF * i) = 1;
		}
		if (nodes[i].constraint & Node::ConstraintType::Y)
		{
			for (int j = 0; j < stiffnessMatrix.cols(); ++j)
			{
				stiffnessMatrix.coeffRef(PlanarElement::DOF * i + 1, j) = 0;
				stiffnessMatrix.coeffRef(j, PlanarElement::DOF * i + 1) = 0;
			}
			stiffnessMatrix.coeffRef(PlanarElement::DOF * i + 1, PlanarElement::DOF * i + 1) = 1;
		}
	}
}

void PlanarSolver::resultsReport(std::ostream &stream, ScaleStatus scale) const
{
	/**
 	* Отображает результаты в выбранный поток
 	*
 	* @param stream Поток вывода результатов
	* @param scale Включить/отключить скалирование результатов
 	*/
	showDisplacements(stream, scale);

	showElementsGlobalStiffnessMatrix(stream, scale);

	showCorrespondenceMatrix(stream);

	showSystemStiffnessMatrixBeforeApplyingBoundaryConditions(stream, scale);

	showSystemStiffnessMatrix(stream, scale);

	showNodalForces(stream);

	showDeformations(stream, scale);

	showStrains(stream, scale);

	showShapeFunctions(stream, scale);
}

void PlanarSolver::showDisplacements(std::ostream &stream, ScaleStatus scale) const
{
	/**
 	* Отображает перемещения в выбранный поток
 	*
 	* @param stream Поток вывода результатов
	* @param scale Включить/отключить скалирование результатов
 	*/
 	double scalar = 1;
	if(scale == ScaleStatus::On) scalar = 1e3;
	stream << "X Displacements: " << std::endl;
	for (int i = 0; i < displacements.size(); i += PlanarElement::DOF) stream << displacements[i] * scalar << std::endl;

	stream << "Y Displacements: " << std::endl;
	for (int i = 1; i < displacements.size(); i += PlanarElement::DOF) stream << displacements[i] * scalar << std::endl;
}

void PlanarSolver::showElementsGlobalStiffnessMatrix(std::ostream &stream, ScaleStatus scale) const
{
	/**
 	* Отображает матрицы жесткости элементов в выбранный поток
 	*
 	* @param stream Поток вывода результатов
 	*/
	for (int i = 0; i < elements.size(); ++i)
	{
		stream << i + 1 << " elements global stiffness matrix:" << std::endl;
		double scalar = 1;
		if(scale == ScaleStatus::On) scalar = materials.at(elements[i].materialId).youngModulus * elements[i].thickness;
		stream << elements[i].stiffnessMatrix / scalar << std::endl;
	}
}

void PlanarSolver::showCorrespondenceMatrix(std::ostream &stream) const
{
	/**
 	* Отображает матрицу соответствия в выбранный поток
 	*
 	* @param stream Поток вывода результатов
 	*/
	stream << "Correspondence matrix: " << std::endl;
	stream << correspondenceMatrix << std::endl;
}

void PlanarSolver::showSystemStiffnessMatrix(std::ostream &stream, ScaleStatus scale) const
{
	/**
 	* Отображает матрицу жесткости системы в выбранный поток
 	*
 	* @param stream Поток вывода результатов
	* @param scale Включить/отключить скалирование результатов
 	*/
 	double scalar = 1;
	if(scale == ScaleStatus::On) scalar = materials.at(0).youngModulus * elements.at(0).thickness;
	stream << "System stiffness matrix: " << std::endl;
	stream << stiffnessMatrix / scalar << std::endl;
}

void PlanarSolver::showSystemStiffnessMatrixBeforeApplyingBoundaryConditions(std::ostream &stream, ScaleStatus scale) const
{
	/**
 	* Отображает матрицу жесткости системы до применения граничных условий в выбранный поток
 	*
 	* @param stream Поток вывода результатов
	* @param scale Включить/отключить скалирование результатов
 	*/
	double scalar = 1;
	if(scale == ScaleStatus::On) scalar = materials.at(0).youngModulus * elements.at(0).thickness;
	stream << "System stiffness matrix before applying boundary conditions" << std::endl;
	stream << stiffnessMatrixBeforeBoundaryConditions / scalar << std::endl;
}

void PlanarSolver::showDeformations(std::ostream &stream, ScaleStatus scale) const
{
	/**
 	* Отображает деформации в выбранный поток
 	*
 	* @param stream Поток вывода результатов
	* @param scale Включить/отключить скалирование результатов
 	*/
	stream << "Elements deformations: " << std::endl;
	for (int i = 0; i < elements.size(); i++)
	{
		auto &lowerLeftNode = nodes.at(elements[i].posToID.at(PlanarElement::NodePosition::LowerLeft));
		auto &lowerRightNode = nodes.at(elements[i].posToID.at(PlanarElement::NodePosition::LowerRight));
		auto &upperLeftNode = nodes.at(elements[i].posToID.at(PlanarElement::NodePosition::UpperLeft));
		auto &upperRightNode = nodes.at(elements[i].posToID.at(PlanarElement::NodePosition::UpperRight));

		auto a = lowerRightNode.x - lowerLeftNode.x;
		auto b = upperLeftNode.y - lowerLeftNode.y;
		auto c = a / (6 * b);
		auto d = b / (6 * a);

		auto m = a / b;
		auto node = nodes[i];
		auto ksi = node.x / a;
		auto eta = node.y / b;

		auto x = 1 / a * (-d * displacements[2 * correspondenceMatrix(i, 0)] +
						  d * displacements[2 * correspondenceMatrix(i, 1)] +
						  eta * displacements[2 * correspondenceMatrix(i, 2)] -
						  eta * displacements[2 * correspondenceMatrix(i, 3)]);
		auto y = m / a * (c * displacements[2 * correspondenceMatrix(i, 0) + 1] -
						  ksi * displacements[2 * correspondenceMatrix(i, 1) + 1] +
						  ksi * displacements[2 * correspondenceMatrix(i, 2) + 1] +
						  c * displacements[2 * correspondenceMatrix(i, 3) + 1]);
		auto angle = m / a * (-c * displacements[2 * correspondenceMatrix(i, 0)] -
							  ksi * displacements[2 * correspondenceMatrix(i, 1)] +
							  ksi * displacements[2 * correspondenceMatrix(i, 2)] +
							  c * displacements[2 * correspondenceMatrix(i, 3)]);

		double scalar = 1;
		if(scale == ScaleStatus::On) scalar = 1e3;
		stream << x * scalar << ' ' << y * scalar << ' ' << angle * scalar << std::endl;
	}
}

void PlanarSolver::showStrains(std::ostream &stream, ScaleStatus scale) const
{
	/**
 	* Отображает напряжений в выбранный поток
 	*
 	* @param stream Поток вывода результатов
	* @param scale Включить/отключить скалирование результатов
 	*/
	stream << "Elements strains: " << std::endl;
	for (int i = 0; i < elements.size(); i++)
	{
		auto &lowerLeftNode = nodes.at(elements[i].posToID.at(PlanarElement::NodePosition::LowerLeft));
		auto &lowerRightNode = nodes.at(elements[i].posToID.at(PlanarElement::NodePosition::LowerRight));
		auto &upperLeftNode = nodes.at(elements[i].posToID.at(PlanarElement::NodePosition::UpperLeft));
		auto &upperRightNode = nodes.at(elements[i].posToID.at(PlanarElement::NodePosition::UpperRight));

		auto &nu = materials.at(elements[i].materialId).poissonRatio;
		auto &E = materials.at(elements[i].materialId).youngModulus;

		auto a = lowerRightNode.x - lowerLeftNode.x;
		auto b = upperLeftNode.y - lowerLeftNode.y;;

		auto mu = (1 - nu) / 2;
		auto m = a / b;
		auto node = nodes[i];
		auto ksi = node.x / a;
		auto eta = node.y / b;
		auto theta = 1 - ksi;
		auto hi = 1 - eta;


		auto x = E / (a * (1 - nu * nu)) * (-hi * displacements[2 * correspondenceMatrix(i, 0)] +
											hi * displacements[2 * correspondenceMatrix(i, 1)] +
											eta * displacements[2 * correspondenceMatrix(i, 2)] -
											eta * displacements[2 * correspondenceMatrix(i, 3)] +
											m * nu * (-theta * displacements[2 * correspondenceMatrix(i, 0) + 1] -
													  ksi * displacements[2 * correspondenceMatrix(i, 1) + 1] +
													  ksi * displacements[2 * correspondenceMatrix(i, 2) + 1] +
													  theta * displacements[2 * correspondenceMatrix(i, 3) + 1]));
		auto y = E / (a * (1 - nu * nu)) * (nu * (-hi * displacements[2 * correspondenceMatrix(i, 0)] +
												  hi * displacements[2 * correspondenceMatrix(i, 1)] +
												  eta * displacements[2 * correspondenceMatrix(i, 2)] -
												  eta * displacements[2 * correspondenceMatrix(i, 3)]) +
											m * (-theta * displacements[2 * correspondenceMatrix(i, 0) + 1] -
												 ksi * displacements[2 * correspondenceMatrix(i, 1) + 1] +
												 ksi * displacements[2 * correspondenceMatrix(i, 2) + 1] +
												 theta * displacements[2 * correspondenceMatrix(i, 3)]));
		auto angle = E / (a * (1 - nu * nu)) * (mu * m * (-theta * displacements[2 * correspondenceMatrix(i, 0)] -
														  ksi * displacements[2 * correspondenceMatrix(i, 1)] +
														  ksi * displacements[2 * correspondenceMatrix(i, 2)] +
														  eta * displacements[2 * correspondenceMatrix(i, 3)]) +
												mu * (-hi * displacements[2 * correspondenceMatrix(i, 0) + 1] +
													 hi * displacements[2 * correspondenceMatrix(i, 1) + 1] +
													 eta * displacements[2 * correspondenceMatrix(i, 2) + 1] -
													 eta * displacements[2 * correspondenceMatrix(i, 3) + 1]));

		double scalar = 1;
		if(scale == ScaleStatus::On) scalar = 1e6;
		stream << x / scalar << ' ' << y / scalar << ' ' << angle / scalar << std::endl;
	}
}

void PlanarSolver::showShapeFunctions(std::ostream &stream, ScaleStatus scale) const
{
	/**
 	* Отображает функции формы в выбранный поток
 	*
 	* @param stream Поток вывода результатов
	* @param scale Включить/отключить скалирование результатов
 	*/
	stream << "Elements shape functions: " << std::endl;
	double scalar = 1;
	if(scale == ScaleStatus::On) scalar = 1e3;
	for (int i = 0; i < nodes.size(); i++)
		stream << nodes[i].x << ' ' << nodes[i].y << ' ' << displacements[i + 1] * scalar << std::endl;
}

void PlanarSolver::showNodalForces(std::ostream &stream) const
{
	/**
 	* Отображает узловые реакции в выбранный поток
 	*
 	* @param stream Поток вывода результатов
 	*/
	stream << "Nodal forces: " << std::endl;
	stream << stiffnessMatrixBeforeBoundaryConditions * displacements << std::endl;
}

void PlanarSolver::setCorrespondenceMatrix()
{
	/**
 	* Определяет матрицу соответствия
 	*/
	for (int i = 0; i < elements.size(); i++)
	{
		correspondenceMatrix(i, 0) = elements[i].posToID[PlanarElement::NodePosition::LowerLeft];
		correspondenceMatrix(i, 1) = elements[i].posToID[PlanarElement::NodePosition::LowerRight];
		correspondenceMatrix(i, 2) = elements[i].posToID[PlanarElement::NodePosition::UpperRight];
		correspondenceMatrix(i, 3) = elements[i].posToID[PlanarElement::NodePosition::UpperLeft];
	}
}
