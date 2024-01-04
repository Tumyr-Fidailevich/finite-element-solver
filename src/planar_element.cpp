#include "planar_element.h"

void PlanarElement::calculateStiffnessMatrix(std::vector<Node> &nodes, std::unordered_map<int, Material> &materials)
{
	auto &material = materials[materialId];

	Node &lowerLeft = nodes[posToID[NodePosition::LowerLeft]];
	Node &lowerRight = nodes[posToID[NodePosition::LowerRight]];
	Node &upperLeft = nodes[posToID[NodePosition::UpperLeft]];
	Node &upperRight = nodes[posToID[NodePosition::UpperRight]];

	auto a = lowerRight.x - lowerLeft.x;
	auto b = upperLeft.y - lowerLeft.y;

	auto mu = (1 - material.poissonRatio) / 2;
	auto hardnessCoefficient = material.youngModulus * thickness / (1 - material.poissonRatio * material.poissonRatio);

	auto c = a / (6 * b);
	auto c_mu = mu * c;

	auto d = b / (6 * a);
	auto d_mu = mu * d;

	auto s = material.poissonRatio / 4;
	auto t = mu / 4;

	stiffnessMatrix << 2 * (d + c_mu), s + t, c_mu - 2 * d, s - t, -d - c_mu, -s - t, d - 2 * c_mu, t - s,
						s + t, 2 * (c + d_mu), t - s, -c - d_mu, -s - t, -c - d_mu, s - t, d_mu - 2 * c,
						c_mu - 2 * d, t - s, 2 * (d + c_mu), -s - t, d_mu - 2 * c, s - t, -d - c_mu, s + t,
						s - t, -c - d_mu, -s - t, 2 * (c + d_mu), t - s, d_mu - 2 * c, s + t, -c - d_mu,
						-d - c_mu, -s - t, d_mu - 2 * c, t - s, 2 * (d + c_mu), s + t, c_mu - 2 * d, s - t,
						-s - t, -c - d_mu, s - t, d_mu - 2 * c, s + t, 2 *(c + d_mu), t - s, c - 2 * d_mu,
						d - 2 * c_mu, s - t, -d - c_mu, s + t, c_mu - 2 * d, t - s, 2 *(d + c_mu), -s - t,
						t - s, d_mu - 2 * c, s + t, -c - d_mu, s - t, c - 2 * d_mu, -s - t, 2 * (c + d_mu);
	stiffnessMatrix = hardnessCoefficient * stiffnessMatrix;
}
