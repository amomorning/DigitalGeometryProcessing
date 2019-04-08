#ifndef SINGLE_PATCH_PARAMETERIZATION_H
#define SINGLE_PATCH_PARAMETERIZATION_H

#include <vector>
#include <Eigen/Dense>

namespace jy_mesh {
	int count_flip_faces(const Eigen::Matrix2Xd &V, const Eigen::Matrix3Xi &F);

	int harmonic_mapping(const Eigen::MatrixXd &V, const Eigen::Matrix3Xi &F,
		const std::vector<std::pair<size_t, Eigen::Vector2d>> &bc,
		Eigen::Matrix3Xd &Vout);

	int lscm_mapping(const Eigen::MatrixXd &V, const Eigen::Matrix3Xi &F,
		const std::vector<std::pair<size_t, Eigen::Vector2d>> &bc,
		Eigen::Matrix2Xd &Vout);

	int arap_mapping(const Eigen::MatrixXd &V, const Eigen::Matrix3Xi &F,
		const std::vector<std::pair<size_t, Eigen::Vector2d>> &bc,
		Eigen::Matrix2Xd &Vout, size_t max_iter = 1000);

	int MVC_mapping(const Eigen::MatrixXd &V, const Eigen::Matrix3Xi &F,
		const std::vector<std::pair<size_t, Eigen::Vector2d>> &bc,
		Eigen::Matrix2Xd &Vout);
}

#endif
