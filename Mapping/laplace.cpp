#include "laplace.h"

namespace jy_mesh {
	void cal_cot_angles(const Eigen::MatrixXd &V, const Eigen::Matrix3Xi &F,
		Eigen::Matrix3Xd &cot_angles)
	{
		cot_angles.resize(3, F.cols());
		for (size_t j = 0; j < F.cols(); ++j) {
			const Eigen::Vector3i &fv = F.col(j);
			for (size_t vi = 0; vi < 3; ++vi) {
				const Eigen::VectorXd &p0 = V.col(fv[vi]);
				const Eigen::VectorXd &p1 = V.col(fv[(vi + 1) % 3]);
				const Eigen::VectorXd &p2 = V.col(fv[(vi + 2) % 3]);
				const double angle = std::acos(std::max(-1.0,
					std::min(1.0, (p1 - p0).normalized().dot((p2 - p0).normalized()))));
				cot_angles(vi, j) = 1.0 / std::tan(angle);
			}
		}
	}

	int cal_cot_laplace(const Eigen::MatrixXd &V, const Eigen::Matrix3Xi &F,
		Eigen::SparseMatrix<double> &L)
	{
		Eigen::Matrix3Xd cot_angles;
		cal_cot_angles(V, F, cot_angles);
		std::vector<Eigen::Triplet<double>> triple;
		triple.reserve(F.cols() * 9);
		for (size_t j = 0; j < F.cols(); ++j) {
			const Eigen::Vector3i &fv = F.col(j);
			const Eigen::Vector3d &ca = cot_angles.col(j);
			for (size_t vi = 0; vi < 3; ++vi) {
				const size_t j1 = (vi + 1) % 3;
				const size_t j2 = (vi + 2) % 3;
				const int fv0 = fv[vi];
				const int fv1 = fv[j1];
				const int fv2 = fv[j2];
				triple.push_back(Eigen::Triplet<double>(fv0, fv0, ca[j1] + ca[j2]));
				triple.push_back(Eigen::Triplet<double>(fv0, fv1, -ca[j2]));
				triple.push_back(Eigen::Triplet<double>(fv0, fv2, -ca[j1]));
			}
		}
		L.resize(V.cols(), V.cols());
		L.setFromTriplets(triple.begin(), triple.end());
		return 1;
	}

	int cal_cot_laplace(const Eigen::MatrixXd &V, const Eigen::Matrix3Xi &F,
		const Eigen::VectorXi &g2l, Eigen::SparseMatrix<double> &L)
	{
		Eigen::Matrix3Xd cot_angles;
		cal_cot_angles(V, F, cot_angles);
		std::vector<Eigen::Triplet<double>> triple;
		triple.reserve(F.cols() * 9);
		Eigen::Vector3i fv;
		for (size_t j = 0; j < F.cols(); ++j) {
			const Eigen::Vector3d &ca = cot_angles.col(j);
			fv << g2l[F(0, j)], g2l[F(1, j)], g2l[F(2, j)];
			for (size_t vi = 0; vi < 3; ++vi) {
				const size_t j1 = (vi + 1) % 3;
				const size_t j2 = (vi + 2) % 3;
				const int fv0 = fv[vi];
				const int fv1 = fv[j1];
				const int fv2 = fv[j2];
				triple.push_back(Eigen::Triplet<double>(fv0, fv0, ca[j1] + ca[j2]));
				triple.push_back(Eigen::Triplet<double>(fv0, fv1, -ca[j2]));
				triple.push_back(Eigen::Triplet<double>(fv0, fv2, -ca[j1]));
			}
		}
		L.resize(V.cols(), V.cols());
		L.setFromTriplets(triple.begin(), triple.end());
		return 1;
	}

	int cal_topo_laplace(const Eigen::MatrixXd &V, const Eigen::Matrix3Xi &F,
		Eigen::SparseMatrix<double> &L)
	{
		const size_t num_faces = F.cols();
		std::vector<Eigen::Triplet<double>> triple;
		triple.reserve(num_faces * 9);
		for (size_t j = 0; j < num_faces; ++j) {
			const Eigen::Vector3i &fv = F.col(j);
			for (size_t vi = 0; vi < 3; ++vi) {
				const int fv0 = fv[vi];
				const int fv1 = fv[(vi + 1) % 3];
				const int fv2 = fv[(vi + 2) % 3];
				triple.push_back(Eigen::Triplet<double>(fv0, fv0, 1));
				triple.push_back(Eigen::Triplet<double>(fv0, fv1, -0.5));
				triple.push_back(Eigen::Triplet<double>(fv0, fv2, -0.5));
			}
		}
		L.resize(V.cols(), V.cols());
		L.setFromTriplets(triple.begin(), triple.end());
		return 1;
	}

	int cal_topo_laplace(const Eigen::MatrixXd &V, const Eigen::Matrix3Xi &F,
		const Eigen::VectorXi &g2l, Eigen::SparseMatrix<double> &L)
	{
		const size_t num_faces = F.cols();
		std::vector<Eigen::Triplet<double>> triple;
		triple.reserve(num_faces * 9);
		Eigen::Vector3i fv;
		for (size_t j = 0; j < num_faces; ++j) {
			fv << g2l[F(0, j)], g2l[F(1, j)], g2l[F(2, j)];
			for (size_t vi = 0; vi < 3; ++vi) {
				const int fv0 = fv[vi];
				const int fv1 = fv[(vi + 1) % 3];
				const int fv2 = fv[(vi + 2) % 3];
				triple.push_back(Eigen::Triplet<double>(fv0, fv0, 1));
				triple.push_back(Eigen::Triplet<double>(fv0, fv1, -0.5));
				triple.push_back(Eigen::Triplet<double>(fv0, fv2, -0.5));
			}
		}
		L.resize(V.cols(), V.cols());
		L.setFromTriplets(triple.begin(), triple.end());
		return 1;
	}

}