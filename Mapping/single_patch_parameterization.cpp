#include <surface_mesh/Surface_mesh.h>
#include "tri_mesh.h"
#include "single_patch_parameterization.h"
#include "laplace.h"
//#include "util.h"

namespace jy_mesh {
	int count_flip_faces(const Eigen::Matrix2Xd &V, const Eigen::Matrix3Xi &F)
	{
		size_t count = 0;
		for (size_t i = 0; i < F.cols(); ++i) {
			//count += util_helper::signed_area(V.col(F(0, i)), V.col(F(1, i)), V.col(F(2, i))) > 0 ? 0 : 1;
		}
		return count;
	}

	void init_para_info(const Eigen::MatrixXd &V, const Eigen::Matrix3Xi &F,
		Eigen::VectorXd &sqrt_areas, std::vector<Eigen::Matrix2d> &inv,
		Eigen::SparseMatrix<double> &Gt)
	{
		Eigen::Matrix2d   m;
		Eigen::Matrix<double, 2, 3>  t;
		sqrt_areas.resize(F.cols());
		inv.reserve(F.cols());
		std::vector<Eigen::Triplet<double>>  triplets;
		triplets.reserve(F.cols() * 6);
		for (size_t i = 0, col = 0; i < F.cols(); ++i) {
			const double A = 0;
			//const double A = util_helper::flatten_triangle(V.col(F(0, i)), V.col(F(1, i)), V.col(F(2, i)), t);
			for (size_t j = 0; j < 2; ++j) {
				m.col(j) = t.col(j + 1) - t.col(0);
			}
			sqrt_areas[i] = sqrt(A);
			inv.push_back(m.inverse());
			for (size_t j = 0; j < 2; ++j, ++col) {
				double sum = 0;
				for (size_t k = 0; k < 2; ++k) {
					const double val = inv[i](k, j) * sqrt_areas[i];
					sum -= val;
					triplets.push_back(Eigen::Triplet<double>(F(k + 1, i), col, val));
				}
				triplets.push_back(Eigen::Triplet<double>(F(0, i), col, sum));
			}

		}
		//build the matrix
		Gt.resize(V.cols(), 2 * F.cols());
		Gt.setFromTriplets(triplets.begin(), triplets.end());
	}

	int harmonic_mapping(const Eigen::MatrixXd &V, const Eigen::Matrix3Xi &F,
		const std::vector<std::pair<size_t, Eigen::Vector2d>> &bc,
		Eigen::Matrix3Xd &Vout)
	{
		Eigen::SparseMatrix<double> L;
		jy_mesh::cal_cot_laplace(V, F, L);
		const double w = 1e6;
		Eigen::MatrixX2d b = Eigen::MatrixX2d::Zero(V.cols(), 3);
		for (size_t i = 0; i < bc.size(); ++i) {
			L.coeffRef(bc[i].first, bc[i].first) += w;
			b.row(bc[i].first) = w * bc[i].second;
		}
		Eigen::SimplicialLLT<Eigen::SparseMatrix<double>>  llt(L);
		if (llt.info() == Eigen::Success) {
			Vout = llt.solve(b).transpose();
			return 1;
		}
		std::cerr << "solve failure!" << std::endl;
		return 0;
	}

	int MVC_mapping(const Eigen::MatrixXd &V, const Eigen::Matrix3Xi &F,
		const std::vector<std::pair<size_t, Eigen::Vector2d>> &bc,
		Eigen::Matrix2Xd &Vout)
	{
		Eigen::Matrix3Xd tan_angles(3, F.cols());
		for (size_t j = 0; j < F.cols(); ++j) {
			const Eigen::Vector3i &fv = F.col(j);
			for (size_t vi = 0; vi < 3; ++vi) {
				const Eigen::VectorXd &p0 = V.col(fv[vi]);
				const Eigen::VectorXd &p1 = V.col(fv[(vi + 1) % 3]);
				const Eigen::VectorXd &p2 = V.col(fv[(vi + 2) % 3]);
				const double angle = std::acos(std::max(-1.0,
					std::min(1.0, (p1 - p0).normalized().dot((p2 - p0).normalized()))));
				tan_angles(vi, j) = std::tan(angle * 0.5);
			}
		}
		Eigen::VectorXi  flags(V.cols());
		flags.setZero();
		for (size_t i = 0; i < bc.size(); ++i) {
			flags[bc[i].first] = 1;
		}
		surface_mesh::Surface_mesh mesh;
		build_mesh(V, F, mesh);
		std::vector<Eigen::Triplet<double>> triple;
		triple.reserve(F.cols() * 6);
		surface_mesh::Surface_mesh::Halfedge_around_vertex_circulator vh_it, vh_end;
		for (auto const &vit : mesh.vertices()) {
			const int vi = vit.idx();
			if (flags[vi])
				continue;
			vh_it = vh_end = mesh.halfedges(vit);
			std::vector<std::pair<int, double>>  wij;
			double sum = 0;
			do {
				const int vj = mesh.to_vertex(*vh_it).idx();
				int fi = mesh.face(*vh_it).idx();
				double val = 0;
				for (size_t i = 0; i < 2; ++i) {
					if (fi >= 0) {
						for (size_t j = 0; j < 3; ++j) {
							if (F(j, fi) == vi) {
								val += tan_angles(j, fi);
								break;
							}
						}
					}
					fi = mesh.face(mesh.opposite_halfedge(*vh_it)).idx();
				}
				val /= (V.col(vi) - V.col(vj)).norm();
				sum += val;
				triple.push_back(Eigen::Triplet<double>(vi, vj, val));
			} while (++vh_it != vh_end);
			triple.push_back(Eigen::Triplet<double>(vi, vi, -sum));
		}
		const double w = 1e0;
		Eigen::MatrixX2d b = Eigen::MatrixX2d::Zero(V.cols(), 2);
		for (size_t i = 0; i < bc.size(); ++i) {
			triple.push_back(Eigen::Triplet<double>(bc[i].first, bc[i].first, w));
			b.row(bc[i].first) = w * bc[i].second;
		}
		Eigen::SparseMatrix<double> L(V.cols(), V.cols());
		L.setFromTriplets(triple.begin(), triple.end());
		Eigen::BiCGSTAB<Eigen::SparseMatrix<double>> cg(L);
		if (cg.info() == Eigen::Success) {
			Vout = cg.solve(b).transpose();
			return 1;
		}
		std::cerr << "solve failure!" << std::endl;
		return 0;
	}

	int lscm_mapping(const Eigen::MatrixXd &V, const Eigen::Matrix3Xi &F,
		const std::vector<std::pair<size_t, Eigen::Vector2d>> &bc,
		Eigen::Matrix2Xd &Vout)
	{
		Eigen::Matrix3Xd E(3, 2);
		E << -1, -1, 1, 0, 0, 1;
		Eigen::Matrix<double, 2, 3> t;
		std::vector<Eigen::Triplet<double>>  triplets;
		triplets.reserve(F.cols() * 12);
		for (int i = 0, col = 0; i < F.cols(); ++i, col += 2) {
			if(V.rows() == 2){
				t << V.col(F(0, i)), V.col(F(1, i)), V.col(F(2, i));
			}else{
				//util_helper::flatten_triangle(V.col(F(0, i)), V.col(F(1, i)), V.col(F(2, i)), t);
			}
			const double area = (t(0, 1) - t(0, 0))* (t(1, 2) - t(1, 0)) - (t(0, 2) - t(0, 0))* (t(1, 1) - t(1,0));
			const Eigen::Matrix3Xd m = sqrt(area) * E * (t * E).inverse();
			for (size_t k = 0; k < 3; ++k) {
				const size_t idx = 2 * F(k, i);
				// a - d = 0
				triplets.push_back(Eigen::Triplet<double>(idx, col, m(k, 0)));
				triplets.push_back(Eigen::Triplet<double>(idx + 1, col, -m(k, 1)));
				// c + b = 0
				triplets.push_back(Eigen::Triplet<double>(idx, col + 1, m(k, 1)));
				triplets.push_back(Eigen::Triplet<double>(idx + 1, col + 1, m(k, 0)));
			}
		}
		Eigen::SparseMatrix<double> Gt(2 * V.cols(), 2 * F.cols());
		Gt.setFromTriplets(triplets.begin(), triplets.end());
		Eigen::SparseMatrix<double> L = Gt * Gt.transpose();
		Eigen::VectorXd b;
		b.setZero(2 * V.cols());
		const double w = 1e8;
		for (size_t i = 0; i < bc.size(); ++i) {
			const size_t idx = 2 * bc[i].first;
			for (size_t j = 0; j < 2; ++j) {
				L.coeffRef(idx + j, idx + j) += w;
				b[idx + j] += bc[i].second[j] * w;
			}
		}

		Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> llt(L);
		if (llt.info() != Eigen::Success) {
			std::cout << "factorization failure!" << std::endl;
			return __LINE__;
		}

		const Eigen::VectorXd &x = llt.solve(b);
		Vout = Eigen::Map<const Eigen::Matrix2Xd>(x.data(), 2, x.size() / 2);
		return 0;
	}

	int arap_mapping(const Eigen::MatrixXd &V, const Eigen::Matrix3Xi &F,
		const std::vector<std::pair<size_t, Eigen::Vector2d>> &bc,
		Eigen::Matrix2Xd &Vout, size_t max_iter)
	{
		Eigen::MatrixX2d x(V.cols(), 2);
		if (Vout.cols() != V.cols()) {
			x.setZero();
		}
		else {
			x = Vout.transpose();
		}
		//initilize
		Eigen::VectorXd  sqrt_areas, areas(F.cols());
		std::vector<Eigen::Matrix2d> R(F.cols()), inv;
		Eigen::SparseMatrix<double>  Gt;
		init_para_info(V, F, sqrt_areas, inv, Gt);
		for (size_t i = 0; i < sqrt_areas.size(); ++i) {
			areas[i] = sqrt_areas[i] * sqrt_areas[i];
		}
		Eigen::SparseMatrix<double> L = Gt * Gt.transpose();
		const double w = 1e8;
		for (size_t i = 0; i < bc.size(); ++i) {
			L.coeffRef(bc[i].first, bc[i].first) += w;
		}
		Eigen::SimplicialLLT<Eigen::SparseMatrix<double>>  llt(L);

		size_t iter = 0;
		double err0 = -1, err1 = 0;
		Eigen::MatrixXd J(2, 2);
		Eigen::Matrix2d m;
		Eigen::MatrixX2d b(2 * F.cols(), 2);
		while (fabs(err1 - err0) > 1e-6 && iter < max_iter) {
			err0 = err1;
			//local_step(nodes);
			for (size_t i = 0; i < F.cols(); ++i) {
				for (size_t j = 0; j < 2; ++j) {
					m.col(j) = x.row(F(j + 1, i)) - x.row(F(0, i));
				}
				//util_helper::optimal_rotation(m * inv[i], J);
				R[i] = J;
			}
			//global_step(nodes);
			for (size_t i = 0; i < F.cols(); ++i) {
				for (size_t j = 0; j < 2; ++j) {
					b.row(2 * i + j) = sqrt_areas[i] * R[i].col(j);
				}
			}
			Eigen::MatrixX2d Gtb = Gt * b;
			for (size_t i = 0; i < bc.size(); ++i) {
				Gtb.row(bc[i].first) += w * bc[i].second;
			}
			x = llt.solve(Gtb); //GtGx = Gtb
			//err1 = energy_value(nodes);
			err1 = 0;
			for (size_t i = 0; i < F.cols(); ++i) {
				for (size_t j = 0; j < 2; ++j) {
					m.col(j) = x.row(F(j + 1, i)) - x.row(F(0, i));
				}
				err1 += areas[i] * (m * inv[i] - R[i]).squaredNorm();
			}
			std::cout << ++iter << ":" << err1 << std::endl;
		}
		Vout = x.transpose();
		return 0;
	}
}