#ifndef LAPLACE_H_H
#define LAPLACE_H_H

#include <Eigen/Dense>
#include <Eigen/Sparse>

namespace jy_mesh
{
int cal_cot_laplace(const Eigen::MatrixXd &V, const Eigen::Matrix3Xi &F,
					Eigen::SparseMatrix<double> &L);

int cal_cot_laplace(const Eigen::MatrixXd &V, const Eigen::Matrix3Xi &F,
					const Eigen::VectorXi &g2l, Eigen::SparseMatrix<double> &L);

int cal_topo_laplace(const Eigen::MatrixXd &V, const Eigen::Matrix3Xi &F,
					 Eigen::SparseMatrix<double> &L);

int cal_topo_laplace(const Eigen::MatrixXd &V, const Eigen::Matrix3Xi &F,
					 const Eigen::VectorXi &g2l, Eigen::SparseMatrix<double> &L);
}

#endif
