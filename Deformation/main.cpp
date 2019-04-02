#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <surface_mesh/Surface_mesh.h>
#include "tri_mesh.h"
#include "mesh_io.h"
#include "binary_io.h"
using namespace std;
using namespace surface_mesh;


const int VERTS = 12500;

void laplacianCotanWeight(const Surface_mesh &mesh,
	Eigen::SparseMatrix<double> &cotan)
{
	Surface_mesh::Face_iterator fit;
	auto points = mesh.get_vertex_property<Point>("v:point");

	fit = mesh.faces_begin();
	std::vector<Eigen::Triplet<double> > tri;
	do {
		Surface_mesh::Vertex_around_face_circulator vf = mesh.vertices(*fit);
		Point p[3];
		int id[3];
		double cot[3];
		for (int i = 0; i < 3; ++i, ++vf) {
			p[i] = points[*vf];
			id[i] = (*vf).idx();
		}

		double sum = 0;
		for (int i = 0; i < 3; ++i) {
			int j = (i + 1) % 3, k = (j + 1) % 3;
			cot[i] = dot(p[j] - p[i], p[k] - p[i]) /
				norm(cross(p[j] - p[i], p[k] - p[i]));

			tri.push_back({ id[j], id[k], -0.5*cot[i] });
			tri.push_back({ id[k], id[j], -0.5*cot[i] });
		}

		for (int i = 0; i < 3; ++i) {
			tri.push_back({ id[i], id[i], 0.5*(cot[(i + 1) % 3], cot[(i + 2) % 3]) });
		}

	} while (++fit != mesh.faces_end());
	cotan.setFromTriplets(tri.begin(), tri.end());
}

void work(Eigen::SparseMatrix<double> &cotan) {
	std::vector<int> fix_idx;
	for (int i = 0; i < 10; ++i) fix_idx.push_back(i);

	std::vector<int> move_idx;
	move_idx.push_back(5792);


	Eigen::SparseMatrix<double> A = cotan, ATA, ATb;
	A.conservativeResize(VERTS + fix_idx.size() + move_idx.size(), VERTS);

	std::vector<Eigen::Triplet<double> > tri;
	for (int i = 0; i < fix_idx.size(); ++i) {
		A.coeffRef(VERTS + i, fix_idx[i]) = 1;
	}
	for (int i = 0; i < move_idx.size(); ++i) {
		A.coeffRef(VERTS + i + fix_idx.size(), move_idx[i]);
	}
	ATA = A.transpose() * A;
}

int main() {
	Eigen::Matrix3Xd V;
	Eigen::Matrix3Xi F;
	common::read_obj("../data/AVE.obj", V, F);
	cout << V.rows() << " " << V.cols() << endl;
	cout << F.rows() << " " << F.cols() << endl;

	Surface_mesh mesh;
	build_mesh(V, F, mesh);

	Eigen::SparseMatrix<double> cotan(VERTS, VERTS);
	laplacianCotanWeight(mesh, cotan);


	work(cotan);

}