#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <surface_mesh/Surface_mesh.h>
#include "tri_mesh.h"
#include "mesh_io.h"
#include "binary_io.h"
#include "vtk.h"
using namespace std;
using namespace surface_mesh;


const int VERTS = 71;
const double pi = acos(-1.0);

void laplacianCotanWeight(const Surface_mesh &mesh,
	Eigen::SparseMatrix<double> &cotan)
{
	std::vector<Eigen::Triplet<double> > tri;
	for (const auto &fit : mesh.faces())
	{
		int i = 0;
		Point p[3];
		int id[3];
		double cot[3];
		for (const auto &vit : mesh.vertices(fit))
		{
			p[i] = mesh.position(vit);
			id[i] = vit.idx();
			++i;
		}

		for (int i = 0; i < 3; ++i) {
			int j = (i + 1) % 3, k = (i + 2) % 3;
			cot[i] = 0.5*dot(p[j] - p[i], p[k] - p[i]) /
				norm(cross(p[j] - p[i], p[k] - p[i]));

			tri.push_back({ id[j], id[k], cot[i] });
			tri.push_back({ id[k], id[j], cot[i] });
		}

		for (int i = 0; i < 3; ++i) {
			tri.push_back({ id[i], id[i], -(cot[(i + 1) % 3] + cot[(i + 2) % 3]) });
		}
	}
	cotan.setFromTriplets(tri.begin(), tri.end());
}

void work(const Eigen::SparseMatrix<double> &L,
	const Eigen::Matrix3Xd &V,
	const std::vector<int> &move_idx,
	const std::vector<Point> &move_coord,
	Eigen::MatrixXd &newV) 
{

	std::vector<int> fix_idx = { 0 };


	Eigen::SparseMatrix<double> A = L;
	A.conservativeResize(VERTS + fix_idx.size() + move_idx.size(), VERTS);

	Eigen::SparseMatrix<double> v = V.transpose().sparseView();
	Eigen::SparseMatrix<double> b = L * v;
	
	b.conservativeResize(VERTS + fix_idx.size() + move_idx.size(), 3);

	for (int i = 0; i < fix_idx.size(); ++i) {
		A.coeffRef(VERTS + i, fix_idx[i]) = 1;

		b.insert(VERTS + i, 0) = V(0, fix_idx[i]);
		b.insert(VERTS + i, 1) = V(1, fix_idx[i]);
		b.insert(VERTS + i, 2) = V(2, fix_idx[i]);
	}
	for (int i = 0; i < move_idx.size(); ++i) {
		A.coeffRef(VERTS + i + fix_idx.size(), move_idx[i]) = 1;

		b.insert(VERTS + i + fix_idx.size(), 0) = move_coord[i][0];
		b.insert(VERTS + i + fix_idx.size(), 1) = move_coord[i][1];
		b.insert(VERTS + i + fix_idx.size(), 2) = move_coord[i][2];
		cout << V.col(move_idx[i]) << endl;
		cout << "move to : \n";
		cout << move_coord[i][0] << " " << move_coord[i][1] << " " << move_coord[i][2] << endl;
	}

	Eigen::SimplicialCholesky<Eigen::SparseMatrix<double> > solver(A.transpose()*A);
	newV.resize(VERTS, 3);
	for (int i = 0; i < 3; ++i) {
		newV.col(i) = solver.solve(A.transpose()*b.col(i)).toDense();
	}
}



int main() {
	Eigen::Matrix3Xd V;
	Eigen::Matrix3Xi F;
	common::read_obj("../data/round.obj", V, F);
	cout << V.rows() << " " << V.cols() << endl;
	cout << F.rows() << " " << F.cols() << endl;

	Surface_mesh mesh;
	build_mesh(V, F, mesh);

	double tot = 0;

	for (const auto &e : mesh.edges()) {
		if (mesh.is_boundary(e)) {
			//cout << he.idx() << endl;
			tot += mesh.edge_length(e);
		}
	}
	cout << "pi = " << pi << endl;
	cout << "tot = " << tot << endl;;
	std::vector<int> move_idx;
	std::vector<Point> move_coord;
	const auto he0 = Surface_mesh::Halfedge(1);
	auto &v = mesh.from_vertex(he0);
	auto &he = mesh.halfedge(v);

	Point p = mesh.position(v);

	while (true) {
		double len = mesh.edge_length(mesh.edge(he));
		v = mesh.to_vertex(he);
		he = mesh.halfedge(v);
		if (he.idx() == 1) break;

		double theta = len / tot * pi*2;
		cout << len << endl;

		p = Point(p[0]*cos(theta) - p[1]*sin(theta), p[0]*sin(theta) + p[1]*cos(theta), 0);

		move_idx.push_back(v.idx());
		move_coord.push_back(p);
		
	}
	
	Eigen::SparseMatrix<double> cotan(VERTS, VERTS);
	laplacianCotanWeight(mesh, cotan);
	
	for (int i = 0; i < move_idx.size(); ++i) {
		cout << move_idx[i] << " : " << norm(move_coord[i]) << endl;;
	}

	
	puts("ok");
	Eigen::MatrixXd newV;


	work(cotan, V, move_idx, move_coord, newV);

	common::save_obj("../data/newRound.obj", newV.transpose(), F);
	puts("saved");

	//work(cotan);

}