#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <surface_mesh/Surface_mesh.h>
#include "tri_mesh.h"
#include "mesh_io.h"
#include "binary_io.h"
#include "vtk.h"
#include "laplace.h"
using namespace std;
using namespace surface_mesh;


const int VERTS = 1415;
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

void work(Eigen::SparseMatrix<double> &L,
	const Eigen::Matrix3Xd &V,
	const std::vector<int> &move_idx,
	const std::vector<Point> &move_coord,
	Eigen::Matrix3Xd &newV) 
{




	Eigen::SparseMatrix<double> v = V.transpose().sparseView();

	Eigen::MatrixXd b(VERTS, 3);
	b.setZero();

	double w = 1e6;
	
	for (int i = 0; i < move_idx.size(); ++i) {
		int id = move_idx[i];
		L.coeffRef(id, id) += w;
		for (int j = 0; j < 3; ++j) 
			b(id, j) = w*move_coord[i][j];
	}

	Eigen::SimplicialCholesky<Eigen::SparseMatrix<double> > solver(L);
	newV = solver.solve(b).transpose();

}


void writeVTK(std::vector<Point> &path) {
	std::vector<double> nodes;
	std::vector<int> lines;
	std::ofstream os("../Mapping/data/line.vtk");

	int n = path.size();
	for (int i = 0; i < n; ++i) {
		nodes.push_back(path[i][0]);
		nodes.push_back(path[i][1]);
		nodes.push_back(path[i][2]);

		if (i) {
			lines.push_back(i - 1);
			lines.push_back(i);
		}
	}

	line2vtk(os, nodes.data(), nodes.size() / 3, lines.data(), lines.size() / 2);
	os.close();
	return;
}


void gao() {
	Eigen::Matrix3Xd V;
	Eigen::Matrix3Xi F;
	common::read_obj("../Mapping/data/nefertiti.obj", V, F);
	cout << V.rows() << " " << V.cols() << endl;
	cout << F.rows() << " " << F.cols() << endl;



	for (int i = 0; i < V.cols(); ++i) {
		V(2, i) = 0;
	}
	cout << V.rows() << " " << V.cols() << endl;
	common::save_obj("../Mapping/data/newface.obj", V, F);
	double tot = 0;
	Surface_mesh mesh;
	build_mesh(V, F, mesh);
	for (const auto &he : mesh.halfedges()) {
		if (mesh.is_boundary(he)) {
			cout << he.idx() << " ";
			cout << mesh.from_vertex(he).idx() << endl;
			tot += mesh.edge_length(mesh.edge(he));
		}
	}
	//cout << "pi = " << pi << endl;
	//cout << "tot = " << tot << endl;;
	std::vector<int> move_idx;
	std::vector<Point> move_coord;
	const auto he0 = Surface_mesh::Halfedge(4425);
	auto &v = mesh.from_vertex(he0);
	auto &he = mesh.halfedge(v);

	Point p = mesh.position(v);

	move_idx.push_back(v.idx());
	move_coord.push_back(p);

	std::vector<std::pair<size_t, Eigen::Vector3d>> bc;

	while (true) {
		double len = mesh.edge_length(mesh.edge(he));
		v = mesh.to_vertex(he);
		he = mesh.halfedge(v);

		double theta = -len / tot * pi * 2;
		//cout << len << endl;

		if (he.idx() == 4425) break;
		p = Point(p[0] * cos(theta) - p[1] * sin(theta), p[0] * sin(theta) + p[1] * cos(theta), 0);

		move_idx.push_back(v.idx());
		move_coord.push_back(p);
	
		Eigen::Vector3d c;
		c[0] = p[0];
		c[1] = p[1];
		c[2] = p[2];
		bc.push_back({ v.idx(), c });

	}

	Eigen::SparseMatrix<double> cotan(VERTS, VERTS);
	//jy_mesh::cal_cot_laplace(V, F, cotan);
	laplacianCotanWeight(mesh, cotan);
	puts("ok");

	for (int i = 0; i < move_idx.size(); ++i) {
		cout << move_idx[i] << " : " << norm(move_coord[i]) << endl;;
	}

	writeVTK(move_coord);

	Eigen::Matrix3Xd newV;

	//harmonic_mapping(V, F, bc, newV );
	work(cotan, V, move_idx, move_coord, newV);

	common::save_obj("../Mapping/data/roundFace.obj", newV, F);
	puts("saved");
}

int main() {
	gao();
}