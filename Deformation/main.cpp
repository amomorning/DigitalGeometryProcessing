#include <Eigen/Dense>
#include "tri_mesh.h"
#include "mesh_io.h"
#include "binary_io.h"
using namespace std;

int main() {
	Eigen::Matrix3Xd V;
	Eigen::Matrix3Xi F;
	common::read_obj("../data/AVE.obj", V, F);

	cout << V.rows() << " " << V.cols() << endl;
	cout << F.rows() << " " << F.cols() << endl;
}