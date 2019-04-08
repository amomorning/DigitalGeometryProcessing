#include <Eigen/Dense>
#include <iostream>
using namespace std;

int main() {
	Eigen::MatrixXd f(1, 9);
	Eigen::Matrix3d ran = Eigen::Matrix3d::Random();

	cout << ran << endl;
	int cnt = 0;
	for (int j = 0; j < 3; ++j) {
		for (int i = 0; i < 3; ++i) {
			f(0, cnt++) = ran(i, j);
		}
	}
	Eigen::MatrixXd t = f.transpose().col(0);
	t.resize(3, 3);
	cout << endl << t << endl;
}