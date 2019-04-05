#include <Eigen/Dense>
#include <iostream>

int main() {

	Eigen::Matrix<double, 6, 4> m = Eigen::Matrix<double, 6, 4>::Random();

	
	std::cout << "m = \n" << m << "\n";

	std::cout << "m = \n" << m.leftCols(1) << "\n";
}