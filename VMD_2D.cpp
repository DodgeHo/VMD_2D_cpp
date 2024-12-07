#define cimg_use_bmp
#include "VMD_2D.h"
#include "CImg\CImg.h"
#include <fstream>
#include <iostream>
using namespace Eigen;
using namespace std;
/*
<VMD_CPP: C++ implementation of Variational Mode Decomposition using Eigen.>
Copyright (C) <2019>  <Lang He: asdsay@gmail.com>
Mozilla Public License v. 2.0.
*/


Eigen::MatrixXd load_image(const std::string& filename) {
	cimg_library::CImg<unsigned char> image(filename.c_str());
	Eigen::MatrixXd matrix(image.height(), image.width());
	for (int i = 0; i < image.height(); i++)
		for (int j = 0; j < image.width(); j++)
			matrix(i, j) = image(j, i);  // Differ between CImg'scoordinate system and Eigen's
	return matrix;
}
void save_image(const std::string& filename, const Eigen::MatrixXd& matrix) {
	cimg_library::CImg<double> image(int(matrix.cols()), int(matrix.rows()));
	for (int i = 0; i < image.height(); i++)
		for (int j = 0; j < image.width(); j++)
			image(j, i) = matrix(i, j);  // ditto
	image.save(filename.c_str());
}


int main() {
	// create a signal2D to simulation the procedure.
	MatrixXd signal2D = load_image("Sample.bmp");

	// initial some input parameters
	const int K = 5, DC = 1, init = 1;
	const double alpha = 5000.0, tau = 0.25, tol = K * 1e-6, eps = 2.2204e-16;
	vector<MatrixXd> u;
	vector<vector<vector<double>>> omega;
	vector<MatrixXcd> u_hat;

	//VMD 2D Process
	VMD_2D(u, u_hat, omega, signal2D, alpha, tau, K, DC, init, tol, eps);

	//Example 2: If you only wants to get sum result of the first n mode of signals.
	for (int k = 0; k < K; k++) {
		MatrixXd eachSignal = u[k];
		string filename = "DecomResult_" + std::to_string(k) + ".bmp";
		save_image(filename, eachSignal);
	}
	
	return 0;
}; 