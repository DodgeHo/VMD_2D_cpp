#pragma once
#include <vector>
#include <cmath>
#include <ctime>
#include "Eigen/Eigen/Eigen"
#include "Eigen/unsupported/Eigen/FFT"
//#include <Eigen/Core>
//#include <unsupported/Eigen/FFT>

#define myPI acos(-1)
using namespace Eigen;
using std::vector;

void VMD_2D
(vector<MatrixXd>& u, vector<MatrixXcd>& u_hat, vector<vector<vector<double>>>& omega,
	MatrixXd& signal, const double alpha, const double tau,
	const int K, const int DC, const int init, const double tol, const double eps);

MatrixXcd fftshift(const MatrixXcd& mat);
MatrixXd sign_plus_one(MatrixXd val);
MatrixXcd sum_u_hat_3d(const std::vector<MatrixXcd>& u_hat);
MatrixXcd fft2D(MatrixXd mat, const int option = 0);
MatrixXcd fft2D(MatrixXcd mat, const int option = 0);