#include "VMD_2D.h"
#include<iostream>
using namespace Eigen;
using namespace std;

void VMD_2D
(vector<MatrixXd>& u, vector<MatrixXcd>& u_hat, vector<vector<vector<double>>>& omega,
	MatrixXd& signal2D, const double alpha, const double tau,
	const int K, const int DC, const int init, const double tol, const double eps) {
	/* ---------------------

	Output:
	-------
	u - the collection of decomposed modes (3D double Matrix in Eigen -vector<MatrixXd>)
	u_hat - spectra of the modes (3D complex<double> Matrix in Eigen -vector<MatrixXcd>)
	omega - estimated mode center - frequencies (2D double Matrix in Eigen -vector<vector<vector<double>>>)
	-------
	Input:
	-------
	signal - the time domain signal2D(2D matrix) to be decomposed
	alpha - the balancing parameter of the data - fidelity constraint
	tau - time - step of the dual ascent(pick 0 for noise - slack)
	K - the number of modes to be recovered
	DC - true if the first mode is putand kept at DC(0 - freq)
	init - 0 = all omegas start at 0
			1 = all omegas start initialized randomly
	tol - tolerance of convergence criterion; typically around 1e-7
	*/

	// ----------Preparations
	// Resolution of image
	int Hy = static_cast<int>(signal2D.rows());
	int Hx = static_cast<int>(signal2D.cols());

	MatrixXd X = MatrixXd::Constant(Hy, Hx, 0);
	X.colwise() = VectorXd::LinSpaced(Hx, 1, Hx) / Hx;

	MatrixXd Y = MatrixXd::Constant(Hy, Hx, 0);
	Y.rowwise() = VectorXd::LinSpaced(Hy, 1, Hy).transpose() / Hy;

	// Spectral Domain discretization
	double fx = 1.0 / Hx;
	double fy = 1.0 / Hy;
	MatrixXd freqs_1 = X.array() - 0.5 - fx;
	MatrixXd freqs_2 = Y.array() - 0.5 - fy;

	// N is the maximum number of iterations
	int N = 3000;

	// For future generalizations: alpha might be individual for each mode
	VectorXd Alpha = alpha * VectorXd::Ones(K);

	// Construct f and f_hat
	MatrixXcd f_hat_pre = fft2D(signal2D);
	MatrixXcd f_hat = fftshift(f_hat_pre);

	// Storage matrices for (Fourier) modes. All iterations are not recorded.
	u_hat.resize(K, MatrixXcd::Zero(Hy, Hx));
	vector<MatrixXcd> u_hat_old = u_hat;
	MatrixXcd sum_uk = MatrixXcd::Zero(Hy, Hx);

	// Storage matrices for (Fourier) Lagrange multiplier.
	MatrixXcd mu_hat = MatrixXcd::Zero(Hy, Hx);

	// N iterations at most, 2 spatial coordinates, K clusters
	omega.resize(N, vector<vector<double>>(2, vector<double>(K, 0.0)));

	int maxK = 0;
	// Initialization of omega_k
	switch (init){
		case 0:
			// spread omegas radially
			// if DC, keep first mode at 0,0
			maxK = DC ? K - 1 : K;
			for (int k = DC; k < maxK; k++)	{
				omega[0][0][k] = 0.25 * cos(myPI * k / maxK);
				omega[0][1][k] = 0.25 * sin(myPI * k / maxK);
			}
			break;

		// Case 1: random on half-plane
		case 1:
			for (int k = 0; k < K; k++)	{
				omega[0][0][k] = ((double)rand() / RAND_MAX) - 0.5;
				omega[0][1][k] = ((double)rand() / RAND_MAX) / 2.0;
			}

			// DC component (if expected)
			if (DC)
				omega[0][0][0] =  omega[0][1][0] = 0;

			break;
	}


	// Main loop for iterative updates
	// Stopping criteria tolerances
	double uDiff = tol + eps;
	double omegaDiff = tol + eps;
	MatrixXcd u_hat_mid_pre;
	// first run
	int n = 0;

	// run until convergence or max number of iterations
	while ((uDiff > tol || omegaDiff > tol) && n < N){
		// first things first
		int k = 0;

		// compute the halfplane mask for the 2D "analytic signal"
		MatrixXd HilbertMask = sign_plus_one(freqs_1 * omega[n][0][k] + freqs_2 * omega[n][1][k]);

		// update first mode accumulator
		sum_uk = u_hat[K-1] + sum_uk - u_hat[k];

		// update first mode's spectrum through wiener filter (on half plane)
		u_hat[k] = ((f_hat.array() - sum_uk.array() - mu_hat.array() / 2.0) * HilbertMask.array()) /
			(1 + Alpha[k] * ((freqs_1.array() - omega[n][0][k])).array().square() +
				(freqs_2.array() - omega[n][1][k]).array().square());

		// update first mode's central frequency as spectral center of gravity
		if (!DC){
			// In Eigen, .sum() computes the sum of all coefficients, .square() computes the square of each coefficient, and .abs() computes the absolute value of each coefficient
			double denominator = (u_hat[k].array().abs().square()).sum();
			omega[n + 1][0][k] = (freqs_1.array() * u_hat[k].array().abs().square()).sum() / denominator;
			omega[n + 1][1][k] = (freqs_2.array() * u_hat[k].array().abs().square()).sum() / denominator;

			// keep omegas on same halfplane
			if ( omega[n + 1][1][k] < 0){
				for (int i = 0; i < omega[0].size(); i++)
					omega[n + 1][i][k] = - omega[n + 1][i][k];
			}
		}

		// recover full spectrum from analytic signal
		u_hat_mid_pre = fft2D(fftshift(u_hat[k]), 1).real();
		u_hat[k] = fftshift(fft2D(u_hat_mid_pre));

		// work on other modes
		for (k = 1; k < K; k++) {
			// recompute Hilbert mask
			HilbertMask = sign_plus_one(freqs_1 * omega[n][0][k] + freqs_2 * omega[n][1][k]) + MatrixXd::Ones(freqs_1.rows(), freqs_1.cols());

			// update accumulator
			sum_uk = u_hat[k - 1] + sum_uk - u_hat[k];

			// update signal spectrum
			u_hat[k] = ((f_hat - sum_uk - mu_hat / 2.0).cwiseProduct(HilbertMask)).cwiseQuotient(
				(1 + Alpha[k] * ((freqs_1.array() - omega[n][0][k])).array().square() +
					(freqs_2.array() - omega[n][1][k]).array().square()).matrix());


			// update signal frequencies
			double denominator = (u_hat[k].cwiseAbs()).array().square().sum();
			omega[n + 1][0][k] = (freqs_1.cwiseProduct((u_hat[k].cwiseAbs()).array().square().matrix())).sum() / denominator;
			omega[n + 1][1][k] = (freqs_2.cwiseProduct((u_hat[k].cwiseAbs()).array().square().matrix())).sum() / denominator;

			// keep omegas on same halfplane
			if (omega[n + 1][1][k] < 0) {
				for (int i = 0; i < omega[0].size(); i++)
					omega[n + 1][i][k] = -omega[n + 1][i][k];
			}

			// recover full spectrum from analytic signal
			u_hat_mid_pre = fft2D(fftshift(u_hat[k]),1).real();
			u_hat[k] = fftshift(fft2D(u_hat_mid_pre));
		}


		// Gradient ascent for augmented Lagrangian
		mu_hat = mu_hat + tau * (sum_u_hat_3d(u_hat) - f_hat);

		// Increment iteration counter
		n++;

		// Convergence?
		
		omegaDiff = 0;
		for (int i = 0; i < 2; i++) {
			for (int k = 0; k < K; k++) {
				double diff = omega[n][i][k] - omega[n - 1][i][k];
				omegaDiff += diff * diff;
			}
		}
		omegaDiff = omegaDiff * K + eps;


		uDiff = eps;
		/*
		for (int k = 0; k < K; k++){
			MatrixXcd diff = u_hat[k] - u_hat_old[k];
			uDiff += (diff.array() * diff.conjugate().array()).sum().real() / (Hx * Hy);
		}*/
		for (int k = 0; k < K; ++k) {
			MatrixXcd pre_diff = u_hat[k] - u_hat_old[k];
			MatrixXcd diff = pre_diff.array() * pre_diff.array().conjugate();
			double sumDiff = diff.real().sum();
			uDiff += sumDiff / (Hx * Hy);
		}

		uDiff = std::abs(uDiff);

		u_hat_old = u_hat;
		std::cout << n << " time; uDiff: " <<uDiff<< " ; omegaDiff: " << omegaDiff << std::endl;
	}

	// Signal Reconstruction
	// Inverse Fourier Transform to compute (spatial) modes
	u.resize(K, MatrixXd::Zero(Hy, Hx));
	for (int k = 0; k < K; k++) {
		u[k] = fft2D(fftshift(u_hat[k]), 1).real();
	}

	// Should the omega-history be returned, or just the final results?
	//vector<vector<vector<double>>> finalOmega = {omega[n]};

	return;
}

#pragma region Ancillary Functions 
MatrixXcd fftshift(const MatrixXcd& mat) {
	MatrixXcd shifted(mat.rows(), mat.cols());
	int halfRows = static_cast<int>(mat.rows() / 2);
	int halfCols = static_cast<int>(mat.cols() / 2);

	for (int i = 0; i < mat.rows(); i++) {
		for (int j = 0; j < mat.cols(); j++) {
			int newRow = (i + halfRows) % mat.rows();
			int newCol = (j + halfCols) % mat.cols();
			shifted(newRow, newCol) = mat(i, j);
		}
	}

	return shifted;
}

MatrixXd sign_plus_one(MatrixXd val) {
	MatrixXd result = val;
	for (int i = 0; i < val.rows(); i++) {
		for (int j = 0; j < val.cols(); j++) {
			result(i, j) = int((0 < val(i, j)) - (val(i, j) < 0) + 1);
		}
	}
	return result;
}

MatrixXcd sum_u_hat_3d(const std::vector<MatrixXcd>& u_hat) {
	Eigen::Index rows = u_hat[0].rows();
	Eigen::Index cols = u_hat[0].cols();

	MatrixXcd sum_matrix = MatrixXcd::Zero(rows, cols); // initialize with zeros

	for (const auto& mat : u_hat) {
		sum_matrix += mat; // element-wise addition
	}

	return sum_matrix;
}



MatrixXcd fft2D(MatrixXcd mat, const int option) {
	//option: 0-fft.fwd; 1-fft.inv
	Eigen::Index rows = mat.rows();
	Eigen::Index cols = mat.cols();
	MatrixXcd output(rows, cols);
	Eigen::FFT<double> fft;
	for (int i = 0; i < rows; ++i) {
		VectorXcd rowResult;
		if (option)
			fft.inv(rowResult, mat.row(i));
		else
			fft.fwd(rowResult, mat.row(i));
		output.row(i) = rowResult;
	}
	for (int i = 0; i < cols; ++i) {
		VectorXcd colResult;
		if (option)
			fft.inv(colResult, output.col(i));
		else
			fft.fwd(colResult, output.col(i));
		output.col(i) = colResult;
	}
	return output;
}

MatrixXcd fft2D(MatrixXd mat, const int option) {
	MatrixXcd newMat = mat;
	return fft2D(newMat, option);
}
#pragma endregion

