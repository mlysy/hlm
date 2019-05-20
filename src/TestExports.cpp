/// @file TestExports.cpp
///
/// @brief C++ test functions.
///
/// These functions are exported by `Rcpp` but not the R package.  Their purpose is to test the C++ code.

#include <Rcpp.h>
#include "hlm/LMFit.h"
#include "hlm/LVLMFit.h"

// [[Rcpp::export("lm_fit")]]
Eigen::VectorXd lm_fit(Eigen::VectorXd y, Eigen::MatrixXd X,
		       Eigen::VectorXd weights) {
  int n = X.rows();
  int p = X.cols();
  hlm::LMFit lm(n, p);
  Eigen::VectorXd beta(p);
  lm.fit(beta, y, X, weights);
  return beta;
}

// [[Rcpp::export("lvlm_fit")]]
Eigen::VectorXd lvlm_fit(Eigen::VectorXd y2, Eigen::MatrixXd Z,
			 Eigen::VectorXd gamma0, int maxit = 25, double epsilon = 1e-8) {
  int n = Z.rows();
  int p = Z.cols();
  int niter;
  double error;
  Eigen::VectorXd gamma(p);
  hlm::LVLMFit lvlm(n, p);
  lvlm.fitControl(maxit, epsilon);
  lvlm.fit(gamma, niter, error, y2, Z, gamma0);
  return gamma;
}

// [[Rcpp::export("lvlm_fitLS")]]
Eigen::VectorXd lvlm_fitLS(Eigen::VectorXd logY2, Eigen::MatrixXd Z) {
  int n = Z.rows();
  int p = Z.cols();
  Eigen::VectorXd gamma(p);
  hlm::LVLMFit lvlm(n, p);
  lvlm.fitLS(gamma, logY2, Z);
  return gamma;
}
