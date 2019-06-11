/// @file TestExports.cpp
///
/// @brief C++ test functions.
///
/// These functions are exported by `Rcpp` but not the R package.  Their purpose is to test the C++ code.

// #undef NDEBUG
#define EIGEN_RUNTIME_NO_MALLOC
#include <Rcpp.h>
#include "hlm/LMFit.h"
#include "hlm/LVLMFit.h"
#include "hlm/HLMFit.h"
using namespace Rcpp;

// [[Rcpp::export("wlm_fit")]]
Eigen::VectorXd wlm_fit(Eigen::VectorXd y, Eigen::MatrixXd X,
			Eigen::VectorXd w) {
  int n = X.rows();
  int p = X.cols();
  hlm::LMFit lm(n, p);
  Eigen::VectorXd beta(p);
  Eigen::internal::set_is_malloc_allowed(false);
  lm.fit(beta, y, X, w);
  Eigen::internal::set_is_malloc_allowed(true);
  return beta;
}

// [[Rcpp::export("lm_fit")]]
Eigen::VectorXd lm_fit(Eigen::VectorXd y, Eigen::MatrixXd X) {
  int n = X.rows();
  int p = X.cols();
  hlm::LMFit lm(n, p);
  Eigen::VectorXd beta(p);
  Eigen::internal::set_is_malloc_allowed(false);
  lm.fit(beta, y, X);
  Eigen::internal::set_is_malloc_allowed(true);
  return beta;
}

// [[Rcpp::export("lvlm_fitFS")]]
Eigen::VectorXd lvlm_fitFS(Eigen::VectorXd y2, Eigen::MatrixXd Z,
			   Eigen::VectorXd gamma0,
			   int maxit = 25, double epsilon = 1e-8,
			   bool initLS = false) {
  int n = Z.rows();
  int p = Z.cols();
  int niter;
  double error;
  double llik;
  Eigen::VectorXd gamma(p);
  gamma = gamma0;
  hlm::LVLMFit lvlm(n, p);
  lvlm.fitControl(maxit, epsilon);
  lvlm.set_ZtZ(Z);
  Eigen::internal::set_is_malloc_allowed(false);
  if(initLS) {
    lvlm.fitFS(gamma, y2, Z, false);
  } else {
    lvlm.fitFS(gamma, y2, Z, gamma, false);
  }
  Eigen::internal::set_is_malloc_allowed(true);	
  return gamma;
}

// [[Rcpp::export("lvlm_fitIRLS")]]
Eigen::VectorXd lvlm_fitIRLS(Eigen::VectorXd y2, Eigen::MatrixXd Z,
			     Eigen::VectorXd gamma0,
			     int maxit = 25, double epsilon = 1e-8,
			     bool initLS = false) {
  int n = Z.rows();
  int p = Z.cols();
  int niter;
  double error;
  double llik;
  Eigen::VectorXd gamma(p);
  hlm::LVLMFit lvlm(n, p);
  lvlm.fitControl(maxit, epsilon);
  lvlm.set_ZtZ(Z);
  Eigen::internal::set_is_malloc_allowed(false);
  if(initLS) {
    lvlm.fitIRLS(gamma, y2, Z);
  } else {
    lvlm.fitIRLS(gamma, y2, Z, gamma0);
  }
  Eigen::internal::set_is_malloc_allowed(true);	
  return gamma;
}

// [[Rcpp::export("lvlm_fitLS")]]
Eigen::VectorXd lvlm_fitLS(Eigen::VectorXd logY2, Eigen::MatrixXd Z) {
  int n = Z.rows();
  int p = Z.cols();
  Eigen::VectorXd gamma(p);
  hlm::LVLMFit lvlm(n, p);
  lvlm.set_ZtZ(Z);
  Eigen::internal::set_is_malloc_allowed(false);
  lvlm.fitLS(gamma, logY2, Z, false);
  Eigen::internal::set_is_malloc_allowed(true);    
  return gamma;
}

// [[Rcpp::export("hlm_fit")]]
Rcpp::List hlm_fit(Eigen::VectorXd y,
		   Eigen::MatrixXd X, Eigen::MatrixXd Z,
		   Eigen::VectorXd beta0, Eigen::VectorXd gamma0,
		   int maxit = 1000, double epsilon = 1e-5,
		   int method = 1) {
  int n = X.rows();
  int p = X.cols();
  int q = Z.cols();
  double llik;
  int niter;
  double error;
  Eigen::VectorXd beta(p);
  Eigen::VectorXd gamma(q);
  hlm::HLMFit hlm(n, p, q);
  hlm.fitControl(maxit, epsilon);
  hlm.fit(beta, gamma, y, X, Z, beta0, gamma0, method);
  hlm.fitStats(llik, niter, error);
  return List::create(_["beta"] = beta,
		      _["gamma"] = gamma,
		      _["loglik"] = llik,
		      _["niter"] = niter,
		      _["error"] = error);
}
