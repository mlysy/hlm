/// @file HLMExports.cpp
///
/// @brief `Rcpp` wrappers for selected functions in the hlm namespace.
///
/// The exported functions are:
///
/// - `LM_Fit`: possibly weighted linear regression.
/// - `LVLM_Fit`: both algorithms for log-variance linear regression
/// - `HLM_Fit`: fit the heteroscedastic linear model 

// #undef NDEBUG
// #define EIGEN_RUNTIME_NO_MALLOC
#include <Rcpp.h>
#include "hlm/LMFit.h"
#include "hlm/LVLMFit.h"
#include "hlm/HLMFit.h"
using namespace Rcpp;

/// Fit a linear model.
///
/// @param[in] y Response vector of length `n`.
/// @param[in] X Covariate matrix of size `n x p`.
///
/// @return The MLE coefficient vector of length `p`.
// [[Rcpp::export("LM_Fit")]]
Eigen::VectorXd lm_fit(Eigen::VectorXd y, Eigen::MatrixXd X) {
  int n = X.rows();
  int p = X.cols();
  hlm::LMFit lm(n, p);
  Eigen::VectorXd beta(p);
  lm.fit(beta, y, X);
  return beta;
}

/// Fit a weighted linear model.
///
/// @param[in] y Response vector of length `n`.
/// @param[in] X Covariate matrix of size `n x p`.
/// @param[in] w Weight vector of length `n`.
///
/// @return The MLE coefficient vector of length `p`.
// [[Rcpp::export("WLM_Fit")]]
Eigen::VectorXd wlm_fit(Eigen::VectorXd y, Eigen::MatrixXd X,
			Eigen::VectorXd w) {
  int n = X.rows();
  int p = X.cols();
  hlm::LMFit lm(n, p);
  Eigen::VectorXd beta(p);
  lm.fit(beta, y, X, w);
  return beta;
}

/// Fit a log-variance linear model using Fisher Scoring algorithm.
///
/// For the definition of an LVLM model, see file `LVLMFit.h`.
///
/// @param[in] y2 Squared response vector of length `n`.
/// @param[in] Z Covariate matrix of size `n x q`.
/// @param[in] gamma0 Initial coefficient vector of length `q`, or scalar.
/// @param[in] initLS If `true` use `gamma0` as initial value.  Otherwise, use least squares estimate calculated by LVLMFit::fitLS.
/// @param[in] maxit Maximum number of algorithm steps.
/// @param[in] epsilon Relative tolerance on the maximum loglikelihood value.
///
/// @return An `Rcpp::List` with elements:
///
/// - `gamma`: The MLE of the coefficient vector of length `q`.
/// - `loglik`: The value of the loglikelihood at the final step of the fitting algorithm.
/// - `iter`: The actual number of algorithm steps.
/// - `error`: The relative error on the loglikelihood at the terminal step.
///
// [[Rcpp::export("LVLM_FitFS")]]
Rcpp::List lvlm_fitFS(Eigen::VectorXd y2, Eigen::MatrixXd Z,
		      Eigen::VectorXd gamma0,
		      int maxit = 25, double epsilon = 1e-8,
		      bool initLS = false) {
  int n = Z.rows();
  int q = Z.cols();
  double llik;
  int niter;
  double error;
  Eigen::VectorXd gamma(q);
  hlm::LVLMFit lvlm(n, q);
  lvlm.fitControl(maxit, epsilon);
  lvlm.set_ZtZ(Z);
  if(initLS) {
    lvlm.fitFS(gamma, y2, Z, false);
  } else {
    lvlm.fitFS(gamma, y2, Z, gamma0, false);
  }
  lvlm.fitStats(llik, niter, error);
  return List::create(_["gamma"] = gamma,
		      _["loglik"] = llik,
		      _["iter"] = niter,
		      _["error"] = error);
}

/// Fit a log-variance linear model using Iteratively Reweighted Least Squares algorithm.
///
/// For the definition of an LVLM model, see file `LVLMFit.h`.
///
/// @param[in] y2 Squared response vector of length `n`.
/// @param[in] Z Covariate matrix of size `n x q`.
/// @param[in] gamma0 Initial coefficient vector of length `q`, or scalar.
/// @param[in] initLS If `true` use `gamma0` as initial value.  Otherwise, use least squares estimate calculated by LVLMFit::fitLS.
/// @param[in] maxit Maximum number of algorithm steps.
/// @param[in] epsilon Relative tolerance on the maximum loglikelihood value.
///
/// @return An `Rcpp::List` with elements:
///
/// - `gamma`: The MLE of the coefficient vector of length `q`.
/// - `loglik`: The value of the loglikelihood at the final step of the fitting algorithm.
/// - `iter`: The actual number of algorithm steps.
/// - `error`: The relative error on the loglikelihood at the terminal step.
///
// [[Rcpp::export("LVLM_FitIRLS")]]
Rcpp::List lvlm_fitIRLS(Eigen::VectorXd y2, Eigen::MatrixXd Z,
			Eigen::VectorXd gamma0,
			int maxit = 25, double epsilon = 1e-8,
			bool initLS = false) {
  int n = Z.rows();
  int q = Z.cols();
  double llik;
  int niter;
  double error;
  Eigen::VectorXd gamma(q);
  hlm::LVLMFit lvlm(n, q);
  lvlm.fitControl(maxit, epsilon);
  if(initLS) {
    lvlm.fitIRLS(gamma, y2, Z);
  } else {
    lvlm.fitIRLS(gamma, y2, Z, gamma0);
  }
  lvlm.fitStats(llik, niter, error);
  return List::create(_["gamma"] = gamma,
		      _["loglik"] = llik,
		      _["iter"] = niter,
		      _["error"] = error);
}

/// Fit a log-variance linear model using the Least Squares algorithm.
///
/// For the definition of an LVLM model, see file `LVLMFit.h`.
///
/// @param[in] logY2 Log of squared response vector of length `n`.
/// @param[in] Z Covariate matrix of size `n x q`.
///
/// @return The least-squares estimated coefficient vector of length `q`.
// [[Rcpp::export("LVLM_FitLS")]]
Eigen::VectorXd lvlm_fitLS(Eigen::VectorXd logY2, Eigen::MatrixXd Z) {
  int n = Z.rows();
  int q = Z.cols();
  Eigen::VectorXd gamma(q);
  hlm::LVLMFit lvlm(n, q);
  lvlm.set_ZtZ(Z);
  lvlm.fitLS(gamma, logY2, Z, false);
  return gamma;
}

/// Fit a heteroscedastic linear model.
///
/// For the definition of the HLM model see file `HLMFit.h`.
///
/// @param[in] y Response vector of length `n`.
/// @param[in] X Mean covariate matrix of size `n x p`.
/// @param[in] Z Variance covariate matrix of size `n x q`.
/// @param[in] beta0 Initial mean coefficient vector of length `p`.
/// @param[in] gamma0 Initial variance coefficient vector of length `q`.
/// @param[in] maxit Maximum number of block ascent cycles.
/// @param[in] epsilon Relative tolerance on the maximum loglikelihood value.
/// @param[in] method Integer describing the conditional LVLM optimization algorithm to use.  0: Fisher Scoring.  1: IRLS.  In both cases the control parameters are the default values in LVLMFit::fitControl.
///
/// @return An `Rcpp::List` with the following elements:
///
/// - `beta`, `gamma`: The MLE mean and variance coefficient vectors of length `p` and `q`.
/// - `loglik`: The value of the loglikelihood at the final step of the fitting algorithm.
/// - `iter`: The actual number of algorithm steps.
/// - `error`: The relative error on the loglikelihood at the terminal step.
///
// [[Rcpp::export("HLM_Fit")]]
Rcpp::List hlm_fit(Eigen::VectorXd y,
		   Eigen::MatrixXd X, Eigen::MatrixXd Z,
		   Eigen::VectorXd beta0, Eigen::VectorXd gamma0,
		   int maxit = 100, double epsilon = 1e-5,
		   int method = 0) {
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
		      _["iter"] = niter,
		      _["error"] = error);
}
