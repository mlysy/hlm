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
			 Eigen::VectorXd gamma0) {
  int n = Z.rows();
  int p = Z.cols();
  hlm::LVLMFit lvlm(n, p);
  Eigen::VectorXd gamma(p);
  lvlm.fit(gamma, y2, Z, gamma0);
  return gamma;
}
