#include <Rcpp.h>
#include "hlm/LMFit.h"

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
