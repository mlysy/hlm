/// @file HLMFit.h

#ifndef HLM_LVLMFIT_H
#define HLM_LVMFIT_H 1

// [[Rcpp::depends("RcppEigen")]]
#include <RcppEigen.h>

namespace hlm {

  using namespace Eigen;
  typedef const Ref<const VectorXd> cRVectorXd;
  typedef const Ref<const MatrixXd> cRMatrixXd;
  typedef Ref<VectorXd> RVectorXd;
  typedef Ref<MatrixXd> RMatrixXd;

  /// Log-Variance Linear Models.
  ///
  /// The log-variance linear model (LVLM) is of the form...
  class LVLMFit {
  private:
    // problem dimensions
    int n_;
    int p_;
    LLT<MatrixXd> ZtZ_;
    VectorXd wy2_;
    VectorXd Zg_;
    // control parameters
    int maxit_;
    int epsilon_;
  public:
    void fitOptions(int maxit = 25, double epsilon = 1e-8);
    LVLMFit();
    LVLMFit(int n, int p);
    void fit(RVectorXd gamma,
	     cRVectorXd& y2, cRMatrixXd& Z, cRVectorXd& gamma0);
  };

  LVLMFit::LVLMFit() {
    fitOptions();
  }

  LVLMFit::LVLMFit(int n, int p) {
    fitOptions();
    n_ = n;
    p_ = p;
    ZtZ_.compute(MatrixXd::Identity(p_,p_));
    Zg_ = VectorXd::Zero(n);
    wy2_ = VectorXd::Zero(n);
  }

  void LVLMFit::fitOptions(int maxit, double epsilon) {
    maxit_ = maxit;
    epsilon_ = epsilon;
    return;
  }

  void LVLMFit::fit(RVectorXd gamma,
		    cRVectorXd& y2, cRMatrixXd& Z, cRVectorXd& gamma0) {
    // precomputations
    ZtZ_.compute(Z.adjoint() * Z);
    gamma = gamma0;
    for(int ii=0; ii<1; ii++) {
      // calculate score function
      Zg_ = Z * gamma;
      wy2_ = y2.array() / Zg_.array().exp() - 1.0;
      gamma += ZtZ_.solve((wy2_.asDiagonal() * Z).colwise().sum().adjoint());
    }
    return;
  }

}

#endif
