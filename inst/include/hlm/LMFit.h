/// @file LMFit.h

#ifndef HLM_LMFIT_H
#define HLM_LMFIT_H 1

// [[Rcpp::depends("RcppEigen")]]
#include <RcppEigen.h>

namespace hlm {

  using namespace Eigen;
  typedef const Ref<const VectorXd> cRVectorXd;
  typedef const Ref<const MatrixXd> cRMatrixXd;
  typedef Ref<VectorXd> RVectorXd;
  typedef Ref<MatrixXd> RMatrixXd;

  class LMFit {
  private:
    int n_;
    int p_;
    LLT<MatrixXd> XtX_;
    MatrixXd Xw_;
  public:
    LMFit() {}
    LMFit(int n, int p);
    void fit(RVectorXd beta, cRVectorXd& y, cRMatrixXd& X, cRVectorXd& w);
  };

  LMFit::LMFit(int n, int p) {
    n_ = n;
    p_ = p;
    XtX_.compute(MatrixXd::Identity(p_,p_));
    Xw_ = MatrixXd::Zero(n_,p_);
  }

  void LMFit::fit(RVectorXd beta, cRVectorXd& y, cRMatrixXd& X, cRVectorXd& w) {
    Xw_ = w.asDiagonal() * X;
    XtX_.compute(Xw_.adjoint() * X);
    beta = XtX_.solve(Xw_.adjoint() * y);
    return;
  }

}

#endif
