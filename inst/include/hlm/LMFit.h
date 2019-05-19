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

  /// Weighted linear regression models.
  class LMFit {
  private:
    int n_;
    int p_;
    LLT<MatrixXd> XtX_;
    MatrixXd Xw_;
  public:
    /// Default constructor.
    LMFit() {}
    /// Size-specific constructor.
    LMFit(int n, int p);
    /// Fit the MLE of an LM.
    void fit(RVectorXd beta, cRVectorXd& y, cRMatrixXd& X);
    /// Fit the MLE of a weighted LM.
    void fit(RVectorXd beta, cRVectorXd& y, cRMatrixXd& X, cRVectorXd& w);
  };

  /// Pre-allocates memory for problems of a specific size.
  ///
  /// @param[in] n Integer number of observations.
  /// @param[in] p Integer number of covariates.
  LMFit::LMFit(int n, int p) {
    n_ = n;
    p_ = p;
    XtX_.compute(MatrixXd::Identity(p_,p_));
    Xw_ = MatrixXd::Zero(n_,p_);
  }

  /// @param[out] beta MLE vector.
  /// @param[in] y Vector of `n` responses.
  /// @param[in] X `n x p` covariate matrix.
  /// @note Should not need to store `X.adjoint() * y`, as solving triangular systems only requires the RHS to be provided sequentially...
  void LMFit::fit(RVectorXd beta, cRVectorXd& y, cRMatrixXd& X) {
    XtX_.compute(X.adjoint() * X);
    beta = XtX_.solve(X.adjoint() * y);
    return;
  }

  /// @param[out] beta MLE vector.
  /// @param[in] y Vector of `n` responses.
  /// @param[in] X `n x p` covariate matrix.
  /// @param[in] w Vector of `n` positive weights.
  void LMFit::fit(RVectorXd beta, cRVectorXd& y, cRMatrixXd& X, cRVectorXd& w) {
    Xw_ = w.asDiagonal() * X;
    XtX_.compute(Xw_.adjoint() * X);
    beta = XtX_.solve(Xw_.adjoint() * y);
    return;
  }

}

#endif
