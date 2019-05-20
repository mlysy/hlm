/// @file LMFit.h

#ifndef HLM_LMFIT_H
#define HLM_LMFIT_H 1

// [[Rcpp::depends("RcppEigen")]]
#include <RcppEigen.h>

namespace hlm {

  using namespace Eigen;

  /// Weighted Linear Models (LM).
  ///
  /// The weighted linear model (LM) is of the form
  ///
  /// \f[
  /// y_i \stackrel{\textrm{ind}}{\sim} \mathcal{N}(\boldsymbol{x}_i'\boldsymbol{\beta}, \sigma^2/w_i).
  /// \f]
  class LMFit {
  private:
    // eigen typedefs
    typedef const Ref<const VectorXd> cRVectorXd;
    typedef const Ref<const MatrixXd> cRMatrixXd;
    typedef Ref<VectorXd> RVectorXd;
    typedef Ref<MatrixXd> RMatrixXd;
    // dimensions of problem
    int n_;
    int p_;
    // memory allocation
    LLT<MatrixXd> XtX_;
    MatrixXd Xw_;
  public:
    /// Default constructor.
    LMFit() {}
    /// Constructor with memory allocation.
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
  inline LMFit::LMFit(int n, int p) {
    n_ = n;
    p_ = p;
    XtX_.compute(MatrixXd::Identity(p_,p_));
    Xw_ = MatrixXd::Zero(n_,p_);
  }

  /// @param[out] beta MLE vector.
  /// @param[in] y Vector of `n` responses.
  /// @param[in] X `n x p` covariate matrix.
  /// @note Should not need to store `X.adjoint() * y`, as solving triangular systems only requires the RHS to be provided sequentially...
  inline void LMFit::fit(RVectorXd beta, cRVectorXd& y, cRMatrixXd& X) {
    XtX_.compute(X.adjoint() * X);
    beta = XtX_.solve(X.adjoint() * y);
    return;
  }

  /// @param[out] beta MLE vector.
  /// @param[in] y Vector of `n` responses.
  /// @param[in] X `n x p` covariate matrix.
  /// @param[in] w Vector of `n` positive weights.
  inline void LMFit::fit(RVectorXd beta, cRVectorXd& y, cRMatrixXd& X, cRVectorXd& w) {
    Xw_ = w.asDiagonal() * X;
    XtX_.compute(Xw_.adjoint() * X);
    beta = XtX_.solve(Xw_.adjoint() * y);
    return;
  }

}

#endif
