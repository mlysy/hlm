/// @file LVLMFit.h
///
/// **TODO:**
///
/// - HLMFit::fit should return loglikelihood as well.
///
/// - HLMFit::fit should be overloaded such that precomputed `ZtZ_` can be supplied.  In fact, should perhaps restructure io arguments:
///
///     - `niter` and `error` (and even `llik`) can be returned separately via e.g. `fitStats`.
///     - boolean argument `calcZtZ` combined with `setZtZ` method.
///
/// - check aliasing!  In fact, am I understanding "temporaries" correctly, i.e., is is more/less efficient to preallocate memory for temporary storage?
///
///    Answer: I think more efficient, because otherwise how would the compiler know the size of the temporary beforehand?  Can it really figure this out at compile-time? Btw, see `Eigen::internal::set_is_malloc_allowed` and discussion [here](https://forum.kde.org/viewtopic.php?f=74&t=124882).
///
///    OK, thanks to `set_is_malloc_allowed` this is now done.
///
/// - stop putting all output arguments at the beginning of the function.

#ifndef HLM_LVLMFIT_H
#define HLM_LVLMFIT_H 1

// [[Rcpp::depends("RcppEigen")]]
#include <RcppEigen.h>
// #include <iostream>

namespace hlm {

  using namespace Eigen;
  static const double RHO = -1.270362845461477707686981375445611775;

  /// Log-Variance Linear Models.
  ///
  /// The log-variance linear model (LVLM) is of the form
  ///
  /// \f[
  /// y_i \stackrel{\textrm{ind}}{\sim} \mathcal{N}\big(0, \exp(\boldsymbol{z}_i'\boldsymbol{\gamma})\big).
  /// \f]
  class LVLMFit {
  private:
    // eigen typedefs
    typedef const Ref<const VectorXd> cRVectorXd;
    typedef const Ref<const MatrixXd> cRMatrixXd;
    typedef Ref<VectorXd> RVectorXd;
    typedef Ref<MatrixXd> RMatrixXd;
    // problem dimensions
    int n_;
    int p_;
    // memory allocation
    LLT<MatrixXd> ZtZ_;
    VectorXd wy2_;
    VectorXd logY2_;
    VectorXd Zg_;
    VectorXd score_;
    // control parameters
    int maxit_;
    double epsilon_;
    // convergence output
    int niter_;
    double error_;
    double llik_;
    /// Loglikelihood with precomputations.
    double loglik();
    /// Relative error.
    double rel_err(double x_new, double x_old) {
      return fabs(x_new - x_old)/(0.1 + fabs(x_new));
    }
  public:
    /// Set fitting control parameters.
    void fitControl(int maxit = 25, double epsilon = 1e-8);
    /// Convergence statistics from the fitting algorithm.
    void fitStats(double& llik, int& niter, double& error);
    /// Initialize the Fisher scoring matrix solver.
    void setZtZ(cRMatrixXd& Z);
    /// Default constructor.
    LVLMFit();
    /// Constructor with memory allocation.
    LVLMFit(int n, int p);
    /// Loglikelihood function.
    double loglik(cRVectorXd& gamma, cRVectorXd& y2, cRMatrixXd& Z);
    /// Fit the MLE of an LVLM.
    void fit(RVectorXd gamma,
	     cRVectorXd& y2, cRMatrixXd& Z,
	     cRVectorXd& gamma0, bool calcZtZ = true);
    /// Fit the MLE of an LVLM.
    void fit(RVectorXd gamma,
    	     cRVectorXd& y2, cRMatrixXd& Z, bool calcZtZ = true);
    /// Least-Squares estimate of the parameter vector
    void fitLS(RVectorXd gamma, cRVectorXd& logY2, cRMatrixXd& Z,
	       bool calcZtZ = true);
  };

  /// Just initializes the control parameters.
  inline LVLMFit::LVLMFit() {
    fitControl();
  }

  /// Pre-allocates memory for problems of a specific size.
  ///
  /// @param[in] n Integer number of observations.
  /// @param[in] p Integer number of covariates.
  inline LVLMFit::LVLMFit(int n, int p) {
    fitControl();
    n_ = n;
    p_ = p;
    ZtZ_.compute(MatrixXd::Identity(p_,p_));
    Zg_ = VectorXd::Zero(n_);
    wy2_ = VectorXd::Zero(n_);
    logY2_ = VectorXd::Zero(n_);
    score_ = VectorXd::Zero(p_);
  }

  /// The Fisher scoring algorithm LVLMFit::fit terminates when either a maximum number of iterations has been reached, or when the relative error reaches a specific bound.  Here, the relative error between step \f$m+1\f$ and step \f$m\f$ is defined as
  /// \f[
  /// \frac{|\ell(\boldsymbol{\gamma}_{m+1}) - \ell(\boldsymbol{\gamma}_{m})|}{0.1 + |\ell(\boldsymbol{\gamma}_{m+1})|},
  /// \f]
  /// where \f$\ell(\boldsymbol{\gamma}_{m})\f$ is the loglikelihoood evaluated at step \f$m\f$.
  ///
  /// @param[in] maxit Maximum number of Fisher scoring interations.
  /// @param[in] epsilon Relative tolerance on the loglikelihood.
  inline void LVLMFit::fitControl(int maxit, double epsilon) {
    maxit_ = maxit;
    epsilon_ = epsilon;
    return;
  }

  /// Calculates `ZtZ_` when it is more efficient to reuse it for Fisher scoring with the same covariate matrix but multiple response vectors.
  ///
  /// @param[in] Z Covariate matrix of size `n x p`.
  inline void LVLMFit::setZtZ(cRMatrixXd& Z) {
    ZtZ_.compute(Z.adjoint() * Z);
    return;
  }
  
  /// @param[in] gamma Parameter vector of length `p` at which to calculate the loglikelihood.
  /// @param[in] y2 Vector of `n` observations.
  /// @param[in] Z Covariate matrix of size `n x p`.
  /// @return Loglikelihood (`double` scalar).
  inline double LVLMFit::loglik(cRVectorXd& gamma, cRVectorXd& y2, cRMatrixXd& Z) {
    Zg_.noalias() = Z * gamma;
    wy2_ = y2.array() / Zg_.array().exp();
    return loglik();
  }

  /// Calculates the loglikelihood based on the internal values of `Zg_ = Z * gamma` and `wy2_ = y2/exp(Zg_)`.
  /// @return Loglikelihood (`double` scalar).
  inline double LVLMFit::loglik() {
    return -.5 * (wy2_ + Zg_).sum();
  }
  
  /// @param[out] llik Loglikelihood at the fitted value of `gamma`.
  /// @param[out] niter Number of iterations of the algorithm.
  /// @param[out] error Relative error between algorithm steps `niter-1` and `niter`.
  inline void LVLMFit::fitStats(double& llik, int& niter, double& error) {
    llik = llik_;
    niter = niter_;
    error = error_;
    return;
  }
  
  /// Applies a Fisher scoring algorithm to compute the MLE.  The algorithm terminates when either `maxit_` iterations has been reached, or when the relative error falls below `epsilon_`.  Here, the relative error between step \f$m+1\f$ and step \f$m\f$ is defined as
  /// \f[
  /// \frac{|\ell(\boldsymbol{\gamma}_{m+1}) - \ell(\boldsymbol{\gamma}_{m})|}{0.1 + |\ell(\boldsymbol{\gamma}_{m+1})|},
  /// \f]
  /// where \f$\ell(\boldsymbol{\gamma}_{m})\f$ is the loglikelihoood evaluated at step \f$m\f$.
  ///
  /// @param[out] gamma MLE vector of length `p` as calculated by the Fisher scoring algorithm.
  /// @param[in] y2 Vector of `n` *squares* of normal responses: `y2 = y^2`.
  /// @param[in] Z Covariate matrix of size `n x p`.
  /// @param[in] gamma0 Parameter vector of length `p` to initialize the algorithm.
  /// @param[in] calcZtZ If `false`, uses the Cholesky solver `ZtZ_` computed from a previous step, or from LVLMFit::setZtZ.
  inline void LVLMFit::fit(RVectorXd gamma,
			   cRVectorXd& y2, cRMatrixXd& Z,
			   cRVectorXd& gamma0, bool calcZtZ) {
    double ll_new, ll_old;
    if(calcZtZ) {
      // precomputations
      ZtZ_.compute(Z.adjoint() * Z);
    }
    gamma = gamma0;
    // Rprintf("epsilon_ = %f\n", epsilon_);
    for(niter_=0; niter_<maxit_; niter_++) {
      // elements of score function
      Zg_.noalias() = Z * gamma;
      wy2_ = y2.array() / Zg_.array().exp();
      // likelihood along the way
      ll_new = loglik();
      // compute and check relative error
      error_ = (niter_ > 0) ? rel_err(ll_new, ll_old) : (epsilon_ + 1.0);
      // Rprintf("ll[%i] = %f\n", niter, ll_new);
      // Rprintf("error[%i] = %f\n", niter, error);
      if(error_ < epsilon_) {
	// Rprintf("tolerance reached.\n");
	break;
      } else {
	ll_old = ll_new;
      }
      // finish score function & update gamma
      wy2_.array() -= 1.0; 
      score_ = (wy2_.asDiagonal() * Z).colwise().sum().adjoint();
      ZtZ_.solveInPlace(score_);
      gamma += score_;
    }
    if(niter_ == maxit_) {
      // final error calculation
      ll_new = loglik(gamma, y2, Z);
      error_ = rel_err(ll_new, ll_old);
      // Rprintf("ll[%i] = %f\n", niter, ll_new);
      // Rprintf("error[%i] = %f\n", niter, error);
    }
    llik_ = ll_new;
    return;
  }

  /// Same as the longer version of LVLMFit::fit, but initializes `gamma0` with the least-squares estimator LVLMFit::fitLS.
  ///
  /// @param[out] gamma MLE vector of length `p` as calculated by the Fisher scoring algorithm.
  /// @param[in] y2 Vector of `n` *squares* of normal responses: `y2 = y^2`.
  /// @param[in] Z Covariate matrix of size `n x p`.
  /// @param[in] calcZtZ If `false`, uses the Cholesky solver `ZtZ_` computed from a previous step, or from LVLMFit::setZtZ.
  inline void LVLMFit::fit(RVectorXd gamma,
  			   cRVectorXd& y2, cRMatrixXd& Z, bool calcZtZ) {
    if(calcZtZ) setZtZ(Z); // precomputations
    logY2_ = y2.array().log();
    fitLS(gamma, logY2_, Z, false);
    // std::cout << "fitLS(gamma) = \n" << gamma << std::endl;
    fit(gamma, y2, Z, gamma, false);
    // std::cout << "fit(gamma) = \n" << gamma << std::endl;
    return;
  }

  /// The least-squares parameter estimation method uses the fact that
  ///
  /// \f[
  /// \log(y_i^2) = \boldsymbol{z}_i'\boldsymbol{\gamma} + \varepsilon_i, \qquad \exp(\varepsilon_i) \stackrel{\textrm{iid}}{\sim} \chi^2_{(1)}.
  /// \f]
  /// Since \f$\rho = E[\varepsilon_i]\f$ can be calculated analytically, the least-squares estimator for \f$\boldsymbol{\gamma}\f$ is the OLS estimator with observation \f$i\f$ having covariates \f$\boldsymbol{z}_i\f$ and response \f$u_i = \log(y_i^2)-\rho\f$.
  ///
  /// @param[out] gamma Fitted value of the LVLM parameters.
  /// @param[in] logY2 Vector of `n` log-square responses, `logY2 = log(y^2)`.
  /// @param[in] Z Covariate matrix of size `n x p`.
  /// @param[in] calcZtZ If `false`, uses the Cholesky solver `ZtZ_` computed from a previous step, or from LVLMFit::setZtZ.
  inline void LVLMFit::fitLS(RVectorXd gamma, cRVectorXd& logY2,
			     cRMatrixXd& Z, bool calcZtZ) {
    if(calcZtZ) {
      // precomputations
      ZtZ_.compute(Z.adjoint() * Z);
    }
    logY2_ = logY2.array() - RHO;
    gamma.noalias() = Z.adjoint() * logY2_;
    ZtZ_.solveInPlace(gamma);
    return;
  }
  
}

#endif
