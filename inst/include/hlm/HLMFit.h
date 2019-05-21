/// @file HLMFit.h

#ifndef HLM_HLMFIT_H
#define HLM_HLMFIT_H 1

// [[Rcpp::depends("RcppEigen")]]
#include <RcppEigen.h>
#include "hlm/LMFit.h"
#include "hlm/LVLMFit.h"
// #include <iostream>

namespace hlm {
  
  using namespace Eigen;

  /// Heteroscedastic Linear Models.
  ///
  /// The heteroscedastic linear model (HLM) is defined as
  ///
  /// \f[
  /// y_i \stackrel{\textrm{ind}}{\sim} \mathcal{N}\big(\boldsymbol{x}_i'\boldsymbol{\beta}, \exp(\boldsymbol{z}_i'\boldsymbol{\gamma})\big).
  /// \f]
  class HLMFit {
  private:
    // eigen typedefs
    typedef const Ref<const VectorXd> cRVectorXd;
    typedef const Ref<const MatrixXd> cRMatrixXd;
    typedef Ref<VectorXd> RVectorXd;
    typedef Ref<MatrixXd> RMatrixXd;
    // problem dimensions
    int n_;
    int p_;
    int q_;
    // convergence output
    int niter_;
    double error_;
    double llik_;
    // memory allocation
    // LLT<MatrixXd> ZtZ_;
    VectorXd w_;
    VectorXd y2_;
    // control parameters
    int maxit_;
    double epsilon_;
    /// Relative error.
    double rel_err(double x_new, double x_old) {
      return fabs(x_new - x_old)/(0.1 + fabs(x_new));
    }
    // fitting components
    LMFit *lm_;
    LVLMFit *lvlm_;
  public:
    /// Default constructor.
    HLMFit();
    /// Constructor with memory allocation.
    HLMFit(int n, int p, int q);
    /// Destructor.
    ~HLMFit();
    /// Set fitting control parameters.
    void fitControl(int maxit = 1000, double epsilon = 1e-5);
    /// Convergence statistics from the fitting algorithm.
    void fitStats(double& llik, int& niter, double& error);
    /// Fit the MLE of the model parameters.
    void fit(RVectorXd beta, RVectorXd gamma,
	     cRVectorXd& y, cRMatrixXd& X, cRMatrixXd& Z,
	     cRVectorXd& beta0, cRVectorXd& gamma0);
    /// Loglikelihood function.
    double loglik(cRVectorXd& beta, cRVectorXd& gamma,
		  cRVectorXd& y, cRMatrixXd& X, cRMatrixXd& Z);
  };

  
  /// Initializes the control parameters and creates pointers to internal `LMFit` and `LVLMFit` objects.
  inline HLMFit::HLMFit() {
    fitControl();
    lm_ = new LMFit();
    lvlm_ = new LVLMFit();
  }

  /// Same as HLMFit::HLMFit, but pre-allocates memory for problems of a specific size.
  ///
  /// @param[in] n Integer number of observations.
  /// @param[in] p Integer number of mean covariates.
  /// @param[in] q Integer number of variance covariates.
  inline HLMFit::HLMFit(int n, int p, int q) {
    fitControl();
    n_ = n;
    p_ = p;
    q_ = q;
    w_ = VectorXd::Zero(n_);
    y2_ = VectorXd::Zero(n_);
    lm_ = new LMFit(n_, p_);
    lvlm_ = new LVLMFit(n_, q_);
  }

  inline HLMFit::~HLMFit() {
    delete lm_;
    delete lvlm_;
  }

  /// The iterative algorithm HLMFit::fit terminates when either a maximum number of iterations has been reached, or when the relative error reaches a specific bound.  Here, the relative error between step \f$m+1\f$ and step \f$m\f$ is defined as
  /// \f[
  /// \frac{|\ell(\boldsymbol{\gamma}_{m+1}) - \ell(\boldsymbol{\gamma}_{m})|}{0.1 + |\ell(\boldsymbol{\gamma}_{m+1})|},
  /// \f]
  /// where \f$\ell(\boldsymbol{\gamma}_{m})\f$ is the loglikelihoood evaluated at step \f$m\f$.
  ///
  /// @param[in] maxit Maximum number of Fisher scoring interations.
  /// @param[in] epsilon Relative tolerance on the loglikelihood.
  inline void HLMFit::fitControl(int maxit, double epsilon) {
    maxit_ = maxit;
    epsilon_ = epsilon;
    return;
  }

  /// @param[out] llik Loglikelihood at the fitted value of `beta` and `gamma`.
  /// @param[out] niter Number of iterations of the algorithm.
  /// @param[out] error Relative error between algorithm steps `niter-1` and `niter`.
  inline void HLMFit::fitStats(double& llik, int& niter, double& error) {
    llik = llik_;
    niter = niter_;
    error = error_;
    return;
  }

  /// @param[in] beta Mean parameter vector of length `p`.
  /// @param[in] gamma Variance parameter vector of length `q`.
  /// @param[in] y Vector of `n` of normal responses.
  /// @param[in] X Mean covariate matrix of size `n x p`.
  /// @param[in] Z Variance covariate matrix of size `n x q`.
  /// @return Loglikelihood (`double` scalar).
  inline double HLMFit::loglik(cRVectorXd& beta, cRVectorXd& gamma,
				cRVectorXd& y,
				cRMatrixXd& X, cRMatrixXd& Z) {
    y2_.noalias() = y - X * beta;
    y2_.array() *= y2_.array();
    return lvlm_->loglik(gamma, y2_, Z);
  }

  
  /// Iterates between fitting the mean parameters \f$\boldsymbol{\beta}\f$ and the variance parameters \f$\boldsymbol{\gamma}\f$ via calls to LMFit::fit and LVLMFit::fit, respectively.  The algorithm terminates when either the maximum number of steps has been reached, or the loglikelihood is not changing by more than a given tolerance.
  ///
  /// @param[out] beta MLE of the mean parameters; a vector of length `p`.
  /// @param[out] gamma MLE of the variance parameters; a vector of length `q`.
  /// @param[out] llik Loglikelihood at the fitted value of `beta` and `gamma` (scalar).
  /// @param[out] niter Number of iterations of the algorithm.
  /// @param[out] error Relative error between steps `niter` and `niter-1`.
  /// @param[in] y Vector of `n` of normal responses.
  /// @param[in] X Mean covariate matrix of size `n x p`.
  /// @param[in] Z Variance covariate matrix of size `n x q`.
  /// @param[in] beta0 Parameter vector of length `p` to initialize the algorithm.
  /// @param[in] gamma0 Parameter vector of length `q` to initialize the algorithm.
  inline void HLMFit::fit(RVectorXd beta, RVectorXd gamma,
			  cRVectorXd& y, cRMatrixXd& X, cRMatrixXd& Z,
			  cRVectorXd& beta0, cRVectorXd& gamma0) {
    double ll_new, ll_old;
    int tmp_niter;
    double tmp_error;
    // precomputations
    lvlm_->setZtZ(Z);
    // initialize
    beta = beta0;
    gamma = gamma0;
    // std::cout << "beta0 = \n" << beta << std::endl;
    // std::cout << "gamma0 = \n" << gamma << std::endl;
    // Rprintf("maxit_ = %i\n", maxit_);
    // Rprintf("epsilon_ = %f\n", epsilon_);
    ll_old = loglik(beta, gamma, y, X, Z);
    for(niter_=0; niter_<maxit_; niter_++) {
      // Rprintf("niter_ = %i\n", niter_);
      // update beta
      w_.noalias() = -Z * gamma;
      w_ = w_.array().exp();
      lm_->fit(beta, y, X, w_);
      // std::cout << "beta = \n" << beta << std::endl;
      // update gamma
      y2_.noalias() = y - X * beta;
      y2_.array() *= y2_.array();
      // std::cout << "y2_ = \n" << y2_ << std::endl;
      // Rprintf("all(y2_ > 0) = %i\n", (y2_.array() > 0).all());
      lvlm_->fit(gamma, y2_, Z, false);
      lvlm_->fitStats(ll_new, tmp_niter, tmp_error);
      // std::cout << "gamma = \n" << gamma << std::endl;
      // check relative error
      error_ = rel_err(ll_new, ll_old);
      // Rprintf("ll_old = %f, ll_new = %f\n", ll_old, ll_new);
      // Rprintf("error_[%i] = %f\n", niter_, error_);
      if(error_ < epsilon_) {
	break;
      } else {
	ll_old = ll_new;
      }
    }
    llik_ = ll_new;
    return;
  }

}

#endif

