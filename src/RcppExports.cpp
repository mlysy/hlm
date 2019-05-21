// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

// wlm_fit
Eigen::VectorXd wlm_fit(Eigen::VectorXd y, Eigen::MatrixXd X, Eigen::VectorXd w);
RcppExport SEXP _hlm_wlm_fit(SEXP ySEXP, SEXP XSEXP, SEXP wSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type y(ySEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type X(XSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type w(wSEXP);
    rcpp_result_gen = Rcpp::wrap(wlm_fit(y, X, w));
    return rcpp_result_gen;
END_RCPP
}
// lm_fit
Eigen::VectorXd lm_fit(Eigen::VectorXd y, Eigen::MatrixXd X);
RcppExport SEXP _hlm_lm_fit(SEXP ySEXP, SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type y(ySEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(lm_fit(y, X));
    return rcpp_result_gen;
END_RCPP
}
// lvlm_fit
Eigen::VectorXd lvlm_fit(Eigen::VectorXd y2, Eigen::MatrixXd Z, Eigen::VectorXd gamma0, int maxit, double epsilon, bool initLS);
RcppExport SEXP _hlm_lvlm_fit(SEXP y2SEXP, SEXP ZSEXP, SEXP gamma0SEXP, SEXP maxitSEXP, SEXP epsilonSEXP, SEXP initLSSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type y2(y2SEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type gamma0(gamma0SEXP);
    Rcpp::traits::input_parameter< int >::type maxit(maxitSEXP);
    Rcpp::traits::input_parameter< double >::type epsilon(epsilonSEXP);
    Rcpp::traits::input_parameter< bool >::type initLS(initLSSEXP);
    rcpp_result_gen = Rcpp::wrap(lvlm_fit(y2, Z, gamma0, maxit, epsilon, initLS));
    return rcpp_result_gen;
END_RCPP
}
// lvlm_fitLS
Eigen::VectorXd lvlm_fitLS(Eigen::VectorXd logY2, Eigen::MatrixXd Z);
RcppExport SEXP _hlm_lvlm_fitLS(SEXP logY2SEXP, SEXP ZSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type logY2(logY2SEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type Z(ZSEXP);
    rcpp_result_gen = Rcpp::wrap(lvlm_fitLS(logY2, Z));
    return rcpp_result_gen;
END_RCPP
}
// hlm_fit
Rcpp::List hlm_fit(Eigen::VectorXd y, Eigen::MatrixXd X, Eigen::MatrixXd Z, Eigen::VectorXd beta0, Eigen::VectorXd gamma0, int maxit, double epsilon);
RcppExport SEXP _hlm_hlm_fit(SEXP ySEXP, SEXP XSEXP, SEXP ZSEXP, SEXP beta0SEXP, SEXP gamma0SEXP, SEXP maxitSEXP, SEXP epsilonSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type y(ySEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type X(XSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type beta0(beta0SEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type gamma0(gamma0SEXP);
    Rcpp::traits::input_parameter< int >::type maxit(maxitSEXP);
    Rcpp::traits::input_parameter< double >::type epsilon(epsilonSEXP);
    rcpp_result_gen = Rcpp::wrap(hlm_fit(y, X, Z, beta0, gamma0, maxit, epsilon));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_hlm_wlm_fit", (DL_FUNC) &_hlm_wlm_fit, 3},
    {"_hlm_lm_fit", (DL_FUNC) &_hlm_lm_fit, 2},
    {"_hlm_lvlm_fit", (DL_FUNC) &_hlm_lvlm_fit, 6},
    {"_hlm_lvlm_fitLS", (DL_FUNC) &_hlm_lvlm_fitLS, 2},
    {"_hlm_hlm_fit", (DL_FUNC) &_hlm_hlm_fit, 7},
    {NULL, NULL, 0}
};

RcppExport void R_init_hlm(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
