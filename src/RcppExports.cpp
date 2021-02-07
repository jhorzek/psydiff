// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include "../inst/include/psydiff.h"
#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// setParameterValues
Rcpp::DataFrame setParameterValues(Rcpp::DataFrame parameterTable, Rcpp::NumericVector parameterValues, Rcpp::StringVector parameterLabels);
RcppExport SEXP _psydiff_setParameterValues(SEXP parameterTableSEXP, SEXP parameterValuesSEXP, SEXP parameterLabelsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::DataFrame >::type parameterTable(parameterTableSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type parameterValues(parameterValuesSEXP);
    Rcpp::traits::input_parameter< Rcpp::StringVector >::type parameterLabels(parameterLabelsSEXP);
    rcpp_result_gen = Rcpp::wrap(setParameterValues(parameterTable, parameterValues, parameterLabels));
    return rcpp_result_gen;
END_RCPP
}
// getParameterValues
Rcpp::NumericVector getParameterValues(Rcpp::List psydiffModel);
RcppExport SEXP _psydiff_getParameterValues(SEXP psydiffModelSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type psydiffModel(psydiffModelSEXP);
    rcpp_result_gen = Rcpp::wrap(getParameterValues(psydiffModel));
    return rcpp_result_gen;
END_RCPP
}
// setParameterList
void setParameterList(const Rcpp::DataFrame& parameterTable, Rcpp::List& parameterList, int person);
RcppExport SEXP _psydiff_setParameterList(SEXP parameterTableSEXP, SEXP parameterListSEXP, SEXP personSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::DataFrame& >::type parameterTable(parameterTableSEXP);
    Rcpp::traits::input_parameter< Rcpp::List& >::type parameterList(parameterListSEXP);
    Rcpp::traits::input_parameter< int >::type person(personSEXP);
    setParameterList(parameterTable, parameterList, person);
    return R_NilValue;
END_RCPP
}
// getSigmaPoints
arma::mat getSigmaPoints(arma::colvec m, arma::mat& A, const bool& covarianceIsRoot, const double& c);
RcppExport SEXP _psydiff_getSigmaPoints(SEXP mSEXP, SEXP ASEXP, SEXP covarianceIsRootSEXP, SEXP cSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::colvec >::type m(mSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type A(ASEXP);
    Rcpp::traits::input_parameter< const bool& >::type covarianceIsRoot(covarianceIsRootSEXP);
    Rcpp::traits::input_parameter< const double& >::type c(cSEXP);
    rcpp_result_gen = Rcpp::wrap(getSigmaPoints(m, A, covarianceIsRoot, c));
    return rcpp_result_gen;
END_RCPP
}
// getMeanWeights
arma::colvec getMeanWeights(const int& n, const double& alpha, const double& beta, const double& kappa);
RcppExport SEXP _psydiff_getMeanWeights(SEXP nSEXP, SEXP alphaSEXP, SEXP betaSEXP, SEXP kappaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int& >::type n(nSEXP);
    Rcpp::traits::input_parameter< const double& >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const double& >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< const double& >::type kappa(kappaSEXP);
    rcpp_result_gen = Rcpp::wrap(getMeanWeights(n, alpha, beta, kappa));
    return rcpp_result_gen;
END_RCPP
}
// getCovWeights
arma::colvec getCovWeights(const int& n, const double& alpha, const double& beta, const double& kappa);
RcppExport SEXP _psydiff_getCovWeights(SEXP nSEXP, SEXP alphaSEXP, SEXP betaSEXP, SEXP kappaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int& >::type n(nSEXP);
    Rcpp::traits::input_parameter< const double& >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const double& >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< const double& >::type kappa(kappaSEXP);
    rcpp_result_gen = Rcpp::wrap(getCovWeights(n, alpha, beta, kappa));
    return rcpp_result_gen;
END_RCPP
}
// getWMatrix
arma::mat getWMatrix(const arma::colvec& meanWeights, const arma::colvec& covWeights);
RcppExport SEXP _psydiff_getWMatrix(SEXP meanWeightsSEXP, SEXP covWeightsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec& >::type meanWeights(meanWeightsSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type covWeights(covWeightsSEXP);
    rcpp_result_gen = Rcpp::wrap(getWMatrix(meanWeights, covWeights));
    return rcpp_result_gen;
END_RCPP
}
// computeMeanFromSigmaPoints
arma::colvec computeMeanFromSigmaPoints(const arma::mat& sigmaPoints, const arma::colvec& meanWeights);
RcppExport SEXP _psydiff_computeMeanFromSigmaPoints(SEXP sigmaPointsSEXP, SEXP meanWeightsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type sigmaPoints(sigmaPointsSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type meanWeights(meanWeightsSEXP);
    rcpp_result_gen = Rcpp::wrap(computeMeanFromSigmaPoints(sigmaPoints, meanWeights));
    return rcpp_result_gen;
END_RCPP
}
// getAFromSigmaPoints
arma::mat getAFromSigmaPoints(const arma::mat& sigmaPoints, const arma::colvec& meanWeights, const double& c);
RcppExport SEXP _psydiff_getAFromSigmaPoints(SEXP sigmaPointsSEXP, SEXP meanWeightsSEXP, SEXP cSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type sigmaPoints(sigmaPointsSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type meanWeights(meanWeightsSEXP);
    Rcpp::traits::input_parameter< const double& >::type c(cSEXP);
    rcpp_result_gen = Rcpp::wrap(getAFromSigmaPoints(sigmaPoints, meanWeights, c));
    return rcpp_result_gen;
END_RCPP
}
// computePFromSigmaPoints
arma::mat computePFromSigmaPoints(const arma::mat& sigmaPoints, const arma::colvec& meanWeights, const double& c);
RcppExport SEXP _psydiff_computePFromSigmaPoints(SEXP sigmaPointsSEXP, SEXP meanWeightsSEXP, SEXP cSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type sigmaPoints(sigmaPointsSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type meanWeights(meanWeightsSEXP);
    Rcpp::traits::input_parameter< const double& >::type c(cSEXP);
    rcpp_result_gen = Rcpp::wrap(computePFromSigmaPoints(sigmaPoints, meanWeights, c));
    return rcpp_result_gen;
END_RCPP
}
// getPhi
arma::mat getPhi(const arma::mat& squareMatrix);
RcppExport SEXP _psydiff_getPhi(SEXP squareMatrixSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type squareMatrix(squareMatrixSEXP);
    rcpp_result_gen = Rcpp::wrap(getPhi(squareMatrix));
    return rcpp_result_gen;
END_RCPP
}
// predictMu
arma::mat predictMu(const arma::mat& Y_, const arma::colvec& meanWeights);
RcppExport SEXP _psydiff_predictMu(SEXP Y_SEXP, SEXP meanWeightsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type Y_(Y_SEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type meanWeights(meanWeightsSEXP);
    rcpp_result_gen = Rcpp::wrap(predictMu(Y_, meanWeights));
    return rcpp_result_gen;
END_RCPP
}
// predictS
arma::mat predictS(const arma::mat& Y_, const arma::mat& W, const arma::mat& R);
RcppExport SEXP _psydiff_predictS(SEXP Y_SEXP, SEXP WSEXP, SEXP RSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type Y_(Y_SEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type W(WSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type R(RSEXP);
    rcpp_result_gen = Rcpp::wrap(predictS(Y_, W, R));
    return rcpp_result_gen;
END_RCPP
}
// predictC
arma::mat predictC(const arma::mat& sigmaPoints, const arma::mat& Y_, const arma::mat& W);
RcppExport SEXP _psydiff_predictC(SEXP sigmaPointsSEXP, SEXP Y_SEXP, SEXP WSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type sigmaPoints(sigmaPointsSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Y_(Y_SEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type W(WSEXP);
    rcpp_result_gen = Rcpp::wrap(predictC(sigmaPoints, Y_, W));
    return rcpp_result_gen;
END_RCPP
}
// computeK
arma::mat computeK(const arma::mat& C, const arma::mat& S);
RcppExport SEXP _psydiff_computeK(SEXP CSEXP, SEXP SSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type C(CSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type S(SSEXP);
    rcpp_result_gen = Rcpp::wrap(computeK(C, S));
    return rcpp_result_gen;
END_RCPP
}
// updateM
arma::mat updateM(const arma::colvec& m_, const arma::mat& K, const arma::colvec& residual);
RcppExport SEXP _psydiff_updateM(SEXP m_SEXP, SEXP KSEXP, SEXP residualSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec& >::type m_(m_SEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type K(KSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type residual(residualSEXP);
    rcpp_result_gen = Rcpp::wrap(updateM(m_, K, residual));
    return rcpp_result_gen;
END_RCPP
}
// updateP
arma::mat updateP(const arma::mat& P_, const arma::mat& K, const arma::mat& S);
RcppExport SEXP _psydiff_updateP(SEXP P_SEXP, SEXP KSEXP, SEXP SSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type P_(P_SEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type K(KSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type S(SSEXP);
    rcpp_result_gen = Rcpp::wrap(updateP(P_, K, S));
    return rcpp_result_gen;
END_RCPP
}
// computeIndividualM2LL
double computeIndividualM2LL(const int& nObservedVariables, const arma::colvec& rawData, const arma::colvec& expectedMeans, const arma::mat& expectedCovariance);
RcppExport SEXP _psydiff_computeIndividualM2LL(SEXP nObservedVariablesSEXP, SEXP rawDataSEXP, SEXP expectedMeansSEXP, SEXP expectedCovarianceSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int& >::type nObservedVariables(nObservedVariablesSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type rawData(rawDataSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type expectedMeans(expectedMeansSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type expectedCovariance(expectedCovarianceSEXP);
    rcpp_result_gen = Rcpp::wrap(computeIndividualM2LL(nObservedVariables, rawData, expectedMeans, expectedCovariance));
    return rcpp_result_gen;
END_RCPP
}
// cholupdate
arma::mat cholupdate(arma::mat L, arma::colvec x, double v, std::string direction);
RcppExport SEXP _psydiff_cholupdate(SEXP LSEXP, SEXP xSEXP, SEXP vSEXP, SEXP directionSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type L(LSEXP);
    Rcpp::traits::input_parameter< arma::colvec >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type v(vSEXP);
    Rcpp::traits::input_parameter< std::string >::type direction(directionSEXP);
    rcpp_result_gen = Rcpp::wrap(cholupdate(L, x, v, direction));
    return rcpp_result_gen;
END_RCPP
}
// qr_
arma::mat qr_(arma::mat X);
RcppExport SEXP _psydiff_qr_(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(qr_(X));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_psydiff_setParameterValues", (DL_FUNC) &_psydiff_setParameterValues, 3},
    {"_psydiff_getParameterValues", (DL_FUNC) &_psydiff_getParameterValues, 1},
    {"_psydiff_setParameterList", (DL_FUNC) &_psydiff_setParameterList, 3},
    {"_psydiff_getSigmaPoints", (DL_FUNC) &_psydiff_getSigmaPoints, 4},
    {"_psydiff_getMeanWeights", (DL_FUNC) &_psydiff_getMeanWeights, 4},
    {"_psydiff_getCovWeights", (DL_FUNC) &_psydiff_getCovWeights, 4},
    {"_psydiff_getWMatrix", (DL_FUNC) &_psydiff_getWMatrix, 2},
    {"_psydiff_computeMeanFromSigmaPoints", (DL_FUNC) &_psydiff_computeMeanFromSigmaPoints, 2},
    {"_psydiff_getAFromSigmaPoints", (DL_FUNC) &_psydiff_getAFromSigmaPoints, 3},
    {"_psydiff_computePFromSigmaPoints", (DL_FUNC) &_psydiff_computePFromSigmaPoints, 3},
    {"_psydiff_getPhi", (DL_FUNC) &_psydiff_getPhi, 1},
    {"_psydiff_predictMu", (DL_FUNC) &_psydiff_predictMu, 2},
    {"_psydiff_predictS", (DL_FUNC) &_psydiff_predictS, 3},
    {"_psydiff_predictC", (DL_FUNC) &_psydiff_predictC, 3},
    {"_psydiff_computeK", (DL_FUNC) &_psydiff_computeK, 2},
    {"_psydiff_updateM", (DL_FUNC) &_psydiff_updateM, 3},
    {"_psydiff_updateP", (DL_FUNC) &_psydiff_updateP, 3},
    {"_psydiff_computeIndividualM2LL", (DL_FUNC) &_psydiff_computeIndividualM2LL, 4},
    {"_psydiff_cholupdate", (DL_FUNC) &_psydiff_cholupdate, 4},
    {"_psydiff_qr_", (DL_FUNC) &_psydiff_qr_, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_psydiff(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
