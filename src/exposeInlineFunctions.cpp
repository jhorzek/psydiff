// Call the package's header function
#include <psydiff.h>

//// [[Rcpp::export]]
//double test2_R(double a){return test2(a);}

// [[Rcpp::export]]
Rcpp::List panelSDE(arma::mat observed, Rcpp::DataFrame parameterTable){
  return(panelSDE_C(observed, parameterTable));
}

// [[Rcpp::export]]
void setParameterValues(Rcpp::DataFrame &parameterTable,
                        Rcpp::NumericVector parameterValues,
                        Rcpp::StringVector parameterLabels){
  setParameterValues_C(parameterTable,
                       parameterValues,
                       parameterLabels);
}

// [[Rcpp::export]]
Rcpp::NumericVector getParameterValues(Rcpp::List panelSDEModel){
  return getParameterValues_C(panelSDEModel);
}

// [[Rcpp::export]]
void setParameterTable(const Rcpp::DataFrame &parameterTable, Rcpp::List &parameterList,
              int person){
  return setParameterTable_C(parameterTable, parameterList, person);
}

// [[Rcpp::export]]
arma::mat getSigmaPoints(arma::colvec m, arma::mat &A,
                         const bool &covarianceIsRoot,
                         const double &c){
  return (getSigmaPoints_C(m, A, covarianceIsRoot, c));
}

// [[Rcpp::export]]
arma::colvec getMeanWeights(const int &n, const double &alpha,
                            const double &beta, const double &kappa){
  return(getMeanWeights_C(n, alpha, beta, kappa));
}

// [[Rcpp::export]]
arma::mat getWMatrix(const arma::colvec &meanWeights, const arma::colvec &covWeights){
  return(getWMatrix_C(meanWeights, covWeights));
}

// [[Rcpp::export]]
arma::colvec computeMeanFromSigmaPoints(const arma::mat &sigmaPoints,
                                        const arma::colvec &meanWeights){
  return(computeMeanFromSigmaPoints_C(sigmaPoints,
                                      meanWeights));
}

// [[Rcpp::export]]
arma::mat getAFromSigmaPoints(const arma::mat &sigmaPoints,
                              const arma::colvec &meanWeights, const double &c){
  return(getAFromSigmaPoints_C(sigmaPoints,
                               meanWeights, c));
}

// [[Rcpp::export]]
arma::mat computePFromSigmaPoints(const arma::mat &sigmaPoints,
                                  const arma::colvec &meanWeights,
                                  const double &c){
  return(computePFromSigmaPoints_C(sigmaPoints,
                                   meanWeights,
                                   c));
}

// [[Rcpp::export]]
arma::mat getPhi(const arma::mat &squareMatrix){
  return(getPhi_C(squareMatrix));
}

// [[Rcpp::export]]
arma::mat predictMu(const arma::mat &Y_, const arma::colvec &meanWeights){
  return(predictMu_C(Y_, meanWeights));
}

// [[Rcpp::export]]
arma::mat predictS(const arma::mat &Y_, const arma::mat &W, const arma::mat &R){
  return(predictS_C(Y_, W, R));
}

// [[Rcpp::export]]
arma::mat predictC(const arma::mat &sigmaPoints, const arma::mat &Y_,
                   const arma::mat &W){
  return(predictC_C(sigmaPoints, Y_,
                    W));
}

// [[Rcpp::export]]
arma::mat computeK(const arma::mat &C, const arma::mat &S){
  return(computeK_C(C, S));
}

// [[Rcpp::export]]
arma::mat updateM(const arma::colvec &m_, const arma::mat &K,
                  const arma::colvec &residual){
  return(updateM_C(m_, K,
                   residual));
}

// [[Rcpp::export]]
arma::mat updateP(const arma::mat &P_, const arma::mat &K, const arma::mat &S){
  return(updateP_C(P_, K, S));
}

// [[Rcpp::export]]
double computeIndividualM2LL(const int &nObservedVariables,
                        const arma::colvec &rawData,
                        const arma::colvec &expectedMeans,
                        const arma::mat &expectedCovariance){
  return(computeIndividualM2LL_C(nObservedVariables,
                                 rawData,
                                 expectedMeans,
                                 expectedCovariance));
}

