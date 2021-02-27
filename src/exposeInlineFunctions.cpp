// Call the package's header function
#include <psydiff.h>

//' setParameterValues
//'
//' Change the parameter values in the parameterTable
//'
//' @param parameterTable parameterTable
//' @param parameterValues values to which the parameters in parameterTable will be changed
//' @param parameterLabels labels of the parameters (must correspond to the label column of parameterTable)
//' @return the parameterTable is changed and returned
// [[Rcpp::export]]
Rcpp::DataFrame setParameterValues(Rcpp::DataFrame parameterTable,
                        Rcpp::NumericVector parameterValues,
                        Rcpp::StringVector parameterLabels){
  setParameterValues_C(parameterTable,
                       parameterValues,
                       parameterLabels);
  return(parameterTable);
}

//' getParameterValues
//'
//' Get the parameter values of a psydiffModel
//'
//' @param psydiffModel psydiffModel
//' @return named vector with parameter values
// [[Rcpp::export]]
Rcpp::NumericVector getParameterValues(Rcpp::List psydiffModel){
  return getParameterValues_C(psydiffModel);
}

//' setParameterList
//'
//' The parameterList holds all parameters for a single person.
//' With setParameterTable these parameters are changed to the
//' values specified in the parameterTable for a person
//'
//' @param parameterTable parameterTable
//' @param parameterList list of parameters
//' @param person integer indicating the person for which the parameters should be written in the parameterList
//' @return the parameterList is changed by reference; nothing is returned explicitly
// [[Rcpp::export]]
void setParameterList(const Rcpp::DataFrame &parameterTable, Rcpp::List &parameterList,
              int person){
  return setParameterList_C(parameterTable, parameterList, person);
}

//' getSigmaPoints
//'
//' Computes the sigma-points of the unscented Kalman filter
//'
//' @param m colvec of latent means
//' @param A lower triangular Cholesky decomposition of latent covariance matrix
//' @param covarianceIsRoot boolean: is A the lower triangular Cholesky decomposition?
//' @param c double: hyperparameter c for computing the sigma points (controls the spread)
//' @return matrix with sigma points
// [[Rcpp::export]]
arma::mat getSigmaPoints(arma::colvec m, arma::mat &A,
                         const bool &covarianceIsRoot,
                         const double &c){
  return (getSigmaPoints_C(m, A, covarianceIsRoot, c));
}

//' getMeanWeights
//'
//' Computes the mean weights of the unscented Kalman filter
//'
//' @param n number of latent variables
//' @param alpha hyperparameter
//' @param beta hyperparameter
//' @param kappa hyperparameter
//' @return matrix with mean weights
// [[Rcpp::export]]
arma::colvec getMeanWeights(const int &n, const double &alpha,
                            const double &beta, const double &kappa){
  return(getMeanWeights_C(n, alpha, beta, kappa));
}

//' getCovWeights
//'
//' Computes the covariance weights of the unscented Kalman filter
//'
//' @param n number of latent variables
//' @param alpha hyperparameter
//' @param beta hyperparameter
//' @param kappa hyperparameter
//' @return matrix with covariance weights
// [[Rcpp::export]]
arma::colvec getCovWeights(const int &n, const double &alpha,
                            const double &beta, const double &kappa){
  return(getCovWeights_C(n, alpha, beta, kappa));
}

//' getWMatrix
//'
//' Computes the weight matrix
//'
//' @param meanWeights matrix with mean weights from getMeanWeights
//' @param covWeights matrix with covariance weights from getCovWeights
//' @return matrix with weights
// [[Rcpp::export]]
arma::mat getWMatrix(const arma::colvec &meanWeights, const arma::colvec &covWeights){
  return(getWMatrix_C(meanWeights, covWeights));
}

//' computeMeanFromSigmaPoints
//'
//' Computes the means given the sigma points and the mean weights
//'
//' @param sigmaPoints sigma points
//' @param meanWeights matrix with mean weights from getMeanWeights
//' @return colvec with means
// [[Rcpp::export]]
arma::colvec computeMeanFromSigmaPoints(const arma::mat &sigmaPoints,
                                        const arma::colvec &meanWeights){
  return(computeMeanFromSigmaPoints_C(sigmaPoints,
                                      meanWeights));
}

//' getAFromSigmaPoints
//'
//' Computes the lower triangular covariance matrix given the sigma points and the covariance weights
//'
//' @param sigmaPoints sigma points
//' @param meanWeights matrix with mean weights from getMeanWeights
//' @param c double: hyperparameter c for computing the sigma points (controls the spread)
//' @return lower triangular matrix
// [[Rcpp::export]]
arma::mat getAFromSigmaPoints(const arma::mat &sigmaPoints,
                              const arma::colvec &meanWeights, const double &c){
  return(getAFromSigmaPoints_C(sigmaPoints,
                               meanWeights, c));
}

//' computePFromSigmaPoints
//'
//' Computes the latent covariance matrix given the sigma points and the covariance weights
//'
//' @param sigmaPoints sigma points
//' @param meanWeights matrix with mean weights from getMeanWeights
//' @param c double: hyperparameter c for computing the sigma points (controls the spread)
//' @return latent covariance matrix
// [[Rcpp::export]]
arma::mat computePFromSigmaPoints(const arma::mat &sigmaPoints,
                                  const arma::colvec &meanWeights,
                                  const double &c){
  return(computePFromSigmaPoints_C(sigmaPoints,
                                   meanWeights,
                                   c));
}

//' getPhi
//'
//' Computes the lower triangular matrix in Formula (33) of Saerkkae (2007)
//'
//' @param squareMatrix a square matrix
//' @return lower triangular matrix
// [[Rcpp::export]]
arma::mat getPhi(const arma::mat &squareMatrix){
  return(getPhi_C(squareMatrix));
}

//' predictMu
//'
//' Computes the predicted manifest means
//'
//' @param Y_ predicted manifest values based on the sigma points
//' @param meanWeights matrix with mean weights from getMeanWeights
//' @return colvec with predicted manifest means
// [[Rcpp::export]]
arma::mat predictMu(const arma::mat &Y_, const arma::colvec &meanWeights){
  return(predictMu_C(Y_, meanWeights));
}

//' predictS
//'
//' Computes the predicted manifest covariance
//'
//' @param Y_ predicted manifest values based on the sigma points
//' @param W matrix with covariance and mean weights computed with getWMatrix
//' @param R manifest covariance matrix
//' @return matrix with predicted manifest covariances
// [[Rcpp::export]]
arma::mat predictS(const arma::mat &Y_, const arma::mat &W, const arma::mat &R){
  return(predictS_C(Y_, W, R));
}

//' predictC
//'
//' Computes the predicted cross-covariance between latent and manifest variables
//'
//' @param sigmaPoints sigma points
//' @param Y_ predicted manifest values based on the sigma points
//' @param W matrix with covariance and mean weights computed with getWMatrix
//' @return matrix with predicted cross-covariances
// [[Rcpp::export]]
arma::mat predictC(const arma::mat &sigmaPoints, const arma::mat &Y_,
                   const arma::mat &W){
  return(predictC_C(sigmaPoints, Y_,
                    W));
}

//' computeK
//'
//' Computes the Kalman Gain
//'
//' @param C matrix with predicted cross-covariances
//' @param S matrix with predicted manifest covariances
//' @return matrix: Kalman Gain
// [[Rcpp::export]]
arma::mat computeK(const arma::mat &C, const arma::mat &S){
  return(computeK_C(C, S));
}

//' updateM
//'
//' Updates the latent means
//'
//' @param m_ colvec with predicted latent means
//' @param K Kalman Gain matrix
//' @param residual colvec with difference between predicted manifest means (mu) and observed manifest values
//' @return colvec with updated latent means
// [[Rcpp::export]]
arma::mat updateM(const arma::colvec &m_, const arma::mat &K,
                  const arma::colvec &residual){
  return(updateM_C(m_, K,
                   residual));
}

//' updateP
//'
//' Updates the latent covariances
//'
//' @param P_ matrix with predicted latent covariances
//' @param K Kalman Gain matrix
//' @param S matrix with predicted manifest covariances
//' @return matrix with updated latent covariances
// [[Rcpp::export]]
arma::mat updateP(const arma::mat &P_, const arma::mat &K, const arma::mat &S){
  return(updateP_C(P_, K, S));
}

//' computeIndividualM2LL
//'
//' Computes the -2 log likelihood for a single subject
//'
//' @param nObservedVariables number of non-NA observations for this person
//' @param rawData colvec with observed data for this person
//' @param expectedMeans colvec with predicted manifest means for this person
//' @param expectedCovariance matrix with predicted manifest covariances for this person
//' @return double: -2 log likelihood
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

//' cholupdate
//'
//' Computes the cholupdate for the square root unscented update (see Van der Merwe and Wan, 2001)
//'
//' @param L cholesky matrix which will be updated
//' @param x colvec which will be used for the update
//' @param v double: additional factor for cholupdate, see Van der Merwe and Wan, 2001
//' @param direction string: either downdate or update
//' @return matrix: downdated / updated matrix
// [[Rcpp::export]]
arma::mat cholupdate(arma::mat L, arma::colvec x, double v, std::string direction){
  return(cholupdate_C(L, x, v, direction));
}

//' qr_
//'
//' qr for the square root unscented update (see van der Merwe and Wan, 2001). Returns R-tilde
//'
//' @param X matrix which will be decomposed
//' @return matrix: R-tilde
// [[Rcpp::export]]
arma::mat qr_(arma::mat X){
  return(qr_C(X));
}

//' logChol2Chol
//'
//' transforms the log-Cholesky decomposition of a matrix to the Cholesky (see Pinheiro, J. C., and Bates, D. M. (1996). Unconstrained parametrizations for variance-covariance matrices. Statistics and Computing, 6(3), 289–296)
//' @param logChol log-Cholesky decomposition of a matrix
//' @return Cholesky decomposition of the matrix
// [[Rcpp::export]]
arma::mat logChol2Chol(arma::mat logChol){
  return logChol2Chol_C(logChol);
}

//' clonePsydiffModel
//'
//' additional function to deep-copy a psydiffModel
//' @param model psydiff model
//' @return clone of the model
// [[Rcpp::export]]
Rcpp::List clonePsydiffModel(Rcpp::List model){
  Rcpp::List modelClone = Rcpp::clone(model);
  return(modelClone);
}

