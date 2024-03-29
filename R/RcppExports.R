# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#' setParameterValues
#'
#' Change the parameter values in the parameterTable
#'
#' @param parameterTable parameterTable
#' @param parameterValues values to which the parameters in parameterTable will be changed
#' @param parameterLabels labels of the parameters (must correspond to the label column of parameterTable)
#' @return the parameterTable is changed and returned
setParameterValues <- function(parameterTable, parameterValues, parameterLabels) {
    .Call(`_psydiff_setParameterValues`, parameterTable, parameterValues, parameterLabels)
}

#' getParameterValues
#'
#' Get the parameter values of a psydiffModel
#'
#' @param psydiffModel psydiffModel
#' @return named vector with parameter values
getParameterValues <- function(psydiffModel) {
    .Call(`_psydiff_getParameterValues`, psydiffModel)
}

#' setParameterList
#'
#' The parameterList holds all parameters for a single person.
#' With setParameterTable these parameters are changed to the
#' values specified in the parameterTable for a person
#'
#' @param parameterTable parameterTable
#' @param parameterList list of parameters
#' @param person integer indicating the person for which the parameters should be written in the parameterList
#' @return the parameterList is changed by reference; nothing is returned explicitly
setParameterList <- function(parameterTable, parameterList, person) {
    invisible(.Call(`_psydiff_setParameterList`, parameterTable, parameterList, person))
}

#' getSigmaPoints
#'
#' Computes the sigma-points of the unscented Kalman filter
#'
#' @param m colvec of latent means
#' @param A lower triangular Cholesky decomposition of latent covariance matrix
#' @param covarianceIsRoot boolean: is A the lower triangular Cholesky decomposition?
#' @param c double: hyperparameter c for computing the sigma points (controls the spread)
#' @return matrix with sigma points
getSigmaPoints <- function(m, A, covarianceIsRoot, c) {
    .Call(`_psydiff_getSigmaPoints`, m, A, covarianceIsRoot, c)
}

#' getMeanWeights
#'
#' Computes the mean weights of the unscented Kalman filter
#'
#' @param n number of latent variables
#' @param alpha hyperparameter
#' @param beta hyperparameter
#' @param kappa hyperparameter
#' @return matrix with mean weights
getMeanWeights <- function(n, alpha, beta, kappa) {
    .Call(`_psydiff_getMeanWeights`, n, alpha, beta, kappa)
}

#' getCovWeights
#'
#' Computes the covariance weights of the unscented Kalman filter
#'
#' @param n number of latent variables
#' @param alpha hyperparameter
#' @param beta hyperparameter
#' @param kappa hyperparameter
#' @return matrix with covariance weights
getCovWeights <- function(n, alpha, beta, kappa) {
    .Call(`_psydiff_getCovWeights`, n, alpha, beta, kappa)
}

#' getWMatrix
#'
#' Computes the weight matrix
#'
#' @param meanWeights matrix with mean weights from getMeanWeights
#' @param covWeights matrix with covariance weights from getCovWeights
#' @return matrix with weights
getWMatrix <- function(meanWeights, covWeights) {
    .Call(`_psydiff_getWMatrix`, meanWeights, covWeights)
}

#' computeMeanFromSigmaPoints
#'
#' Computes the means given the sigma points and the mean weights
#'
#' @param sigmaPoints sigma points
#' @param meanWeights matrix with mean weights from getMeanWeights
#' @return colvec with means
computeMeanFromSigmaPoints <- function(sigmaPoints, meanWeights) {
    .Call(`_psydiff_computeMeanFromSigmaPoints`, sigmaPoints, meanWeights)
}

#' getAFromSigmaPoints
#'
#' Computes the lower triangular covariance matrix given the sigma points and the covariance weights
#'
#' @param sigmaPoints sigma points
#' @param meanWeights matrix with mean weights from getMeanWeights
#' @param c double: hyperparameter c for computing the sigma points (controls the spread)
#' @return lower triangular matrix
getAFromSigmaPoints <- function(sigmaPoints, meanWeights, c) {
    .Call(`_psydiff_getAFromSigmaPoints`, sigmaPoints, meanWeights, c)
}

#' computePFromSigmaPoints
#'
#' Computes the latent covariance matrix given the sigma points and the covariance weights
#'
#' @param sigmaPoints sigma points
#' @param meanWeights matrix with mean weights from getMeanWeights
#' @param c double: hyperparameter c for computing the sigma points (controls the spread)
#' @return latent covariance matrix
computePFromSigmaPoints <- function(sigmaPoints, meanWeights, c) {
    .Call(`_psydiff_computePFromSigmaPoints`, sigmaPoints, meanWeights, c)
}

#' getPhi
#'
#' Computes the lower triangular matrix in Formula (33) of Saerkkae (2007)
#'
#' @param squareMatrix a square matrix
#' @return lower triangular matrix
getPhi <- function(squareMatrix) {
    .Call(`_psydiff_getPhi`, squareMatrix)
}

#' predictMu
#'
#' Computes the predicted manifest means
#'
#' @param Y_ predicted manifest values based on the sigma points
#' @param meanWeights matrix with mean weights from getMeanWeights
#' @return colvec with predicted manifest means
predictMu <- function(Y_, meanWeights) {
    .Call(`_psydiff_predictMu`, Y_, meanWeights)
}

#' predictS
#'
#' Computes the predicted manifest covariance
#'
#' @param Y_ predicted manifest values based on the sigma points
#' @param W matrix with covariance and mean weights computed with getWMatrix
#' @param R manifest covariance matrix
#' @return matrix with predicted manifest covariances
predictS <- function(Y_, W, R) {
    .Call(`_psydiff_predictS`, Y_, W, R)
}

#' predictC
#'
#' Computes the predicted cross-covariance between latent and manifest variables
#'
#' @param sigmaPoints sigma points
#' @param Y_ predicted manifest values based on the sigma points
#' @param W matrix with covariance and mean weights computed with getWMatrix
#' @return matrix with predicted cross-covariances
predictC <- function(sigmaPoints, Y_, W) {
    .Call(`_psydiff_predictC`, sigmaPoints, Y_, W)
}

#' computeK
#'
#' Computes the Kalman Gain
#'
#' @param C matrix with predicted cross-covariances
#' @param S matrix with predicted manifest covariances
#' @return matrix: Kalman Gain
computeK <- function(C, S) {
    .Call(`_psydiff_computeK`, C, S)
}

#' updateM
#'
#' Updates the latent means
#'
#' @param m_ colvec with predicted latent means
#' @param K Kalman Gain matrix
#' @param residual colvec with difference between predicted manifest means (mu) and observed manifest values
#' @return colvec with updated latent means
updateM <- function(m_, K, residual) {
    .Call(`_psydiff_updateM`, m_, K, residual)
}

#' updateP
#'
#' Updates the latent covariances
#'
#' @param P_ matrix with predicted latent covariances
#' @param K Kalman Gain matrix
#' @param S matrix with predicted manifest covariances
#' @return matrix with updated latent covariances
updateP <- function(P_, K, S) {
    .Call(`_psydiff_updateP`, P_, K, S)
}

#' computeIndividualM2LL
#'
#' Computes the -2 log likelihood for a single subject
#'
#' @param nObservedVariables number of non-NA observations for this person
#' @param rawData colvec with observed data for this person
#' @param expectedMeans colvec with predicted manifest means for this person
#' @param expectedCovariance matrix with predicted manifest covariances for this person
#' @return double: -2 log likelihood
computeIndividualM2LL <- function(nObservedVariables, rawData, expectedMeans, expectedCovariance) {
    .Call(`_psydiff_computeIndividualM2LL`, nObservedVariables, rawData, expectedMeans, expectedCovariance)
}

#' cholupdate
#'
#' Computes the cholupdate for the square root unscented update (see Van der Merwe and Wan, 2001)
#'
#' @param L cholesky matrix which will be updated
#' @param x colvec which will be used for the update
#' @param v double: additional factor for cholupdate, see Van der Merwe and Wan, 2001
#' @param direction string: either downdate or update
#' @return matrix: downdated / updated matrix
cholupdate <- function(L, x, v, direction) {
    .Call(`_psydiff_cholupdate`, L, x, v, direction)
}

#' qr_
#'
#' qr for the square root unscented update (see van der Merwe and Wan, 2001). Returns R-tilde
#'
#' @param X matrix which will be decomposed
#' @return matrix: R-tilde
qr_ <- function(X) {
    .Call(`_psydiff_qr_`, X)
}

#' logChol2Chol
#'
#' transforms the log-Cholesky decomposition of a matrix to the Cholesky (see Pinheiro, J. C., and Bates, D. M. (1996). Unconstrained parametrizations for variance-covariance matrices. Statistics and Computing, 6(3), 289–296)
#' @param logChol log-Cholesky decomposition of a matrix
#' @return Cholesky decomposition of the matrix
logChol2Chol <- function(logChol) {
    .Call(`_psydiff_logChol2Chol`, logChol)
}

#' clonePsydiffModel
#'
#' additional function to deep-copy a psydiffModel
#' @param model psydiff model
#' @return clone of the model
clonePsydiffModel <- function(model) {
    .Call(`_psydiff_clonePsydiffModel`, model)
}

