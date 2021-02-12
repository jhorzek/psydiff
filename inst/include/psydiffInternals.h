#ifndef PSYDIFFINTERNALS_H
#define PSYDIFFINTERNALS_H

inline arma::mat getSigmaPoints_C(arma::colvec m, arma::mat &A, const bool &covarianceIsRoot,
                                  const double &c){
  int n = m.n_elem;
  arma::mat sigmaPoints(n,2*n+1, arma::fill::zeros);
  if(!covarianceIsRoot){
    A = chol(A, "lower");
  }

  sigmaPoints.submat(0,1,n-1,n) = sqrt(c)*A;
  sigmaPoints.submat(0,n+1,n-1,2*n) = -1*sqrt(c)*A;
  sigmaPoints.each_col() += m;
  return(sigmaPoints);
}

inline arma::colvec getMeanWeights_C(int n, double alpha, double beta, double kappa){
  double lambda = pow(alpha, 2)*(n+kappa)-n;
  arma::colvec meanWeights(2*n+1);
  meanWeights.fill(1/(2*(n+lambda)));
  meanWeights(0) = lambda/(n+lambda);
  return(meanWeights);
}

inline arma::colvec getCovWeights_C(const int &n, const double &alpha,
                                  const double &beta, const double &kappa){
  double lambda = pow(alpha, 2)*(n+kappa)-n;
  arma::colvec covWeights(2*n+1);
  covWeights.fill(1/(2*(n+lambda)));
  covWeights(0) = lambda/(n+lambda + 1-pow(alpha,2) + beta);
  return(covWeights);
}

inline arma::mat getWMatrix_C(const arma::colvec &meanWeights,
                              const arma::colvec &covWeights){
  int nelem = meanWeights.n_elem;
  arma::mat W(nelem, nelem);
  arma::mat diagMinusMeanWeights = arma::eye(nelem, nelem);
  diagMinusMeanWeights.each_col() -= meanWeights;
  arma::mat element2(nelem, nelem);
  W = diagMinusMeanWeights*diagmat(covWeights)*arma:: trans(diagMinusMeanWeights);
  return(W);
}

inline arma::colvec computeMeanFromSigmaPoints_C(const arma::mat &sigmaPoints,
                                               const arma::colvec &meanWeights){
  return(sigmaPoints*meanWeights);
}

inline arma::mat getAFromSigmaPoints_C(const arma::mat &sigmaPoints,
                                     const arma::colvec &meanWeights, const double &c){
  // compute means from sigma points
  arma::colvec means = computeMeanFromSigmaPoints_C(sigmaPoints, meanWeights);
  int n = means.n_elem;
  arma::mat covarianceRoot = sigmaPoints.submat(0,1,n-1,n);
  covarianceRoot.each_col() -= means;
  covarianceRoot = covarianceRoot/sqrt(c);
  return(covarianceRoot);
}

inline arma::mat computePFromSigmaPoints_C(const arma::mat &sigmaPoints,
                                         const arma::colvec &meanWeights,
                                         const double &c){
  return (getAFromSigmaPoints_C(sigmaPoints, meanWeights, c) *
    arma::trans(getAFromSigmaPoints_C(sigmaPoints, meanWeights, c)));
}

inline arma::mat getPhi_C(const arma::mat &squareMatrix){
  // computes the lower triangular matrix in Formula (33) of Särkkä (2007)
  arma::mat phi;
  phi = trimatl(squareMatrix);
  phi.diag() = .5*squareMatrix.diag();
  return(phi);
}

inline arma::colvec predictMu_C(const arma::mat &Y_, const arma::colvec &meanWeights){
  return(Y_ * meanWeights);
}

inline arma::mat predictS_C(const arma::mat &Y_, const arma::mat &W, const arma::mat &R){
  return Y_ * W * arma::trans(Y_) + R;
}

inline arma::mat predictC_C(const arma::mat &sigmaPoints, const arma::mat &Y_,
                            const arma::mat &W){
  return sigmaPoints * W * arma::trans(Y_);
}

inline arma::mat computeK_C(const arma::mat &C, const arma::mat &S){
  return C * arma::pinv(S);
}

inline arma::mat updateM_C(const arma::colvec &m_, const arma::mat &K,
                           const arma::colvec &residual){
  return m_ + K*residual;
}

inline arma::mat updateP_C(const arma::mat &P_, const arma::mat &K, const arma::mat &S){
  return P_ - K*S*arma::trans(K);
}

inline double computeIndividualM2LL_C(const int &nObservedVariables,
                                      const arma::colvec &rawData,
                                      const arma::colvec &expectedMeans,
                                      const arma::mat &expectedCovariance){
  double m2LL;
  double klog2pi = nObservedVariables*std::log(2*M_PI);
  double logDetExpCov = std::log(arma::det(expectedCovariance));
  arma::mat dist = arma::trans(rawData - expectedMeans)*arma::pinv(expectedCovariance)*(rawData - expectedMeans);
  m2LL = klog2pi +
    logDetExpCov +
    dist(0,0); // note: dist is a 1x1 matrix; extraction is necessary for the data type to be compatible
  return(m2LL);
}

inline double computeIndividualM2LLChol_C(const int &nObservedVariables,
                                      const arma::colvec &rawData,
                                      const arma::colvec &expectedMeans,
                                      const arma::mat &expectedCovarianceChol){
  double m2LL;
  double klog2pi = nObservedVariables*std::log(2*M_PI);
  double logDetExpCov = 2*sum(arma::log(expectedCovarianceChol.diag()));
  arma::colvec residual = rawData - expectedMeans;
  arma::mat choInv = arma::pinv(arma::trimatl(expectedCovarianceChol));
  arma::mat dist = arma::trans(residual)*arma::trans(choInv)*choInv*(residual);
  m2LL = klog2pi +
    logDetExpCov +
    dist(0,0); // note: dist is a 1x1 matrix; extraction is necessary for the data type to be compatible
  return(m2LL);
}

inline arma::mat cholupdate_C(arma::mat L, arma::colvec x, double v, std::string direction){
  // Computes the cholupdate for the square root unscented update (see van der Merwe & Wan, 2001)
  // follows the code example found here: https://en.wikipedia.org/wiki/Cholesky_decomposition
  // and the adaption for the cholupdate in Van der Merwe and Wan (2001) found here:
  // https://github.com/piotrek-orlowski/ukfRcpp/blob/master/src/cholupdate.cpp
  int n = x.n_elem, dir;
  double r = 0, c, s;
  if(direction == "update"){
    dir = 1;
  }else if(direction == "downdate"){
    dir = -1;
  }else{
    Rcpp::stop("Unknown direction in cholupdate. Possible are downdate and update.");
  }
  if(v<0){
    Rcpp::warning("v in cholupdate is negative. Using absolute value");
    v = -1*v;
  }
  v = pow(v,.25);
  x = v*x;

  for(int k = 0; k < n; k++){
    r = sqrt(pow(L(k,k),2.0) + dir*pow(x(k),2.0));
    c = r/L(k,k);
    s = x(k)/L(k,k);
    L(k,k) = r;
    if(k < n-1){
      L(arma::span(k+1, n-1), k ) = (L(arma::span(k+1, n-1), k ) + dir*s*x(arma::span(k+1, n-1)))/c;
      x(arma::span(k+1, n-1)) = c*x(arma::span(k+1, n-1)) - s*L(arma::span(k+1, n-1), k );
    }
  }
  return(L);
}

inline arma::mat qr_C(arma::mat X){
  arma::mat Q,R;
  arma:: qr_econ(Q, R, X);
  return(R);
}


inline arma::mat predictSquareRootOfS_C(arma::mat Y_, arma::colvec mu, arma::mat Rchol,
                                        double covWeight_0, double covWeight_1){
  arma::mat srS;
  double abscovWeight_0;
  std::string direction;
  arma::mat tempMat(Y_.n_rows, Y_.n_cols + Rchol.n_cols, arma::fill::zeros);

  for(int co = 0; co < Y_.n_cols; co++){
    tempMat.col(co) = sqrt(covWeight_1)*(Y_.col(co) - mu);
  }
  tempMat.cols(arma::span(Y_.n_cols, Y_.n_cols+Rchol.n_cols-1)) = Rchol;

  srS = arma::trans(qr_C(arma::trans(tempMat)));
  if(covWeight_0 < 0){
    direction = "downdate";
    abscovWeight_0 = -1*covWeight_0;
  }else{
    direction = "update";
    abscovWeight_0 = covWeight_0;
  }

  srS = cholupdate_C(srS, Y_.col(0) - mu, abscovWeight_0, direction);

  return(arma::trimatl(srS));
}


inline arma::mat computeK_SR_C(const arma::mat &C, const arma::mat &srS){
  arma::mat srSInv = arma::pinv(arma::trimatl(srS));
  return C * arma::trans(srSInv)*srSInv;
}


#endif
