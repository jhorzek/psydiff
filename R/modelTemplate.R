#' modelTemplate
#'
#' @param srUpdate boolean: Should the square root version be used for the updates?
#' @return String with C++ code for the actual model. LATENTEQUATIONPLACEHOLDER has to be replaced with the latent equation and MEASUREMENTEQUATIONPLACEHOLDER with the manifest equation.
modelTemplate <- function(srUpdate = TRUE){
  if(srUpdate){
    mod <- '// [[Rcpp::depends(RcppArmadillo, psydiff)]]

#include <psydiff.h>
#include <RcppArmadillo.h>

/* The type of container used to hold the state vector */
  typedef arma::mat state_type;

// By default, odeint will not work with armadillo
// The following code is copied from headmyshouder
// and allows the type of the states to be arma::mat
// https://stackoverflow.com/questions/41753823/armadillo-conflicts-with-boost-odeint-odeint-resizes-the-state-vector-to-zero-d
namespace boost { namespace numeric { namespace odeint {
  // template for arma::mat - define as resizable
  template <>
    struct is_resizeable<arma::mat>
    {
      typedef boost::true_type type;
      const static bool value = type::value;
    };

  // define size comparisons for two arma::mat
  template <>
    struct same_size_impl<arma::mat, arma::mat>
    {
      static bool same_size(const arma::mat& X, const arma::mat& Y)
      {
        return (X.n_cols == Y.n_cols && X.n_rows == Y.n_rows);
      }
    };

  // define resize command for arma::mat
  template<>
    struct resize_impl<arma::mat, arma::mat>
    {
      static void resize(arma::mat& X, const arma::mat& Y)
      {
        X.reshape(Y.n_rows, Y.n_cols);
      }
    };

} } } // namespace boost::numeric::odeint

// Next, we set up the structs in which the parameters
// and internal information for odeint are saved

struct odeintpars{
  Rcpp::List parameterList, additional;
  arma::mat meanWeights, covWeights, W;
  int nlatent, nmanifest, person;
  double c, alpha, beta, kappa;
};

MEASUREMENTEQUATIONPLACEHOLDER

arma::mat getY_(const arma::mat sigmaPoints, const odeintpars &pars,
                const int &person,
                const double &t){
  arma::mat Y_(pars.nmanifest, sigmaPoints.n_cols, arma::fill::zeros);
  arma::colvec currentSigma;
  for(int co = 0; co < sigmaPoints.n_cols; co++){
    currentSigma = sigmaPoints.col(co);
    Y_.col(co) = measurementEquation(currentSigma, pars,
                                     person, t);
  }
  return Y_;
}

LATENTEQUATIONPLACEHOLDER

arma::mat computeL(const odeintpars &pars,
                   const int &person,
                   const double &t){
  Rcpp::List modelPars = pars.parameterList;
  arma::mat L = modelPars["L"];
  return(L);
}

arma::mat getMMatrix(const arma::mat &sigmaPoints, const odeintpars &pars,
                     const int &person,
                     const double &t){
  // Note: This function has to be compiled later on as the users are supposed to
  // provide their own models. The current implementation is a placeholder.
  // sigmaPoints are called x in Sarkka (2007)

  arma::mat M, Qc, invCov;
  arma::mat sigmaPredict(sigmaPoints.n_rows, sigmaPoints.n_cols, arma::fill::zeros);
  arma::mat A =  getAFromSigmaPoints_C(sigmaPoints,
                                       pars.meanWeights,
                                       pars.c);
  invCov = arma::inv(arma::trimatl(A));


  // Qc
  // Qc is the diffusion matrix of the Brownian Motion.
  // In our case this corresponds to the identity matrix (standard Brownian motion).
  // Modulation of the Brownian Motion is "outsourced" to Matrix L
  Qc = arma::eye(pars.nlatent, pars.nlatent);

  for(int co = 0; co < sigmaPoints.n_cols; co++){
    sigmaPredict.col(co) = latentEquation(sigmaPoints.col(co), pars,
                                          person, t);
  }

  arma::colvec m_x(pars.nlatent, 1, arma::fill::zeros);
  arma::colvec m_f(pars.nlatent, 1, arma::fill::zeros);
  for(int co = 0; co < sigmaPoints.n_cols; co++){
    m_x += pars.meanWeights(co) * sigmaPoints.col(co);
    m_f += pars.meanWeights(co) * sigmaPredict.col(co);
  }

  // evaluate summations
  arma::mat sum1(pars.nlatent, pars.nlatent, arma::fill::zeros);
  arma::mat sum2(pars.nlatent, pars.nlatent, arma::fill::zeros);
  for(int co = 0; co < sigmaPoints.n_cols; co++){
    sum1 += pars.covWeights(co) * (sigmaPoints.col(co) - m_x) * arma::trans(sigmaPredict.col(co) - m_f);
    sum2 += pars.covWeights(co) * (sigmaPredict.col(co) - m_f)*arma::trans(sigmaPoints.col(co) - m_x);
  }

  arma::mat L = computeL(pars, person, t);

  // compute M
  M = invCov*(
    sum1+
      sum2+
      L*Qc*arma::trans(L)
  )*arma::trans(invCov);
  return(M);
}

// Model for integration
class odeintModel{
  struct odeintpars modelParameters;
  public:
    odeintModel(struct odeintpars newParameters) : modelParameters(newParameters){}
  void operator()(const state_type &x , state_type &dxdt , const double  t ){
    arma::mat A = getAFromSigmaPoints_C(x, modelParameters.meanWeights, modelParameters.c);
    // compute M matrix
    arma::mat M = getMMatrix(x, modelParameters, modelParameters.person, t);

    arma::mat phi_M = getPhi_C(M);

    arma::mat augmentedMatrix(x.n_rows, x.n_cols, arma::fill::zeros);
    augmentedMatrix.submat(0,1,modelParameters.nlatent-1,modelParameters.nlatent) = A*phi_M;
    augmentedMatrix.submat(0,modelParameters.nlatent+1,modelParameters.nlatent-1,2*modelParameters.nlatent) = -1*A*phi_M;

    arma::mat sigmaPredict(x.n_rows, x.n_cols, arma::fill::zeros);

    for(int co = 0; co < x.n_cols; co++){
      sigmaPredict.col(co) = latentEquation(x.col(co), modelParameters, modelParameters.person, t);
    }

    for(int co = 0; co < dxdt.n_cols; co++){
      dxdt.col(co) = sigmaPredict * modelParameters.meanWeights + sqrt(modelParameters.c)*augmentedMatrix.col(co);
    }
  }
};

// [[Rcpp::export]]
Rcpp::List fitModel(Rcpp::List psydiffModel, double alpha = 0.6, double beta = 2.0, double kappa = -1.0,
                    double timeStep = 0.01, std::string integrateFunction = "default", bool breakEarly = true, int verbose = 0){

  // extract parameters from model
  Rcpp::List pars = psydiffModel["pars"];
  Rcpp::List parameterList = pars["parameterList"]; // has all starting values
  Rcpp::DataFrame parameterTable = Rcpp::as<Rcpp::DataFrame>(pars["parameterTable"]);
  Rcpp::NumericVector personInParameterTable = parameterTable["person"]; // persons
  Rcpp::List additional = psydiffModel["additional"];

  odeintpars odeintparam;
  odeintparam.parameterList = parameterList;
  odeintparam.nlatent = psydiffModel["nlatent"];
  odeintparam.nmanifest = psydiffModel["nmanifest"];
  odeintparam.alpha = alpha;
  odeintparam.beta = beta;
  odeintparam.kappa = kappa;
  if(odeintparam.nlatent + odeintparam.kappa == 0){
    Rcpp::warning("Setting of kappa results in division by 0. Kappa was changed to kappa - 1");
    kappa -= 1;
    odeintparam.kappa = kappa;
  }
  double c = pow(odeintparam.alpha,2)*(odeintparam.nlatent + odeintparam.kappa);
  odeintparam.c = c;

  // data
  Rcpp::List data = psydiffModel["data"];
  Rcpp::NumericVector personsInData = data["person"]; // persons
  Rcpp::NumericVector uniquePersons = unique(personsInData);
  int sampleSize = uniquePersons.length();
  arma::mat observations = data["observations"]; // observations
  arma::colvec dt = data["dt"]; // dt

  // changed
  Rcpp::LogicalVector changed = parameterTable["changed"];

  // initialize Unscented matrices
  int numsteps, nObservedVariables, timePoints;
  arma::mat individualObservations, // data for one individual
  Y_(odeintparam.nmanifest, 2*odeintparam.nlatent+1), // predicted observations from sigma points
  P(odeintparam.nlatent, odeintparam.nlatent), // updated latent covariance
  P_(odeintparam.nlatent, odeintparam.nlatent), // latent covariance
  A(odeintparam.nlatent, odeintparam.nlatent), // root of (updated) latent covariance
  S(odeintparam.nmanifest,odeintparam.nmanifest), // predicted observed covariance
  C(odeintparam.nlatent, odeintparam.nmanifest), // predicted cross-covariance
  K(odeintparam.nlatent, odeintparam.nmanifest), // Kalman Gain
  latentScores(observations.n_rows, odeintparam.nlatent, arma::fill::zeros), // collects predicted latent scores
  predictedManifest(observations.n_rows, odeintparam.nmanifest, arma::fill::zeros); // collects predicted observations
  //individualLatentScores,
  //individualPredictedManifest;
  arma::colvec individualDts, // person specific dt
  m_(odeintparam.nlatent), // predicted latent means
  m(odeintparam.nlatent), // updated latent means
  mu(odeintparam.nmanifest), // predicted manifest means
  individualXTimeObservations, // for the observed data of a person at a specific time point
  observedNotNA, // will be used to store observations without missings
  residual, // difference between observed and predicted
  m2LL(sampleSize, arma::fill::zeros); // individual likelihoods
  arma::uvec nonmissing, // vector for indices of nonmissing data
  missing; // vector for indices of missing data

  // define type of state
  state_type x;

  // initialize start and end time of integration
  double startTime, endTime;

  // iterate over all persons
  for(int person = 0; person < sampleSize; person++){
    int selectedPerson = uniquePersons(person);
    odeintparam.person = selectedPerson;
    // check for this person if any parameters changed
    Rcpp::LogicalVector parschanged = changed[personInParameterTable == selectedPerson];

    if(is_false(any(parschanged))){
      // if nothing changed, we can skip this person
      if(verbose == 1){Rcpp::Rcout << "Skipping person " << selectedPerson << std::endl;}
      continue;
    }
    // if something changed, set parameters
    // Step 1: Set parameters
    // (will pass by reference and change the parameters directly in the parameterList)
    setParameterList_C(parameterTable, parameterList, selectedPerson);

    // extract initial parameters
    m = Rcpp::as<arma::colvec>(parameterList["m0"]);
    P = Rcpp::as<arma::mat>(parameterList["A0"])*arma::trans(Rcpp::as<arma::mat>(parameterList["A0"]));
    A = chol(P, "lower");
    arma::mat Rchol = odeintparam.parameterList["Rchol"];
    arma::mat R = Rchol*arma::trans(Rchol); // manifest covariance

    // extract individual observations
    arma::uvec rowInd = arma::find(Rcpp::as<arma::rowvec>(personsInData) == selectedPerson);
    individualObservations = observations.rows(rowInd);
    individualDts = dt.rows(rowInd);

    // iterate over all time points
    double timeSum = 0.0;
    timePoints = individualDts.n_elem;
    for(int timePoint = 0; timePoint < timePoints; timePoint++){
      if(verbose == 1){Rcpp::Rcout << "Fitting Person" <<  selectedPerson << " time " << timePoint << std::endl;}
      timeSum += individualDts(timePoint);

      // extract individual- and time-specific data as colvec
      individualXTimeObservations = arma::trans(individualObservations.row(timePoint));

      // find non-missing
      nonmissing = arma::find_finite(individualXTimeObservations);
      missing = arma::find_nonfinite(individualXTimeObservations);
      nObservedVariables = nonmissing.size();

      observedNotNA = individualXTimeObservations(nonmissing);

      // set latent mean
      m_ = m;

      // PREDICTION

      // Compute sigma-points
      x = getSigmaPoints_C(m_, A, true, odeintparam.c);

      // Compute mean weights, cov weights and weight matrix M
      arma::colvec meanWeights = getMeanWeights_C(odeintparam.nlatent, odeintparam.alpha, odeintparam.beta, odeintparam.kappa);
      odeintparam.meanWeights = meanWeights;
      arma::colvec covWeights = getCovWeights_C(odeintparam.nlatent, odeintparam.alpha, odeintparam.beta, odeintparam.kappa);
      odeintparam.covWeights = covWeights;
      arma::mat W = getWMatrix_C(meanWeights, covWeights);
      odeintparam.W = W;

      // Set up model
      odeintModel individualOdeintModel(odeintparam);

      // Integrate
      startTime = 0.0;
      endTime = individualDts(timePoint);

      if(abs(startTime - endTime) > 0){
        double currentTimeStep = timeStep;
        for(int integrateLoop = 0; integrateLoop < 20; integrateLoop++){
          if(integrateFunction == "rk4"){
            boost::numeric::odeint::runge_kutta4< state_type > stepper;
            numsteps = boost::numeric::odeint::integrate_const(stepper, individualOdeintModel, x,
                                                               startTime, endTime,
                                                               currentTimeStep);//, push_back_state_and_time( x_vec , times ));
          }else{
            numsteps = boost::numeric::odeint::integrate(individualOdeintModel, x,
                                                         startTime, endTime,
                                                         currentTimeStep);
          }
          if(arma::is_finite(x)){break;}
          currentTimeStep += timeStep*.1;
        }
        if(!arma::is_finite(x) &&  verbose > 0){Rcpp::warning("Non-finite value in integration.");}
        m_ = x.col(0);
        A = arma::trimatl(getAFromSigmaPoints_C(x, meanWeights, c));
      }

      // MANIFEST PREDICTION
      Y_ = getY_(x, odeintparam,
                 selectedPerson,
                 timeSum);

      // expected manifest mean
      mu = predictMu_C(Y_, meanWeights);

      // save prediction
      predictedManifest.row(min(rowInd) + timePoint) = arma::trans(mu);

      if(nObservedVariables > 0){
        // root of expected manifest covariance
        arma::mat srS = predictSquareRootOfS_C(Y_.rows(nonmissing), mu(nonmissing), Rchol.submat(nonmissing, nonmissing),
                                               covWeights(0), covWeights(1));

        // expected manifest cross-covariance
        C = predictC_C(x, Y_.rows(nonmissing), odeintparam.W);

        K = computeK_SR_C(C, srS);

        // compute residuals
        residual = individualXTimeObservations - mu;

        // UPDATE
        // Update latent mean
        m = updateM_C(m_, K, residual(nonmissing));
        latentScores.row(min(rowInd) + timePoint) = arma::trans(m);

        // Update latent covariance root
        arma::mat U = K*srS;

        for(int co = 0; co < U.n_cols; co++){
          A = arma::trimatl(cholupdate_C(A, U.col(co), 1, "downdate"));
        }

        // compute Likelihood --nonmissing

        arma::uvec m2LLInd = arma::find(Rcpp::as<arma::rowvec>(uniquePersons) == selectedPerson);
        m2LL.elem(m2LLInd) += computeIndividualM2LLChol_C(nObservedVariables,
                                                      individualXTimeObservations(nonmissing),
                                                      mu(nonmissing), srS);

       if(breakEarly && !arma::is_finite(m2LL.elem(m2LLInd))){
          if(verbose > 0){
            Rcpp::warning("Non-finite m2LL");
          }
          Rcpp::List ret  = Rcpp::List::create(Rcpp::Named("latentScores") = latentScores,
                                               Rcpp::Named("predictedManifest") = predictedManifest,
                                               Rcpp::Named("m2LL") = m2LL);
          return(ret);
        }

      }else{
        // if all data is missing / no update is requested
        m = m_;
        latentScores.row(min(rowInd) + timePoint) = arma::trans(m);
      }
    }
  }

  Rcpp::List ret  = Rcpp::List::create(Rcpp::Named("latentScores") = latentScores,
                                       Rcpp::Named("predictedManifest") = predictedManifest,
                                       Rcpp::Named("m2LL") = m2LL);
  return(ret);
}

// [[Rcpp::export]]
Rcpp::NumericVector getGradients(Rcpp::List psydiffModel, double eps, double alpha = 0.6, double beta = 2.0, double kappa = -1.0,
                                 double timeStep = 0.01, std::string integrateFunction = "default", bool breakEarly = true, int verbose = 0){
  Rcpp::List psydiffModelClone = Rcpp::clone(psydiffModel);
  Rcpp::List pars = psydiffModelClone["pars"];
  Rcpp::DataFrame parameterTable = Rcpp::as<Rcpp::DataFrame>(pars["parameterTable"]);
  Rcpp::NumericVector currentParameterValues = getParameterValues_C(psydiffModelClone);
  Rcpp::NumericMatrix m2LLs(currentParameterValues.length() , 2);
  Rcpp::NumericVector gradients(currentParameterValues.length());
  arma::colvec m2LL;
  Rcpp::List fittedModel;

  for(int par = 0; par < currentParameterValues.length(); par++){

    currentParameterValues(par) += eps;
    setParameterValues_C(parameterTable, currentParameterValues, currentParameterValues.names());
    fittedModel = fitModel(psydiffModelClone, alpha, beta, kappa, timeStep, integrateFunction, breakEarly, verbose);
    m2LL = Rcpp::as<arma::colvec>(fittedModel["m2LL"]);
    m2LLs(par,0) = sum(m2LL);

    currentParameterValues(par) -= 2*eps;
    setParameterValues_C(parameterTable, currentParameterValues, currentParameterValues.names());
    fittedModel = fitModel(psydiffModelClone, alpha, beta, kappa, timeStep, integrateFunction, breakEarly, verbose);
    m2LL = Rcpp::as<arma::colvec>(fittedModel["m2LL"]);
    m2LLs(par,1) = sum(m2LL);

    currentParameterValues(par) += eps;

  }

  gradients = (m2LLs(Rcpp::_,0) - m2LLs(Rcpp::_,1))/(2*eps);
  gradients.names() = currentParameterValues.names();
  return(gradients);
}

'
  }else{
    mod <- '
  // [[Rcpp::depends(RcppArmadillo, psydiff)]]

#include <psydiff.h>
#include <RcppArmadillo.h>

/* The type of container used to hold the state vector */
typedef arma::mat state_type;

// By default, odeint will not work with armadillo
// The following code is copied from headmyshouder
// and allows the type of the states to be arma::mat
// https://stackoverflow.com/questions/41753823/armadillo-conflicts-with-boost-odeint-odeint-resizes-the-state-vector-to-zero-d
namespace boost { namespace numeric { namespace odeint {
// template for arma::mat - define as resizable
template <>
struct is_resizeable<arma::mat>
{
  typedef boost::true_type type;
  const static bool value = type::value;
};

// define size comparisons for two arma::mat
template <>
struct same_size_impl<arma::mat, arma::mat>
{
  static bool same_size(const arma::mat& X, const arma::mat& Y)
  {
    return (X.n_cols == Y.n_cols && X.n_rows == Y.n_rows);
  }
};

// define resize command for arma::mat
template<>
struct resize_impl<arma::mat, arma::mat>
{
  static void resize(arma::mat& X, const arma::mat& Y)
  {
    X.reshape(Y.n_rows, Y.n_cols);
  }
};

} } } // namespace boost::numeric::odeint

// Next, we set up the structs in which the parameters
// and internal information for odeint are saved

struct odeintpars{
  Rcpp::List parameterList, additional;
  arma::mat meanWeights, covWeights, W;
  int nlatent, nmanifest, person;
  double c, alpha, beta, kappa;
};

MEASUREMENTEQUATIONPLACEHOLDER

arma::mat getY_(const arma::mat sigmaPoints, const odeintpars &pars,
                const int &person,
                const double &t){
  arma::mat Y_(pars.nmanifest, sigmaPoints.n_cols, arma::fill::zeros);
  arma::colvec currentSigma;
  for(int co = 0; co < sigmaPoints.n_cols; co++){
    currentSigma = sigmaPoints.col(co);
    Y_.col(co) = measurementEquation(currentSigma, pars,
           person, t);
  }
  return Y_;
}

LATENTEQUATIONPLACEHOLDER

arma::mat computeL(const odeintpars &pars,
                   const int &person,
                   const double &t){
  Rcpp::List modelPars = pars.parameterList;
  arma::mat L = modelPars["L"];
  return(L);
}

arma::mat getMMatrix(const arma::mat &sigmaPoints, const odeintpars &pars,
                     const int &person,
                     const double &t){
  // Note: This function has to be compiled later on as the users are supposed to
  // provide their own models. The current implementation is a placeholder.
  // sigmaPoints are called x in Sarkka (2007)

  arma::mat M, Qc, invCov;
  arma::mat sigmaPredict(sigmaPoints.n_rows, sigmaPoints.n_cols, arma::fill::zeros);
  arma::mat A =  getAFromSigmaPoints_C(sigmaPoints,
                                       pars.meanWeights,
                                       pars.c);
  invCov = arma::inv(arma::trimatl(A));


  // Qc
  // Qc is the diffusion matrix of the Brownian Motion.
  // In our case this corresponds to the identity matrix (standard Brownian motion).
  // Modulation of the Brownian Motion is "outsourced" to Matrix L
  Qc = arma::eye(pars.nlatent, pars.nlatent);

  for(int co = 0; co < sigmaPoints.n_cols; co++){
    sigmaPredict.col(co) = latentEquation(sigmaPoints.col(co), pars,
                     person, t);
  }

  arma::colvec m_x(pars.nlatent, 1, arma::fill::zeros);
  arma::colvec m_f(pars.nlatent, 1, arma::fill::zeros);
  for(int co = 0; co < sigmaPoints.n_cols; co++){
    m_x += pars.meanWeights(co) * sigmaPoints.col(co);
    m_f += pars.meanWeights(co) * sigmaPredict.col(co);
  }

  // evaluate summations
  arma::mat sum1(pars.nlatent, pars.nlatent, arma::fill::zeros);
  arma::mat sum2(pars.nlatent, pars.nlatent, arma::fill::zeros);
  for(int co = 0; co < sigmaPoints.n_cols; co++){
    sum1 += pars.covWeights(co) * (sigmaPoints.col(co) - m_x) * arma::trans(sigmaPredict.col(co) - m_f);
    sum2 += pars.covWeights(co) * (sigmaPredict.col(co) - m_f)*arma::trans(sigmaPoints.col(co) - m_x);
  }

  arma::mat L = computeL(pars, person, t);

  // compute M
  M = invCov*(
    sum1+
    sum2+
    L*Qc*arma::trans(L)
  )*arma::trans(invCov);
  return(M);
}

// Model for integration
class odeintModel{
  struct odeintpars modelParameters;
public:
  odeintModel(struct odeintpars newParameters) : modelParameters(newParameters){}
  void operator()(const state_type &x , state_type &dxdt , const double  t ){
    arma::mat A = getAFromSigmaPoints_C(x, modelParameters.meanWeights, modelParameters.c);
    // compute M matrix
    arma::mat M = getMMatrix(x, modelParameters, modelParameters.person, t);

    arma::mat phi_M = getPhi_C(M);

    arma::mat augmentedMatrix(x.n_rows, x.n_cols, arma::fill::zeros);
    augmentedMatrix.submat(0,1,modelParameters.nlatent-1,modelParameters.nlatent) = A*phi_M;
    augmentedMatrix.submat(0,modelParameters.nlatent+1,modelParameters.nlatent-1,2*modelParameters.nlatent) = -1*A*phi_M;

    arma::mat sigmaPredict(x.n_rows, x.n_cols, arma::fill::zeros);

    for(int co = 0; co < x.n_cols; co++){
      sigmaPredict.col(co) = latentEquation(x.col(co), modelParameters, modelParameters.person, t);
    }

    for(int co = 0; co < dxdt.n_cols; co++){
      dxdt.col(co) = sigmaPredict * modelParameters.meanWeights + sqrt(modelParameters.c)*augmentedMatrix.col(co);
    }
  }
};

// [[Rcpp::export]]
Rcpp::List fitModel(Rcpp::List psydiffModel, double alpha = 0.6, double beta = 2.0, double kappa = -1.0,
                    double timeStep = 0.01, std::string integrateFunction = "default", bool breakEarly = true, int verbose = 0){

  // extract parameters from model
  Rcpp::List pars = psydiffModel["pars"];
  Rcpp::List parameterList = pars["parameterList"]; // has all starting values
  Rcpp::DataFrame parameterTable = Rcpp::as<Rcpp::DataFrame>(pars["parameterTable"]);
  Rcpp::NumericVector personInParameterTable = parameterTable["person"]; // persons
  Rcpp::List additional = psydiffModel["additional"];

  odeintpars odeintparam;
  odeintparam.parameterList = parameterList;
  odeintparam.nlatent = psydiffModel["nlatent"];
  odeintparam.nmanifest = psydiffModel["nmanifest"];
  odeintparam.alpha = alpha;
  odeintparam.beta = beta;
  odeintparam.kappa = kappa;
  if(odeintparam.nlatent + odeintparam.kappa == 0){
    Rcpp::warning("Setting of kappa results in division by 0. Kappa was changed to kappa - 1");
    kappa -= 1;
    odeintparam.kappa = kappa;
  }
  double c = pow(odeintparam.alpha,2)*(odeintparam.nlatent + odeintparam.kappa);
  odeintparam.c = c;

  // data
  Rcpp::List data = psydiffModel["data"];
  Rcpp::NumericVector personsInData = data["person"]; // persons
  Rcpp::NumericVector uniquePersons = unique(personsInData);
  int sampleSize = uniquePersons.length();
  arma::mat observations = data["observations"]; // observations
  arma::colvec dt = data["dt"]; // dt

  // changed
  Rcpp::LogicalVector changed = parameterTable["changed"];

  // initialize Unscented matrices
  int numsteps, nObservedVariables, timePoints;
  arma::mat individualObservations, // data for one individual
  Y_(odeintparam.nmanifest, 2*odeintparam.nlatent+1), // predicted observations from sigma points
  P(odeintparam.nlatent, odeintparam.nlatent), // updated latent covariance
  P_(odeintparam.nlatent, odeintparam.nlatent), // latent covariance
  A(odeintparam.nlatent, odeintparam.nlatent), // root of (updated) latent covariance
  S(odeintparam.nmanifest,odeintparam.nmanifest), // predicted observed covariance
  C(odeintparam.nlatent, odeintparam.nmanifest), // predicted cross-covariance
  K(odeintparam.nlatent, odeintparam.nmanifest), // Kalman Gain
  latentScores(observations.n_rows, odeintparam.nlatent, arma::fill::zeros), // collects predicted latent scores
  predictedManifest(observations.n_rows, odeintparam.nmanifest, arma::fill::zeros); // collects predicted observations
  //individualLatentScores,
  //individualPredictedManifest;
  arma::colvec individualDts, // person specific dt
  m_(odeintparam.nlatent), // predicted latent means
  m(odeintparam.nlatent), // updated latent means
  mu(odeintparam.nmanifest), // predicted manifest means
  individualXTimeObservations, // for the observed data of a person at a specific time point
  observedNotNA, // will be used to store observations without missings
  residual, // difference between observed and predicted
  m2LL(sampleSize, arma::fill::zeros); // individual likelihoods
  arma::uvec nonmissing, // vector for indices of nonmissing data
  missing; // vector for indices of missing data

  // define type of state
  state_type x;

  // initialize start and end time of integration
  double startTime, endTime;

  // iterate over all persons
  for(int person = 0; person < sampleSize; person++){
    int selectedPerson = uniquePersons(person);
    odeintparam.person = selectedPerson;
    // check for this person if any parameters changed
    Rcpp::LogicalVector parschanged = changed[personInParameterTable == selectedPerson];

    if(is_false(any(parschanged))){
      // if nothing changed, we can skip this person
      if(verbose == 1){Rcpp::Rcout << "Skipping person " << selectedPerson << std::endl;}
      continue;
    }
    // if something changed, set parameters
    // Step 1: Set parameters
    // (will pass by reference and change the parameters directly in the parameterList)
    setParameterList_C(parameterTable, parameterList, selectedPerson);

    // extract initial parameters
    m = Rcpp::as<arma::colvec>(parameterList["m0"]);
    P = Rcpp::as<arma::mat>(parameterList["A0"])*arma::trans(Rcpp::as<arma::mat>(parameterList["A0"]));
    arma::mat Rchol = odeintparam.parameterList["Rchol"];
    arma::mat R = Rchol*arma::trans(Rchol); // manifest covariance

    // extract individual observations
    arma::uvec rowInd = arma::find(Rcpp::as<arma::rowvec>(personsInData) == selectedPerson);
    individualObservations = observations.rows(rowInd);
    individualDts = dt.rows(rowInd);
    //individualPredictedManifest = predictedManifest.rows(rowInd);
    //individualLatentScores = latentScores.rows(rowInd);

    // iterate over all time points
    double timeSum = 0.0;
    timePoints = individualDts.n_elem;
    for(int timePoint = 0; timePoint < timePoints; timePoint++){
      if(verbose == 1){Rcpp::Rcout << "Fitting Person" <<  selectedPerson << " time " << timePoint << std::endl;}
      timeSum += individualDts(timePoint);

      // extract individual- and time-specific data as colvec
      individualXTimeObservations = arma::trans(individualObservations.row(timePoint));

      // find non-missing
      nonmissing = arma::find_finite(individualXTimeObservations);
      missing = arma::find_nonfinite(individualXTimeObservations);
      nObservedVariables = nonmissing.size();

      observedNotNA = individualXTimeObservations(nonmissing);

      // set latent mean and covariance as well as root of covariance
      m_ = m;
      P_ = P;
      A = chol(P_, "lower");

      // PREDICTION

      // Compute sigma-points
      x = getSigmaPoints_C(m_, A, true, odeintparam.c);

      // Compute mean weights, cov weights and weight matrix M
      arma::colvec meanWeights = getMeanWeights_C(odeintparam.nlatent, odeintparam.alpha, odeintparam.beta, odeintparam.kappa);
      odeintparam.meanWeights = meanWeights;
      arma::colvec covWeights = getCovWeights_C(odeintparam.nlatent, odeintparam.alpha, odeintparam.beta, odeintparam.kappa);
      odeintparam.covWeights = covWeights;
      arma::mat W = getWMatrix_C(meanWeights, covWeights);
      odeintparam.W = W;

      // Set up model
      odeintModel individualOdeintModel(odeintparam);

      // Integrate
      startTime = 0.0;
      endTime = individualDts(timePoint);

      if(abs(startTime - endTime) > 0){
double currentTimeStep = timeStep;
        for(int integrateLoop = 0; integrateLoop < 20; integrateLoop++){
          if(integrateFunction == "rk4"){
            boost::numeric::odeint::runge_kutta4< state_type > stepper;
            numsteps = boost::numeric::odeint::integrate_const(stepper, individualOdeintModel, x,
                                                               startTime, endTime,
                                                               currentTimeStep);//, push_back_state_and_time( x_vec , times ));
          }else{
            numsteps = boost::numeric::odeint::integrate(individualOdeintModel, x,
                                                         startTime, endTime,
                                                         currentTimeStep);
          }
          if(arma::is_finite(x)){break;}
          currentTimeStep += timeStep*.1;
        }
        if(!arma::is_finite(x) &&  verbose > 0){Rcpp::warning("Non-finite value in integration.");}
        m_ = x.col(0);
        P_ = computePFromSigmaPoints_C(x, meanWeights, c);
      }

      // MANIFEST PREDICTION
      Y_ = getY_(x, odeintparam,
                 selectedPerson,
                 timeSum);

      // expected manifest mean
      mu = predictMu_C(Y_, meanWeights);

      // save prediction
      predictedManifest.row(min(rowInd) + timePoint) = arma::trans(mu);

      if(nObservedVariables > 0){
        // expected manifest covariance
        S = predictS_C(Y_.rows(nonmissing), odeintparam.W, R.submat(nonmissing, nonmissing));

        // make symmetric if rounding errors lead to non-symmetry
        S = (S + arma::trans(S))/2;

        // expected manifest cross-covariance
        C = predictC_C(x, Y_.rows(nonmissing), odeintparam.W);

        // Kalman Gain
        K = computeK_C(C, S);

        // compute residuals
        residual = individualXTimeObservations - mu;

        // UPDATE
        // Update latent mean
        m = updateM_C(m_, K, residual(nonmissing));
        latentScores.row(min(rowInd) + timePoint) = arma::trans(m);

        // Update latent covariance
        P = updateP_C(P_, K, S);
        // make symmetric if rounding errors lead to non-symmetry
        P = (P + arma::trans(P))/2;

        // compute Likelihood --nonmissing

        arma::uvec m2LLInd = arma::find(Rcpp::as<arma::rowvec>(uniquePersons) == selectedPerson);
        m2LL.elem(m2LLInd) += computeIndividualM2LL_C(nObservedVariables,
                  individualXTimeObservations(nonmissing),
                  mu(nonmissing), S);

        if(breakEarly && !arma::is_finite(m2LL.elem(m2LLInd))){
          if(verbose > 0){
            Rcpp::warning("Non-finite m2LL");
          }
          Rcpp::List ret  = Rcpp::List::create(Rcpp::Named("latentScores") = latentScores,
                                               Rcpp::Named("predictedManifest") = predictedManifest,
                                               Rcpp::Named("m2LL") = m2LL);
          return(ret);
        }

      }else{
        // if all data is missing / no update is requested
        m = m_;
        latentScores.row(min(rowInd) + timePoint) = arma::trans(m);

        P = P_;
        // make symmetric if rounding errors lead to non-symmetry
        P = (P + arma::trans(P))/2;
      }
    }
  }

  Rcpp::List ret  = Rcpp::List::create(Rcpp::Named("latentScores") = latentScores,
                                       Rcpp::Named("predictedManifest") = predictedManifest,
                                       Rcpp::Named("m2LL") = m2LL);
  return(ret);
}

// [[Rcpp::export]]
Rcpp::NumericVector getGradients(Rcpp::List psydiffModel, double eps, double alpha = 0.6, double beta = 2.0, double kappa = -1.0,
                                 double timeStep = 0.01, std::string integrateFunction = "default", int verbose = 0){
  Rcpp::List psydiffModelClone = Rcpp::clone(psydiffModel);
  Rcpp::List pars = psydiffModelClone["pars"];
  Rcpp::DataFrame parameterTable = Rcpp::as<Rcpp::DataFrame>(pars["parameterTable"]);
  Rcpp::NumericVector currentParameterValues = getParameterValues_C(psydiffModelClone);
  Rcpp::NumericMatrix m2LLs(currentParameterValues.length() , 2);
  Rcpp::NumericVector gradients(currentParameterValues.length());
  arma::colvec m2LL;
  Rcpp::List fittedModel;

  for(int par = 0; par < currentParameterValues.length(); par++){

    currentParameterValues(par) += eps;
    setParameterValues_C(parameterTable, currentParameterValues, currentParameterValues.names());
    fittedModel = fitModel(psydiffModelClone, alpha, beta, kappa, timeStep, integrateFunction, breakEarly, verbose);
    m2LL = Rcpp::as<arma::colvec>(fittedModel["m2LL"]);
    m2LLs(par,0) = sum(m2LL);

    currentParameterValues(par) -= 2*eps;
    setParameterValues_C(parameterTable, currentParameterValues, currentParameterValues.names());
    fittedModel = fitModel(psydiffModelClone, alpha, beta, kappa, timeStep, integrateFunction, breakEarly, verbose);
    m2LL = Rcpp::as<arma::colvec>(fittedModel["m2LL"]);
    m2LLs(par,1) = sum(m2LL);

    currentParameterValues(par) += eps;

  }

  gradients = (m2LLs(Rcpp::_,0) - m2LLs(Rcpp::_,1))/(2*eps);
  gradients.names() = currentParameterValues.names();
  return(gradients);
}
'
  }

return(mod)
}
