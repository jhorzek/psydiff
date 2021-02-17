#' getDataTemplate
#'
#' Generate a template for data simulation
#' @param nSubjects number of persons
#' @param nVariables number of simulatedObservation variables
#' @param nTimpoints number of time points per person
#' @param dts either a single dt, a vector with dts for an individual or a vector with dts for all individuals
#' @export
getDataTemplate <- function(nSubjects, nVariables, nTimpoints, dts){
  if(length(dts) == 1){
    dts <- rep(c(0,rep(dts, nTimpoints-1)),nSubjects)
  }else if(length(dts) == nTimpoints){
    dts <- rep(dts,nSubjects)
  }else if(length(dts) == nTimpoints*nSubjects){}else{
    stop("The length of dts has to be 1, nTimepoints, or nTimpoints*nSubjects")
  }
  dat <- list("person" = rep(seq_len(nSubjects), each = nTimpoints),
              "observations" = matrix(NA,
                                      nrow = nSubjects*nTimpoints ,
                                      ncol = nVariables),
              "dt" = dts)
  return(dat)
}

#' simulateData
#'
#' Simulate data from a given model. The dataset in the model should match the pattern of the
#' data you are interested in. A new data set can for instance be generate with getDataTemplate.
#'
#' @param model psydiff model
#' @examples
#' ### Simulation of a second order stochastic differential equation model
#'
#' manifestEquations <- "
#' y = x(0);
#' "
#'
#' latentEquations <- "
#' dx(0) = x(1);
#' dx(1) = cint + a*x(0) +b*x(1);
#' "
#'
#' A0 <- diag(1,2)
#' m0 <- c(0,1)
#'
#' Rchol <- matrix(10,1,1)
#' L <- matrix(0,2,2)
#'
#' parameters <- list("cint" = 1, "a" = -.005, "b" = -.04)
#'
#' dataTemplate <- getDataTemplate(nSubjects = 1,
#'                                 nVariables = 1,
#'                                 nTimpoints = 200,
#'                                 dts = 1)
#'
#' # set up model
#' model <- newPsydiff(dataset = dataTemplate,
#'                     latentEquations = latentEquations,
#'                     manifestEquations = manifestEquations,
#'                     L = L, Rchol = Rchol, A0 = A0, m0 = m0,
#'                     parameters = parameters)
#'
#' compileModel(model)
#'
#' datSim <- simulateData(model)
#'
#' plot(datSim$predictedManifest, type = "l")
#' points(datSim$simulatedObservation)
#' @export
simulateData <- function(model){

  modelClone <- clonePsydiffModel(model)
  modelClone$pars$parameterTable$changed <- T # force refit for all individuals
  # generate predictions
  f <- fitModel(modelClone, skipUpdate = TRUE)

  ## extract manifest residual covariance matrix
  R <- model$pars$parameterList$Rchol%*%t(model$pars$parameterList$Rchol)

  # add noise to observations
  simulatedObservation <- f$predictedManifest
  for(ro in seq_len(nrow(simulatedObservation))){
    simulatedObservation[ro, ] <- simulatedObservation[ro, ] + mvtnorm::rmvnorm(n = 1,
                     mean = rep(0, model$nmanifest),
                     sigma = R)
  }
  return(list("latentScores" = f$latentScores,
              "predictedManifest" = f$predictedManifest,
              "simulatedObservation" = simulatedObservation)
  )
}
