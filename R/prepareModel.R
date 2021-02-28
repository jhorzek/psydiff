#' newPsydiff
#'
#' @param dataset list with fields person, observations, and dt
#' @param latentEquations string with latent equations
#' @param manifestEquations string with manifest equations
#' @param L lower triangular Cholesky decomposition of the diffusion matrix
#' @param Rchol lower triangular Cholesky decomposition of the manifest variance
#' @param A0 lower triangular Cholesky decomposition of the initial latent variance
#' @param m0 vector of initial latent means
#' @param grouping string specifying the groupings
#' @param parameters list with named parameters
#' @param groupingvariables list with variables used for grouping
#' @param additional list for anything additional that should be passed to the latent or manifest equation
#' @param srUpdate boolean: Should the square root version be used for the updates?
#' @param alpha controls the alpha parameter of the unscented transform. Should be relatively small.
#' @param beta controls the beta parameter of the unscented transform. 2 should work fine for normal distributions
#' @param kappa controls the kappa parameter of the unscented transform. 0 should work fine for normal distributions
#' @param timeStep timeStep of the numerical integration. You should pass a vector as the integration might run into problems for a specific value (e.g., .01) but work fine for another (e.g., .005). The function will stop once one of the integrations resulted in no errors
#' @param integrateFunction which function should be used for integration? Possible are rk4, runge_kutta_dopri5, and default. default will first try rk4 and - if this fails - runge_kutta_dopri5. rk4 is often a lot slower on average but runge_kutta_dopri5 tends to get stuck from time to time.
#' @param breakEarly boolean: Should the integration be stopped if the prediction for at least one time point did not work? Setting to FALSE can be useful for debugging
#' @param verbose Values > 0 will print additional information
#' @param eps controls the step size in the numerical approximation of the gradients. You should pass a vector, as this will allow psydiff to try computing the gradients for an alternative step size if one of them fails
#' @param direction direction of the steps for gradient approximation. Possible are right, left, and central.
#' @param sanityChecks boolean: Should the parameters be checked and adjusted if they might cause problems?
#'
#' @return psydiffModel that can be compiled with compileModel()
#'
#' @examples
#' library(psydiff)
#' library(ctsemOMX)
#'
#' ## Example 3: Kalman Filter
#' set.seed(175446)
#'
#' ## define the population model:
#' n.subjects = 10
#' # set the drift matrix. Note that drift eta_1_eta2 is set to equal 0 in the population.
#' ct_drift <- matrix(c(-.3,.2,0,-.5), ncol = 2)
#'
#' generatingModel<-ctsem::ctModel(Tpoints=5,n.latent=2,n.TDpred=0,n.TIpred=0,n.manifest=2,
#'                                 MANIFESTVAR=diag(0.5,2),
#'                                 LAMBDA=diag(1,2),
#'                                 DRIFT=ct_drift,
#'                                 DIFFUSION=matrix(c(.5,0,0,.5),2),
#'                                 CINT=matrix(0,nrow = 2, ncol = 1),
#'                                 T0MEANS=matrix(0,ncol=1,nrow=2),
#'                                 T0VAR=diag(1,2), type = "omx")
#'
#' # simulate a training data and testing data set
#' traindata <- ctsem::ctGenerate(generatingModel,n.subjects = n.subjects, wide = TRUE)
#' # introduce some missings:
#' traindata[1:4,2] <- NA
#' traindata[5,2:5] <- NA
#' traindata[6,7] <- NA
#' traindata[19,4] <- NA
#'
#' ## Build the analysis model. Note that drift eta1_eta2 is freely estimated
#' # although it is 0 in the population.
#' myModel <- ctsem::ctModel(Tpoints=5,n.latent=2,n.TDpred=0,n.TIpred=0,n.manifest=2,
#'                           LAMBDA=diag(1,2),
#'                           MANIFESTVAR=diag(.5,2), MANIFESTMEANS = "auto",
#'                           CINT=0,
#'                           DIFFUSION=matrix(c('eta1_eta1',0,0,'eta2_eta2'),2),
#'                           T0MEANS=matrix(c(1,2),ncol=1,nrow=2),
#' T0VAR="auto", type = "omx")
#'
#' myModel <- ctFit(myModel, dat = traindata, objective = "Kalman", useOptimizer = T,
#'                  stationary = c('T0TRAITEFFECT','T0TIPREDEFFECT'))
#' myModel$mxobj$fitfunction$result[[1]]
#'
#'
#' ## with psydiff
#' # prepare data
#' longdata <- ctWideToLong(traindata, n.manifest = 2, Tpoints =  5)
#' dat <- list("person" = longdata[,"id"], "observations" = longdata[,c("Y1", "Y2")], "dt" = longdata[,"dT"])
#'
#' ## prepare model
#'
#' latentEquations <- "
#' dx = DRIFT*x;
#' "
#'
#' manifestEquations <- "
#' y = MANIFESTMEANS + LAMBDA*x;
#' "
#'
#' LAMBDA <- fromMxMatrix(myModel$mxobj$LAMBDA)
#' DRIFT <- fromMxMatrix(myModel$mxobj$DRIFT)
#' MANIFESTMEANS <- fromMxMatrix(myModel$mxobj$MANIFESTMEANS)
#'
#' parameters <- list("LAMBDA" = LAMBDA, "DRIFT" = DRIFT,
#'                    "MANIFESTMEANS" = MANIFESTMEANS)
#'
#' m0  <- fromMxMatrix(myModel$mxobj$T0MEANS)
#' A0 <- sdeModelMatrix(values = t(chol(myModel$mxobj$T0VAR$result)),
#'                      labels = matrix(c("T0var_eta1", "",
#'                                        "T0var_eta2_eta1", "T0var_eta2"),2,2,T))
#' Rchol = sdeModelMatrix(values = t(chol(myModel$mxobj$MANIFESTVAR$result)),
#'                        labels = matrix(c("", "",
#'                                          "", ""),2,2,T))
#' L = sdeModelMatrix(values = t(chol(myModel$mxobj$DIFFUSION$result)),
#'                    labels = matrix(c("eta1_eta1", "",
#'                                      "", "eta2_eta2"),2,2,T))
#'
#' # set up model
#' model <- newPsydiff(dataset = dat, latentEquations = latentEquations,
#'                   manifestEquations = manifestEquations,
#'                   L = L, Rchol = Rchol, A0 = A0, m0 = m0,
#'                   parameters = parameters)
#'
#' # compile model
#' compileModel(model)
#'
#'
#' # fit model
#' out <- fitModel(psydiffModel = model)
#' sum(out$m2LL)
#'
#' # change parameter values
#' # inspect the model
#' newValues <- inspectModel(model)
#' psydiff::setParameterValues(parameterTable = model$pars$parameterTable,
#'                             parameterValues = newValues, parameterLabels = names(newValues))
#'
#' # fit model
#' out <- fitModel(psydiffModel = model)
#' sum(out$m2LL)
#'
#' ## optimize model with optimx
#'
#' optimized <- psydiffOptimx(model)
#'
#' ## additional grouping
#' # we will make the initial mean mm_Y1 person specific and mm_Y2 depend on a grouping parameter
#' grouping <- "
#' mm_Y1 | person;
#' mm_Y2 | group1;
#' "
#' groupinglist <- list("group1" = c(rep(1,5), rep(2,5)))
#'
#' # set up model
#' model <- newPsydiff(dataset = dat, latentEquations = latentEquations,
#'                   manifestEquations = manifestEquations, grouping = grouping,
#'                   L = L, Rchol = Rchol, A0 = A0, m0 = m0,
#'                   parameters = parameters, groupingvariables = groupinglist, compile = TRUE)
#' parval <- psydiff::getParameterValues(model)
#'
#' optimized <- psydiffOptimx(model, control = list(trace = 1))
#'
#' ## The following example is taken from ctsem and also demonstrates the use of the GIST optimizer
#' library(ctsemOMX)
#' data('ctExample3')
#' model <- ctModel(n.latent = 1, n.manifest = 3, Tpoints = 100,
#'                  LAMBDA = matrix(c(1, 'lambda2', 'lambda3'), nrow = 3, ncol = 1),
#'                  CINT= matrix('cint'), T0VAR = diag(1),
#'                  MANIFESTMEANS = matrix(c(0, 'manifestmean2', 'manifestmean3'), nrow = 3,
#'                                         ncol = 1))
#' fit <- ctFit(dat = ctExample3, ctmodelobj = model, objective = 'Kalman',
#'              stationary = c("T0TRAITEFFECT", "T0TIPREDEFFECT"), useOptimizer = F)
#' omxGetParameters(fit$mxobj)
#' fit$mxobj$fitfunction$result[[1]]
#'
#' latentEquations <- "
#' dx(0) = cint + driftEta1*x(0);
#' "
#' manifestEquations <- "
#' y(0) = x(0);
#' y(1) = manifestmean2 + lambda2*x(0);
#' y(2) = manifestmean3 + lambda3*x(0);
#' "
#' # The optimization depends highly on the starting values. The values blow are taken from ctsem
#' parameters <- list("driftEta1" = -1.150227824,
#'                    "cint" = 11.214338733,
#'                    "lambda2" = 0.480098989,
#'                    "lambda3" = 0.959200513,
#'                    "manifestmean2" = 2.824753235,
#'                    "manifestmean3" = 5.606085485)
#'
#' m0  <- fromMxMatrix(fit$mxobj$T0MEANS)
#' A0 <- sdeModelMatrix(values = fit$mxobj$T0VARchol$result,
#'                      labels = matrix("",1,1))
#' Rchol = sdeModelMatrix(values = fit$mxobj$MANIFESTVARchol$result,
#'                        labels = matrix(c("mvar1", "", "",
#'                                          "", "mvar2", "",
#'                                          "", "", "mvar3"),3,3,T))
#' L = sdeModelMatrix(values = fit$mxobj$DIFFUSIONchol$result,
#'                    labels = matrix(c("lvar"),1,1,T))
#'
#' longdata <- ctWideToLong(ctExample3, n.manifest = 3, Tpoints =  100)
#' dat <- list("person" = longdata[,"id"], "observations" = longdata[,c("Y1", "Y2", "Y3")], "dt" = longdata[,"dT"])
#'
#' # set up model
#' model <- newPsydiff(dataset = dat, latentEquations = latentEquations,
#'                     manifestEquations = manifestEquations,
#'                     L = L, Rchol = Rchol, A0 = A0, m0 = m0,
#'                     parameters = parameters, verbose = 0, kappa = 0, alpha = .81, beta = 2)
#'
#' compileModel(model)
#' out <- fitModel(model)
#' out$m2LL
#' getParameterValues(model)
#'
#' opt <- psydiffOptimx(model, method = c('Nelder-Mead', 'BFGS', 'nlm', 'nlminb'), control = list(trace = 1))
#'
#' startingValues <- psydiff::getParameterValues(model)
#' adaptiveLassoWeights <- rep(1, length(startingValues))
#' names(adaptiveLassoWeights) <- names(startingValues)
#' regularizedParameters <- "lambda2"
#' lambda <- 10
#'
#' opt <- GIST(model = model, startingValues = startingValues, lambda = lambda,
#'             adaptiveLassoWeights = adaptiveLassoWeights,
#'             regularizedParameters = regularizedParameters,
#'             verbose = 1, maxIter_out = 200, sig = .4, break_outer = 10e-20)
#'
#' getParameterValues(opt$model)
#' out <- fitModel(opt$model)
#' matplot(out$predictedManifest, type = "l")
#' points(dat$observations[,1])
#' points(dat$observations[,2])
#' points(dat$observations[,3])
#'
#' ## second order model
#' set.seed(12391)
#'
#' library(ctsemOMX)
#' DRIFT <- matrix(c(0, 1,
#'                   -.005, -.04),2,2,T)
#' LAMBDA <- matrix(c(1,0),1,2,T)
#' MANIFESTVAR <- matrix(1)
#' LATENTVAR <- matrix(c(0.001, 0,
#'                       0, 0.001),2,2,T)
#' ctMod <- ctModel(LAMBDA = LAMBDA, n.manifest = 1, n.latent = 2,
#'                  Tpoints = 300,
#'                  T0MEANS = c(0,1),
#'                  T0VAR = diag(2), CINT = matrix(c(0, 1),nrow = 2),
#'                  MANIFESTVAR = MANIFESTVAR, MANIFESTMEANS = 0,
#'                  DRIFT = DRIFT, DIFFUSION = LATENTVAR)
#' ctDat <- ctGenerate(ctMod, n.subjects = 1)
#' plot(ctDat$Y1, type = "p")
#'
#' ## Define model in psydiff
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
#' m0 <- matrix(c("0","1"), nrow = 2)
#'
#' Rchol <- matrix("1",1,1)
#' L <- ctMod$DIFFUSION
#'
#' parameters <- list("cint" = 1, "a" = -.5, "b" = -.5)
#'
#' dat <- list("person" = ctDat$id,
#'             "observations" = matrix(ctDat$Y1, ncol = 1),
#'             "dt" = c(0,rep(1,length(ctDat$time)-1)))
#'
#' ## optimize model
#'
#' model <- newPsydiff(dataset = dat, latentEquations = latentEquations,
#'                     manifestEquations = manifestEquations,
#'                     L = L, Rchol = Rchol, A0 = A0, m0 = m0,
#'                     parameters = parameters)
#' compileModel(model)
#'
#' (startingValues <- psydiff::getParameterValues(model))
#'
#' ## Of the optimizers I tried, Nelder-Mead results in  the best fit for this model
#' nm <- psydiffOptimNelderMead(model, control = list("trace" = 1))
#' lines(model$predictedManifest, col = "#008080", lwd = 3) # Note: model object was changed by reference
#' @export

newPsydiff <- function(dataset, latentEquations, manifestEquations, L, Rchol, A0, m0,
                       grouping = NULL, parameters, groupingvariables = NULL,
                       additional = NULL,
                       # settings for the integration
                       srUpdate = TRUE, alpha = .2, beta = 2.0, kappa = 0.0,
                       timeStep = NULL, integrateFunction = "default", breakEarly = TRUE, verbose = 0,
                       eps = rev(c(1e-4, 1e-5, 1e-6, 1e-8)), direction = "central",
                       sanityChecks = TRUE){
  nlatent <- ncol(L)
  nmanifest <- ncol(Rchol)

  # checks
  checkEquation(latentEquations)
  checkEquation(manifestEquations)
  if((!is.list(dataset)) | is.data.frame(dataset)){
    stop("dataset must be a list")
  }
  datasetNames <- names(dataset)
  for(datasetName in datasetNames){
    if(!(datasetName %in% c("person", "observations", "dt"))){
      stop("Unknown field in dataset. Dataset must have the fields 'person',
           'observations', 'dt'")
    }
  }
  if(!ncol(dataset[["observations"]]) == nmanifest){
    stop("Dataset and R differ in their number of manifest variables.")
  }
  if(!is.list(parameters)){
    stop("parameters must be of class data.frame")
  }
  if(!(is.null(additional) || is.list(additional))){
    stop("additional must be of class data.frame or NULL")
  }

  # set time step defaults
  if(is.null(timeStep)){
    # minimum of 10 steps
    timeStep <- seq(min(dataset$dt[dataset$dt > 0])/30, min(dataset$dt[dataset$dt > 0])/10, length.out = 5)
  }

  ## Extract persons
  persons <- dataset$person

  ## Define matrices, vectors, etc.
  ## redefine A0, L, and Rchol
  A0Sep <- sepNameFromValue(as.matrix(A0))
  LSep <- sepNameFromValue(as.matrix(L))
  RcholSep <- sepNameFromValue(as.matrix(Rchol))


  # build log-Cholesky
  A0LogChol <- chol2LogChol(matValues = A0Sep$values, matLabels = A0Sep$labels, sanityChecks = sanityChecks, matrixName = "A0")
  LLogChol <- chol2LogChol(matValues = LSep$values, matLabels = LSep$labels, sanityChecks = sanityChecks, matrixName = "L")
  RLogChol <- chol2LogChol(matValues = RcholSep$values, matLabels = RcholSep$labels, sanityChecks = sanityChecks, matrixName = "Rchol")

  ## Element refers to containers in which the parameters are stored
  parameters$A0LogChol <- A0LogChol
  parameters$LLogChol <- LLogChol
  parameters$RLogChol <- RLogChol
  parameters$m0 <- m0

  elementNames <- names(parameters)

  ## initialize parameterTable and starting values. parameterTable will hold the
  # parameter values for each individual and tell other functions,
  # where those parameters can be found (e.g., a4 is in matrix A in row
  # 1, col 1)
  parTabAndStart <- extractParameters(elementNames = elementNames, parameters = parameters)
  parameterTableInit <- parTabAndStart$parameterTableInit
  parameterList <- parTabAndStart$startValues

  parameterTable <- do.call("rbind", replicate(length(unique(persons)), parameterTableInit, simplify = FALSE))
  parameterTable$person <- rep(unique(persons), each = nrow(parameterTableInit))

  ## check groupings
  parameterTable <- makeGroups(grouping = grouping, groupingvariables = groupingvariables,
                               parameterTable = parameterTable)

  ## prepare equations
  manifestEquationsStart <- "
  arma::colvec measurementEquation(const arma::colvec x, const odeintpars &pars,
                                                             const int &person,
                                                             const double &t){
    arma::colvec y(pars.nmanifest, arma::fill::zeros);
  "
  manifestEquationsMiddle <- prepareEquations(equations = manifestEquations,
                                              parameters = parameters)
  manifestEquationsEnd <- "
    return(y);
  }
  "
  manifestEquationsCpp <- paste0(manifestEquationsStart, manifestEquationsMiddle, manifestEquationsEnd)

  latentEquationsStart <- "
  arma::colvec latentEquation(const arma::colvec x, const odeintpars &pars,
                            const int &person,
                            const double &t){
                            arma::colvec dx(pars.nlatent, arma::fill::zeros);
  "
  latentEquationsMiddle <- prepareEquations(equations = latentEquations,
                                            parameters = parameters)
  latentEquationsEnd <- "
    return(dx);
  }
  "
  latentEquationsCpp <- paste0(latentEquationsStart, latentEquationsMiddle, latentEquationsEnd)

  ## set up model

  modelCpp <- modelTemplate(srUpdate)

  modelCpp <- stringr::str_replace_all(modelCpp, "MEASUREMENTEQUATIONPLACEHOLDER", manifestEquationsCpp)
  modelCpp <- stringr::str_replace_all(modelCpp, "LATENTEQUATIONPLACEHOLDER", latentEquationsCpp)

  ## return model
  pars <- list("parameterList" = parameterList,
               "parameterTable" = parameterTable
  )
  psydiffModel <- list("pars" = pars,
                       "nlatent" = nlatent,
                       "nmanifest" = nmanifest,
                       "data" = dataset,
                       "additional" = additional,
                       "m2LL" = rep(0.0, length(unique(persons))),
                       "predictedManifest" = matrix(-99999.99,
                                                    ncol = ncol(dataset[["observations"]]),
                                                    nrow = nrow(dataset[["observations"]])),
                       "latentScores" = matrix(-99999.99,
                                               ncol = nlatent,
                                               nrow = nrow(dataset[["observations"]])),
                       "alpha" = alpha,
                       "beta" = beta,
                       "kappa" = kappa,
                       "timeStep" = timeStep,
                       "integrateFunction" = integrateFunction,
                       "breakEarly" = breakEarly,
                       "verbose" = verbose,
                       "eps" = eps,
                       "direction" = direction,
                       "cppCode" = modelCpp
  )

  message("Use compileModel to compile the model.")
  return("psydiffModel" = psydiffModel)

}

#' compileModel
#'
#' compile the model from newPsydiff
#'
#' @param psydiffModel model from newPsydiff
#' @export
compileModel <- function(psydiffModel){
  filename <- tempfile(pattern = "psydiff_", tmpdir = tempdir(), fileext = ".cpp")
  fileConn<-file(filename)
  writeLines(psydiffModel$cppCode, fileConn)
  close(fileConn)
  cat("Compiling model...")
  Rcpp::sourceCpp(file = filename)
  cat("Done.\n")
  message("The model can be fitted with the fitModel()-function and the returned object. Use getGradients() to compute the central gradients of the model. If you want to compute gradients in parallel, use compileParallelGradients to set up the cluster.")
}

#' compileParallelGradients
#'
#' compile the model from newPsydiff for use with parallel gradient computation. This function will invoke parallel::makeCluster(nCores, type = "PSOCK"). Stop the cluster with stopParallelGradients().
#'
#' @param psydiffModel model from newPsydiff
#' @param nCores number of cores to use
#' @return psydiffModel with additional cl field for workers
#' @export
compileParallelGradients <- function(psydiffModel, nCores){
  availableCores <- parallel::detectCores()
  if(nCores > availableCores){
    stop("More cores requested than available on your machine. Use parallel::detectCores() to detect the number of cores you can maximally use.")
  }
  cl <- parallel::makeCluster(nCores, type = "PSOCK")

  parallel::clusterExport(cl, c("psydiffModel"), envir = environment())
  parallel::clusterEvalQ(cl, library(psydiff))

  cat("Compiling models...")
  parallel::clusterEvalQ(cl, compileModel(psydiffModel))
  cat("Done.\n")
  message("Gradients can now be computed with getGradientsParallel(). Stop the cluster with stopParallelGradients()")

  psydiffModel$cl <- cl
  return(psydiffModel)
}

#' stopParallelGradients
#'
#' Stops the workers.
#'
#' @param psydiffModel model from newPsydiff
#' @return
#' @export
stopParallelGradients <- function(psydiffModel){
  availableCores <- parallel::stopCluster(psydiffModel$cl)
}

#' getGradientsParallel
#'
#' Computes gradients in parallel. Requires compileParallelGradients to be used first!
#'
#' @param psydiffModel model from newPsydiff
#' @export
getGradientsParallel <- function(psydiffModel){
  pars <- getParameterValues(psydiffModel)
  parallel::clusterExport(psydiffModel$cl, c("psydiffModel"), envir = environment())
  parGrad <- parallel::parLapplyLB(cl = psydiffModel$cl, X = seq_len(length(getParameterValues(psydiffModel))), fun = function(X,psydiffModel,gradFun, pars){
    getGradient(psydiffModel = psydiffModel, par = X-1, currentParameterValues = pars)
  }, psydiffModel = psydiffModel, gradFun = getGradient, pars = pars)
  parGrad <- unlist(parGrad)
  return(parGrad[names(pars)])
}

checkEquation <- function(equation){
  if(!stringr::str_detect(equation, ";", negate = FALSE)){
    stop(paste0("Expected each command in ...\n", equation, "\n to end with ; "))
  }
}

extractParameters <- function(elementNames, parameters){
  parameterLabels <- c()
  parameterValues <- c()
  parameterTargets <- c()
  parameterTargetClass <- c()
  parameterRow <- c()
  parameterCol <- c()
  parameterChanged <- c()

  startValues <- parameters
  for(elementName in elementNames){
    currentElement <- parameters[[elementName]]
    if(is.matrix(currentElement)){
      elementValues <- matrix(NA, nrow(currentElement), ncol(currentElement))
      for(ro in 1:nrow(currentElement)){
        for(co in 1:ncol(currentElement)){
          matrixelement <- currentElement[ro,co]
          if(!stringr::str_detect(matrixelement, "=", negate = FALSE)){
            elementValues[ro,co] <- as.numeric(matrixelement)
            next
          }
          matrixelementSplit <- stringr::str_split(matrixelement, "=")[[1]]
          matrixelementSplit <- stringr::str_remove_all(matrixelementSplit, " ")
          elementValues[ro,co] <- as.numeric(matrixelementSplit[2])
          parameterLabels <- c(parameterLabels, matrixelementSplit[1])
          parameterValues <- c(parameterValues, as.numeric(matrixelementSplit[2]))
          parameterTargets <- c(parameterTargets, elementName)
          parameterTargetClass <- c(parameterTargetClass, "matrix")
          parameterRow <- c(parameterRow, ro-1)
          parameterCol <- c(parameterCol, co-1)
          parameterChanged <- c(parameterChanged, TRUE)
        }
      }
    }else if(is.numeric(currentElement) && length(currentElement)==1){
      elementValues <- currentElement
      parameterLabels <- c(parameterLabels, elementName)
      parameterValues <- c(parameterValues, currentElement)
      parameterTargets <- c(parameterTargets, elementName)
      parameterTargetClass <- c(parameterTargetClass, "double")
      parameterRow <- c(parameterRow, 0)
      parameterCol <- c(parameterCol, 0)
      parameterChanged <- c(parameterChanged, TRUE)
    }else if(is.vector(currentElement)){
      warning(paste0(elementName, " is of type vector. It will be interpreted as column vector. Please make sure that this is considered in your model definition."))
      elementValues <- rep(NA, length(currentElement))
      for(ro in 1:length(currentElement)){
        matrixelement <- currentElement[ro]
        if(!stringr::str_detect(matrixelement, "=", negate = FALSE)){
          elementValues[ro] <- as.numeric(matrixelement)
          next
        }
        matrixelementSplit <- stringr::str_split(matrixelement, "=")[[1]]
        matrixelementSplit <- stringr::str_remove_all(matrixelementSplit, " ")
        elementValues[ro] <- as.numeric(matrixelementSplit[2])
        parameterLabels <- c(parameterLabels, matrixelementSplit[1])
        parameterValues <- c(parameterValues, as.numeric(matrixelementSplit[2]))
        parameterTargets <- c(parameterTargets, elementName)
        parameterTargetClass <- c(parameterTargetClass, "colvec")
        parameterRow <- c(parameterRow, ro-1)
        parameterCol <- c(parameterCol, 0)
        parameterChanged <- c(parameterChanged, TRUE)
      }
    }else{
      stop(paste0("Unknown data type for parameter ", currentElement, ". Allowed are numeric, matrix, and vector."))
    }
    startValues[[elementName]] <- elementValues
  }

  parameterTableInit <- data.frame("label" = parameterLabels,
                                   "value" = parameterValues,
                                   "target" = parameterTargets,
                                   "targetClass" = parameterTargetClass,
                                   "row" = parameterRow,
                                   "col" = parameterCol,
                                   "changed" = parameterChanged)
  return(list("parameterTableInit" = parameterTableInit,
              "startValues" = startValues))

}

makeGroups <- function(grouping, groupingvariables, parameterTable){
  if(!is.null(grouping)){
    persons <- parameterTable$person
    checkEquation(grouping)
    ## split groupings
    grouping <- stringr::str_remove_all(grouping, "\n")
    grouping <- stringr::str_remove_all(grouping, " ")
    splitgrouping <- stringr::str_split(grouping, ";")[[1]]
    splitgrouping <- splitgrouping[!splitgrouping == ""]

    for(grp in splitgrouping){
      grpsplit <- stringr::str_split(grp, "\\|")[[1]]
      parname <- grpsplit[1]
      if(grpsplit[2] %in% c("person", "persons", "individual", "individuals")){
        for(p in unique(persons)){
          parameterTable[parameterTable$label == parname & parameterTable$person == p,"label"] <- paste0(parname, "_p", p)
        }
      }else{
        ## user-specified splitting
        groupingvariableNames <- names(groupingvariables)
        if(!grpsplit[2] %in% groupingvariableNames){
          stop(paste0("Grouping variable ", grpsplit[2], " not found in groupingvariables"))}
        groupingIndicator <- groupingvariables[[grpsplit[2]]]
        for(p in unique(persons)){
          parameterTable[parameterTable$label == parname & parameterTable$person == p,"label"] <- paste0(parname, "_p", groupingIndicator[p])
        }
      }
    }
  }
  return(parameterTable)
}


prepareEquations <- function(equations, parameters){
  ## split at any special symbols (+, -, ...) except for _
  elementsInEquations <- stringr::str_split(equations, '[^[:alnum:]|^_]')[[1]]
  elementsInEquations <- elementsInEquations[!elementsInEquations== ""]
  equationsCombined <- "Rcpp::List modelPars = pars.parameterList;
  "

  targets <- names(parameters)

  for(elementInEquation in elementsInEquations){
    if(elementInEquation %in% targets){
      targetClass <- parameters[[elementInEquation]]
      if(is.matrix(targetClass)){
        equationsCombined <- paste0(equationsCombined, '
               arma::mat ', elementInEquation, ' = modelPars["',elementInEquation,'"]; \n')
      }else if(is.vector(targetClass)){
        if(length(targetClass) == 1){
          equationsCombined <- paste0(equationsCombined, '
               double ', elementInEquation, ' = modelPars["',elementInEquation,'"]; \n')
        }else{
          equationsCombined <- paste0(equationsCombined, '
               arma::colvec ', elementInEquation, ' = modelPars["',elementInEquation,'"]; \n')
        }
      }else{
        stop("Error in prepareEquations: targetClass unknown.")
      }
      # remove element from target to prevent multiple entries
      targets <- targets[!targets == elementInEquation]
    }
  }
  equationsCombined <- paste0(equationsCombined, "
                              ", equations)
  return(equationsCombined)
}

fromMxMatrix <- function(mxMatrixObject){
  values <- as.character(mxMatrixObject$values)
  labels <- as.character(mxMatrixObject$labels)
  values[!is.na(labels)] <- paste0(labels[!is.na(labels)], " = ", values[!is.na(labels)])
  values <- matrix(values,
                   nrow = nrow(mxMatrixObject$values),
                   ncol = ncol(mxMatrixObject$values))
  return(values)

}

sdeModelMatrix <- function(values, labels){
  if(!(nrow(values) == nrow(labels) && ncol(values) == ncol(labels))){
    stop("values and labels have different dimensions")
  }
  ro <- nrow(values)
  co <- ncol(values)
  values <- as.character(values)
  labels <- as.character(labels)
  naOrEmpty <- is.na(labels) | labels == ""
  values[!naOrEmpty] <- paste0(labels[!naOrEmpty], " = ", values[!naOrEmpty])
  values <- matrix(values,
                   nrow = ro,
                   ncol = co)
  return(values)

}

sepNameFromValue <- function(mat){
  if(is.numeric(mat)){
    return(list("labels" = NULL, "values" = mat))
  }
  labelMat <- matrix("", ncol = ncol(mat), nrow = nrow(mat))
  valueMat <- matrix(NA, ncol = ncol(mat), nrow = nrow(mat))
  for(i in seq_len(nrow(mat))){
    for(ii in seq_len(ncol(mat))){
      if(grepl("=", mat[i,ii])){
        splitted <- stringr::str_split(mat[i,ii], "=")[[1]]
        labelMat[i,ii] <- stringr::str_remove_all(splitted[1], " ")
        valueMat[i,ii] <- as.numeric(splitted[2])
      }else{
        valueMat[i,ii] <- as.numeric(mat[i,ii])
      }
    }
  }
  return(list("labels" = labelMat, "values" = valueMat))
}

chol2LogChol <- function(matValues, matLabels, sanityChecks, matrixName = ""){
  if(sanityChecks){
    if(any(abs(diag(matValues)) < .01)){
      warning(paste0("Setting the diagonal elements of ", matrixName, " to very small values will often result in errors. The small values will be replaced by .01. Set sanityChecks = FALSE to prevent this."))
      matValues[(abs(matValues) < .01) & diag(nrow(matValues))] <- .01
    }
  }
  diag(matValues) <- log(diag(matValues))
  diagLabels <- diag(matLabels)
  diag(matLabels)[!diagLabels == ""] <- paste0("ln", diagLabels[!diagLabels == ""])

  combinedMat <- matLabels
  for(i in seq_len(nrow(combinedMat))){
    for(ii in seq_len(ncol(combinedMat))){
      if(combinedMat[i,ii] == ""){
        combinedMat[i,ii] <- matValues[i,ii]
      }else{
        combinedMat[i,ii] <- paste0(combinedMat[i,ii], " = ", matValues[i,ii])
      }
    }
  }
  return(combinedMat)
}

