#' panelsde
#'
#' @example
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
#' ## with panelsde
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
#' model <- panelsde(dataset = dat, latentEquations = latentEquations,
#'                   manifestEquations = manifestEquations,
#'                   L = L, Rchol = Rchol, A0 = A0, m0 = m0,
#'                   parameters = parameters, compile = TRUE)
#'
#' # fit model
#' out <- fitModel(panelSDEModel = model)
#' sum(out$m2LL)
#'
#' # change parameter values
#' parval <- psydiff::getParameterValues(model)+1
#' psydiff::setParameterValues(parameterTable = model$pars$parameterTable,
#'                             parameterValues = parval, parameterLabels = names(parval))
#'
#' # fit model
#' out <- fitModel(panelSDEModel = model)
#' sum(out$m2LL)
#'
#' ## optimize model with optim
#' fitfun <- function(parval, model){
#'   psydiff::setParameterValues(parameterTable = model$pars$parameterTable,
#'                               parameterValues = parval, parameterLabels = names(parval))
#'   out <- try(fitModel(panelSDEModel = model))
#'   if(any(class(out) == "try-error")){
#'     return(NA)
#'   }
#'   return(sum(out$m2LL))
#' }
#'
#' optimized <- optim(par = parval, fn = fitfun, gr = NULL, model, method = "BFGS")
#' optimized$value
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
#' model <- panelsde(dataset = dat, latentEquations = latentEquations,
#'                   manifestEquations = manifestEquations, grouping = grouping,
#'                   L = L, Rchol = Rchol, A0 = A0, m0 = m0,
#'                   parameters = parameters, groupingvariables = groupinglist, compile = TRUE)
#' parval <- psydiff::getParameterValues(model)
#'
#' optimized <- optim(par = parval, fn = fitfun, gr = NULL, model, method = "BFGS")
#' optimized$par
#' optimized$value

panelsde <- function(dataset, latentEquations, manifestEquations, L, Rchol, A0, m0,
                     grouping = NULL, parameters, groupingvariables = NULL,
                     additional = NULL, compile = TRUE){
  nlatent <- ncol(L)
  nmanifest <- ncol(Rchol)

  # checks
  checkEquation(latentEquations)
  checkEquation(manifestEquations)
  if(!is.list(dataset)){
    stop("dataset must be of class list")
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

  ## Extract persons
  persons <- dataset$person

  ## Define matrices, vectors, etc.
  ## Element refers to containers in which the parameters are stored
  parameters$L <- L
  parameters$A0 <- A0
  parameters$m0 <- m0
  parameters$Rchol <- Rchol

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
  arma::colvec measurementEquation(const arma::colvec &x, const odeintpars &pars,
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
  arma::colvec latentEquation(const arma::colvec &x, const odeintpars &pars,
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
  modelCpp <- modelTemplate()
  modelCpp <- stringr::str_replace_all(modelCpp, "MEASUREMENTEQUATIONPLACEHOLDER", manifestEquationsCpp)
  modelCpp <- stringr::str_replace_all(modelCpp, "LATENTEQUATIONPLACEHOLDER", latentEquationsCpp)

  ## compile
  if(compile){
    filename <- tempfile(pattern = "panelsde_", tmpdir = tempdir(), fileext = ".cpp")
    fileConn<-file(filename)
    writeLines(modelCpp, fileConn)
    close(fileConn)
    cat("Compiling model...")
    Rcpp::sourceCpp(file = filename)
    cat("Done.\n The model can be fitted with the fitModel()-function and the returned object.")
  }

  ## return model
  pars <- list("parameterList" = parameterList,
               "parameterTable" = parameterTable
  )
  panelSDEModel <- list("pars" = pars,
                        "nlatent" = nlatent,
                        "nmanifest" = nmanifest,
                        "data" = dataset,
                        "additional"= additional
  )

  if(compile){
    return(panelSDEModel)
  }

  return(list("panelSDEModel" = panelSDEModel, "cppCode" = modelCpp))

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
    }else if(is.numeric(currentElement)){
      elementValues <- currentElement
      parameterLabels <- c(parameterLabels, elementName)
      parameterValues <- c(parameterValues, currentElement)
      parameterTargets <- c(parameterTargets, elementName)
      parameterTargetClass <- c(parameterTargetClass, "double")
      parameterRow <- c(parameterRow, 0)
      parameterCol <- c(parameterCol, 0)
      parameterChanged <- c(parameterChanged, TRUE)
    }else if(is.vector(currentElement)){
      warning(paste0(elementName, " is of type vector. It will be interpreted as column vector.
                     Please make sure that this is considered in your model definition."))
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
  ## remove any special symbols (+, +, ...)
  elementsInEquations <- stringr::str_split(equations, '[^[:alnum:]]')[[1]]
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
      }else if(is.vec(targetClass)){
        equationsCombined <- paste0(equationsCombined, '
               arma::colvec ', elementInEquation, ' = modelPars["',elementInEquation,'"]; \n')
      }else if(is.numeric(targetClass)){
        equationsCombined <- paste0(equationsCombined, '
               double ', elementInEquation, ' = modelPars["',elementInEquation,'"]; \n')
      }else{
        stop("Error in prepareEquations: targetClass unknown.")
      }
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
