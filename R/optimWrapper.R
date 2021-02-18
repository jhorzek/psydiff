#' psydiffOptimNelderMead
#'
#' wrapper for Nelder-Mead optimizer from optim
#' @param model psydiff model object (must have been compiled before calling psydiffOptimNelderMead)
#' @param ... additional parameters to pass to optim
#' @export
psydiffOptimNelderMead <- function(model, ...){
  startingValues <- psydiff::getParameterValues(model)

  out <- optim(par = startingValues, fn = psydiffFitInternal, method = "Nelder-Mead",
               model = model,
               parsnames = names(startingValues), ...)
  return(out)
}

#' psydiffOptimBFGS
#'
#' wrapper for BFGS optimizer from optim
#' @param model psydiff model object (must have been compiled before calling psydiffOptimBFGS)
#' @param ... additional parameters to pass to optim
#' @export
#'
psydiffOptimBFGS <- function(model, ...){
  startingValues <- psydiff::getParameterValues(model)

  out <- optim(par = startingValues, fn = psydiffFitInternal, method = "BFGS",
               model = model,
               parsnames = names(startingValues), ...)
  return(out)
}

#' psydiffOptimx
#'
#' wrapper for optimx optimization
#' @param model psydiff model object (must have been compiled before calling psydiffOptimx)
#' @param ... additional parameters to pass to optimx
#' @import optimx
#' @export
#'
psydiffOptimx <- function(model, method=c("Nelder-Mead","BFGS", "nlm", "nlminb"), ...){
  startingValues <- psydiff::getParameterValues(model)

  out <- optimx::optimx(par = startingValues, fn = psydiffFitInternal, method = method,
                        model = model,
                        parsnames = names(startingValues), ...)
  return(out)
}

#' extractOptimx
#'
#' sets the model parameters to the best values obtained from optimx
#' @param model psydiff model
#' @param opt result from calling psydiffOptimx
#' @export
#'
extractOptimx <- function(model, opt){
  if(!any(class(opt) == "optimx")){
    stop("opt has to be of class optimx")
  }
  pars <- getParameterValues(model)
  values <- opt$value
  bestValue <- which(values == min(values))[1] # if multiple optimizers find the same optimum, the first will be used
  optimizer <- rownames(opt)[bestValue]
  optimizedPars <- data.matrix(opt[optimizer,names(pars)])
  setParameterValues(parameterTable = model$pars$parameterTable,
                     parameterValues = optimizedPars,
                     parameterLabels = colnames(optimizedPars))
  message(paste0("The lowest minimum of ", min(values), " was found with ", optimizer, ". Returning the model with updated parameter values. WARNING: The model has not been fitted yet! Use fitModel on the returned object."))
  return(model)
}

#' psydiffOptimxMultiStart
#'
#' wrapper for optimx optimization
#' @param model psydiff model object (must have been compiled before calling psydiffOptimx)
#' @param ... additional parameters to pass to optimx
#' @import optimx
#' @export
#'
psydiffOptimxMultiStart <- function(model, method=c("Nelder-Mead","BFGS"), startMatrix = NULL, nstart = 20, sd = 1, ...){

  startingValues <- psydiff::getParameterValues(model)
  if(is.null(startMatrix)){
    startMatrix <- matrix(rep(startingValues, nstart), nrow = nstart, byrow = TRUE)
    colnames(startMatrix) <- names(startingValues)
    startMatrix <- startMatrix + rnorm(nstart * length(startingValues), mean = 0, sd = sd)
  }

  out <- optimx::multistart(parmat = startMatrix, fn = psydiffFitInternal, method = method,
                            model = model,
                            parsnames = names(startingValues), ...)
  return(out)
}

#' psydiffDEoptim
#'
#' wrapper for DEoptim optimization
#' @param model psydiff model object (must have been compiled before calling psydiffDEoptim)
#' @param lower vector specifying a lower bound for each parameter in the model
#' @param upper vector specifying an upper bound for each parameter in the model
#' @param ... additional parameters to pass to DEoptim
#' @import DEoptim
#' @export
#'
psydiffDEoptim <- function(model, lower, upper, ...){
  startingValues <- psydiff::getParameterValues(model)

  out <- DEoptim::DEoptim(fn = ffit, lower = lower, upper = upper,
                          model = model,
                          parsnames = names(startingValues),
                          ...)
  return(out)
}

#' psydiffFitInternal
#'
#' internal function for optimization
#' @param pars vector with parameter values
#' @param parsnames vector with parameter names
#' @param model psydiff model object (must have been compiled before calling psydiffDEoptim)
#' @export
#'
psydiffFitInternal <- function(pars, parsnames, model){
  setParameterValues(model$pars$parameterTable, parameterLabels = parsnames, parameterValues = pars)

  invisible(capture.output(fit <- try(fitModel(model),
                                      silent = TRUE),
                           type = "message"))

  if(any(class(fit) == "try-error") || anyNA(fit$m2LL)){
    return(99999999999)
  }
  return(sum(fit$m2LL))
}

#' getStandardErrors
#'
#' compute standard errors for a fitted model using the diagonal elements in the inverse of the hessian approximation from optimHess
#' @param model psydiff model
#' @param ... additional arguments to pass to optimHess
#' @return vector with standard errors
#' @export
getStandardErrors <- function(model, ...){
  pars <- getParameterValues(model)
  hess <- optimHess(fn = psydiffFitInternal, par = pars,
                            model = model,
                            parsnames = names(pars),
                    ...)
  # Note: We are minimizing the negative log likelihood.
  # The Hessian is therefore the "observed Fisher Information"
  # and it's inverse is the covariance matrix of the parameters
  standardErrors <- sqrt(diag(solve(hess)))
  names(standardErrors) <- names(pars)
  return(standardErrors)
}

#' getStandardErrorsHW
#'
#' compute Huber White robust standard errors for a fitted model
#' @param model psydiff model
#' @param ... additional arguments to pass to optimHess
#' @return vector with standard errors
#' @export
getStandardErrorsHW <- function(model, ...){
  warning("Huber White Standard Errors are not working correctly!")
  pars <- getParameterValues(model)
  grad <- getGradients(model)
  hess <- optimHess(fn = psydiffFitInternal, par = pars,
                    model = model,
                    parsnames = names(pars),
                    ...)
  # Note: We are minimizing the negative log likelihood.
  # The Hessian is therefore the "observed Fisher Information"
  InformationInverse <- solve(hess)
  covMat <- InformationInverse%*%matrix(grad, ncol = 1)%*%matrix(grad, nrow = 1)%*%InformationInverse
  standardErrorsHW <- sqrt(diag(covMat))
  names(standardErrorsHW) <- names(pars)
  return(standardErrorsHW)
}
