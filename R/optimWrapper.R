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
