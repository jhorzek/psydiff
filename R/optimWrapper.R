#' psydiffOptimNelderMead
#'
#' wrapper for Nelder-Mead optimizer from optim
#' @param model psydiff model object (must have been compiled before calling psydiffOptimNelderMead)
#' @param ... additional parameters to pass to optim
#' @export
psydiffOptimNelderMead <- function(model, ...){
  startingValues <- psydiff::getParameterValues(model)

  ffit <- function(pars, parsnames, model){
    setParameterValues(model$pars$parameterTable, parameterLabels = parsnames, parameterValues = pars)
    fit <- fitModel(model)
    return(sum(fit$m2LL))
  }

  out <- optim(par = startingValues, fn = ffit, method = "Nelder-Mead",
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

  ffit <- function(pars, parsnames, model){
    setParameterValues(model$pars$parameterTable, parameterLabels = parsnames, parameterValues = pars)
    fit <- fitModel(model)
    return(sum(fit$m2LL))
  }

  out <- optim(par = startingValues, fn = ffit, method = "BFGS",
               model = model,
               parsnames = names(startingValues), ...)
  return(out)
}

#' psydiffOptimx
#'
#' wrapper for optimx optimization
#' @param model psydiff model object (must have been compiled before calling psydiffOptimx)
#' @param ... additional parameters to pass to optimx
#' @export
#'
psydiffOptimx <- function(model, method=c("Nelder-Mead","BFGS"), ...){
  startingValues <- psydiff::getParameterValues(model)

  ffit <- function(pars, parsnames, model){
    setParameterValues(model$pars$parameterTable, parameterLabels = parsnames, parameterValues = pars)
    fit <- fitModel(model)
    return(sum(fit$m2LL))
  }

  out <- optimx::optimx(par = startingValues, fn = ffit, method = method,
                        model = model,
                        parsnames = names(startingValues), ...)
  return(out)
}

