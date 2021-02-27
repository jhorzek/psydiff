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
psydiffOptimBFGS <- function(model, fn = psydiffFitInternal, gr = psydiffGradientsInternal, ...){
  startingValues <- psydiff::getParameterValues(model)

  out <- optim(par = startingValues, fn = fn, gr = gr,
               method = "BFGS",
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
psydiffOptimx <- function(model, method=c("Nelder-Mead","BFGS", "nlm", "nlminb"), fn = psydiffFitInternal, gr = psydiffGradientsInternal,...){
  startingValues <- psydiff::getParameterValues(model)

  out <- optimx::optimx(par = startingValues, fn = fn, gr = gr, method = method,
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

  out <- DEoptim::DEoptim(fn = psydiffFitInternal, lower = lower, upper = upper,
                          model = model,
                          parsnames = names(startingValues),
                          ...)
  return(out)
}

#' psydiffSubplex
#'
#' wrapper for subplex optimization
#' @param model psydiff model object (must have been compiled before calling psydiffSubplex)
#' @param ... additional parameters to pass to subplex
#' @import subplex
#' @export
#'
psydiffSubplex <- function(model, ...){
  startingValues <- psydiff::getParameterValues(model)

  out <- subplex::subplex(fn = psydiffFitInternal,
                          par = startingValues,
                          model = model,
                          parsnames = names(startingValues),
                          ...)
  return(out)
}

#' psydiffNLopt
#'
#' wrapper for NLopt optimization
#' @param model psydiff model object (must have been compiled before calling psydiffNLopt)
#' @param ... additional parameters to pass to NLopt
#' @import nloptr
#' @export
#'
psydiffNLopt <- function(model, eval_f = psydiffFitInternal,...){
  startingValues <- psydiff::getParameterValues(model)
  local_opts <-list("algorithm"="NLOPT_LD_MMA",
                    "xtol_rel"=1.0e-15)
  opts <-list("algorithm"="NLOPT_LD_LBFGS_NOCEDAL",
              "xtol_rel"=1.0e-15,
              "maxeval"=500,
              "local_opts"= local_opts,
              "print_level"=1)
  out <- nloptr::nloptr(eval_f = eval_f,
                        x0 = startingValues,
                        opts = opts,
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

#' psydiffGradientsInternal
#'
#' internal function for optimization
#' @param pars vector with parameter values
#' @param parsnames vector with parameter names
#' @param model psydiff model object
#' @export
#'
psydiffGradientsInternal <- function(pars, parsnames, model){
  setParameterValues(model$pars$parameterTable, parameterLabels = parsnames, parameterValues = pars)
  if(!is.null(model$cl)){
    invisible(capture.output(gradients <- try(getGradientsParallel(model),
                                        silent = TRUE),
                             type = "message"))
  }else{
    invisible(capture.output(gradients <- try(getGradients(model),
                                        silent = TRUE),
                             type = "message"))
  }


  if(any(class(gradients) == "try-error")){
    gradients <- rep(1, length(pars))
    names(gradients) <- parsnames
    return(gradients)
  }
  if(anyNA(gradients)){
    gradients[is.na(gradients)] <- 1
    return(gradients)
  }
  return(gradients)
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
  # The Hessian is therefore .5*"observed Fisher Information"
  # and 2 times it's inverse is the covariance matrix of the parameters
  standardErrors <- sqrt(diag(2*solve(hess)))
  names(standardErrors) <- names(pars)
  return(standardErrors)
}


