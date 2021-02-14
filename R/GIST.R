#' GIST
#'
#' General Iterative Shrinkage and Thresholding Algorithm based on Gong, P., Zhang, C., Lu, Z., Huang, J., & Ye, J. (2013). A General Iterative Shrinkage and Thresholding Algorithm for Non-convex Regularized Optimization Problems. In S. Dasgupta & D. McAllester (Eds.), Proceedings of Machine Learning Research (PMLR; Vol. 28, Issue 2, pp. 37--45). PMLR. http://proceedings.mlr.press
#'
#' GIST minimizes a function of form f(theta) = l(theta) + g(theta), where l is the likelihood and g is a penalty function. Various penalties are supported, however currently only lasso and adaptive lasso are implemented.
#'
#' @param model model
#' @param startingValues named vector with starting values
#' @param lambda penalty value
#' @param adaptiveLassoWeights named vector with adaptive lasso weights
#' @param regularizedParameters named vector of regularized parameters
#' @param eta if the current step size fails, eta will decrease the step size. Must be > 1
#' @param sig GIST: sigma value in Gong et al. (2013). Sigma controls the inner stopping criterion and must be in (0,1). Generally, a larger sigma enforce a steeper decrease in the regularized likelihood while a smaller sigma will result in faster acceptance of the inner iteration.
#' @param initialStepsize initial stepsize to be tried in the outer iteration
#' @param stepsizeMin Minimal acceptable step size. Must be > 0. A larger number corresponds to a smaller step from one to the next iteration. All step sizes will be computed as described by Gong et al. (2013)
#' @param stepsizeMax Maximal acceptable step size. Must be > stepsizeMin. A larger number corresponds to a smaller step from one to the next iteration. All step sizes will be computed as described by Gong et al. (2013)
#' @param GISTLinesearchCriterion criterion for accepting a step. Possible are 'monotone' which enforces a monotone decrease in the objective function or 'non-monotone' which also accepts some increase.
#' @param GISTNonMonotoneNBack in case of non-monotone line search: Number of preceding regM2LL values to consider
#' @param maxIter_out maximal number of outer iterations
#' @param maxIter_in maximal number of inner iterations
#' @param break_outer stopping criterion for the outer iteration.
#' @param numDeriv boolean should numDeriv package be used for derivatives?
#' @param verbose set to 1 to print additional information and plot the convergence
#' @export
GIST <- function(model, startingValues, lambda, adaptiveLassoWeights, regularizedParameters,
                 eta = 1.5, sig = .2, initialStepsize = 1, stepsizeMin = 0, stepsizeMax = 999999999,
                 GISTLinesearchCriterion = "monotone", GISTNonMonotoneNBack = 5,
                 maxIter_out = 100, maxIter_in = 100,
                 break_outer = .00000001,
                 numDeriv = FALSE,
                 verbose = 0, silent = FALSE){
  # iteration counter
  k_out <- 1
  convergence <- TRUE

  parameterNames <- names(startingValues)
  adaptiveLassoWeightsMatrix <- diag(adaptiveLassoWeights[parameterNames])
  rownames(adaptiveLassoWeightsMatrix) <- parameterNames
  colnames(adaptiveLassoWeightsMatrix) <- parameterNames

  # set parameters
  psydiff::setParameterValues(parameterTable = model$pars$parameterTable,
                              parameterValues = startingValues,
                              parameterLabels = names(startingValues))
  invisible(capture.output(out1 <- try(fitModel(model), silent = T), type = "message"))
  if(any(class(out1) == "try-error")  ||
     anyNA(out1$m2LL)){
    stop("Infeasible starting values in GIST.")
  }

  parameters_km1 <- NULL
  parameters_k <- getParameterValues(model)
  parameterNames <- names(parameters_k)
  gradients_km1 <- NULL
  if(numDeriv){
    gradients_k <- try(numDeriv::grad(fitf, x = parameters_k, method = "Richardson", model = model))
  }else{
    gradients_k <- try(getGradients(model))
  }

  m2LL_k <- sum(out1$m2LL)

  regM2LL_k <- m2LL_k + regCtsem::exact_getRegValue(lambda = lambda,
                                                    theta = parameters_k,
                                                    regIndicators = regularizedParameters,
                                                    adaptiveLassoWeights = adaptiveLassoWeights)
  regM2LL <- rep(NA, maxIter_out)
  regM2LL[1] <- regM2LL_k

  if(verbose == 2){
    convergencePlotValues <- matrix(NA, nrow = length(parameters_k), ncol = maxIter_out, dimnames = list(parameterNames, 1:maxIter_out))
  }
  resetStepSize <- TRUE
  resetIteration <- -1
  while(k_out < maxIter_out){
    k <- 0
    # set initial step size
    if(is.null(parameters_km1)){
      stepsize <- initialStepsize
    }else{
      x_k <- parameters_k - parameters_km1
      y_k <- gradients_k - gradients_km1

      stepsize <- (t(x_k)%*%y_k)/(t(x_k)%*%x_k)
      # sometimes this step size is extremely large and the algorithm converges very slowly
      # we found that in these cases it can help to reset the stepsize
      # this will be done randomly here:
      if(runif(1,0,1) > .9){
        stepsize <- initialStepsize
      }

      if(is.na(stepsize)){
        if(!silent){
          warning(paste0("Outer iteration ", k_out, ": NA or infinite step size..."))}
        break
      }
      if(is.infinite(stepsize)){
        if(!silent){
          warning(paste0("Outer iteration ", k_out, ": NA or infinite step size..."))}
        break
      }
      if(stepsize < stepsizeMin){

        if(!resetStepSize){
          if(!silent){warning(paste0("Outer iteration ", k_out, ": Stepsize below specified minimum..."))}
          break
        }
        stepsize <- initialStepsize
      }
      if(stepsize > stepsizeMax){
        if(!resetStepSize){
          if(!silent){warning(paste0("Outer iteration ", k_out, ": Stepsize above specified maximum..."))}
          break
        }
        stepsize <- initialStepsize
      }

    }


    # inner iteration
    while(k < maxIter_in){
      u_k <- parameters_k - gradients_k/as.vector(stepsize)
      parameters_kp1 <- rep(NA, length(parameters_k))
      names(parameters_kp1) <- parameterNames
      for(i in seq_len(length(parameters_kp1))){
        parameterName <- parameterNames[i]
        lambda_i <- lambda*adaptiveLassoWeights[parameterName]
        if(parameterName %in% regularizedParameters){
          # update parameter i with lasso
          parameters_kp1[parameterName] <- sign(u_k[parameterName])*max(c(0,abs(u_k[parameterName]) - lambda_i/stepsize))
        }else{
          parameters_kp1[parameterName] <- u_k[parameterName]
        }
      }

      psydiff::setParameterValues(parameterTable = model$pars$parameterTable,
                                  parameterValues = parameters_kp1,
                                  parameterLabels = names(parameters_kp1))
      invisible(capture.output(out1 <- try(fitModel(model), silent = T), type = "message"))
      if(any(class(out1) == "try-error") || anyNA(out1$m2LL)){

        # update step size
        stepsize <- eta*stepsize

        # update iteration counter
        k <- k+1

        # skip rest
        next
      }

      m2LL_kp1 <- sum(out1$m2LL)

      regM2LL_kp1 <- m2LL_kp1 + regCtsem::exact_getRegValue(lambda = lambda,
                                                            theta = parameters_kp1,
                                                            regIndicators = regularizedParameters,
                                                            adaptiveLassoWeights = adaptiveLassoWeights)

      if(verbose == 2){
        cat(paste0("\r",
                   "###### [", sprintf("%*d", 3, k_out), "-", sprintf("%*d", 3, k),
                   " | regM2LL:  ", sprintf('%.3f',regM2LL_kp1),
                   " | zeroed: ", sprintf("%*d", 3, sum(parameters_kp1[regularizedParameters] == 0)),
                   " | stepsize: ", sprintf("%.3f", stepsize),
                   " ######"
        )
        )
      }

      # break if line search condition is satisfied
      if(GISTLinesearchCriterion == "monotone"){
        breakCriterion <- regM2LL_kp1 <= regM2LL_k - (sig/2) * stepsize * sum((parameters_k - parameters_kp1)^2)
      }else if(GISTLinesearchCriterion == "non-monotone"){
        nBack <- max(1,k_out-GISTNonMonotoneNBack)
        breakCriterion <- regM2LL_kp1 <= max(regM2LL[nBack:k_out]) - (sig/2) * stepsize * sum((parameters_k - parameters_kp1)^2)
      }else{
        stop("Unknown GISTLinesearchCriterion. Possible are monotone and non-monotone.")
      }
      if(breakCriterion){
        #print("breaking inner")
        break
      }

      # update step size
      stepsize <- eta*stepsize

      # update iteration counter
      k <- k+1
    }

    if(k == maxIter_in){
      if(!silent){
        warning("Maximal number of inner iterations used by GIST. Consider increasing the number of inner iterations.")}
    }

    if(numDeriv){
      out3 <- try(numDeriv::grad(fitf, x = parameters_k, method = "Richardson", model = model))
    }else{
      out3 <- try(getGradients(model))
    }

    if(any(class(out3) == "try-error") ||
       anyNA(out3)){
      if(!silent){
        convergence <- FALSE
        stop("No gradients in GIST")}

    }

    gradients_kp1 <- out3

    # update parameters for next iteration
    parameters_km1 <- parameters_k
    parameters_k <- parameters_kp1
    gradients_km1 <- gradients_k
    gradients_k <- gradients_kp1
    m2LL_k <- m2LL_kp1
    regM2LL_k <- regM2LL_kp1

    # break outer loop if stopping criterion is satisfied
    k_out <- k_out + 1
    regM2LL[k_out] <- regM2LL_k
    if(verbose == 1){
      plot(x=1:maxIter_out, y = regM2LL, xlab = "iteration", ylab = "f(theta)", type = "l", main = "Convergence Plot")
      cat(paste0("\r",
                 "## [", sprintf("%*d", 3, k_out),
                 "] m2LL: ", sprintf('%.3f',m2LL_k),
                 " | regM2LL:  ", sprintf('%.3f',regM2LL_k),
                 " | zeroed: ", sprintf("%*d", 3, sum(parameters_k[regularizedParameters] == 0)),
                 " ##"
      )
      )
    }

    if(verbose == 2){
      convergencePlotValues[,k_out] <- parameters_k
      matplot(x=1:maxIter_out, y = t(convergencePlotValues), xlab = "iteration", ylab = "value", type = "l", main = "Convergence Plot")
    }


    breakOuter <- (sqrt(sum((parameters_k - parameters_km1)^2))/sqrt(sum((parameters_km1)^2))) < break_outer

    if(breakOuter){
      break
    }

  }
  if(is.na(regM2LL_k) | is.infinite(regM2LL_k) |
     anyNA(parameters_k) | any(is.infinite(parameters_k)) |
     anyNA(gradients_k) | any(is.infinite(gradients_k))){
    convergence <- FALSE
  }

  if(k_out == maxIter_out){
    if(!silent){
      warning("Maximal number of outer iterations used by GIST. Consider increasing the number of outer iterations.")}
  }

  adaptiveLassoWeightsMatrix <- diag(adaptiveLassoWeights[parameterNames])
  rownames(adaptiveLassoWeightsMatrix) <- parameterNames
  colnames(adaptiveLassoWeightsMatrix) <- parameterNames
  parameterMatrix <- matrix(parameters_k[parameterNames], ncol = 1)
  rownames(parameterMatrix) <- parameterNames
  gradientMatrix <- matrix(gradients_k[parameterNames], ncol = 1)
  rownames(gradientMatrix) <- parameterNames
  subgradients <- regCtsem::exact_getSubgradients(theta = parameterMatrix, jacobian = gradientMatrix,
                                                  regIndicators = regularizedParameters, lambda = lambda,
                                                  lineSearch = NULL, adaptiveLassoWeightsMatrix = adaptiveLassoWeightsMatrix)

  if(any(abs(subgradients)>1)){
    # Problem: gradients extremely dependent on the epsilon in the approximation
    if(!silent){
      #warning("Model did not reach convergence")
    }
  }

  return(list("type" = "psydiff",
              "model" = model,
              "regM2LL" = regM2LL,
              "convergence" = convergence))

}

fitf <- function(pars, model){
  setParameterValues(model$pars$parameterTable, parameterValues = pars, parameterLabels = names(pars))
  fit <- fitModel(model)
  return(sum(fit$m2LL))
}

