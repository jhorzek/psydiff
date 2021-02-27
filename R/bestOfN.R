bestOfN <- function(model, minVal = -5, maxVal = 5, sd = .5, nstart = 100){
  pars <- getParameterValues(model)
  model2 <- model
  if(length(minVal) == 1){
    minVal <- rep(minVal, length(pars))
  }
  if(length(maxVal) == 1){
    maxVal <- rep(maxVal, length(pars))
  }
  if(length(sd) == 1){
    sd <- rep(sd, length(pars))
  }
  if(!length(minVal) == length(pars)){
    stop("Length of minVal does not correspond to the number of parameters in the model")
  }
  if(!length(maxVal) == length(pars)){
    stop("Length of maxVal does not correspond to the number of parameters in the model")
  }
  if(!length(sd) == length(pars)){
    stop("Length of sd does not correspond to the number of parameters in the model")
  }

  bestpars <- c()
  bestfit <- Inf
  for(n in seq_len(nstart)){
    currentPars <- pars
    for(p in seq_len(length(pars))){
      sampleFrom <- seq(minVal[p], maxVal[p], length.out = 1000)
      currentPars[p] <- sample(x = sampleFrom, size = 1, prob = dnorm(sampleFrom, mean = pars[p], sd = sd[p]))
    }
    setParameterValues(parameterTable = model2$pars$parameterTable,
                       parameterValues = currentPars,
                       parameterLabels = names(currentPars))
    f <- fitModel(model2)
    if(is.null(f$m2LL)){next}
    if(is.finite(sum(f$m2LL)) && (sum(f$m2LL) < bestfit)){
      bestfit <- sum(f$m2LL)
      bestpars <- currentPars
    }
  }
  return(bestpars)
}
