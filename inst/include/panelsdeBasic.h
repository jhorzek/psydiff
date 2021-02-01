#ifndef PANELSDEBASIC_H
#define PANELSDEBASIC_H

inline Rcpp::List panelSDE_C(arma::mat observed, Rcpp::DataFrame parameterTable){
  Rcpp::List panelSDEModel = Rcpp::List::create(Rcpp::_["observed"] = observed ,
                                                Rcpp::_["parameterTable"] = parameterTable);
  return(panelSDEModel);
}

inline void setParameterValues_C(Rcpp::DataFrame &parameterTable,
                                 Rcpp::NumericVector parameterValues,
                                 Rcpp::StringVector parameterLabels){
  if(parameterValues.length() != parameterLabels.length()){
    Rcpp::stop("Length of parameterValues and parameterLabels does not match");
  }
  // extract relevant elements from model
  Rcpp::NumericVector values = parameterTable["value"];
  Rcpp::StringVector labels = parameterTable["label"];
  Rcpp::LogicalVector changed = parameterTable["changed"];

  bool wasfound;
  for(int i = 0; i < parameterLabels.length(); i++){
    wasfound = false;
    for(int j = 0; j < labels.length(); j++){
      if(labels(j) == parameterLabels(i)){
        wasfound = true;
        // only change if the value changed
        if(values(j) != parameterValues(i)){
          values(j) = parameterValues(i);
          // set flag to changed
          changed(j) = true;
        }
      }
    }
    if(!wasfound){
      Rcpp::warning("The following parameter was not found in the parameterTable: " + parameterLabels(i));
    }
  }
}

// returns the parameter values of the model
inline Rcpp::NumericVector getParameterValues_C(const Rcpp::List panelSDEModel) {
  // extract relevant elements from model
  Rcpp::List pars = panelSDEModel["pars"];
  Rcpp::DataFrame parameterTable = Rcpp::as<Rcpp::DataFrame>(pars["parameterTable"]);

  Rcpp::StringVector parameterLabels = parameterTable["label"];
  Rcpp::StringVector uniqueParameterLabels = unique(parameterLabels);
  Rcpp::NumericVector parameterValuesRep = parameterTable["value"]; // values with repeated elements

  Rcpp::NumericVector paramterValues (uniqueParameterLabels.length(),-9999.99);

  for(int i = 0; i < uniqueParameterLabels.size(); i++){

    for(int j = 0; j < parameterLabels.length(); j++)
      if(uniqueParameterLabels(i) == parameterLabels(j)){
        paramterValues(i) = parameterValuesRep(j);
        break;
      }
  }
  paramterValues.names() = uniqueParameterLabels;

  return(Rcpp::clone(paramterValues));
}

// setParameterValues only defines the values in parameterTable
// when fitting the model, the drift, etc. have to be filled the
// the individual parameter values
inline void setParameterTable_C(const Rcpp::DataFrame &parameterTable, Rcpp::List &parameterList,
                                int person){
  // extract relevant elements from model
  Rcpp::NumericVector persons = parameterTable["person"];
  // select the rows of the specified person
  Rcpp::LogicalVector personSelector = (persons==person);

  if(is_false(any(personSelector))){
    Rcpp::stop("Person " + std::to_string(person) + " not found in data set.");
  }

  Rcpp::StringVector target = parameterTable["target"];
  Rcpp::StringVector personTarget = target[personSelector];

  Rcpp::NumericVector value = parameterTable["value"];
  Rcpp::NumericVector personValue = value[personSelector];

  Rcpp::StringVector targetClass = parameterTable["targetClass"];
  Rcpp::StringVector personTargetClass = targetClass[personSelector];

  Rcpp::NumericVector rw = parameterTable["row"];
  Rcpp::NumericVector personRow = rw[personSelector];

  Rcpp::NumericVector cl = parameterTable["col"];
  Rcpp::NumericVector personCol = cl[personSelector];

  Rcpp::LogicalVector changed = parameterTable["changed"];
  Rcpp::LogicalVector personChanged = changed[personSelector];

  Rcpp::StringVector parameterListNames = parameterList.names();

  for(int i = 0; i < personTarget.length(); i++){
    if(!personChanged(i)){continue;}
    bool targetFound= false;
      for(int j = 0; j < parameterListNames.length(); j++){
        if(parameterListNames(j) == personTarget(i)){
          targetFound = true;
          // change parameter value
          Rcpp::String sel = parameterListNames(j);
          if(personTargetClass(i) == "matrix"){
            Rcpp::NumericMatrix targetMatrix = parameterList[sel];
            targetMatrix(personRow(i), personCol(i)) = personValue(i);
          }else if(personTargetClass(i) == "double"){
            parameterList[sel] = personValue(i);
          }else if(personTargetClass(i) == "vec"){
            Rcpp::NumericVector targetVec = parameterList[sel];
            targetVec(personRow(i)) = personValue(i);
          }else{
            Rcpp::stop("Parameter class not found");
          }
          break;
        }

      }
      if(!targetFound){
        Rcpp::warning("Target " + personTarget(i) + " not found");
      }
  }
}


#endif
