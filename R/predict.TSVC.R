#' Prediction from Varying Coefficient Trees
#' 
#' @description 
#' Obtains predictions from a fitted TSVC object. 
#' 
#' @param object a fitted object of class \code{\link[TSVC]{TSVC}}.
#' @param X_new optionally, data frame of class \code{\link{data.frame}} which contains the variables with which to predict. If \code{NULL}, the fitted 
#' linear predictors are use. 
#' @param ... further arguments passed to \code{\link[stats]{predict.glm}}. 
#' 
#' @details 
#' \code{predict.TSVC} is a wrapper function of \code{predict.glm}, which obtains predictions for objects of class \code{\link[stats]{glm}}. Further
#' arguments can be passed to \code{predict.glm} via the '...'-argument. 
#' 
#' @author Moritz Berger <moritz.berger@imbie.uni-bonn.de> \cr \url{https://www.imbie.uni-bonn.de/people/dr-moritz-berger/}
#' 
#' @references 
#' Berger, M., G. Tutz and M. Schmid (2019). Tree-Structured Modelling of Varying Coefficients. Statistics and Computing 29, 217-229,
#' https://doi.org/10.1007/s11222-018-9804-8. 
#' 
#' @seealso \code{\link[TSVC]{TSVC}}, \code{\link[TSVC]{plot.TSVC}}, \code{\link[TSVC]{summary.TSVC}}
#' 
#' @examples 
#' # Swiss Labour Market 
#' library(AER)
#' data("SwissLabor")
#' 
#' # recode factors
#' sl <- SwissLabor
#' sl$participation <- as.numeric(sl$participation)-1
#' sl$foreign       <- as.numeric(sl$foreign)-1
#' 
#' X_new <- data.frame("income"=c(10,12), "age"=c(4.5,5.8))
#' 
#' \dontrun{
#' fit1 <- TSVC(participation~income+age, data=sl, family=binomial(link="logit"), 
#'              nperm=1000, trace=TRUE)
#' predict(fit1, X_new, type="response")
#' }
#' 
#' 
#' @method predict TSVC 
#' @export
#' @importFrom stats predict 

predict.TSVC  <- function(object, 
                          X_new=NULL, 
                          ...){

  DM_kov <- object$X
  sb     <- object$sb
  nvar   <- ncol(DM_kov)
  
  if(is.null(X_new)){
    pred <- predict(object$model, ...)
  } else{
    
    if("Intercept"%in%names(DM_kov)){
      X_new <- as.data.frame(cbind(X_new,"Intercept"=1))
    }
    if(!is.null(sb)){
      for(s in 1:length(sb)){
        svar <- paste0(sb[s],"_em")
        X_new[,svar] <- X_new[,sb[s]]
      }
    }
    
    for(i in 1:ncol(X_new)){
      if(is.factor(X_new[,i])){
        if(nlevels(X_new[,i])==2){
          X_new[,i] <- as.numeric(X_new[,i])-1
        } 
      }
    }
    
    if(!is.null(names(DM_kov))){
      var_names <- names(DM_kov)
    } else{
      var_names <- paste0("x",1:nvar)
    }
    n_quantile     <- ifelse(is.null(object$call$n_quantile), 20, object$call$n_quantile)
    ordered_values <- lapply(DM_kov, ord_values, n_quantile = n_quantile)
    n_levels       <- sapply(ordered_values,length)
    thresholds     <- lapply(ordered_values,thresh)
    var_list       <- attributes(object)$vl
    X_new          <- X_new[,var_names]
    
    to_build <- c(1:nvar)[sapply(1:nvar, function(j) length(var_list[[j]])>0)]
    
    design_upper <- lapply(to_build,function(j) designlist(X_new,var_list[[j]],j,thresholds,var_names))
    design_upper <- lapply(1:length(design_upper), function(j) do.call(cbind, design_upper[[j]]))
    design_upper <- do.call(cbind,design_upper)
    design_lower <- lapply(to_build,function(j) designlist(X_new,var_list[[j]],j,thresholds,var_names, upper=FALSE))
    design_lower <- lapply(1:length(design_lower), function(j) do.call(cbind, design_lower[[j]]))
    design_lower <- do.call(cbind,design_lower)
    
    X_new <- as.data.frame(cbind(X_new, design_upper, design_lower))[,all.vars(formula(object$model))[-1], drop=FALSE]
    
    pred  <- predict(object$model, newdata=X_new, ...)
  }
  
  return(pred)
}
