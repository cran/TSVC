#' @export 

TSVC.default <- function(formula, 
                         data,
                         family=gaussian,
                         alpha=0.05, 
                         nperm=1000, 
                         effmod=NULL,
                         notmod=NULL, 
                         only_effmod=NULL,
                         smooth=NULL,
                         split_intercept=FALSE,
                         trace=FALSE,
                         ...){
  
  # check input 
  if(missing(formula)){
    stop("Argument formula is missing with no default.")
  }
  if(missing(data)){
    stop("Argument data is missing with no default.")
  }
  
  # check family 
  if (is.character(family)) 
    family <- get(family, mode = "function", envir = parent.frame())
  if (is.function(family)) 
    family <- family()
  if (is.null(family$family)) {
    print(family)
    stop("'family' not recognized")
  }
  
  # predefinition 
  comp   <- specification(formula, data, effmod, notmod, only_effmod, smooth, split_intercept)
  y      <- comp$y
  DM_kov <- comp$DM_kov
  effmod <- comp$effmod
  notmod <- comp$notmod
  exclude<- comp$exclude
  smooth <- comp$smooth
  
  # model fit 
  output <- effmodTree(y, DM_kov, family, alpha, nperm, effmod, notmod, exclude, smooth, split_intercept, trace, ...)
  coefficients <- list("beta_constant"= output$beta_noeffmod,
                       "beta_varying" = output$beta_effmod)
  
  # return 
  to_return <- list("splits"=output$splits,
                    "coefficients"=coefficients,
                    "pvalues"=output$pvalues,
                    "devs"=output$devs,
                    "crits"=output$crits,
                    "y"=y,
                    "X"=DM_kov,
                    "model"=output$model,
                    "call"=match.call())
  
  class(to_return) <- "TSVC"
  attr(to_return, "vl") <- output$vl 
  return(to_return)

}
