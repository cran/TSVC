#' @export 

TSVC.default <- function(formula, 
                         data,
                         family=gaussian,
                         alpha=0.05, 
                         nperm=1000, 
                         nodesize_min=5, 
                         bucket_min=1,
                         depth_max=NULL,
                         splits_max=NULL,
                         perm_test=TRUE,
                         test_linear=FALSE,
                         gpd_approx=FALSE,
                         effmod=NULL,
                         notmod=NULL, 
                         only_effmod=NULL,
                         smooth=NULL,
                         split_intercept=FALSE,
                         sb_slope=NULL,
                         sb_slope_c=FALSE,
                         n_quantile=20,
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
  comp   <- specification(formula, data, effmod, notmod, only_effmod, smooth, split_intercept, sb_slope)
  y      <- comp$y
  DM_kov <- comp$DM_kov
  effmod <- comp$effmod
  notmod <- comp$notmod
  exclude<- comp$exclude
  smooth <- comp$smooth
  
  # model fit 
  output <- effmodTree(y, DM_kov, family, alpha, nperm, nodesize_min, bucket_min, depth_max, splits_max, perm_test, test_linear, 
                       gpd_approx, effmod, notmod, exclude, smooth, split_intercept, sb_slope, sb_slope_c, n_quantile, trace, ...)
  coefficients <- list("beta_constant"= output$beta_noeffmod,
                       "beta_varying" = output$beta_effmod)
  
  # return 
  to_return <- list("splits"=output$splits,
                    "coefficients"=coefficients,
                    "pvalues"=output$pvalues,
                    "pvalues_linear"=output$pvalues_linear,
                    "devs"=output$devs,
                    "crits"=output$crits,
                    "y"=y,
                    "X"=DM_kov,
                    "sb"=sb_slope,
                    "model"=output$model,
                    "all_models"=output$all_models,
                    "call"=match.call())
  
  class(to_return) <- "TSVC"
  attr(to_return, "vl") <- output$vl 
  attr(to_return, "family") <- family
  return(to_return)

}
