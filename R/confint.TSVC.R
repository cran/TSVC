#' Confidence intervals for Varying Coefficient Trees
#'
#' @param object a fitted object of class \code{\link[TSVC]{TSVC}}.
#' @param ... further arguments passed to or from other methods.
confint <- function (object, ...) {
  UseMethod("confint", object)
}
#' Confidence intervals for Varying Coefficient Trees
#' 
#'@description
#'Construct parametric bootstrap percentile confidence intervals of effects of covariates that vary with the values of one or several effect modifiers. The basic method 
#'is described in Spuck et al. (2025). 
#'
#'@param object a fitted object of class \code{\link[TSVC]{TSVC}}. 
#'@param bootstrap_n the number of bootstrap samples to be drawn. 
#'@param alpha significance level \eqn{alpha} of the confidence interval.
#'@param post_pruning method to select the maximal number of splits; can be \code{"AIC"} or \code{"BIC"}. If \code{NULL} (default), no post-pruning is performed. 
#'@param splits_max maximal number of splits to be considered. If \code{post_pruning} is \code{NULL} it is ignored.
#'@param trace if \code{TRUE}, information about the estimation progress is printed.    
#'@param ... further arguments passed to or from other methods.  
#'
#'
#'@details 
#'The method is so far mainly put to the test for gaussian (family=gaussian) and binary (family=binary(link="logit")) outcome. It should be taken with care for 
#'differently scaled outcomes. 
#'
#'
#' @author 
#' Moritz Berger <Moritz.Berger@imbie.uni-bonn.de> \cr \url{https://www.imbie.uni-bonn.de/people/dr-moritz-berger/} 
#' 
#' Nikolai Spuck <Spuck@imbie.uni-bonn.de> \cr \url{https://www.imbie.uni-bonn.de/people/nikolai-spuck/}
#' 
#' @references 
#' Berger, M., G. Tutz and M. Schmid (2019). Tree-Structured Modelling of Varying Coefficients. Statistics and Computing 29, 217-229,
#' https://doi.org/10.1007/s11222-018-9804-8. 
#' 
#' Spuck, N., M. Schmid, M. Monin and M. Berger (2025). Confidence intervals for tree-structured varying coefficients. Computational Statistics and Data Analysis. 
#' 
#' @seealso 
#' \code{\link[TSVC]{TSVC}}, \code{\link[TSVC]{plot.TSVC}}, \code{\link[TSVC]{predict.TSVC}}, \code{\link[TSVC]{summary.TSVC}}
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
#' \dontrun{
#' fit1 <- TSVC(participation~income+age, data=sl, family=binomial(link="logit"), 
#' perm_test=FALSE, test_linear=FALSE, splits_max=3)
#' confint(fit1, bootstrap_n=500, alpha=0.05, trace=TRUE)
#' }
#' 
#' @method confint TSVC
#' @export
#' @importFrom stats as.formula na.omit rbinom rnorm sd
#' @importFrom tibble add_column
#' 
confint.TSVC <- function(object, 
                         bootstrap_n = 500, 
                         alpha = 0.05, 
                         post_pruning = NULL, 
                         splits_max = 5,
                         trace=FALSE,
                         ...){
  
  if(is.null(object$call$test_linear)){
    stop("test_linear needs to be FALSE in order to construct confidence intervals.")
  }
  if(object$call$test_linear){
    stop("test_linear needs to be FALSE in order to construct confidence intervals.")    
  }
  if(!is.null(post_pruning) && !post_pruning%in%c("AIC","BIC")){
    stop("post_pruning needs to be 'AIC' or 'BIC'")
  }
  data <- object$X
  response <- all.vars(object$call$formula)[1]
  data[,response] <- object$y
  call <- object$call
  n <- dim(data)[1]
  family <- attributes(object)$family
  
  origVars<- NULL
  only_effmod <- object$call$only_effmod
  if(is.null(only_effmod)){
    origVars <- names(data)[names(data) != response]
  }else{
    origVars <- names(data)[names(data) != response & names(data) != only_effmod]
  }
  vars <- names(data)[names(data) != response]
  form_base <- object$call$formula
  data_hi <- constructTSVCdata(object = object, data = data, out = response)
  data_h <- data_hi$data_h
  data_i <- data_hi$data_i
  data_h[, response] <- data[, response]
  if(is.null(only_effmod)){
    vars <- names(data_h)[names(data_h) != response]
  }else{
    vars <- names(data_h)[names(data_h) != response & names(data_h) != only_effmod]
  }
  form_h <- paste(vars, collapse = "+")
  form_h <- as.formula(paste(response, form_h, sep = "~"))
  model_h <- glm(form_h, data = data_h, family = family)
  bootstrappedEstimates <- matrix(rep(NA, length(model_h$coefficients) * bootstrap_n), ncol = length(model_h$coefficients))
  bootstrappedEstimates <- data.frame(bootstrappedEstimates)
  names(bootstrappedEstimates) <- names(model_h$coefficients)
  CI1 <- data.frame(Variable = names(model_h$coefficients), Coefficient = model_h$coefficients,
                    CI_L = rep(NA, length(model_h$coefficients)), CI_U = rep(NA, length(model_h$coefficients)))
  for(i in 1:bootstrap_n){
    if(trace){
      cat(".")
    }
    data_b <- data
    response_mu <- predict(object, X_new = data_b, type = "response")
    if(family$family == "gaussian"){
      data_b[,response] <- rnorm(n, response_mu, sd(object$model$residuals))
    }else{
      if(family$family == "binomial")
        data_b[,response] <- rbinom(n, 1, prob = response_mu) 
    }
    data_effects <- constructTSVCdata(object = object, data = data_b, out = response)$data_i
    if(is.null(post_pruning)){
      call[["data"]] <- data_b
      model_b <- eval(call)
      data_bhi <- constructTSVCdata(object = model_b, data = data_b, out = response)
    }else{
      models <- optimize_nSplits(object, data_b, post_pruning, splits_max)
      if(models$splits == 0){
        data_bhi <- list(data_h = data_b, data_i = data_b)
      }else{
        if(models$splits == splits_max){
          data_bhi <- constructTSVCdata(object = models$model_full, data = data_b, out = response)
        }else{
          data_bhi <- constructTSVCdata_s(object = models$model_full, data = data_b, splits = models$splits, out = response)
        }
      }
    }
    data_bh <- data_bhi$data_h
    data_bi <- data_bhi$data_i
    data_bh[, response] <- data_b[, response]
    if(is.null(only_effmod)){
      vars <- names(data_bh)[names(data_bh) != response]
    }else{
      vars <- names(data_bh)[names(data_bh) != response & names(data_bh) != only_effmod]
    }
    form_bh <- paste(vars, collapse = "+")
    form_bh <- as.formula(paste(response, form_bh, sep = "~"))
    model_bh <- glm(form_bh, data = data_bh, family = family)
    for(var in origVars){
      data_effects[, var] <- rep(NA, n)
    }
    for(nd in names(model_bh$coefficients)[-1]){
      var <- nd
      Effect <- NA
      if(!(nd %in% origVars)){
        var <- origVars[sapply(origVars, function(ov){endsWith(nd, ov)})]
        data_effects[data_bi[, nd] == 1, var] <- model_bh$coefficients[nd]
      }else{
        data_effects[, nd] <- model_bh$coefficients[nd]
      }
    }
    for(nd in names(model_h$coefficients)[-1]){
      if(!(nd %in% origVars)){
        bootstrappedEstimates[i, nd] <- mean(data_effects[data_effects[, nd] == 1, origVars[sapply(origVars, function(ov){endsWith(nd, ov)})]])
      }else{
        bootstrappedEstimates[i, nd] <- mean(data_effects[, nd])
      }
    }
  }
  for(nd in names(bootstrappedEstimates)[-1]){
    CI1[CI1$Variable == nd, "CI_L"] <- quantile(bootstrappedEstimates[, nd], probs = alpha/2)
    CI1[CI1$Variable == nd, "CI_U"] <- quantile(bootstrappedEstimates[, nd], probs = 1 - alpha/2)
  }
  orig_names   <- names(coef(object$model))
  where_names  <- match(round(coef(model_h), 4), round(coef(object$model), 4))
  CI1$Variable <- orig_names[where_names]
  
  rownames(CI1) <- NULL
  class(CI1) <- append(class(CI1),"confint.TSVC")
  return(CI1)
}



