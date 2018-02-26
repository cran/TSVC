#' Tree-Structured Modelling of Varying Coefficients
#' 
#' @description 
#' A function to fit tree-structured varying coefficient (TSVC) models. By recursive splitting the method allows to simultaneously detect covariates 
#' with varying coefficients and the effect modifiers that induce varying coefficients if they are present. The basic method is described in 
#' Berger, Tutz and Schmid (2018). 
#' 
#' @param formula object of class \code{\link{formula}}: a symbolic description of the (linear) model to be fit. See also details. 
#' @param data data frame of class \code{\link{data.frame}} containing the variables in the model.
#' @param family a description of the error distribution and link function to be used in the model (as for \code{\link{glm}}). 
#' This can be a character string naming a family function, a family function or the result of a call to a family function. 
#' See \code{\link{family}} for details of family functions. 
#' @param alpha significance level \eqn{alpha} for the permutation tests. 
#' @param nperm number of permutations used for the permutation tests. 
#' @param effmod optional vector of covariates that serve as effect modifier. If \code{NULL} (default), all covariates are considered as potential effect modifiers. 
#' @param notmod optional list of class \code{\link{list}} containing pairs of covariate/effect modifier that are not considered as candidates for splitting during iteration.
#' If \code{NULL} (default), all combinations of covariates and potential effect modifiers are considered for splitting. 
#' @param trace if \code{TRUE}, information about the estimation progress is printed. 
#' @param x object of class \code{TSVC}.
#' @param ... further arguments passed to or from other methods. 
#' 
#' @details 
#' A typical \link{formula} has the form \code{response ~ covariates}, where \code{response} is the name of the response variable and \code{covariates} is a
#' series of variables that are incorporated in the model. 
#' 
#' With p covariates, \code{TSVC} expects a formula of the form \eqn{y ~ x_1+...+x_p}. If no further specifications are made (\code{effmod=NULL} and \code{notmod=NULL}) it is assumed that each covariate \eqn{x_j, j = {1,...,p}} 
#' can be modified by all the other variables \eqn{x_m, m = {1,...,p} \ j}. 
#' 
#' Remark: Significance of each split is verified by permutation tests. The result of the permutation tests 
#' can strongly depend on the number of permutations \code{nperm}.
#' 
#' @return 
#' Object of class \code{"TSVC"}. An object of class \code{"TSVC"} is a list containing the following components: 
#' 
#' \item{splits}{matrix with detailed information about all executed splits during the fitting process.}
#' \item{coefficients}{list of estimated coefficients for covariates with and without varying coefficients (including a non-varying intercept).}
#' \item{pvalues}{p-values of each permuation test during the fitting process.}
#' \item{devs}{maximal value statistics \eqn{T_m} of the selected effect modifier in each iteration during the fitting process.}
#' \item{crit}{critical values of each permutation test during the fitting process.}
#' \item{y}{response vector.}
#' \item{X}{matrix of all the variables (covariates and effect modifiers) for model fitting.}
#' \item{model}{internaly fitted model in the last iteration of class \code{\link{glm}}.}
#' 
#' @author 
#' Moritz Berger <Moritz.Berger@imbie.uni-bonn.de> \cr \url{http://www.imbie.uni-bonn.de/personen/dr-moritz-berger/}
#' 
#' @references 
#' Berger, M., G. Tutz and M. Schmid (2018). Tree-Structured Modelling of Varying Coefficients. Statistics and Computing. Under review. 
#' 
#' Hastie, T. and R. Tibshirani (1993). Varying-coefficient models. Journal of the Royal Statistical Society B 55, 757-796.
#' 
#' Hothorn T., K. Hornik and A. Zeileis (2006). Unbiased recursive partitioning: A conditional inference framework. Journal of Computational 
#' and Graphical Statistics 15(3), 651-674. 
#' 
#' @seealso 
#' \code{\link[TSVC]{plot.TSVC}}, \code{\link[TSVC]{predict.TSVC}}, \code{\link[TSVC]{summary.TSVC}}
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
#' \dontshow{
#' fit0 <- TSVC(participation~foreign+oldkids, data=sl, family=binomial(link="logit"), 
#'              nperm=50, trace=TRUE)
#' print(fit0)
#' }
#' \donttest{
#' fit1 <- TSVC(participation~income+age, data=sl, family=binomial(link="logit"), 
#'              nperm=300, trace=TRUE)
#' print(fit1)
#' 
#' # In fit2, variable 'foreign' does not serve as effect modifier 
#' # and the effect of 'foreign' is not modified by the other variables.
#' # That means 'foreign' is assumed to only have simple linear effect on the response. 
#' fit2 <- TSVC(participation~income+age+foreign, data=sl, family=binomial(link="logit"), 
#'              nperm=300, trace=TRUE, effmod=c("income","age"), 
#'              notmod=list(c("foreign","income"),c("foreign","age")))
#' print(fit2)
#' }
#' 
#' @exportClass TSVC
#' @export
#' @importFrom stats coef deviance formula gaussian glm quantile 
#' @importFrom utils combn 
#' 

TSVC        <- function(formula, 
                        data,
                        family=gaussian,
                        alpha=0.05, 
                        nperm=1000, 
                        effmod=NULL,
                        notmod=NULL, 
                        trace=FALSE,
                        ...){
  UseMethod("TSVC")
}

#' @rdname TSVC
#' @method print TSVC 
#' @export

print.TSVC <-
  function(x, # object of class TSVC 
           ...){
    
    nobs   <- length(x$y)
    nvar   <- ncol(x$X)
    
    vars_varying  <- unique(x$splits[,"variable"])
    splits        <- x$splits[,c("variable","effect_modifier","threshold")]
    splits$variable <- names(x$X)[splits$variable]
   
    cat("\n")
    cat("Tree-structured varying coefficient model:\n")
    cat("\n")
    cat("Call:\n",paste(deparse(x$call)),"\n")
    cat("\n")
    cat("n=",nobs, "\n")
    cat("p=",nvar, "\n")
    if(!is.null(vars_varying)){
      cat("Variables with varying coefficients:", vars_varying, "\n")
    } else{
      cat("No varying coefficients.\n")
    }
    cat("\n")
    cat("Overview of executed splits:\n")
    if(length(vars_varying)>0){
      print(splits)
    } else{
      cat("no split performed")
    }
    invisible(x)
}  