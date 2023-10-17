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
#' @param nodesize_min minimum number of observations that must exist in a node in order for a split to be attempted. 
#' @param bucket_min the minimum number of observations in any terminal node. 
#' @param depth_max maximum depth of any node in each tree, with the root node counted as depth 0. If \code{NULL} (default), the size of the trees is not restricted.
#' @param splits_max maximum number of splits performed. If \code{NULL} (default), the number of splits is not restricted.
#' @param perm_test if \code{FALSE}, no permutation tests are performed, but each tree is grown until the minimum node size constraint is reached.
#' @param test_linear should linear effects that were not modified during iteration tested for significance? 
#' @param effmod optional vector of covariates that serve as effect modifier. If \code{NULL} (default), all covariates are considered as potential effect modifiers. 
#' @param notmod optional list of class \code{\link{list}} containing pairs of covariate/effect modifier that are not considered as candidates for splitting during iteration.
#' If \code{NULL} (default), all combinations of covariates and potential effect modifiers are considered for splitting. 
#' @param only_effmod optional vector of covariates that serve as effect modifier, only. If \code{NULL} (default), all effect modifiers are included in the predictor of the model and 
#' are allowed to be modified. 
#' @param smooth optional vector of covariates with a smooth effect on the response. The (smooth) effects fo these variables are not allowed to be modified.
#' @param split_intercept if \code{TRUE}, the intercept is allowed to be modified by the covariates. If \code{FALSE} (default), the intercept is set constant.  
#' @param sb_slope optional vector of covariates that are allowed to be modified by itself. Such an effect corresponds to a structural break in the slope. 
#' @param trace if \code{TRUE}, information about the estimation progress is printed. 
#' @param x object of class \code{TSVC}.
#' @param ... further arguments passed to or from other methods. 
#' 
#' @details 
#' A typical \link{formula} has the form \code{response ~ covariates}, where \code{response} is the name of the response variable and \code{covariates} is a
#' series of variables that are incorporated in the model. 
#' 
#' With p covariates, \code{TSVC} expects a formula of the form \eqn{y ~ x_1+...+x_p}. If no further specifications are made (\code{effmod=NULL}, \code{notmod=NULL}, \code{only_effmod=NULL}) it is assumed that each covariate \eqn{x_j, j = {1,...,p}} 
#' can be modified by all the other variables \eqn{x_m, m = {1,...,p} \ j}. 
#' 
#' Remark: Significance of each split is verified by permutation tests. The result of the permutation tests 
#' can strongly depend on the number of permutations \code{nperm}.
#' 
#' Note: The algorithm currently does not support splitting of/by factor variables. If a factor variable is included in the \link{formula} of the model, the variable will not serve as 
#' effect modifier and its effect will not be modified. 
#' 
#' @return 
#' Object of class \code{"TSVC"}. An object of class \code{"TSVC"} is a list containing the following components: 
#' 
#' \item{splits}{matrix with detailed information about all executed splits during the fitting process.}
#' \item{coefficients}{list of estimated coefficients for covariates with and without varying coefficients (including a non-varying intercept).}
#' \item{pvalues}{p-values of each permuation test during the fitting process.}
#' \item{pvalues_linear}{p-values of the permutation tests on the linear effects in the last step of the algorithm.}
#' \item{devs}{maximal value statistics \eqn{T_m} of the selected effect modifier in each iteration during the fitting process.}
#' \item{crit}{critical values of each permutation test during the fitting process.}
#' \item{y}{response vector.}
#' \item{X}{matrix of all the variables (covariates and effect modifiers) for model fitting.}
#' \item{sb}{variables for which a structural break in the slope was allowed.}
#' \item{model}{internally fitted model in the last iteration of class \code{\link{glm}} or \code{\link[mgcv]{gam}}.}
#' \item{all_models}{list of internally fitted models of class \code{\link{glm}} or \code{\link[mgcv]{gam}}.}
#' 
#' @author 
#' Moritz Berger <Moritz.Berger@imbie.uni-bonn.de> \cr \url{https://www.imbie.uni-bonn.de/people/dr-moritz-berger/}
#' 
#' @references 
#' Berger, M., G. Tutz and M. Schmid (2019). Tree-Structured Modelling of Varying Coefficients. Statistics and Computing 29, 217-229,
#' https://doi.org/10.1007/s11222-018-9804-8.
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
#' \dontrun{
#' fit1 <- TSVC(participation~income+age, data=sl, family=binomial(link="logit"), 
#'              nperm=300, trace=TRUE)
#' print(fit1)
#' class(fit$model) # glm 
#' 
#' # In fit2, variable 'foreign' does not serve as effect modifier 
#' # and the effect of 'foreign' is not modified by the other variables.
#' # That means 'foreign' is assumed to only have simple linear effect on the response. 
#' fit2 <- TSVC(participation~income+age+foreign, data=sl, family=binomial(link="logit"), 
#'              nperm=300, trace=TRUE, effmod=c("income","age"), 
#'              notmod=list(c("foreign","income"),c("foreign","age")))
#' print(fit2)
#'
#' # In fit3, variable 'age' does only serve as effect modifier. That means the effect of 'age' 
#' # is not included in the predictor of the model.    
#' fit3 <- TSVC(participation~income+age+foreign, data=sl, family=binomial(link="logit"),
#'              nperm=300, trace=TRUE, only_effmod="age")
#' print(fit3)  
#' 
#' # In fit4, the intercept is allowed to be modified by 'age' and 'income'. 
#' # The two covariates, however, are not allowed to modify each other. 
#' fit4 <- TSVC(participation~income+age, data=sl, family=binomial(link="logit"),
#'              nperm=300, trace=TRUE, split_intercept=TRUE,
#'              notmod=list(c("income","age"), c("age", "income")))
#' print(fit4)    
#' 
#' 
#' # In fit5, variable 'age' has a smooth effect on the response. 
#' # Hence, the (smooth) effect of 'age' will not be modified by the other variables.  
#' fit5 <- TSVC(participation~income+age+foreign, data=sl, family=binomial(link="logit"),
#'              nperm=300, trace=TRUE, smooth="age")                     
#' print(fit5)
#' class(fit5$model) # gam      
#' 
#' # In fit6, the intercept is allowed to be modified by 'age' and 'income', but the two variables are 
#' # not included in the predictor of the model. Here, no permutation tests are performed, but the 
#' # tree is pruned by a minimum node size constraint. 
#' fit6 <- TSVC(participation~income+age, data=sl, family=binomial(link="logit"),
#'              perm_test=FALSE, nodesize_min=100, bucket_min=100, trace=TRUE, split_intercept=TRUE,
#'              effmod=c("income","age"), only_effmod = c("income", "age"))   
#' print(fit6)                          
#' 
#' }
#' 
#' @export
#' @importFrom stats coef deviance formula gaussian glm quantile 
#' @importFrom utils combn 
#' @importFrom mgcv gam 
#' 

TSVC        <- function(formula, 
                        data,
                        family=gaussian,
                        alpha=0.05, 
                        nperm=1000, 
                        nodesize_min=5,
                        bucket_min=1,
                        depth_max=NULL,
                        splits_max=NULL,
                        perm_test=TRUE,
                        test_linear=TRUE,
                        effmod=NULL,
                        notmod=NULL,
                        only_effmod=NULL,
                        smooth=NULL,
                        split_intercept=FALSE,
                        sb_slope=NULL,
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
    cat("n =",nobs, "\n")
    cat("p =",nvar, "\n")
    if(!is.null(vars_varying)){
      cat("Variables with varying coefficients:", vars_varying, "\n")
    } else{
      cat("No varying coefficients.\n")
    }
    cat("\n")
    cat("Overview of executed splits:\n")
    if(!is.null(vars_varying)){
      print(splits)
    } else{
      cat("no split performed")
    }
    invisible(x)
}  
