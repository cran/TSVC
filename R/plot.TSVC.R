#' Plotting of Varying Coefficient Trees 
#' 
#' @description 
#' Visualization of trees of effects of covariates that vary with the values of one or several effect modifiers. 
#' 
#' @param x a fitted object of class \code{\link[TSVC]{TSVC}}.
#' @param variable name of the variable, for which the tree shall be plotted.
#' @param ellipse_a controls width of ellipses containing coefficient estimates.
#' @param ellipse_b controls height of ellipses containing coefficient estimates.
#' @param ellipse_x controls location on x-axis of ellipses containing coefficient estimates.
#' @param ellipse_y controls location on y-axis of ellipses containing coefficient estimates.
#' @param branch_adj vertical adjustment of branch labels.
#' @param cex.lines width of branches of the tree. 
#' @param cex.branches size of the labels of the tree.
#' @param cex.coefs size of the coefficients in the terminal nodes of the tree.
#' @param cex.main size of the title of the tree. 
#' @param cex.numbers size of the internally used node number. 
#' @param draw_numbers if \code{true}, internally used node numbers are displayed. 
#' @param title  optional title, which is addded to the tree; if \code{title=NULL} the title is the name of the variable in the data. 
#' @param decimals number of decimals of coefficient estimates. Per default the coefficient estimates are displayed with three decimals. 
#' @param confint optional fitted object of class \code{\link[TSVC]{confint.TSVC}} with confidence intervals to be plotted 
#' in the terminal nodes of the tree; if \code{confint=NULL} (default) only the coefficient estimates will be plotted.  
#' @param ... further arguments passed to or from other methods. 
#' 
#' @author 
#' Moritz Berger <Moritz.Berger@imbie.uni-bonn.de> \cr \url{https://www.imbie.uni-bonn.de/people/dr-moritz-berger/}
#' 
#' @references 
#' Berger, M., G. Tutz and M. Schmid (2019). Tree-Structured Modelling of Varying Coefficients. Statistics and Computing 29, 217-229,
#' https://doi.org/10.1007/s11222-018-9804-8. 
#' 
#' @seealso 
#' \code{\link[TSVC]{TSVC}}, \code{\link[TSVC]{predict.TSVC}}, \code{\link[TSVC]{summary.TSVC}}
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
#'              nperm=1000, trace=TRUE)
#' plot(fit1, "income")
#' 
#' fit2 <- TSVC(participation~income+age, data=sl, family=binomial(link="logit"), 
#' perm_test=FALSE, test_linear=FALSE, splits_max=3)
#' set.seed(20012025)
#' ci2  <- confint(fit2, bootstrap_n=500, alpha=0.05, trace=TRUE)
#' plot(fit2, variable="income", confint=ci2, ellipse_y=0, draw_numbers=FALSE)
#' }
#' 
#' @method plot TSVC
#' @export
#' @importFrom plotrix draw.ellipse 
#' @importFrom grDevices grey
#' @importFrom graphics lines plot.new plot.window points rect text 
#' @importFrom methods is

plot.TSVC <- function(x, 
                      variable, 
                      ellipse_a=0.8, 
                      ellipse_b=0.2, 
                      ellipse_x=0,
                      ellipse_y=0,
                      branch_adj=0,
                      cex.lines=2,
                      cex.branches=1,
                      cex.coefs=1,
                      cex.main=1,
                      cex.numbers=1, 
                      draw_numbers=TRUE, 
                      title=NULL,
                      decimals=3, 
                      confint=NULL,
                      ...){
  
  if(is.null(x$splits)){
    cat("There is no plot available in the case of no varying coefficients.")
  } else{
    X        <- x$X 
    variable <- which(names(X)==variable)
    var_name <- colnames(X)[variable]
    info     <- x$splits[which(x$splits[,"variable"]==variable),]
    if(nrow(info)==0){
      beta_variable <- x$coefficients$beta_constant[var_name]
      cat("Variable", var_name, "does not have varying coefficients. There is no tree to plot.\n")
      cat("Estimated parameter:", beta_variable)
    } else{
      
      if(!is.null(confint) && !is(confint,"confint.TSVC")){
        stop("confint must be of class 'confint.TSVC'")
      } 

      params <- x$coefficients$beta_varying[[var_name]]
      
      if(!is.null(confint)){
        confint <- confint[confint$Variable%in%names(params),]
        confint <- confint[match(confint$Variable, names(params)),]
      }
      
      ptree(info, var_name, params, X, ellipse_a, ellipse_b, ellipse_x, ellipse_y, branch_adj, 
            cex.lines, cex.branches, cex.coefs, cex.main, cex.numbers, draw_numbers, title, decimals, confint)
    }
  }
  
}