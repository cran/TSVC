#' Plotting of Varying Coefficient Trees 
#' 
#' @description 
#' Visualization of trees of effects of covariates that vary with the values of one or several effect modifiers. 
#' 
#' @param x object of class \code{\link[TSVC]{TSVC}}.
#' @param variable name of the variable, for which the tree shall be plotted.
#' @param cex.lines width of branches of the tree. 
#' @param cex.branches size of the labels of the tree.
#' @param cex.coefs size of the coefficients in the terminal nodes of the tree.
#' @param cex.main size of the title of the tree. 
#' @param title  optional title, which is addded to the tree; if \code{title=NULL} the title is the name of the variable in the data. 
#' @param ... further arguments passed to or from other methods. 
#' 
#' @author 
#' Moritz Berger <Moritz.Berger@imbie.uni-bonn.de> \cr \url{http://www.imbie.uni-bonn.de/personen/dr-moritz-berger/}
#' 
#' @references 
#' Berger, M., G. Tutz and M. Schmid (2018). Tree-Structured Modelling of Varying Coefficients. Statistics and Computing, published online,
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
#' }
#' 
#' @method plot TSVC
#' @export
#' @importFrom plotrix draw.ellipse 
#' @importFrom grDevices grey
#' @importFrom graphics lines plot.new plot.window points rect text 

plot.TSVC <- function(x, 
                            variable, 
                            cex.lines=2,
                            cex.branches=1,
                            cex.coefs=1,
                            cex.main=1,
                            title=NULL, 
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
      params <- x$coefficients$beta_varying[[var_name]]
      ptree(info, var_name, params, X, cex.lines, cex.branches, cex.coefs, cex.main, title)
    }
  }
  
}