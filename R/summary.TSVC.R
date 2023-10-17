#' Summary of Tree-Structured Varying Coefficient Models 
#' 
#' @description 
#' Summary for an object of class \code{TSVC}, with an overview of all executed splits during the fitting process.
#' 
#' @param object object of class \code{\link[TSVC]{TSVC}}.
#' @param x object of class \code{\link[TSVC]{summary.TSVC}}.
#' @param ... further arguments passed to or from other methods.
#' 
#' @return object of class \code{"summary.TSVC"}. An object of class \code{"summary.TSVC"} is a list containing the 
#' following components: 
#' 
#' \item{stats}{overview of detected varying coefficients, responsible effect modifiers and executed splits.}
#' \item{nosplits}{total number of executed splits during the fitting process.} 
#' 
#' @author 
#' Moritz Berger <Moritz.Berger@imbie.uni-bonn.de> \cr \url{https://www.imbie.uni-bonn.de/people/dr-moritz-berger/}
#' 
#' 
#' @references 
#' Berger, M., G. Tutz and M. Schmid (2019). Tree-Structured Modelling of Varying Coefficients. Statistics and Computing 29, 217-229,
#' https://doi.org/10.1007/s11222-018-9804-8. 
#' 
#' @seealso 
#' \code{\link[TSVC]{TSVC}}, \code{\link[TSVC]{plot.TSVC}}, \code{\link[TSVC]{predict.TSVC}}
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
#' summary(fit1)
#' }
#' 
#' 
#' @method summary TSVC
#' @export

summary.TSVC <-
  function(object, # object of class TSVC
           ...){
    
    to_return <- list(call=object$call)
    
    nvar     <- ncol(object$X)
    var_names<- colnames(object$X)
    overview <- infos_summary(object,1:nvar)
    overview$variable <- var_names 
    nos      <- nrow(object$splits)
    
    to_return$stats    <- overview 
    to_return$nosplits <- nos  
    
    class(to_return) <- "summary.TSVC"
    return(to_return)
    
  }

#' @rdname summary.TSVC
#' @method print summary.TSVC
#' @export

print.summary.TSVC <-
  function(x, # object of class summary.TSVC 
           ...){
    
    cat("\n")
    cat("Tree-structured varying coefficient model:\n")
    cat("\n")
    cat("Call:\n",paste(deparse(x$call)),"\n")
    cat("\n")
    cat("----------\n")
    cat("\n")
    cat("Overview:\n")
    cat("\n")
    print(x$stats)
    cat("\n")
    cat("Total number of Splits:", x$nosplits)
  }
