prepare_exclude <- function(var, notmod, name_x){
  
  if(is.null(notmod)){
    notmod <- list() 
  }
  
  for(i in name_x[-var]){
    notmod[[length(notmod)+1]] <- c(name_x[var],i)
  }
 
  return(notmod) 
}