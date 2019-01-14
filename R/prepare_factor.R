prepare_factor <- function(fac, effmod, notmod, name_x){
  
  if(is.null(effmod)){
    effmod <- name_x[-fac]
  } else{
    if(name_x[fac]%in%effmod){
      effmod <- effmod[!which(effmod==name_x[fac])]
    }
  }
  
  if(is.null(notmod)){
    notmod <- list() 
  }
  
  for(i in name_x[-fac]){
    notmod[[length(notmod)+1]] <- c(name_x[fac],i)
  }
 
  return(list(effmod, notmod)) 
}