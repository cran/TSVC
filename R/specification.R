specification <- function(formula, data, effmod, notmod){
  
  formula <- paste(formula)
  
  ### y ###
  name_y  <- formula[2]
  y       <- data[, name_y]
  
  ### DM_kov ### 
  x      <- gsub(" ","",formula[3])
  name_x <- unlist(strsplit(x,"\\+"))
  DM_kov <- data[, name_x, drop=FALSE]
  
  if(any(sapply(1:ncol(DM_kov),function(j) class(DM_kov[,j]))=="character")){
    stop("variable of class 'character' is not suitable")
  }
  if(any(sapply(1:ncol(DM_kov),function(j) class(DM_kov[,j]))=="logical")){
    stop("variable of class 'logical' is not suitable")
  }
  
  for(i in 1:ncol(DM_kov)){
    if(is.factor(DM_kov[,i])){
      if(is.ordered(DM_kov[,i])){
        DM_kov[,i] <- as.numeric(DM_kov[,i])
      } else{
        if(nlevels(DM_kov[,i])==2){
          DM_kov[,i] <- as.numeric(DM_kov[,i])-1
        } else{
          stop("the use of nominal factors is currently not implemented")
        }
      }
    }
  }
  
  if(!is.null(effmod)){
    effmod <- which(name_x%in%effmod)
  }
  if(!is.null(notmod)){
    notmod <- lapply(1:length(notmod), function(j) c(which(name_x==notmod[[j]][1]), which(name_x==notmod[[j]][2])))
  }
  
  return(list("y"=y,
              "DM_kov"=DM_kov,
              "effmod"=effmod,
              "notmod"=notmod))
}