specification <- function(formula, data, effmod, notmod, only_effmod, smooth, split_intercept){
  
  formula <- paste(formula)
  
  ### y ###
  name_y  <- formula[2]
  y       <- data[, name_y]
  
  ### DM_kov ### 
  x      <- gsub(" ","",formula[3])
  name_x <- unlist(strsplit(x,"\\+"))
  DM_kov <- data[, name_x, drop=FALSE]
  if(split_intercept){
    DM_kov <- as.data.frame(cbind(DM_kov, "Intercept"=1))
  }
  
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
          preparation <- prepare_factor(i, effmod, notmod, name_x)
          effmod      <- preparation[[1]]
          notmod      <- preparation[[2]]
        }
      }
    }
  }
  
  if(!is.null(only_effmod)){
    for(i in only_effmod){
      var    <- which(name_x==i)
      notmod <- prepare_exclude(var, notmod, name_x)
    }
    only_effmod <- which(name_x%in%only_effmod)
  }
  
  if(!is.null(smooth)){
    for(i in smooth){
      var <- which(name_x==i)
      notmod <- prepare_exclude(var, notmod, name_x)
    }
    smooth <- which(name_x%in%smooth)
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
              "notmod"=notmod,
              "exclude"=c(only_effmod,smooth),
              "smooth"=smooth))
}