info_summary <-
  function(object,
           variable){
    
    splits <- object$splits
    coefficients <- object$coefficients 
    vars_varying  <- unique(splits[,"variable"])
    
    type <- ifelse(variable %in% vars_varying,"varying","constant")
    if(type=="varying"){
        effmod <- paste(unique(splits[splits[,"variable"]==variable,"effect_modifier"]),collapse=",")
        nos    <- nrow(splits[splits[,"variable"]==variable,])
    } else{
        effmod <- nos <- "---"
        var_name <- names(object$X)[variable]
        if(!(var_name %in% names(coefficients$beta_constant))){
          type <- "---"
        }
        if(class(object$model)[1] == "gam" && var_name %in% sapply(object$model$smooth,"[[","vn")){
          type <- "smooth"
        }
    }
    
    output <- list("variable"=variable,
                   "type"=type,
                   "effmod"=effmod,
                   "nos"=nos)
    return(output)
  }

infos_summary <-
  function(object,
           variables){
    
    nvar   <- length(variables)
    output <- data.frame("variable"=numeric(nvar),
                         "type"=character(nvar),
                         "effect_modifier"=character(nvar),
                         "nosplits"=numeric(nvar),
                         stringsAsFactors=FALSE)
    
    for(i in 1:nvar){
      info <- info_summary(object,variables[i])
      output[i,"variable"]        <- info$variable
      output[i,"type"]            <- info$type
      output[i,"effect_modifier"] <- info$effmod
      output[i,"nosplits"]        <- info$nos
    }
    
    return(output)
  }