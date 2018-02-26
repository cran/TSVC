info_summary <-
  function(splits,
           variable){
    
    vars_varying  <- unique(splits[,"variable"])
    
    type <- ifelse(variable %in% vars_varying,"varying","constant")
    if(type=="varying"){
        effmod <- paste(unique(splits[splits[,"variable"]==variable,"effect_modifier"]),collapse=",")
        nos    <- nrow(splits[splits[,"variable"]==variable,])
    } else{
        effmod <- nos <- "---"
    }
    
    
    output <- list("variable"=variable,
                   "type"=type,
                   "effmod"=effmod,
                   "nos"=nos)
    return(output)
  }

infos_summary <-
  function(splits,
           variables){
    
    nvar   <- length(variables)
    output <- data.frame("variable"=numeric(nvar),
                         "type"=character(nvar),
                         "effect_modifier"=character(nvar),
                         "nosplits"=numeric(nvar),
                         stringsAsFactors=FALSE)
    
    for(i in 1:nvar){
      info <- info_summary(splits,variables[i])
      output[i,"variable"]        <- info$variable
      output[i,"type"]            <- info$type
      output[i,"effect_modifier"] <- info$effmod
      output[i,"nosplits"]        <- info$nos
    }
    
    return(output)
  }