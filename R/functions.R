
# adjust node 
lu <- function(last=last, cd=cd, d=d, erg){
  
  s <- c("l","u")
  
  for(i in 1:length(s)){
    s_new <- paste0(last,s[i])
    if(cd==d){
      erg[length(erg)+1] <- s_new
    } else{
      erg <- lu(s_new,cd+1,d,erg)
    }
  }
  return(erg)
}


# compute ordered values 
ord_values <- function(x, n_quantile){
  if(is.factor(x)){
    ret <- unique(x)
  } else{
      if(!all((x - round(x)) == 0)){
        ret <- quantile(x,seq(0.05,1,length=n_quantile))
      } else{
        ret <- unique(sort(x))
      }
  }
  return(ret)
}

thresh <- function(ordered_values){
  ret <- ordered_values[-length(ordered_values)]
  return(ret)
}

# functions to build design 
design_one  <- function(x,threshold,upper){
  if(upper){
    ret <- ifelse(x > threshold,1,0)
  } else{
    ret <- ifelse(x > threshold,0,1)
  }
  return(ret)
}

design <- function(x,thresholds,upper){
  ret <- sapply(thresholds, function(j) design_one(x,j,upper))
  return(ret)
}

designlist <- function(X,vars,label,thresholds,var_names,upper=TRUE){
  ret <- lapply(vars, function(j) {
                  ret <- design(X[,j],thresholds[[j]],upper)
                  if(!is.matrix(ret)){
                    ret <- t(as.matrix(ret))
                  }
                  colnames(ret) <- paste0("s",which(var_names==j),1:length(thresholds[[j]]),ifelse(upper,"_u","_l"),label)
                  return(ret)})
  names(ret) <- vars 
  return(ret)
}

# function to specify variables
vars <- function(X){
  ret <- rev(combn(colnames(X),ncol(X)-1,simplify=FALSE))
  return(ret)
}

# function to compute all models 
one_model <- function(var,split_var,kn,count,j,design_lower,design_upper,params,dat0,var_names,smooth,split_intercept,family,start,mcs){
  
  # var: integer
  # split_var: name as character 
  dat   <- data.frame(dat0,design_lower[[var]][[split_var]][,j,drop=FALSE],design_upper[[var]][[split_var]][,j,drop=FALSE])
  
  help1 <- params[[count]][[var]][kn]
  help2 <- paste0(unlist(params[[count]])[-which(unlist(params[[count]])==help1)],collapse="+")
  help3 <- paste(help1,c(colnames(design_lower[[var]][[split_var]])[j],colnames(design_upper[[var]][[split_var]])[j]),sep=":")
  help4 <- paste(c(help2,help3),collapse="+")
  help5 <- formula(paste("y~",help4, ifelse(split_intercept,"-1","")))
  
  if(is.null(smooth)){
    mod   <- glm(help5,family=family,data=dat,etastart=start)
  } else{
    help00 <- paste0("s(", var_names[smooth], ", bs='ps'", mcs, ")", collapse="+")
    help5  <- formula(paste("y~", help4, "+", help00, ifelse(split_intercept,"-1","")))
    mod    <- gam(help5,family=family,data=dat,etastart=start)
  }
  attributes(mod)$formula_orig <- help5
  return(mod)
  
}


allmodels <- function(var,split_var,kn,count,design_lower,design_upper,splits_evtl,params,dat0,mod0,n_s,var_names,smooth,split_intercept,family,start,mcs){
  
  # var: integer
  # split_var: name as character 
  
  deviances      <- rep(0,n_s[split_var])
  splits_aktuell <- splits_evtl[[count]][[var]][[split_var]][kn,]
  splits_aktuell <- splits_aktuell[!is.na(splits_aktuell)]
  
  if(length(splits_aktuell)>0){
    for(j in splits_aktuell){
      suppressWarnings(
        mod <- one_model(var,split_var,kn,count,j,design_lower,design_upper,params,dat0,var_names,smooth,split_intercept,family,start,mcs)
      )
      deviances[j] <-  deviance(mod0)-deviance(mod)
    }
  }
  return(deviances)
}

# function to compute permutation test 
one_permutation <- function(var,split_var,kn,count,design_lower,design_upper,
                            DM_kov,thresholds,which_obs,splits_evtl,params,dat0,mod0,n_s,var_names,smooth,split_intercept,family,start,mcs){
  
  dv_perm <- rep(0,n_s[var])
  obs_aktuell <- which_obs[[count]][[var]][kn,]
  obs_aktuell <- obs_aktuell[!is.na(obs_aktuell)]
  DM_kov_perm <- DM_kov
  DM_kov_perm[obs_aktuell,split_var] <- sample(DM_kov_perm[obs_aktuell,split_var],length(obs_aktuell))
  
  design_upper_perm <- design_upper
  design_upper_perm[[var]][[split_var]] <- design(DM_kov_perm[,split_var],thresholds[[split_var]], upper=TRUE)
  colnames(design_upper_perm[[var]][[split_var]]) <- colnames(design_upper[[var]][[split_var]])
  
  design_lower_perm <- design_lower
  design_lower_perm[[var]][[split_var]] <- design(DM_kov_perm[,split_var],thresholds[[split_var]], upper=FALSE)
  colnames(design_lower_perm[[var]][[split_var]]) <- colnames(design_lower[[var]][[split_var]])
  
  dv_perm <- allmodels(var,split_var,kn,count,design_lower_perm,design_upper_perm,splits_evtl,params,dat0,mod0,n_s,var_names,smooth,split_intercept,family,start,mcs)
  
  return(max(dv_perm))
  
}

# function to compute final tests on linear terms 
final_test <- function(var,var_seq,params_opt,dat,var_names,smooth,split_intercept,family,mcs){
  
  help01 <- var_seq[!(var==var_seq)]
  help02 <- unlist(params_opt[help01])
  help03 <- paste(help02,collapse="+")
  help04 <- formula(paste("y~",help03,"-1"))
  
  help12 <- unlist(params_opt[var_seq])
  help13 <- paste(help12,collapse="+")
  help14 <- formula(paste("y~",help13,"-1"))
  
  if(is.null(smooth)){
    mod0   <- glm(help04,family=family,data=dat)
  } else{
    help00 <- paste0("s(", var_names[smooth], ", bs='ps'", mcs, ")", collapse="+")
    help04  <- formula(paste("y~", help03, "+", help00, ifelse(split_intercept,"-1","")))
    mod0    <- gam(help04,family=family,data=dat)
  }
  if(is.null(smooth)){
    mod1   <- glm(help14,family=family,data=dat)
  } else{
    help00 <- paste0("s(", var_names[smooth], ", bs='ps'", mcs, ")", collapse="+")
    help14  <- formula(paste("y~", help13, "+", help00, ifelse(split_intercept,"-1","")))
    mod1    <- gam(help14,family=family,data=dat)
  }

  dev <-  deviance(mod0)-deviance(mod1)
  return(dev)
}

# function to compute one step of final backward elimination 
one_step <- function(vars,var_seq,var_names,params_opt,dat,smooth,split_intercept,family,nperm,mcs){
  
  pvalues <- numeric(length(vars))
  for(v in seq_along(vars)){
    T1 <- final_test(vars[v],var_seq,params_opt,dat,var_names,smooth,split_intercept,family=family,mcs)
    T0 <- numeric(nperm)
    for(perm in 1:nperm){
      var <- var_names[vars[v]]
      dat_a <- dat
      dat_a[,var] <- sample(dat_a[,var],nrow(dat))
      T0[perm] <- final_test(vars[v],var_seq,params_opt,dat_a,var_names,smooth,split_intercept,family=family,mcs)
    }
    pvalues[v] <- sum(T0>=T1)/nperm
  }
  
  return(pvalues)
}

# functions to fix names 
check_names <- function(model, params){
  
  params_c  <- unlist(params)
  names_mod <- names(coef(model))
  params_split <- strsplit(params_c,":")
  g1_c         <- which(sapply(params_split,length)>1)
  params_split <- params_split[g1_c]
  params_c     <- params_c[g1_c]
  names_split  <- strsplit(names_mod,":")
  g1_mod       <- which(sapply(names_split,length)>1)
  names_split  <- names_split[g1_mod]
  names_mod    <- names_mod[g1_mod]
  if(length(params_c)>0){
    if(!all(params_c %in% names_mod)){
      for(k in 1:length(params_split)){
        found <- FALSE
        kk   <- 0
        while(!found){
          kk <- kk+1
          found <- all(params_split[[k]] %in% names_split[[kk]])
        } 
        names(model$coefficients)[g1_mod[kk]] <- params_c[k] 
      }
    }
  }
  return(model)
}

check_names_list <- function(model_list, params){
  for(j in 1:length(model_list)){
    model_list[[j]] <- check_names(model_list[[j]],params[[j]])
  }
  return(model_list)
}

# functions for confidence intervals 
optimize_nSplits <- function(object, data, method = "BIC", splits_max = 10){
  n <- dim(data)[1]
  call <- object$call
  call[["data"]] <- data
  call[["splits_max"]] <- splits_max
  model_TSVC <- eval(call)
  if(method == "AIC"){
    AIC <- sapply(model_TSVC$all_models, AIC)
    splits <- which(AIC == min(na.omit(AIC))) - 1
    return(list(splits = splits, AIC = AIC, model_full = model_TSVC))
  }
  if(method == "BIC"){
    BIC <- sapply(model_TSVC$all_models, BIC)
    splits <- which(BIC == min(na.omit(BIC))) - 1
    return(list(splits = splits, BIC = BIC, model_full = model_TSVC))
  }
}

constructTSVCdata <- function(object, data, out){
  data_h <- data
  model_b <- object
  n_varyingCoeff <- length(model_b$coefficients$beta_varying)
  varyingVars <- unique(model_b$splits$variable)
  if(length(varyingVars) == 0){
    return(list(data_h = data, data_i = data))
  }
  vars <- names(model_b$X)[-varyingVars]
  allVars <- names(model_b$X)
  data_i <- data_h
  for(j in 1:n_varyingCoeff){
    n_nodes <- length(unlist(model_b$coefficients$beta_varying[[j]]))
    splits <- model_b$splits[model_b$splits$variable == varyingVars[j],]
    leaf_nodes <- c()
    for(k in 1:(n_nodes-1)){
      H <- 1*I(data_h[, splits$effect_modifier[k]] <= splits$threshold[k])
      H2 <- 1*I(data_h[, splits$effect_modifier[k]] > splits$threshold[k])
      data_h <- tibble::add_column(data_h, H)
      data_h <- tibble::add_column(data_h, H2)
      data_i <- tibble::add_column(data_i, H)
      data_i <- tibble::add_column(data_i, H2)
      names(data_h)[(length(names(data_h))-1):length(names(data_h))] <- c(paste0("L", as.character(splits$number[k])), paste0("R", as.character(splits$number[k])))
      names(data_i)[(length(names(data_i))-1):length(names(data_i))] <- c(paste0("L", as.character(splits$number[k])), paste0("R", as.character(splits$number[k])))
      if(k > 1){
        node <- splits[k, "number"]
        parent_node <- splits[splits$left == node | splits$right == node, "number"]
        parent_part <- ifelse(splits[splits$number == parent_node, "left"] == node, paste0("L", as.character(parent_node)),
                              paste0("R", as.character(parent_node)))
        data_h[, paste0("L", as.character(splits$number[k]))] <- data_h[, paste0("L", as.character(splits$number[k]))]*
          data_h[, parent_part]
        data_h[, paste0("R", as.character(splits$number[k]))] <- data_h[, paste0("R", as.character(splits$number[k]))]*
          data_h[, parent_part]
        data_i[, paste0("L", as.character(splits$number[k]))] <- data_i[, paste0("L", as.character(splits$number[k]))]*
          data_i[, parent_part]
        data_i[, paste0("R", as.character(splits$number[k]))] <- data_i[, paste0("R", as.character(splits$number[k]))]*
          data_i[, parent_part]
        leaf_nodes <- leaf_nodes[leaf_nodes != parent_part]
        leaf_nodes <- c(leaf_nodes, paste0("L", as.character(splits$number[k])), 
                        paste0("R", as.character(splits$number[k])))
      }else{
        leaf_nodes <- c(leaf_nodes, paste0("L", as.character(splits$number[k])), paste0("R", as.character(splits$number[k])))
      }
    }
    for(ln in leaf_nodes){
      data_h[, ln] <- data_h[, ln] * data_h[, names(model_b$coefficients$beta_varying)[j]]
      names(data_h)[names(data_h) == ln] <- paste0(ln, names(model_b$coefficients$beta_varying)[j])
      names(data_i)[names(data_i) == ln] <- paste0(ln, names(model_b$coefficients$beta_varying)[j])
      vars <- c(vars, paste0(ln, names(model_b$coefficients$beta_varying)[j]))
    }
    data_h <- data_h[, unique(c(allVars, vars))]
    data_i <- data_i[, unique(c(allVars, vars))]
  }
  data_h <- data_h[, vars]
  data_i <- data_i[, vars]
  data_h[, out] <- data[, out]
  return(list(data_h = data_h, data_i = data_i))
}

constructTSVCdata_s <- function(object, data, splits, out){
  data_h <- data
  model_b <- object$all_models[splits + 1][[1]]
  allVars <- names(object$X)
  varyingVars <- unique(object$splits[1:splits,]$variable)
  varyingVars_names <- allVars[varyingVars]
  n_varyingCoeff <- length(varyingVars)
  if(length(varyingVars) == 0){
    return(list(data_h = data, data_i = data))
  }
  vars <- names(object$X)[-varyingVars]
  model_b$splits <- object$splits[1:splits,]
  data_i <- data_h
  for(j in 1:n_varyingCoeff){
    n_nodes <- sum(1*grepl(varyingVars_names[j], names(model_b$coefficients)))
    splits <- model_b$splits[model_b$splits$variable == varyingVars[j],]
    leaf_nodes <- c()
    for(k in 1:(n_nodes-1)){
      H <- 1*I(data_h[, splits$effect_modifier[k]] <= splits$threshold[k])
      H2 <- 1*I(data_h[, splits$effect_modifier[k]] > splits$threshold[k])
      data_h <- tibble::add_column(data_h, H)
      data_h <- tibble::add_column(data_h, H2)
      data_i <- tibble::add_column(data_i, H)
      data_i <- tibble::add_column(data_i, H2)
      names(data_h)[(length(names(data_h))-1):length(names(data_h))] <- c(paste0("L", as.character(splits$number[k])), paste0("R", as.character(splits$number[k])))
      names(data_i)[(length(names(data_i))-1):length(names(data_i))] <- c(paste0("L", as.character(splits$number[k])), paste0("R", as.character(splits$number[k])))
      if(k > 1){
        node <- splits[k, "number"]
        parent_node <- splits[splits$left == node | splits$right == node, "number"]
        parent_part <- ifelse(splits[splits$number == parent_node, "left"] == node, paste0("L", as.character(parent_node)),
                              paste0("R", as.character(parent_node)))
        data_h[, paste0("L", as.character(splits$number[k]))] <- data_h[, paste0("L", as.character(splits$number[k]))]*
          data_h[, parent_part]
        data_h[, paste0("R", as.character(splits$number[k]))] <- data_h[, paste0("R", as.character(splits$number[k]))]*
          data_h[, parent_part]
        data_i[, paste0("L", as.character(splits$number[k]))] <- data_i[, paste0("L", as.character(splits$number[k]))]*
          data_i[, parent_part]
        data_i[, paste0("R", as.character(splits$number[k]))] <- data_i[, paste0("R", as.character(splits$number[k]))]*
          data_i[, parent_part]
        leaf_nodes <- leaf_nodes[leaf_nodes != parent_part]
        leaf_nodes <- c(leaf_nodes, paste0("L", as.character(splits$number[k])), 
                        paste0("R", as.character(splits$number[k])))
      }else{
        leaf_nodes <- c(leaf_nodes, paste0("L", as.character(splits$number[k])), paste0("R", as.character(splits$number[k])))
      }
    }
    for(ln in leaf_nodes){
      data_h[, ln] <- data_h[, ln] * data_h[, varyingVars_names[j]]
      names(data_h)[names(data_h) == ln] <- paste0(ln, varyingVars_names[j])
      names(data_i)[names(data_i) == ln] <- paste0(ln, varyingVars_names[j])
      vars <- c(vars, paste0(ln, varyingVars_names[j]))
    }
    data_h <- data_h[, unique(c(allVars, vars))]
    data_i <- data_i[, unique(c(allVars, vars))]
  }
  data_h <- data_h[, vars]
  data_i <- data_i[, vars]
  data_h[, out] <- data[, out]
  return(list(data_h = data_h, data_i = data_i))
}

