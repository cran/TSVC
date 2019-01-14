
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
ord_values <- function(x){
  if(is.factor(x)){
    ret <- unique(x)
  } else{
      if(!all((x - round(x)) == 0)){
        ret <- quantile(x,seq(0.05,1,by=0.05))
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

