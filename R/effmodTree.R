effmodTree  <- function(y, 
                        DM_kov,
                        family,
                        alpha, 
                        nperm, 
                        nodesize_min,
                        bucket_min,
                        depth_max,
                        perm_test,
                        effmod,
                        notmod, 
                        exclude,
                        smooth,
                        split_intercept,
                        trace, 
                        ...){
  
  n    <- length(y)
  nvar <- ncol(DM_kov)
  if(length(nodesize_min)==1){
    nodesize_min <- rep(nodesize_min,nvar)
  }
  if(!is.null(depth_max) & length(depth_max)==1){
    depth_max <- rep(depth_max, nvar)
  }
  
  mc   <- as.list(match.call(expand.dots=TRUE))
  mcs  <- paste0(names(mc),"=",mc)[-c(1:16)]
  mcs  <- ifelse(length(mcs)==0, "", paste0(",", paste0(mcs, collapse=",")))
    
  if(!is.null(names(DM_kov))){
    var_names <- names(DM_kov)
  } else{
    var_names <- paste0("x",1:nvar)
  }
  
  ordered_values <- lapply(DM_kov,ord_values)
  n_levels       <- sapply(ordered_values,length)
  thresholds     <- lapply(ordered_values,thresh)
  n_s            <- n_levels-1

  mod_potential <- list()
  if(perm_test){
    devs          <- c()
    crits         <- c()
    pvalues       <- c()
  } else{
    devs          <- NULL
    crits         <- NULL
    pvalues       <- NULL
  }
  splits        <- data.frame("variable"=numeric(),
                              "effect_modifier"=character(),
                              "split"=numeric(),
                              "level"=numeric(),
                              "node"=numeric(),
                              "threshold"=numeric(),
                              "number"=numeric(),
                              "left"=numeric(),
                              "right"=numeric(),stringsAsFactors=FALSE)
  params        <- list()
  vars_evtl     <- list()
  splits_evtl   <- list()
  which_obs     <- list()
  numbers       <- list()
  count         <- 1
  var_list      <- vars(DM_kov)
  
  # delete variables that do not modify 
  if(!is.null(effmod)){
    for(i in 1:nvar){
      var_list[[i]] <- var_list[[i]][var_list[[i]]%in%var_names[effmod]]
    }
  }
  # delete single pairs that can not vary/modify (variable/effect modifier)
  if(!is.null(notmod)){
    for(i in 1:length(notmod)){
      var_list[[notmod[[i]][1]]] <- var_list[[notmod[[i]][1]]][!(var_list[[notmod[[i]][1]]]==var_names[notmod[[i]][2]])]
    }
  }
  # delete intercept
  if(split_intercept){
    if(is.null(effmod)){
      for(i in 1:(nvar-1)){
        var_list[[i]] <- var_list[[i]][-length(var_list[[i]])]
      }
    }
  }
  
  params[[1]]      <- as.list(var_names)
  if(!is.null(exclude)){
    for(e in exclude){
      params[[1]][[e]] <- params[[1]][[e]][-1]
    }
  }
  
  which_obs[[1]]   <- replicate(nvar,matrix(1:n,nrow=1),simplify=FALSE)
  vars_evtl[[1]]   <- lapply(1:nvar,function(j) matrix(var_list[[j]], nrow=1))
  splits_evtl[[1]] <- lapply(1:nvar,function(j) {
                        ret <- lapply(seq_along(var_list[[j]]), function(v) matrix(1:n_s[var_list[[j]][v]],nrow=1))
                        names(ret) <- var_list[[j]]
                        return(ret)})
  numbers[[1]]     <- replicate(nvar,1,simplify=FALSE)
  
  help0  <- paste0(if(is.null(exclude)){var_names}else{var_names[-exclude]},collapse="+")
  form0  <- formula(paste("y~",help0,ifelse(split_intercept,"-1","")))
  dat0   <- as.data.frame(cbind(y,DM_kov))
  if(is.null(smooth)){
    mod0   <- mod_potential[[1]] <- glm(form0,family=family,data=dat0)
  } else{
    help00 <- paste0("s(", var_names[smooth], ", bs='ps'", mcs, ")", collapse="+")
    form0  <- formula(paste("y~", help0, "+", help00,ifelse(split_intercept,"-1","")))
    mod0   <- mod_potential[[1]] <- gam(form0,family=family,data=dat0)
  }
  start <- predict(mod0)
  
  design_upper <- lapply(1:nvar,function(j) designlist(DM_kov,var_list[[j]],j,thresholds,var_names))
  design_lower <- lapply(1:nvar,function(j) designlist(DM_kov,var_list[[j]],j,thresholds,var_names, upper=FALSE))
  
  # delete binary vars from intercept (until first split in intercept)
  int_del <- c() 
  if(split_intercept){
    for(i in 1:(nvar-1)){
      if(n_levels[i]==2 && !(i%in%exclude)){
        int_del <- c(int_del,var_names[i])
      }
    }
  }
  split_in_int  <- FALSE
  
  sig      <- TRUE
  anysplit <- TRUE
  
  # minimal node size constraint
  for(v in 1:nvar){
    if(length(var_list[[v]])>0){
      if(length(which(!is.na(which_obs[[1]][[v]][1,])))<nodesize_min[v]){
        for(var in var_list[[v]]){splits_evtl[[1]][[v]][[var]][1,] <- rep(NA, n_s[var])}
      }
    }
  }
  # depth_max==0? 
  if(!is.null(depth_max)){
    for(v in 1:nvar){
      if(length(var_list[[v]])>0){
        if(depth_max[v]==0){
          for(var in var_list[[v]]){splits_evtl[[1]][[v]][[var]][1,] <- rep(NA, n_s[var])}
        }
      }
    }
  }
  anysplit <- !all(is.na(unlist(splits_evtl[[1]])))
  
  
  while(sig & anysplit){
    
    # minbucket
    for(v in 1:nvar){
      if(length(var_list[[v]])>0){
        for(i in var_list[[v]]){
          n_knots <- length(params[[count]][[v]])
          for(kn in 1:n_knots){
            splits_aktuell <- splits_evtl[[count]][[v]][[i]][kn,]
            if(length(splits_aktuell)>0){
              for(j in splits_aktuell){
                thresh <- ordered_values[[i]][1:n_s[i]][j]
                obs0   <- which_obs[[count]][[v]][kn,]
                obs1   <- obs2 <- obs0
                obs1[DM_kov[,i]>thresh] <- NA
                obs2[DM_kov[,i]<=thresh] <- NA
                if(length(which(!is.na(obs1)))<bucket_min | length(which(!is.na(obs2)))<bucket_min){
                  splits_evtl[[count]][[v]][[i]][kn,j] <- NA
                }
              }
            }
          }
        }
      }
    }
    anysplit <- !all(is.na(unlist(splits_evtl[[count]])))

    if(anysplit){
      # estimate all models 
      dv <- lapply(1:nvar,function(var) {
              if(length(var_list[[var]])>0){
                ret <- lapply(var_list[[var]],function(i) {
                        n_knots   <- length(params[[count]][[var]])
                        deviances <- matrix(rep(0,n_s[i]*n_knots),ncol=n_knots)
                        if(!(var==nvar & i%in%int_del & !split_in_int)){ 
                          for(kn in 1:n_knots){
                            deviances[,kn] <- allmodels(var,i,kn,count,design_lower,design_upper,splits_evtl,params,dat0,mod0,n_s,var_names,smooth,split_intercept,family,start,mcs)
                          }
                        }
                        return(deviances)
                })
                names(ret) <- var_list[[var]]
              } else{
                ret <- 0 
              }
              return(ret)
            })
      
      if(max(unlist(dv))>0){
      
        # select optimum 
        variable          <- which.max(lapply(1:nvar,function(j) max(unlist(dv[[j]]))))
        split_variable    <- which.max(lapply(seq_along(var_list[[variable]]), function(j) max(dv[[variable]][[j]])))
        split_variable    <- names(dv[[variable]])[split_variable]
        split    <- as.numeric(which(dv[[variable]][[split_variable]]==max(dv[[variable]][[split_variable]]),arr.ind=TRUE)[,1])
        knoten   <- as.numeric(which(dv[[variable]][[split_variable]]==max(dv[[variable]][[split_variable]]),arr.ind=TRUE)[,2])
        if(length(split)>1){
          split  <- split[1]
          knoten <- knoten[1]
          warning(paste("Maximum in iteration ",count," not uniquely defined"))
        }
        param_old   <- params[[count]][[variable]][knoten]
        level       <- length(strsplit(param_old,":")[[1]])
        number      <- numbers[[count]][[variable]][knoten]
        left        <- max(numbers[[count]][[variable]])+1
        right       <- max(numbers[[count]][[variable]])+2
        
        param_new   <- paste(param_old,c(colnames(design_lower[[variable]][[split_variable]])[split],colnames(design_upper[[variable]][[split_variable]])[split]),sep=":")
        
        
        # compute permutation test 
        if(perm_test){
          dev <- rep(NA,nperm)
          
          for(perm in 1:nperm){
            dev[perm] <- one_permutation(variable,split_variable,knoten,count,design_lower,design_upper,
                                         DM_kov,thresholds,which_obs,splits_evtl,params,dat0,mod0,n_s,var_names,smooth,split_intercept,family,start,mcs)
            if(trace){
              cat(".")
            }
          }
          
          # test decision 
          adaption <- sum(!is.na(vars_evtl[[count]][[variable]][knoten,]))
          crit_val <- quantile(dev,1-(alpha/adaption))
          Tj       <- max(dv[[variable]][[split_variable]])
          proof    <- Tj > crit_val
          devs[count]    <- Tj
          crits[count]    <- crit_val
          pvalues[count] <- sum(dev>Tj)/nperm
        } else{
          proof <- TRUE
        }
        
        if(proof){
          
          # fit new model 
          mod0  <- mod_potential[[count+1]] <- one_model(variable,split_variable,knoten,count,split,design_lower,design_upper,params,dat0,var_names,smooth,split_intercept,family,start,mcs)
          dat0  <- data.frame(dat0,design_lower[[variable]][[split_variable]][,split,drop=FALSE],design_upper[[variable]][[split_variable]][,split,drop=FALSE])
          start <- predict(mod0)
          
          # adjust knoten 
          if(level>1){
            help_kn4 <- lu(c(),1,level-1,c())
            help_kn5 <- unlist(strsplit(param_old,""))
            help_kn6 <- paste0(help_kn5[which(help_kn5=="_")+1],collapse="")
            knoten2  <- which(help_kn4==help_kn6)
          } else{
            knoten2 <- knoten
          }
          
          splits[count,"variable"] <- variable
          splits[count,"effect_modifier"] <- split_variable
          splits[count,"split"] <- split
          splits[count,"level"] <- level
          splits[count,"node"]  <- knoten2
          splits[count,"threshold"] <- thresholds[[split_variable]][[split]]
          splits[count,c(7,8,9)] <- c(number,left,right)
          
          # generiere neue parameter 
          params[[count+1]]                                 <- params[[count]]
          params[[count+1]][[variable]]                     <- rep("",length(params[[count]][[variable]])+1)
          params[[count+1]][[variable]][c(knoten,knoten+1)] <- param_new
          params[[count+1]][[variable]][-c(knoten,knoten+1)]<- params[[count]][[variable]][-knoten]
          
          # passe splits_evtl an
          n_knots                                                       <- length(params[[count+1]][[variable]])
          splits_evtl[[count+1]]                                        <- splits_evtl[[count]]
          for(var in var_list[[variable]]){
            splits_evtl[[count+1]][[variable]][[var]]                       <- matrix(0,nrow=n_knots,ncol=n_s[var])
            splits_evtl[[count+1]][[variable]][[var]][c(knoten,knoten+1),]  <- matrix(rep(splits_evtl[[count]][[variable]][[var]][knoten,],2),nrow=2,byrow=T)
            splits_evtl[[count+1]][[variable]][[var]][-c(knoten,knoten+1),] <- splits_evtl[[count]][[variable]][[var]][-knoten,]
          }
          splits_evtl[[count+1]][[variable]][[split_variable]][knoten,splits_evtl[[count+1]][[variable]][[split_variable]][knoten,]>=split] <- NA 
          splits_evtl[[count+1]][[variable]][[split_variable]][(knoten+1),splits_evtl[[count+1]][[variable]][[split_variable]][(knoten+1),]<=split] <- NA
          
          # passe which_obs an 
          which_obs[[count+1]]                                   <- which_obs[[count]]
          which_obs[[count+1]][[variable]]                       <- matrix(0,nrow=n_knots,ncol=n)
          which_obs[[count+1]][[variable]][c(knoten,knoten+1),]  <- matrix(rep(which_obs[[count]][[variable]][knoten,],2),nrow=2,byrow=T)
          which_obs[[count+1]][[variable]][-c(knoten,knoten+1),] <- which_obs[[count]][[variable]][-knoten,]
          thresh <- ordered_values[[split_variable]][1:n_s[split_variable]][split]
          which_obs[[count+1]][[variable]][knoten,DM_kov[,split_variable]>thresh] <- NA
          which_obs[[count+1]][[variable]][(knoten+1),DM_kov[,split_variable]<=thresh] <- NA
          
          # minimal node size constraint
          if(length(which(!is.na(which_obs[[count+1]][[variable]][knoten,])))<nodesize_min[variable]){
            for(var in var_list[[variable]]){splits_evtl[[count+1]][[variable]][[var]][knoten,] <- rep(NA, n_s[var])}
          }
          if(length(which(!is.na(which_obs[[count+1]][[variable]][knoten+1,])))<nodesize_min[variable]){
            for(var in var_list[[variable]]){splits_evtl[[count+1]][[variable]][[var]][knoten+1,] <- rep(NA, n_s[var])}
          }
          
          # depth_max 
          if(!is.null(depth_max)){
            if(level==depth_max[variable]){
              for(var in var_list[[variable]]){splits_evtl[[count+1]][[variable]][[var]][c(knoten,knoten+1),] <- rep(NA, n_s[var])}
            }
          }
          
          # passe vars_evtl an 
          vars_evtl[[count+1]]                       <- vars_evtl[[count]]
          vars_evtl[[count+1]][[variable]]           <- matrix(0,nrow=n_knots,ncol=length(var_list[[variable]]))
          vars_evtl[[count+1]][[variable]][c(knoten,knoten+1),]  <- matrix(rep(vars_evtl[[count]][[variable]][knoten,],2),nrow=2,byrow=T)
          vars_evtl[[count+1]][[variable]][-c(knoten,knoten+1),] <- vars_evtl[[count]][[variable]][-knoten,]
          
          vars_evtl[[count+1]][[variable]][knoten,sapply(var_list[[variable]], function(var) all(is.na(unlist(splits_evtl[[count+1]][[variable]][[var]][knoten,]))))] <- NA
          vars_evtl[[count+1]][[variable]][knoten+1,sapply(var_list[[variable]], function(var) all(is.na(unlist(splits_evtl[[count+1]][[variable]][[var]][knoten+1,]))))] <- NA
          
          # any split? 
          anysplit <- !all(is.na(unlist(splits_evtl[[count+1]])))
          
          # add binary variables 
          if(!split_in_int){
            if(variable==nvar){
              split_in_int <- TRUE
            }
          }
          
          # passe numbers an 
          numbers[[count+1]]                                  <- numbers[[count]]
          numbers[[count+1]][[variable]]                      <- numeric(length=n_knots)
          numbers[[count+1]][[variable]][c(knoten,knoten+1)]  <- c(left,right)
          numbers[[count+1]][[variable]][-c(knoten,knoten+1)] <- numbers[[count]][[variable]][-knoten] 
          
          # trace
          if(trace){
            cat(paste0("\n Split"," ",count,";"," ","Variable"," ",variable,"\n"))
          }
          
          # erhoehe counter
          count <- count+1 
          
        } else{
          sig <- FALSE
        }
      } else{
        sig <- FALSE
      }
    }
  }
  mod_potential <- check_names_list(mod_potential, params)

  ################################################################################### 
  mod_opt     <- mod_potential[[count]]
  params_opt  <- params[[count]]
  
  if(count>1){
    effmod_vars   <- unique(splits[,1])
    noeffmod_vars <- c(1:nvar)[-effmod_vars]
    
    inmodel <- unique(c(unique(splits[,1]),which(var_names %in% splits[,2])))
    not_inmodel <- c(1:nvar)[-inmodel]
  } else{
    noeffmod_vars <- 1:nvar
    not_inmodel   <- 1:nvar 
  }

  # delete intercept
  if(split_intercept){
    not_inmodel <- not_inmodel[!not_inmodel==nvar]
  }
  
  if(perm_test){
    # permutations test for linear terms 
    not_sig <- logical(length(not_inmodel))
    done    <- FALSE
    if(is.null(exclude)){
      toTest <- not_inmodel
    } else{
      if(any(exclude%in%not_inmodel)){
        toTest <- not_inmodel[!(not_inmodel%in%exclude)]
      } else{
        toTest <- not_inmodel
      }
    }
    h1      <- 1:nvar
    
    pvalues_linear <- numeric(length(toTest))
    names(pvalues_linear) <- var_names[toTest]
    while(!done & length(toTest)>0){
      ps <- one_step(toTest,h1,var_names,params_opt,dat0,smooth,split_intercept,family,nperm,mcs)
      if(max(ps)<alpha){
        done <- TRUE
        pvalues_linear[toTest] <- ps
      } else{
        delvar  <- toTest[which.max(ps)]
        toTest  <- toTest[!(delvar==toTest)]
        not_sig[delvar==not_inmodel] <- TRUE
        h1 <- c(1:nvar)[-not_inmodel[not_sig]]
        pvalues_linear[which.max(ps)] <- max(ps)
      }
    }
  
    # determine final model 
    if(any(not_sig)){
      if(length(which(not_sig))<nvar){
        h1 <- c(1:nvar)[-not_inmodel[not_sig]]
        h2 <- unlist(params_opt[h1])
        h3 <- paste(h2,collapse="+")
        if(h3==""){
          h3 <- 1 
        }
      } else{
        h3 <- 1 
      }
      h4 <- formula(paste("y~",h3,ifelse(split_intercept,"-1","")))
      if(is.null(smooth)){
        mod_opt <- glm(h4,family=family,data=dat0)
      } else{
        h4  <- formula(paste("y~", h3, "+", help00,ifelse(split_intercept,"-1","")))
        mod_opt <- gam(h4,family=family,data=dat0)
      }
      if(h3!=1){
        mod_opt <- check_names(mod_opt, params_opt[h1])
      }
      beta_hat <- coef(mod_opt)
      noeffmod_vars <- noeffmod_vars[!(noeffmod_vars%in%not_inmodel[not_sig])]
    } else{
      beta_hat    <- coef(mod_opt)
      beta_hat[is.na(beta_hat)] <- 0 
    }
  } else{
    beta_hat    <- coef(mod_opt)
    beta_hat[is.na(beta_hat)] <- 0
    pvalues_linear <- NULL
  }
  
  # save results 
  if(count>1){
    beta_effmod          <- lapply(effmod_vars,function(j) beta_hat[params_opt[[j]]])
    names(beta_effmod)   <- var_names[effmod_vars]
    param_names   <- unlist(lapply(effmod_vars, function(j) params_opt[[j]]))
    beta_noeffmod <- beta_hat[!names(beta_hat)%in%param_names]
  } else{
    beta_noeffmod <- beta_hat
    beta_effmod   <- c()
  }
  
  to_return <- (list("model"=mod_opt,
                     "beta_noeffmod"=beta_noeffmod,
                     "beta_effmod"=beta_effmod,
                     "splits"=splits,
                     "pvalues"=pvalues,
                     "pvalues_linear"=pvalues_linear, 
                     "devs"=devs,
                     "crits"=crits,
                     "vl"=var_list))
  
  return(to_return)
  
}