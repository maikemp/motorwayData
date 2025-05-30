




mse = function(Y, Y_hat){
  return(mean((Y - Y_hat)**2, na.rm=TRUE))
}


kappa <- function(m) {
  N <- sum(m)
  No <- sum(diag(m))
  Ne <- 1 / N * sum(colSums(m) * rowSums(m))
  return( (No - Ne) / (N - Ne) )
}

confusion_matrix <- function(Y.hat, Y, a){
  
  cm = table(Y.hat, Y)
  
  if (!"0" %in% rownames(cm)){
    new_row = as.data.frame(list(0, 0), row.names=c("0"), col.names=c("0","1"),
                            check.names=FALSE)
    cm = rbind(cm, as.matrix(new_row))
  }
  if (!"1" %in% rownames(cm)){
    new_row = as.data.frame(list(0, 0), row.names=c("1"), col.names=c("0","1"),
                            check.names=FALSE)
    cm = rbind(cm, as.matrix(new_row))
  }
  dm_list = list()
  dm_list[[paste0("Y.hat.",a)]] = c("0", "1")
  dm_list[[paste0("Y.",a)]] = c("0", "1")
  dimnames(cm) = dm_list
  
  return(cm)
}

get_errors_and_var_imp = function(forest, Y, Y_val, X_val, outcome="Y", rep=NA,
                                  r=4){
  
  # Calculate in-sample error. 
  predictions.oob = predict(forest)
  Y.hat.oob = predictions.oob$predictions
  
  mse.oob = round(mse(Y, Y.hat.oob),r)
  
  # Calculate error in hold-out set. 
  predictions.val = predict(forest, X_val)
  Y.hat.val = predictions.val$predictions
  
  mse.val = round(mse(Y_val, Y.hat.val),r)
  
  # Get variable importance
  var_imp <- variable_importance(forest)
  var_imp_df = as.data.frame(cbind(var=colnames(X_val), imp=round(var_imp,5)))
  var_imp_df = var_imp_df[order(as.numeric(var_imp_df$V2)*(-1)),]
  
  if (outcome=="Y"){
    return_arg = list("mse" = mse.oob, "mse.val" = mse.val, "var_imp_df" = var_imp_df)
  }else if (outcome=="W"){
    
    metrics = c("metric"="mse", "oob"=mse.oob, "val"=mse.val)
    
    # Assigned class labels when rounding. 
    Y.class.hat.oob = round(Y.hat.oob) # weirdly round down 0.5... may change this.
    Y.class.hat.val = round(Y.hat.val)
    
    # Confusion matrix for oob predictions.
    cm.oob = confusion_matrix(Y.class.hat.oob, Y, "oob")
    cm_df.oob = as.data.frame(as.table(cm.oob))
    # confusion matrix for validation set. 
    cm.val = confusion_matrix (Y.class.hat.val, Y_val, "val")
    cm_df.val = as.data.frame(as.table(cm.val))
    
    acc.oob = sum(diag(cm.oob))/sum(cm.oob)
    acc.val = sum(diag(cm.val))/sum(cm.val)
    metrics = rbind(metrics, c("metric"="accuracy", "oob"=round(acc.oob,r), 
                               "val"=round(acc.val,r)))
    
    prec.oob = (cm_df.oob[(cm_df.oob$Y.oob==1)&(cm_df.oob$Y.hat.oob==1),]$Freq / 
                   sum(cm_df.oob[(cm_df.oob$Y.hat.oob==1),]$Freq))
    prec.val = (cm_df.val[(cm_df.val$Y.val==1)&(cm_df.val$Y.hat.val==1),]$Freq / 
                   sum(cm_df.val[(cm_df.val$Y.hat.val==1),]$Freq))
    metrics = rbind(metrics, c("metric"="precision", "oob"=round(prec.oob,r), 
                               "val"=round(prec.val,r))) 
    
    recall.oob = (cm_df.oob[(cm_df.oob$Y.oob==1)&(cm_df.oob$Y.hat.oob==1),]$Freq /
                sum(cm_df.oob[(cm_df.oob$Y.oob==1),]$Freq))
    recall.val = (cm_df.val[(cm_df.val$Y.val==1)&(cm_df.val$Y.hat.val==1),]$Freq /
                sum(cm_df.val[(cm_df.val$Y.val==1),]$Freq))
    metrics = rbind(metrics, c("metric"="recall", "oob"=round(recall.oob,r), 
                               "val"=round(recall.val,r)))
    
    f1.oob = 2*((prec.oob*recall.oob)/(prec.oob+recall.oob))
    f1.val = 2*((prec.val*recall.val)/(prec.val+recall.val))
    metrics = rbind(metrics, c("metric"="f1_score", "oob"=round(f1.oob,r), 
                               "val"=round(f1.val,r)))
    
    kappa.oob = kappa(cm.oob)
    kappa.val = kappa(cm.val)
    metrics = rbind(metrics, c("metric"="kappa", "oob"=round(kappa.oob,r), 
                               "val"=round(kappa.val,r)))
    
    colnames(metrics) = c("metric", paste0("oob_", rep), paste0("val_",rep))
    
    return_arg = list("var_imp_df" = var_imp_df, "metrics"=metrics)
  }
  
  return(return_arg)
}


avrg_and_save_results = function(n,mse.oob, mse_val, var_imp_df, path_out, exp, 
                                 reps, outcome="Y", metrics=NA){
  
  if (outcome=="Y"){
    mse_fin = round(mean(mse.oob), 6)
    mse_val_fin = round(mean(mse_val), 6)
    write.csv(c("avrg. oob MSE"=mse_fin, "avrg. Validation MSE"=mse_val_fin, 
                "oob MSE"=mse.oob, "val MSE" =mse_val),
              file=paste0(path_out,"/mse_", exp,"_me_",n,".csv"))
    e = "_me"
    
  }else if (outcome=="W"){
    cn = colnames(metrics)
    metrics_fin = cbind(metrics,
                        "oob_mean"=round(rowMeans(sapply(metrics[,cn[ifelse(startsWith(cn,"oob"), TRUE, FALSE)]], as.numeric)),5),
                        "val_mean"=round(rowMeans(sapply(metrics[,cn[ifelse(startsWith(cn,"val"), TRUE, FALSE)]], as.numeric)),5))
    write.csv(metrics_fin,
              file=paste0(path_out,"/metrics_", exp,"_ps_",n,".csv"))
    e = "_ps"
    
  }

  
  # If there is any non-NA value in a row, replace NAs by zeros to calculate 
  # average variable importance.
  var_imp_df_2 = var_imp_df[rowSums(sapply(var_imp_df,is.na)*1)<reps,]
  var_imp_df_2[is.na(var_imp_df_2)] = 0
  var_imp_avrg = cbind("var"=var_imp_df_2[,1],
                       "avrg_var_imp"=rowMeans(
                         sapply(var_imp_df_2[,2:length(var_imp_df_2)], 
                                as.numeric)))
  # Merge average with original variable importance to keep overview of 
  # non-selected variables. 
  var_imp_fin = merge(var_imp_avrg, var_imp_df, by="var", all=TRUE)
  var_imp_fin$avrg_var_imp = round(as.numeric(var_imp_fin$avrg_var_imp),4)
  var_imp_fin = var_imp_fin[order(as.numeric(var_imp_fin$avrg_var_imp)*(-1)),]
  
  
  write.csv(var_imp_fin, file=paste0(path_out,"/var_imp_", exp, e, "_",n,".csv"))
  
  print(paste("Saved out results of Experiment", exp, "To", path_out))
}



run_ffs_regression_forest = function(X, Y, x_cols, n_trees_ffs, 
                                     clusters, rep=NULL){
  # Get all possible combinations.
  combs = list()
  for (i in 1:length(x_cols)){
    for (j in 1:length(x_cols)){
      if (i<j){
        combs = c(combs, list(list(x_cols[i], x_cols[j])))
      }
    }
  }
  
  ffs_res = c()
  for (i in combs){
    comb = c(i[[1]], i[[2]])
    forest_c = regression_forest(subset(X, select=comb), Y, num.trees=n_trees_ffs,
                                 clusters=clusters)
    Y.hat.oob = predict(forest_c)$predictions
    mse_f = mse(Y, Y.hat.oob)
    ffs_res = rbind(ffs_res, c(comb, mse_f))
  }
  
  # Take very large starting value for the error. 
  current_error = 99999999
  # Find best combination.
  min_error = min(ffs_res[,3])
  best_comb = ffs_res[ffs_res[,3]==min_error][1:2]
  # Fit forest with best combination. 
  forest_c = regression_forest(subset(X, select=c(best_comb)), 
                               Y, num.trees=n_trees_ffs, clusters=clusters)
  # Get MSE with the oob estimates. 
  mse_f = mse(Y, predict(forest_c)$predictions)
  
  while (min_error<current_error){
    
    # Refit the best model. 
    forest_c = regression_forest(subset(X, select=c(best_comb)), 
                                 Y, num.trees=n_trees_ffs, clusters=clusters)
    current_error = mse(Y, predict(forest_c)$predictions)
    print(paste(c(best_comb)))
    print(paste(round(current_error,5)))
    # Prepare new iteration of combinations. 
    combs = list()
    for (i in x_cols){
      if (! i %in% best_comb){
        combs = c(combs, list(c(best_comb, i)))
      }
    }
    
    # Run new iteration. 
    if (length(combs)>0){
      
      
      ffs_res = c()
      for (comb in combs){
        forest_c = regression_forest(subset(X, select=comb), Y, num.trees=n_trees_ffs,
                                     clusters=clusters)
        Y.hat.oob = predict(forest_c)$predictions
        mse_f = mse(Y, Y.hat.oob)
        ffs_res = rbind(ffs_res, c(comb, mse_f))
      }
      
      # Find best new combination 
      min_error = min(ffs_res[,ncol(ffs_res)])
      best_comb = ffs_res[ffs_res[,ncol(ffs_res)]==min_error][1:ncol(ffs_res)-1]
      final_comb = best_comb[1:length(best_comb)-1]
    }else{
      final_comb = best_comb
    }
  }
  
  return(final_comb)
  
}



run_ffs_regression_forest_r = function(X, Y, x_cols, n_trees_ffs, 
                                       clusters, rep){
  
  combs = list()
  for (i in 1:length(x_cols)){
    for (j in 1:length(x_cols)){
      if (i<j){
        combs = c(combs, list(list(x_cols[i], x_cols[j])))
      }
    }
  }
  
  ffs_res = c()
  for (i in combs){
    comb = c(i[[1]], i[[2]])
    forest_c = regression_forest(subset(X, select=comb), Y, num.trees=n_trees_ffs,
                                 clusters=clusters)
    Y.hat.oob = predict(forest_c)$predictions
    mse_f = mse(Y, Y.hat.oob)
    ffs_res = rbind(ffs_res, c(comb, mse_f))
  }
  
  # Take very large starting value for the error. 
  current_error = 99999999
  # Find best combination.
  min_error = min(ffs_res[,3])
  best_comb = ffs_res[ffs_res[,3]==min_error][1:2]
  # Fit forest with best combination and one additional random variable. 
  forest_c = regression_forest(subset(X, select=c(best_comb, paste0("aux_random_", rep))), 
                               Y, num.trees=n_trees_ffs, clusters=clusters)
  # Get MSE with the oob estimates. 
  mse_f = mse(Y, predict(forest_c)$predictions)
  
  while (min_error<current_error){
    
    # Refit the best model with one additional random variable to 
    forest_c = regression_forest(subset(X, select=c(best_comb, paste0("aux_random_", rep))), 
                                 Y, num.trees=n_trees_ffs, clusters=clusters)
    current_error = mse(Y, predict(forest_c)$predictions)
    print(paste(c(best_comb, paste0("aux_random_", rep))))
    print(paste(round(current_error,5)))
    # Prepare new iteration of combinations. 
    combs = list()
    for (i in x_cols){
      if (! i %in% best_comb){
        combs = c(combs, list(c(best_comb, i)))
      }
    }
    
    # Run new iteration. 
    if (length(combs)>0){
      

      ffs_res = c()
      for (comb in combs){
        forest_c = regression_forest(subset(X, select=comb), Y, num.trees=n_trees_ffs,
                                     clusters=clusters)
        Y.hat.oob = predict(forest_c)$predictions
        mse_f = mse(Y, Y.hat.oob)
        ffs_res = rbind(ffs_res, c(comb, mse_f))
      }
      
      # Find best new combination 
      min_error = min(ffs_res[,ncol(ffs_res)])
      best_comb = ffs_res[ffs_res[,ncol(ffs_res)]==min_error][1:ncol(ffs_res)-1]
      final_comb = best_comb[1:length(best_comb)-1]
    }else{
      final_comb = best_comb
    }
  }
  
  
  return(final_comb)
  
}

