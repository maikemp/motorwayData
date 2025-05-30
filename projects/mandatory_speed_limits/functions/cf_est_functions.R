

# Copied from: https://grf-labs.github.io/grf/articles/ate_transport.html
average_treatment_effect_test <- function(forest, tau.hat.train, X.test, p.hat, 
                                          Y.orig, W.orig){
  # TODO: switch all these print messages to real 
  Y.hat <- forest$Y.hat
  W.hat <- forest$W.hat
  X.orig <- forest$X.orig
  n.rct <- length(Y.orig)
  n.obs <- nrow(X.test)
  
  # Checkup that all rownames make sense. 
  Y.hat <- Y.hat[rownames(p.hat), , drop=F]
  Y.orig <- Y.orig[rownames(p.hat), , drop=F]
  W.hat <- W.hat[rownames(p.hat), , drop=F]
  W.orig <- W.orig[rownames(p.hat), , drop=F]
  X.orig <- X.orig[rownames(p.hat), , drop=F]
  

  if (!(sum(rownames(Y.hat)==rownames(Y.orig))==length(Y.orig) & 
        sum(rownames(W.hat)==rownames(Y.orig))==length(Y.orig) &
        sum(rownames(W.orig)==rownames(Y.orig))==length(Y.orig) &
        sum(rownames(X.orig)==rownames(Y.orig))==length(Y.orig) &
        sum(rownames(p.hat)==rownames(Y.orig))==length(Y.orig))){
    print(paste("Rownames of the Data in estimation of average_treatment_effect_test don't match"))
    
  }
  
  tau.hat.test <- predict(forest, X.test)$predictions
  debiasing.weights <- (1 - p.hat) / p.hat * (W.orig - W.hat) / (W.hat * (1 - W.hat))
  Y.residual <- Y.orig - (Y.hat + tau.hat.train * (W.orig - W.hat))
  
  # Save all tau.hat for possible re-calculations later. 
  tau.hat.all <- as.data.frame(c(tau.hat.train, tau.hat.test))
  rownames(tau.hat.all) <- rownames(rbind(X.orig, X.test))
  
  tau.hat <- mean(tau.hat.test) + sum(debiasing.weights * Y.residual) / n.obs
  sigma2.hat <- var(tau.hat.test) / n.obs + sum(Y.residual^2 * (debiasing.weights / n.obs)^2)
  
  return(list("estimate" = tau.hat, 
              "std_err" = sqrt(sigma2.hat),
              "tau_hat_all" = tau.hat.all))
}



short_outcome <- function(outcome){
  
  if (outcome == 'fatal_rate'){
    out_short = 'fr'
  }else if (outcome == 'severe_rate'){
    out_short = 'sr'
  }else if (outcome == 'light_rate'){
    out_short = 'lr'
  }else{
    out_short = 't'
  }
  
  return(out_short)
}



exp_str <- function(exp){
  if (exp %in% c('b','ab1', 'ab2', 'ba', 'b_res')){
    exp_str = ''
  }else if (exp == 'ab'){
    exp_str="a_"
  }else{
    exp_str = paste0(exp, '_')
  }
  return(exp_str)
}


