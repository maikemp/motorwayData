library(grf)
library(rjson)
library(groupdata2)

library("rstudioapi")
this_file = getSourceEditorContext()$path 
source(paste0(this_file,"/../functions/rate_functions.R"))
source(paste0(this_file,"/../functions/cf_est_functions.R"))

# Cleans up space but slows down code..
gc()

path_in = ".../out/data/analysis_df.csv"
path_out = ".../data_out/analysis/rate/"

df_f <- read.csv(path_in, row.names=1)


cols <- get_colnames(df_f)

n_trees = 10000 # Set to 10000
n_trees_tune = 50 # Set to 50

# Only for testing:
# limit = "maxspeed_100"
# outcome = "total_rate"
# exp = "a"



# a:    Estimation is not cluster robust. 
# b:    Individual estimations are cluster robust but splitting is not. 
# b_cr: Estimations 

for (exp in c("ba")){ #, 
  print(paste("! !   - - - - - - - - - > - - - - > Start with:", exp))
  for (limit in c("maxspeed_100", "maxspeed_120", "maxspeed_130")){ 
    
    print(paste("! !    - - - - > Start with:", limit))
    
    for (outcome in c('light_rate', 'severe_rate', 'fatal_rate', 'total_rate')){
      
      set.seed(845213)
      
      lim <- substring(limit,nchar(limit)-2,nchar(limit))
      e = exp_str(exp)
      out_short <- short_outcome(outcome)
      
      print(paste("! !    - - > Start with:", outcome))
      
      df = df_f[(df_f['maxspeed_none']==1)|(df_f[limit]==1),]
      
      # # # #   # # #   # #   #       #  # #  #        #   # #   # # #   # # # # 
      # Activate following 5 lines to restrict according to prop.score 
      # restriction, conditional speed limit & no change in sl
      # path_out = "d:/speedlimit/data_out/analysis/rate/r_"
      # e = ""
      # df = df[(df[paste0("W_hat_", e, lim)]>0.05) &
      #           (df[paste0("W_hat_", e, lim)]<0.95) & (df["maxspeed_cond"]==0) &
      #           (df["maxspeed_change"]==0) & (df["overtaking_ht"]==0),]

      
      n <- dim(df)[1]
      df[,cols$cluster] <- factor(df[,cols$cluster])  

      
      # Only for b_cr, also separate training and test sample by cluster. 
      if (exp %in% c("a", "b", "ab", "ba")){
        train <- sample(1:n, n / 2)
      }else if(exp %in% c("a_cr", "b_cr")){
        train_c <- sort(sample(unique(df[,cols$cluster]), length(unique(df[,cols$cluster]))/2))
        train <- which(unlist(df[,cols$cluster] %in% train_c))
      }
      
      
      ###########
      # Draw folds & automatically decide if cluster robust. 
      # folds <- fold(df, k=5, id_col=colnames(C))[".folds"]
      
      # Cross fit stuff. 
      # for (k in 1:K){
        
      # Choose training and test set.
      W_train = data.matrix(subset(df, select=c(limit)))[train,]
      X_train = data.matrix(subset(df, select=cols$x_cols))[train,]
      Y_train = data.matrix(subset(df, select=c(outcome)))[train,]
      C_train = data.matrix(subset(df, select=c(cols$cluster)))[train,]
      
      W_test = data.matrix(subset(df, select=c(limit)))[-train,]
      X_test = data.matrix(subset(df, select=cols$x_cols))[-train,]
      Y_test = data.matrix(subset(df, select=c(outcome)))[-train,]
      C_test = data.matrix(subset(df, select=c(cols$cluster)))[-train,]
      rm(df)
      
  
      if (exp=="ab"){ # For "ab", nusiance functions are estimated without cluster robustness.
        C_nuis_train = NULL
        C_nuis_test = NULL}
      if (exp=="ba"){# For "ba", nusiance functions are estimated with cluster robustness.
        C_nuis_train = C_train
        C_nuis_test = C_test
        }
      
      # If a different clustering scheme is used for nuisance functions and cf,
      # fit the nusiance functions separately.
      if (exp %in% c("ba", "ab")){
        # Nuisance functions for the prioritisation forest (train)
        forest_y_train = regression_forest(X_train, Y_train, num.trees = n_trees,
                                           clusters = C_nuis_train, 
                                           tune.parameters = c("mtry", "min.node.size"),
                                           tune.num.trees = n_trees_tune)
        Y.hat.priority = forest_y_train$predictions
        forest_w_train = regression_forest(X_train, W_train, num.trees = n_trees,
                                           clusters = C_nuis_train, 
                                           tune.parameters = c("mtry", "min.node.size"),
                                           tune.num.trees = n_trees_tune)
        W.hat.priority = forest_w_train$predictions   
        
        # Nuisance functions for the evaluation forest (test)
        forest_y_test = regression_forest(X_test, Y_test, num.trees = n_trees,
                                          clusters = C_nuis_test, 
                                          tune.parameters = c("mtry", "min.node.size"),
                                          tune.num.trees = n_trees_tune)
        Y.hat.eval = forest_y_test$predictions
        forest_w_test = regression_forest(X_test, W_test, num.trees = n_trees,
                                          clusters = C_nuis_test, 
                                          tune.parameters = c("mtry", "min.node.size"),
                                          tune.num.trees = n_trees_tune)
        W.hat.eval = forest_w_test$predictions    
        
      }else{# If the nuisance functions have same clustering as cf, no need to pre-compute.
        Y.hat.priority = NULL
        W.hat.priority = NULL
        Y.hat.eval = NULL
        W.hat.eval = NULL
      }
        
      # For these experiments, cf are estimated without cluster robustness.
      # Remove clustering by setting C to NULL.
      if (exp %in% c("a", "a_cr", "ba")){ 
        C_train = NULL; C_test = NULL
      } 
      
      # Estimate priority ranking on the training sample. 
      cf.priority <- causal_forest(X_train, Y_train, W_train, clusters=C_train, 
                                   num.trees=n_trees, 
                                   tune.parameters=c("mtry", "min.node.size"),
                                   Y.hat=Y.hat.priority, W.hat=W.hat.priority)
        
      # Compute prioritization based on the estimated treatment effects.
      # -1: in this example the treatment should reduce the risk of crash.
      priority.cate <- -1 * predict(cf.priority, X_test)$predictions
      # rownames(priority.cate) <- rownames(X_test)



      #
      cf.eval <- causal_forest(X_test, Y_test, W_test, clusters=C_test,
                               num.trees=n_trees,
                               tune.parameters=c("mtry", "min.node.size"),
                               Y.hat=Y.hat.eval, W.hat=W.hat.eval)
      
      q <- seq(0.005, 1, by = 0.005)
      rate_qini <- rank_average_treatment_effect(cf.eval, priority.cate, 
                                            target="QINI", q=q)
      rate_autoc <- rank_average_treatment_effect(cf.eval, priority.cate, 
                                                  target="AUTOC", q=q)
      print(rate_qini)
      print(rate_autoc)

      # Plot the Targeting Operator Characteristic curve.
      plot(rate_qini, 
           main=paste("TOC - ",outcome," - ",limit," - ",exp))
      
      # Save out RATE results. 
      write.csv(c("estimate"=rate_qini$estimate, "std.err"=rate_qini$std.err, 
                  "target"=rate_qini$target), 
                paste0(path_out, "rate_qini_",exp,"_",lim,"_",out_short,".csv"))
      write.csv(c("estimate"=rate_autoc$estimate, "std.err"=rate_autoc$std.err, 
                  "target"=rate_autoc$target), 
                paste0(path_out, "rate_autoc_",exp,"_",lim,"_",out_short,".csv"))
      
      # Also save out TOC values (rather than plots themselves.)
      write.csv(rate_qini$TOC,
                paste0(path_out, "TOC_",exp,"_",lim,"_",out_short,".csv"))
      
        
    }
  }
}





