library(grf)
library(rstudioapi)


this_file = getSourceEditorContext()$path 

source(paste0(this_file,"/../functions/spatial_prediction_functions.R"))

path_out = ".../data_out/analysis/exp_results"
path_in = ".../out/data/df_prep.csv"

reps = 10


df = read.csv(path_in)

# Get relevant column names.
id_col = "X"
aux_cols = colnames(df)[ifelse(startsWith(colnames(df), 'aux_'), TRUE, FALSE)]
y_cols = colnames(df)[ifelse(endsWith(colnames(df), '_rate'), TRUE, FALSE)]
sl_cols = colnames(df)[ifelse(startsWith(colnames(df), 'maxspeed'), TRUE, FALSE)]
w_cols = c("maxspeed_100", "maxspeed_120", "maxspeed_130" )
c = 'location'
x_cols = colnames(df)[ifelse(!colnames(df) %in% c(y_cols, sl_cols, aux_cols, 
                                                  id_col, c, 'overtaking_ht'), TRUE, FALSE)]
# As we can only evaluate the actual outcome, not the outcome, marginalizing 
# over anything, evaluate by including treatment as an ordered variable. 
df["maxspeed"] <- data.frame(as.matrix(df[w_cols]) %*% c(-3,-2,-1)+4)
x_cols <- c(x_cols, "maxspeed")

# Increase for the final run!! # - ! - # - ! - # - ! - # - ! - # - ! - # - ! - # 
n_trees = 4000
n_trees_ffs = 150
n_trees_tune = 50


#outcome = 'fatal_rate'
#out_short = 'fr'

for (outcome in c('fatal_rate', 'severe_rate', 'light_rate')){
  
  set.seed(845213)
  
  print(paste("-----------",outcome, "-----------"))
  print(paste("----------------------",outcome, "----------------------"))
  print(paste("-----------",outcome, "-----------"))
  
  if (outcome == 'fatal_rate'){
    out_short = 'fr'
  }else if (outcome == 'severe_rate'){
    out_short = 'sr'
  }else if (outcome == 'light_rate'){
    out_short = 'lr'
  }else if (outcome == 'total'){
    out_short = 't'
  }else{
    out_short = "unknown_out_short"
  }
  
    
  mse_a = c()
  mse_val_a = c()
  var_imp_df_a = data.frame("var"=c(x_cols))
  mse_b = c()
  mse_val_b = c()
  var_imp_df_b = data.frame("var"=c(x_cols))
  mse_c = c()
  mse_val_c = c()
  var_imp_df_c = data.frame("var"=c(x_cols))
  
  cvres = list()
  
  
  ### a #
  exp = "a"
  for (rep in 0:(reps-1)){
    print(paste(exp, rep))
    # Set holdout_set aside.
    df_val = df[df[paste0('aux_holdout_',rep)]==1,]
    df_t = df[df[paste0('aux_holdout_',rep)]==0,]
    # Training data
    Y = data.matrix(subset(df_t, select=c(outcome)))
    X = data.matrix(subset(df_t, select=x_cols))
    C = data.matrix(subset(df_t, select=c(c)))
    # Validation data.
    Y_val = data.matrix(subset(df_val, select=c(outcome)))
    X_val = data.matrix(subset(df_val, select=x_cols))
    
    
    ### Exp (a) Random CV, all vars  
    # Start with random CV, all vars
    forest_a = regression_forest(X, Y, num.trees=n_trees, 
                                 sample.weights=NULL,
                                 clusters=NULL, equalize.cluster.weights=FALSE, 
                                 sample.fraction=0.5,
                                 tune.parameters=c("mtry", "min.node.size"),
                                 tune.num.trees=n_trees_tune, tune.num.reps=100, 
                                 tune.num.draws=1000)
    
    cvres = rbind(cvres, c("Exp"="a", "mtry"=forest_a$tunable.params$mtry,
                           "min.node.size"=forest_a$tunable.params$min.node.size))
    model_val_a = get_errors_and_var_imp(forest_a, Y, Y_val, X_val)
    mse_a = c(mse_a, model_val_a$mse)
    mse_val_a = c(mse_val_a, model_val_a$mse.val)
    var_imp_df_a = merge(var_imp_df_a, model_val_a$var_imp_df, by="var", all=TRUE)
    
    if (rep == (reps-1)){
      avrg_and_save_results(out_short, mse_a, mse_val_a, var_imp_df_a, path_out, 
                            "a", reps)
    }
    rm(forest_a)
  }
  
  # Reset seed for more independent way of re-running parts.
  set.seed(845213)
  ### b # 
  exp = "b"
  for (rep in 0:(reps-1)){
    print(paste(exp, rep))
    df_val = df[df[paste0('aux_holdout_',rep)]==1,]
    df_t = df[df[paste0('aux_holdout_',rep)]==0,]
    # Training data
    Y = data.matrix(subset(df_t, select=c(outcome)))
    X = data.matrix(subset(df_t, select=x_cols))
    C = data.matrix(subset(df_t, select=c(c)))
    # Validation data.
    Y_val = data.matrix(subset(df_val, select=c(outcome)))
    X_val = data.matrix(subset(df_val, select=x_cols))
    
    ### Exp (b) Spatial CV, all vars  
    forest_b = regression_forest(X, 
                                 Y, num.trees=n_trees, sample.weights=NULL,
                                 clusters=C, equalize.cluster.weights=FALSE, 
                                 sample.fraction=0.5,
                                 tune.parameters=c("mtry", "min.node.size"),
                                 tune.num.trees=n_trees_tune, tune.num.reps=100, 
                                 tune.num.draws=1000)
    cvres = rbind(cvres, c("Exp"="b", "mtry"=forest_b$tunable.params$mtry,
                           "min.node.size"=forest_b$tunable.params$min.node.size))
    
    model_val_b = get_errors_and_var_imp(forest_b, Y, Y_val, X_val)
    mse_b = c(mse_b, model_val_b$mse)
    mse_val_b = c(mse_val_b, model_val_b$mse.val)
    var_imp_df_b = merge(var_imp_df_b, model_val_b$var_imp_df, by="var", all=TRUE)  
    
    if (rep == (reps-1)){
      avrg_and_save_results(out_short, mse_b, mse_val_b, var_imp_df_b, path_out, 
                            "b", reps)
    }
    rm(forest_b)
  }
  
  # Reset seed for more independent way of re-running parts.
  set.seed(845213)
  ### c # 
  exp = "c"
  for (rep in 0:(reps-1)){
    print(paste(exp, rep))
    df_val = df[df[paste0('aux_holdout_',rep)]==1,]
    df_t = df[df[paste0('aux_holdout_',rep)]==0,]
    # Training data
    Y = data.matrix(subset(df_t, select=c(outcome)))
    X = data.matrix(subset(df_t, select=x_cols))
    C = data.matrix(subset(df_t, select=c(c)))
    # Validation data.
    Y_val = data.matrix(subset(df_val, select=c(outcome)))
    X_val = data.matrix(subset(df_val, select=x_cols))
    
    ### Exp (c) Spatial CV, FFS
    ffs_vars = run_ffs_regression_forest(X, Y, x_cols, n_trees_ffs, 
                                         clusters=C, rep=rep)
    
    # Fit final model and validate on holdout set. 
    forest_c = regression_forest(subset(X, select=ffs_vars), 
                                 Y, num.trees=n_trees, clusters=C, 
                                 tune.parameters=c("mtry", "min.node.size"), 
                                 tune.num.trees=n_trees_tune, tune.num.reps=100, 
                                 tune.num.draws=1000)
    cvres = rbind(cvres, c("Exp"="c", "mtry"=forest_c$tunable.params$mtry,
                           "min.node.size"=forest_c$tunable.params$min.node.size))
    
    model_val_c = get_errors_and_var_imp(forest_c, Y, Y_val, 
                                         subset(X_val, select=ffs_vars))
    mse_c = c(mse_c, model_val_c$mse)
    mse_val_c = c(mse_val_c, model_val_c$mse.val)
    var_imp_df_c = merge(var_imp_df_c, model_val_c$var_imp_df, by="var", all=TRUE)
    
    if (rep == (reps-1)){
      avrg_and_save_results(out_short, mse_c, mse_val_c, var_imp_df_c, path_out, 
                            "c", reps)
    }
    rm(forest_c)
  }
  
  # Save out tuning parameters.
  write.csv(cvres, file=paste0(path_out,"/CV_results_me_",out_short,".csv"))
  
}
  







