library(grf)
library(rstudioapi)


this_file = getSourceEditorContext()$path 

source(paste0(this_file,"/../functions/spatial_prediction_functions.R"))

path_out = ".../data_out/analysis/exp_results"
path_in = ".../out/data/df_prep.csv"

reps = 10


df_f = read.csv(path_in)

# Get relevant columns names.
id_col = "X"
aux_cols = colnames(df_f)[ifelse(startsWith(colnames(df_f), 'aux_'), TRUE, FALSE)]
y_cols = colnames(df_f)[ifelse(endsWith(colnames(df_f), '_rate'), TRUE, FALSE)]
sl_cols = colnames(df_f)[ifelse(startsWith(colnames(df_f), 'maxspeed'), TRUE, FALSE)]
w_cols = c("maxspeed_100", "maxspeed_120", "maxspeed_130" )
cl = 'location'
x_cols = colnames(df_f)[ifelse(!colnames(df_f) %in% c(y_cols, sl_cols, aux_cols, 
                                                      id_col, cl, 'overtaking_ht'), TRUE, FALSE)]

n_trees = 4000
n_trees_ffs = 150
n_trees_tune = 50


metric_names = c("mse", "accuracy", "precision", "recall", "f1_score", "kappa")

# limit = 'maxspeed_120'


### a # All variables, random CV
for (limit in c('maxspeed_130','maxspeed_120','maxspeed_100')){
  
  set.seed(845213)
  
  metrics_a = data.frame("metric"=metric_names) 
  var_imp_df_a = data.frame("var"=c(x_cols))
  metrics_b = data.frame("metric"=metric_names) 
  var_imp_df_b = data.frame("var"=c(x_cols))
  metrics_c = data.frame("metric"=metric_names) 
  var_imp_df_c = data.frame("var"=c(x_cols))

  cvres = list()
  
  lim = substring(limit,nchar(limit)-2,nchar(limit))
  print(paste("----------------------------------------------------------------------"))
  print(paste(limit))
  print(paste("----------------------------------------------------------------------"))

  exp = "a"
  for (rep in 0:(reps-1)){
    print(paste(exp, rep))
    # Set holdout_set aside.
    df = df_f[(!is.na(df_f[paste0('aux_holdout_', lim,'_',rep)])),]
    df_val = df[(df[paste0('aux_holdout_', lim,'_',rep)]==1),]
    df_t = df[(df[paste0('aux_holdout_', lim,'_',rep)]==0),]
    print(paste(limit,'-> validation size:',dim(df_val), 'test size:', dim(df_t)))
    rm(df)
    
    # Training data
    W = data.matrix(subset(df_t, select=c(limit)))
    X = data.matrix(subset(df_t, select=x_cols))
    C = NULL
    # Validation data.
    W_val = data.matrix(subset(df_val, select=c(limit)))
    X_val = data.matrix(subset(df_val, select=x_cols))
    rm(df_t, df_val)


    ### Exp (a) Random CV, all vars
    # Start with random CV, all vars
    forest_a = regression_forest(X, W, num.trees = n_trees,
                                 sample.weights = NULL,
                                 clusters = NULL, equalize.cluster.weights = FALSE,
                                 sample.fraction = 0.5,
                                 tune.parameters = c("mtry", "min.node.size"),
                                 tune.num.trees = n_trees_tune, tune.num.reps = 100,
                                 tune.num.draws = 1000)

    cvres = rbind(cvres, c("Exp"=exp, "mtry"=forest_a$tunable.params$mtry,
                           "min.node.size"=forest_a$tunable.params$min.node.size))

    model_val_a = get_errors_and_var_imp(forest_a, W, W_val, X_val,
                                         outcome="W", rep=rep)
    metrics_a = merge(metrics_a, model_val_a$metrics, by="metric", all=TRUE)
    var_imp_df_a = merge(var_imp_df_a, model_val_a$var_imp_df, by="var", all=TRUE)

    if (rep == (reps-1)){
      avrg_and_save_results(lim, NA, NA, var_imp_df_a, path_out, exp, reps,
                            outcome="W", metrics=metrics_a)
    }
    rm(forest_a)
  }

  # Reset seed for more independent way of re-running parts.
  set.seed(845213)
  ### b # Spatial CV, all variables
  exp = "b"
  for (rep in 0:(reps-1)){
    print(paste(exp, rep))
    # Set holdout_set aside.
    df = df_f[(!is.na(df_f[paste0('aux_holdout_', lim,'_',rep)])),]
    df_val = df[(df[paste0('aux_holdout_', lim,'_',rep)]==1),]
    df_t = df[(df[paste0('aux_holdout_', lim,'_',rep)]==0),]
    print(paste(limit,'-> validation size:',dim(df_val), 'test size:', dim(df_t)))
    rm(df)
    
    # Training data
    W = data.matrix(subset(df_t, select=c(limit)))
    X = data.matrix(subset(df_t, select=x_cols))
    C = data.matrix(subset(df_t, select=c(cl)))
    # Validation data.
    W_val = data.matrix(subset(df_val, select=c(limit)))
    X_val = data.matrix(subset(df_val, select=x_cols))
    rm(df_t, df_val)

    ### Exp (b) Spatial CV, all vars
    forest_b = regression_forest(X, W, num.trees = n_trees, sample.weights = NULL,
                                 clusters = C, equalize.cluster.weights = FALSE,
                                 sample.fraction = 0.5,
                                 tune.parameters = c("mtry", "min.node.size"),
                                 tune.num.trees = n_trees_tune, tune.num.reps = 100,
                                 tune.num.draws = 1000)
    cvres = rbind(cvres, c("Exp"=exp, "mtry"=forest_b$tunable.params$mtry,
                           "min.node.size"=forest_b$tunable.params$min.node.size))

    model_val_b = get_errors_and_var_imp(forest_b, W, W_val, X_val,
                                         outcome="W", rep=rep)
    metrics_b = merge(metrics_b, model_val_b$metrics, by="metric", all=TRUE)
    var_imp_df_b = merge(var_imp_df_b, model_val_b$var_imp_df, by="var", all=TRUE)

    if (rep == (reps-1)){
      avrg_and_save_results(lim, NA, NA, var_imp_df_b, path_out, exp, reps,
                            outcome="W", metrics=metrics_b)
    }
    rm(forest_b)
  }
  
  # Reset seed for more independent way of re-running parts.
  set.seed(845213)  
  ## c # Spatial CV, FFS
  exp = "c"
  for (rep in 0:(reps-1)){
    print(paste(exp, rep))
    # Set holdout_set aside.
    df = df_f[(!is.na(df_f[paste0('aux_holdout_', lim,'_',rep)])),]
    df_val = df[(df[paste0('aux_holdout_', lim,'_',rep)]==1),]
    df_t = df[(df[paste0('aux_holdout_', lim,'_',rep)]==0),]
    rm(df)
    print(paste(limit,'-> validation size:',dim(df_val), 'test size:', dim(df_t)))

    # Training data
    W = data.matrix(subset(df_t, select=c(limit)))
    X = data.matrix(subset(df_t, select=x_cols))
    C = data.matrix(subset(df_t, select=c(cl)))
    # Validation data.
    W_val = data.matrix(subset(df_val, select=c(limit)))
    X_val = data.matrix(subset(df_val, select=x_cols))
    rm(df_t, df_val)

    ### Exp (c) Spatial CV, FFS
    ffs_vars = run_ffs_regression_forest(X, W, x_cols, n_trees_ffs,
                                         clusters=C, rep=rep)

    # Fit final model and validate on holdout set.
    forest_c = regression_forest(subset(X, select=ffs_vars), W, num.trees=n_trees,
                                 clusters=C,
                                 tune.parameters = c("mtry", "min.node.size"),
                                 tune.num.trees = n_trees_tune, tune.num.reps = 100,
                                 tune.num.draws = 1000)
    cvres = rbind(cvres, c("Exp"=exp, "mtry"=forest_c$tunable.params$mtry,
                           "min.node.size"=forest_c$tunable.params$min.node.size))

    model_val_c = get_errors_and_var_imp(forest_c, W, W_val,
                                         subset(X_val, select=ffs_vars),
                                         outcome="W", rep=rep)
    metrics_c = merge(metrics_c, model_val_c$metrics, by="metric", all=TRUE)
    var_imp_df_c = merge(var_imp_df_c, model_val_c$var_imp_df, by="var", all=TRUE)

    if (rep == (reps-1)){
      avrg_and_save_results(lim, NA, NA, var_imp_df_c, path_out, exp, reps,
                            outcome="W", metrics=metrics_c)
    }
    rm(forest_c)
  }


  write.csv(cvres, file=paste0(path_out,"/CV_results_ps_",lim,".csv"))
}

