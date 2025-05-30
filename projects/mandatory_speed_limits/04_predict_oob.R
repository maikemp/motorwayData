## Implement estimation of Y.hat for all outcomes. 
# 1st implement Y.hat without clustering. Y_hat_a
# 2st implement Y.hat only with clustering. Y_hat
# No holdout set required. 
library(grf)
library(rstudioapi)


this_file = getSourceEditorContext()$path 

source(paste0(this_file,"/../functions/spatial_prediction_functions.R"))
source(paste0(this_file,"/../functions/cf_est_functions.R"))


path_in = ".../out/data/df_prep.csv"

path_out = ".../data_out/analysis/oob_predictions/"
df_f = read.csv(path_in)
rownames(df_f) = df_f$X
gc()


path_max_ext = '.../out/data/df_max_ext_imputed.csv'


# Get relevant columns names.
id_col = "X"
aux_cols = colnames(df_f)[ifelse(startsWith(colnames(df_f), 'aux_'), TRUE, FALSE)]
y_cols = colnames(df_f)[ifelse(endsWith(colnames(df_f), '_rate'), TRUE, FALSE)]
sl_cols = colnames(df_f)[ifelse(startsWith(colnames(df_f), 'maxspeed'), TRUE, FALSE)]
w_cols = c("maxspeed_100", "maxspeed_120", "maxspeed_130" )
cl = 'location'
x_cols = colnames(df_f)[ifelse(!colnames(df_f) %in% c(y_cols, sl_cols, aux_cols, 
                                                      id_col, cl, 'overtaking_ht'), TRUE, FALSE)]


n_trees = 10000 # Set to 10000
n_trees_ffs = 200 # Set to 200
n_trees_tune = 50 # Set to 100.

# Make dataframes
X = data.matrix(subset(df_f, select=c(x_cols)))
C = data.matrix(subset(df_f, select=c(cl)))

df_maxext = read.csv(path_max_ext, row.names=1)
df_maxext = df_maxext[,c(x_cols, 'location', 'maxspeed_none')]

exps <- c("b")

# predict Y.hat. 
for (outcome in c('total_rate','light_rate', 'severe_rate', 'fatal_rate')){

  set.seed(845213)
  tun_p <- data.frame("target"=outcome, "exp"=exps, "mtry"=NA, "min.node.size"=NA)
  
  out_short = short_outcome(outcome)
  
  # Set up dataframes.
  out_data = data.frame("id" = df_f$X)
  rownames(out_data) = df_f$X
  # ffs_vars_all = data.frame("var"=c(x_cols))
  
  Y = data.matrix(subset(df_f, select=c(outcome)))
  
  # Breaks loop but may allow for process to go through?
  # rm(df_f, x_cols, aux_cols, y_cols, sl_cols, w_cols)
  
  #_________________________________________________________________________####
  #   a                                                                     ####
  # Make naive forest without clustering. 
  forest_a_y = regression_forest(X, Y, num.trees = n_trees,
                                 clusters = NULL, equalize.cluster.weights = FALSE,
                                 tune.parameters = c("mtry", "min.node.size"),
                                 tune.num.trees = n_trees_tune)
  Y_hat_a = forest_a_y$predictions
  rownames(Y_hat_a) = rownames(Y)
  colnames(Y_hat_a) = paste0("Y_hat_a_",out_short)
  
  out_data = merge(out_data, Y_hat_a, by="row.names", all=TRUE)
  rownames(out_data) = out_data$Row.names
  out_data <- subset(out_data,select=-c(Row.names))
  # Save tuning parameteres.
  tun_p[(tun_p["target"]==outcome)&(tun_p["exp"]=="a"),  
        c("mtry", "min.node.size")] <- c(forest_a_y$tunable.params$mtry, 
                                         forest_a_y$tunable.params$min.node.size)
  
  # Use trained model to make predictions on the maximum extent data, as well. 
  y_hat_a_max_ext = predict(forest_a_y, subset(df_maxext, select=x_cols))
  rownames(y_hat_a_max_ext) = rownames(df_maxext)
  write.csv(y_hat_a_max_ext, 
            paste0(path_out,"max_ext/","max_ext_predictions_a_",out_short,".csv"),
            row.names=T)
  
  rm(forest_a_y, Y_hat_a, y_hat_a_max_ext)
  
  print(paste("Done with:",outcome, ", a"))

  #_________________________________________________________________________####
  #   b                                                                     ####
  # use clustering on the location level.
  forest_y = regression_forest(X, Y, num.trees = n_trees,
                               clusters = C, equalize.cluster.weights = FALSE,
                               tune.parameters = c("mtry", "min.node.size"),
                               tune.num.trees = n_trees_tune)
  Y_hat = forest_y$predictions
  rownames(Y_hat) = rownames(Y)
  colnames(Y_hat) = paste0("Y_hat_",out_short)
  
  out_data = merge(out_data, Y_hat, by="row.names", all=TRUE)
  rownames(out_data) = out_data$Row.names
  out_data <- subset(out_data,select=-c(Row.names))
  
  # Save out.
  write.csv(out_data, file=paste0(path_out,"oob_predictions_",out_short,".csv"),
            row.names=T)
  # Add tuning parameters.
  tun_p[(tun_p["target"]==outcome)&(tun_p["exp"]=="b"),  
        c("mtry", "min.node.size")] <- c(forest_y$tunable.params$mtry, 
                                         forest_y$tunable.params$min.node.size)
  
  # Use trained model to make predictions on the maximum extent data, as well. 
  y_hat_b_max_ext = predict(forest_y, subset(df_maxext, select=x_cols))
  rownames(y_hat_b_max_ext) = rownames(df_maxext)
  write.csv(y_hat_b_max_ext, 
            paste0(path_out,"max_ext/","max_ext_predictions_b_",out_short,".csv"),
            row.names=T)
  
  # No FFS for Y as it has proven to have no effect in the experiment. 
  # Save tuning parameteres.
  write.csv(tun_p, file=paste0(path_out,"tuning_params_",out_short,".csv"))
  
  rm(forest_y, Y_hat, Y, tun_p, out_data)
  
  print(paste("Done with:",outcome, ", b"))
}

df_f = read.csv(path_in)
rownames(df_f) = df_f$X
# Predict W.hat.
for (limit in c('maxspeed_100', 'maxspeed_120', 'maxspeed_130')){

  set.seed(845213)
  tun_p <- data.frame("target"=limit, "exp"=exps, "mtry"=NA, "min.node.size"=NA)

  # Set up dataframes.
  out_data = data.frame("id" = df_f$X)
  rownames(out_data) = df_f$X

  lim = substring(limit,nchar(limit)-2,nchar(limit))

  # Select only relevant stuff.
  df = df_f[(df_f['maxspeed_none']==1)|(df_f[limit]==1),]
  W = data.matrix(subset(df, select=c(limit)))
  X = data.matrix(subset(df, select=c(x_cols)))
  C = data.matrix(subset(df, select=c(cl)))

  #_________________________________________________________________________####
  #   a                                                                     ####
  # Make naive forest without clustering.
  forest_a_w = regression_forest(X, W, num.trees = n_trees,
                                 clusters = NULL, equalize.cluster.weights = FALSE,
                                 tune.parameters = c("mtry", "min.node.size"),
                                 tune.num.trees = n_trees_tune)
  W_hat_a = forest_a_w$predictions
  rownames(W_hat_a) = rownames(W)
  colnames(W_hat_a) = paste0("W_hat_a_",lim)
  out_data = merge(out_data, W_hat_a, by=0, all=TRUE)
  rownames(out_data) = out_data$Row.names
  out_data <- subset(out_data,select=-c(Row.names))

  tun_p[(tun_p["target"]==limit)&(tun_p["exp"]=="a"),
        c("mtry", "min.node.size")] <- c(forest_a_w$tunable.params$mtry,
                                         forest_a_w$tunable.params$min.node.size)

  # Use trained model to make predictions on the maximum extent data, as well.
  w_hat_a_max_ext = predict(forest_a_w, subset(df_maxext, select=x_cols))
  rownames(w_hat_a_max_ext) = rownames(df_maxext)
  write.csv(w_hat_a_max_ext,
            paste0(path_out,"max_ext/","max_ext_predictions_a_",lim,".csv"),
            row.names=T)

  rm(forest_a_w, W_hat_a, w_hat_a_max_ext)

  print(paste("Done with:", limit, "a"))

  #_________________________________________________________________________####
  #   b                                                                     ####
  # Use clustering on the location level.
  forest_w = regression_forest(X, W, num.trees = n_trees,
                               clusters = C, equalize.cluster.weights = FALSE,
                               tune.parameters = c("mtry", "min.node.size"),
                               tune.num.trees = n_trees_tune)
  W_hat = forest_w$predictions
  rownames(W_hat) = rownames(W)
  colnames(W_hat) = paste0("W_hat_",lim)
  out_data = merge(out_data, W_hat, by=0, all=TRUE)
  rownames(out_data) = out_data$Row.names
  out_data <- subset(out_data,select=-c(Row.names))

  # Save out results and remove from memory.
  write.csv(out_data, file=paste0(path_out,"oob_predictions_",lim,".csv"))
  tun_p[(tun_p["target"]==limit)&(tun_p["exp"]=="b"),
        c("mtry", "min.node.size")] <- c(forest_w$tunable.params$mtry,
                                         forest_w$tunable.params$min.node.size)

  # Use trained model to make predictions on the maximum extent data, as well.
  w_hat_b_max_ext = predict(forest_w, subset(df_maxext, select=x_cols))
  rownames(w_hat_b_max_ext) = rownames(df_maxext)
  write.csv(w_hat_b_max_ext,
            paste0(path_out,"max_ext/","max_ext_predictions_b_",lim,".csv"),
            row.names=T)

  # Save out tuning parameters and remove from memory
  write.csv(tun_p, file=paste0(path_out,"tuning_params_",lim,".csv"))

  rm(forest_w, W_hat, out_data, tun_p, w_hat_b_max_ext)

  print(paste("Done with:", limit, "b"))
}
# 

#___________________________________________________________________________####
#   p.hat (probability of being in analysis sample S                      ) ####
# For max_extent data, also predict probability of being in the analysis.
# Differs by experiment (cluster robust?) and speed limit (different sample)
for (exp in exps){
  for (limit in c('maxspeed_100', 'maxspeed_120', 'maxspeed_130')){
    for (s in c("","r_","ra_")){

      set.seed(845213)
      lim = substring(limit,nchar(limit)-2,nchar(limit))

      df = df_f[(df_f['maxspeed_none']==1)|(df_f[limit]==1),]

      # For the restricted sample, make sure that only observations really used
      # in the analysis are in RCT, rest goes to OBS.
      if (s=="r_"){
        # load propensity score.
        prop_score <- read.csv(file=paste0(path_out,"oob_predictions_",lim,".csv"),
                                row.names=1)
        if (exp=="b"){e="_"}else{e=paste0("_",exp,"_")}
        prop_score <- prop_score[paste0("W_hat",e,lim)]
        df = transform(merge(df, prop_score, by=0, all.x=F, all.y=F),
                       row.names=Row.names, Row.names=NULL)
        df = df[(df[paste0("W_hat", e, lim)]>0.05) &
                  (df[paste0("W_hat", e, lim)]<0.95) & (df["maxspeed_cond"]==0) &
                  (df["maxspeed_change"]==0) & (df["overtaking_ht"]==0),]
      }
      if (s=="ra_"){
        # load propensity score.
        prop_score <- read.csv(file=paste0(path_out,"oob_predictions_",lim,".csv"),
                               row.names=1)
        if (exp=="b"){e="_"}else{e=paste0("_",exp,"_")}
        prop_score <- prop_score[paste0("W_hat",e,lim)]
        df = transform(merge(df, prop_score, by=0, all.x=F, all.y=F),
                       row.names=Row.names, Row.names=NULL)
        df = df[(df[paste0("W_hat", e, lim)]<0.95),]
      }      

      # Select only relevant training data and test data.
      df_obs <- as.data.frame(subset(df_maxext, select=c(x_cols, 'location',
                                                         'maxspeed_none')))
      # Select only the currently unrestricted samples that are not in df.
      df_obs <- df_obs[!(row.names(df_obs) %in% rownames(df)) &
                         (df_obs$maxspeed_none==1),]

      # Select required things.
      X = data.matrix(subset(df, select=c(x_cols)))

      # Select X for test data, as well.
      X.test = subset(df_obs, select=colnames(X))
      # Cluster robust?
      if (exp %in% c('a', 'ab1')){
        C.test = NULL
        C = NULL
      }else{
        C.test = round(data.matrix(subset(df_obs, select=c(cl))))
        C = data.matrix(subset(df, select=c(cl)))
      }
      # Some cleaning
      rm(df_obs, df)

      n.rct = nrow(X)
      n.obs = nrow(X.test)

      # Train forest that predicts probability of being in the training sample.
      S.forest <- regression_forest(X = rbind(X, X.test),
                                    Y = c(rep(1, n.rct), rep(0, n.obs)),
                                    clusters = c(C, C.test),
                                    num.trees = n_trees,
                                    tune.parameters = c("mtry", "min.node.size"),
                                    tune.num.trees = n_trees_tune)
      rm(C, C.test)

      # Get oob predictions for p.hat.all
      p.hat.all <- predict(S.forest)$predictions
      p.hat.all <- as.data.frame(p.hat.all)
      rownames(p.hat.all) = rownames(rbind(X, X.test))

      # Subset only for analysis sample. This is whats needed later on.
      p.hat <- head(p.hat.all, n.rct)


      # Save all p.hat for possible needs to re-calculate things.
      write.csv(p.hat,
                file=paste0(path_out, s,'p_hat_ate_trans/',
                            "p_hat_",exp,'_',lim,".csv"),
                row.names=T)
      # Save p.hat for predictions of ATE_transport.
      write.csv(p.hat.all,
                file=paste0(path_out, s,'p_hat_ate_trans/',
                            "p_hat_all_",exp,'_',lim,".csv"),
                row.names=T)

      rm(S.forest, p.hat.all, p.hat)
      print(paste("Done with exp", exp, ", limit", limit, s, "."))
    }
  }
}








