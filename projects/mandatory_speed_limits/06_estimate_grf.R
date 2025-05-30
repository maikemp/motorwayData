library(grf)
library(rstudioapi, deplyr)
# library(logr)
library(rjson)
library("rstudioapi")
this_file = getSourceEditorContext()$path 
source(paste0(this_file,"/../functions/cf_est_functions.R"))

# Cleans up space.
gcinfo(FALSE)



path_out = ".../data_out/analysis/grf_results/"
path_in = ".../out/data/analysis_df.csv"
#path_out = "C:/Users/Maike/sciebo/speedLimit/out/regression_results/"

path_max_ext = '.../out/data/df_max_ext_imputed.csv'
path_ate_transp=".../data_out/analysis/oob_predictions/p_hat_ate_trans/"

# log_open(paste0(path_out,"res_log.log"))
# Load data and fix index.
df_f <- read.csv(path_in, row.names=1)

# Load cutoff values for CATEs.
cates_cutoffs <- fromJSON(file=paste0(this_file, 
                                      "/../../../Libraries/cates_cutoffs.json"))

# 
y_cols = colnames(df_f)[ifelse(endsWith(colnames(df_f), '_rate'), T, F)]
y_hat_cols = colnames(df_f)[ifelse(startsWith(colnames(df_f), 'Y_hat'), T, F)]
w_hat_cols = colnames(df_f)[ifelse(startsWith(colnames(df_f), 'W_hat'), T, F)]
sl_cols = colnames(df_f)[ifelse(startsWith(colnames(df_f), 'maxspeed'), T, F)]
w_cols = c("maxspeed_100", "maxspeed_120", "maxspeed_130" )
cl = 'location'
x_cols = colnames(df_f)[ifelse(!colnames(df_f) %in% c(y_cols, sl_cols, cl,  
                              y_hat_cols, w_hat_cols, 'overtaking_ht'), T, F)]

print(x_cols)

n_trees = 15000 # Set to 15000
n_trees_tune = 50 # Set to 50

df_maxext = read.csv(path_max_ext, row.names=1)
df_maxext = df_maxext[,c(x_cols,"maxspeed_cond", "maxspeed_change", "location", 
                         "maxspeed_none")]



for (exp in c("b")){ #, 
  print(paste("! ! ! -  - - - - - - - - - - - - - - >  - - - - - - - - - > - - - - - - > Start with:", exp))
  for (limit in c("maxspeed_130","maxspeed_120", "maxspeed_100")){ 
    
    print(paste("! ! ! -  - - - - - - - - - - - - - - - - >  - - - - - - > Start with:", limit))
    
    
    for (outcome in c('light_rate', 'severe_rate', 'fatal_rate', 'total_rate')){
      # Set seed here for better reproducability.
      set.seed(845213)
      
      te_res = c()
      
      print(paste("! ! ! -  - - - - - - - - - - - - - - - - - > Start with:", outcome))
      
      e = exp_str(exp)
      out_short = short_outcome(outcome)
      lim = substring(limit,nchar(limit)-2,nchar(limit))
      
      df = df_f[(df_f['maxspeed_none']==1)|(df_f[limit]==1),]
      
      # # # #   # # #   # #   #       #  # #  #        #   # #   # # #   # # # # 
      # Activate following 5 lines to restrict according to prop.score 
      # restriction, conditional speed limit & no change in sl
      # path_out = "D:/speedLimit/data_out/analysis/grf_results/r_"
      # df = df[(df[paste0("W_hat_", e, lim)]>0.05) &
      #           (df[paste0("W_hat_", e, lim)]<0.95) & (df["maxspeed_cond"]==0) &
      #            (df["maxspeed_change"]==0) & (df["overtaking_ht"]==0),]      
      # path_ate_transp="D:/speedLimit/data_out/analysis/oob_predictions/r_p_hat_ate_trans/"
      # # # #   # # #   # #   #       #  # #  #        #   # #   # # #   # # # # 
      
      # # # #   # # #   # #   #       #  # #  #        #   # #   # # #   # # # # 
      # Activate following 3 lines to restrict by excluding only most dangerous segments
      # restriction, conditional speed limit & no change in sl
      # path_out = "D:/speedLimit/data_out/analysis/grf_results/ra_"
      # df = df[(df[paste0("W_hat_", e, lim)]<0.95),]      
      # path_ate_transp="D:/speedLimit/data_out/analysis/oob_predictions/ra_p_hat_ate_trans/"
      # # # #   # # #   # #   #       #  # #  #        #   # #   # # #   # # # #       
      
      print(paste("Dimension df:",dim(df)))  

      
      Y = data.matrix(subset(df, select=c(outcome)))
      W = data.matrix(subset(df, select=c(limit)))
      if (exp == 'c'){
        path_in_ffs = paste0("D:/speedLimit/data_out/analysis/oob_predictions/ffs_vars_fin_",lim,".csv")
        ffs_vars = read.csv(path_in_ffs, row.names=1)
        ffs_var_cols = ffs_vars[ffs_vars[,paste0("W_hat_",lim)] == 1,]$var
        # Select only variables that were relevant in FFS. 
        X = data.matrix(subset(df, select=ffs_var_cols))
        rm(path_in_ffs, ffs_vars, ffs_var_cols)
        ey = ''
      }else if (exp %in% c('ab1','ab2')){
        X = data.matrix(subset(df, select=x_cols))
        ey='a_'
      }else{ # includes a, b, ab, ba, that use the setup defined in exp_str(). (a, ab: non-robust nuisance functions; b, ba: robust nuisance functions.)
        X = data.matrix(subset(df, select=x_cols))
        ey = e
      }
      
      if (exp %in% c('a', 'ab1', 'ba')){
        C = NULL
      }else{ # includes 'b', 'ab', 'ab2'
        C = data.matrix(subset(df, select=c(cl)))
      }
      print(paste("cluster on", colnames(C)))
      
      Y.hat = data.matrix(subset(df, select=c(paste0("Y_hat_", ey, out_short))))
      W.hat = data.matrix(subset(df, select=c(paste0("W_hat_", e, lim)))) 
      print(paste("Y_hat: ",paste0("Y_hat_", ey, out_short)))
      print(paste("W_hat: ",paste0("W_hat_", e, lim)))
      
      rm(df)
      
      # print(paste("Mean of Y:", round(mean(Y),3)))
      # print(paste("Share of treated:", round(sum(W)/length(W),3)))
      
      #_____________________________________________________________________####
      # Get causal forest.                                                  ####
      
      cf <- causal_forest(X, Y, W, Y.hat=Y.hat, W.hat = W.hat, 
                          num.trees=n_trees, clusters=C, 
                          tune.parameters=c("mtry", "min.node.size"),
                          tune.num.trees=n_trees_tune)
      rm(C)
      
      #_____________________________________________________________________####
      #   Get tuning parameters and save out.                               ####
      
      tun_p <- data.frame("target"=outcome, "limit"=limit, "exp"=exp, 
                          "mtry"=cf$tunable.params$mtry, 
                          "min.node.size"=cf$tunable.params$min.node.size)
      write.csv(tun_p, 
                paste0(path_out, 
                       "tuning_params_", exp,'_', lim, "_",out_short,".csv"), 
                row.names = F)
      
      #_____________________________________________________________________####
      #   Get variable importance and save them out.                        ####   
      
      var_imp = variable_importance(cf)
      var_imp_df = as.data.frame(cbind(var=colnames(X), imp=var_imp))
      var_imp_df = var_imp_df[order(var_imp_df$V2,decreasing=T),]
      write.csv(var_imp_df, paste0(path_out, "var_imp_", exp,'_', lim, "_",out_short,".csv"), 
                row.names = F)
      rm(var_imp, var_imp_df)
      

      #_____________________________________________________________________####
      #   Get individual treatment effects plus std.err. and save it out.   ####
      
      # cf_pred <- predict(cf, X, estimate.variance = T)
      # tau_hat <- cf_pred$predictions
      # sigma_hat <- sqrt(cf_pred$variance.estimates)
      # indiv_res = cbind(tau_hat, sigma_hat)
      # rownames(indiv_res) = rownames(Y)
      # write.csv(indiv_res, paste0(path_out,"ind_res_",exp,'_',lim,"_",out_short,".csv"), 
      #           row.names = T)
      
      cf_pred <- predict(cf, estimate.variance=T)
      tau_hat_oob <- cf_pred$predictions
      sigma_hat <- sqrt(cf_pred$variance.estimates)
      debiased_err <- cf_pred$debiased.error
      ex_err <- cf_pred$excess.error
      # E[Y | X, W = 0]
      mu_hat_0 <- Y.hat - W.hat * tau_hat_oob
      # E[Y | X, W = 1]
      mu_hat_1 <- Y.hat + (1 - W.hat) * tau_hat_oob
      colnames(mu_hat_0) <- c("mu_hat_0")
      colnames(mu_hat_1) <- c("mu_hat_1")
      indiv_res = cbind(tau_hat_oob, sigma_hat, debiased_err, ex_err, 
                        mu_hat_0, mu_hat_1)
      rownames(indiv_res) = rownames(Y)
      write.csv(indiv_res, paste0(path_out,"ind_res_",exp,'_',lim,"_",out_short,".csv"),
                row.names = T)
      
      #hist(tau_hat_oob, main=paste("Treatment effects for:", limit), breaks=75)
      rm(cf_pred, sigma_hat, debiased_err, ex_err, indiv_res, Y.hat, 
         W.hat, mu_hat_0, mu_hat_1)

      #_____________________________________________________________________####
      #   Predict treatment effect on the full extent data and save it out. ####    
      
      max_ext_pred <- predict(cf, 
                              data.matrix(subset(df_maxext, select=colnames(X))), 
                              estimate.variance = T)
      tau_hat_me <- cbind(max_ext_pred$predictions, 
                          sqrt(max_ext_pred$variance.estimates), 
                          df_maxext$maxspeed_cond,
                          df_maxext$maxspeed_change,
                          df_maxext$AADT)
      colnames(tau_hat_me) = c("tau_hat_i", "simga_hat_i", "maxspeed_cond", 
                               "maxspeed_change", "AADT")
      rownames(tau_hat_me) = rownames(df_maxext)
      write.csv(tau_hat_me, paste0(path_out,"ind_res_max_ext_",exp,'_',lim,"_",out_short,".csv"), 
                row.names = T)
      #hist(tau_hat_me[,"tau_hat_i"], main=paste("Treatment effects for:", limit), breaks=75)
      rm(tau_hat_me)
      
      #_____________________________________________________________________####
      #   Estimate ATE.                                                     ####      

      ate = average_treatment_effect(cf, "all", method="AIPW")
      ATE = ate[1]
      StdErr = ate[2]
      
      print(paste("ATE:", round(ATE,3)))
      print(paste("std.err:", round(StdErr,3)))
      print(paste("t:", round(ATE/StdErr,3)))
      
      # Add Treatment effect on the treated and treatment effect on the untreated
      atet = average_treatment_effect(cf, "treated", method="AIPW")
      ateut = average_treatment_effect(cf, "control", method="AIPW")
      # Also derive for overlap.
      ate_ov = average_treatment_effect(cf, "overlap", method="AIPW")

      te_res = cbind(te_res, 
                     "ate"=c(ATE, StdErr, "W_share"=mean(W), "Y_mean"=mean(Y)),
                     "atet"=c(atet[1], atet[2], "W_share"=1, "Y_mean"=mean(Y[W==1])),
                     "ateut"=c(ateut[1], ateut[2], "W_share"=0, "Y_mean"=mean(Y[W==0])),
                     "ate_ov"=c(ate_ov[1], ate_ov[2], "W_share"=NULL, "Y_mean"=NULL))
      
      #_____________________________________________________________________####
      #   Run omnibus test and save out results.                            ####
      
      ot = test_calibration(cf)
      print(ot)
      write.csv(ot, paste0(path_out,"omn_test_",exp,'_',lim,"_",out_short,".csv"))
      rm(ate, ATE, StdErr, ot)
      

      
      if (exp == 'c'){
        X = data.matrix(subset(X, select=x_cols))
      }
      
      #_____________________________________________________________________####
      #   Get CATEs for high and low values of above defined values.        ####      
        
      for (k in names(cates_cutoffs)){
        cate_h = average_treatment_effect(cf, "all", method="AIPW", subset=X[,k]> cates_cutoffs[k])
        cate_l = average_treatment_effect(cf, "all", method="AIPW", subset=X[,k]<=cates_cutoffs[k])
        cn = c(colnames(te_res), paste0("cate_h_", k), paste0("cate_l_", k))
        te_res = cbind(te_res, c(cate_h[1], cate_h[2], "W_share"=mean(W[X[,k] > cates_cutoffs[k]]), 
                                 "Y_mean"=mean(Y[X[,k] > cates_cutoffs[k]])), 
                               c(cate_l[1], cate_l[2], "W_share"=mean(W[X[,k]<= cates_cutoffs[k]]), 
                                 "Y_mean"=mean(Y[X[,k]<= cates_cutoffs[k]])))
        colnames(te_res) = cn
      }
      
      write.csv(te_res, paste0(path_out, "cates_",exp,"_",lim,"_",out_short,".csv"), 
                row.names = T)
      rm(cate_h, cate_l, cn, te_res)
      
      #_____________________________________________________________________####
      #   Derive ATE_transport for unrestricted roads not in analysis data  ####  
      # See: https://grf-labs.github.io/grf/articles/ate_transport.html
      test_df <- as.data.frame(subset(df_maxext, select=c(x_cols, 
                                                          'maxspeed_none')))
      test_df <- test_df[!(row.names(test_df) %in% rownames(X)) & 
                         (test_df$maxspeed_none>0),]

      X.test = subset(test_df, select=colnames(X))
      rm(test_df, X)
      
      if (exp == 'ba'){exp_2 = "b"}else if (exp == "ab"){exp_2 = "a"}else{exp_2 = exp}
      
      p.hat = read.csv(paste0(path_ate_transp, "p_hat_",exp_2,"_", lim,'.csv'), 
                       row.names=1)
      
      ate_transport <- average_treatment_effect_test(forest = cf, 
                                                     tau.hat.train= tau_hat_oob,
                                                     X.test = X.test,
                                                     p.hat = p.hat,
                                                     Y.orig = Y, 
                                                     W.orig = W)
      
      rm(X.test, cf, p.hat, Y, W)
      # Write out. 
      write.csv(c("ATE_test"=ate_transport$estimate, 
                  "std_err"=ate_transport$std_err), 
                paste0(path_out, "ate_transport/", 
                       "ate_t_",exp,"_",lim,"_",out_short,".csv"))
      # Save individual results for easier re-calculation later on.
      ate_transport_i_res = ate_transport$tau_hat_all
      colnames(ate_transport_i_res)=c("tau_hat_all")
      write.csv(ate_transport_i_res, 
                paste0(path_out, "ate_transport/", 
                       "ate_t_i_res_",exp,"_",lim,"_",out_short,".csv"), 
                row.names = T)
      rm(ate_transport, ate_transport_i_res)

    }
  }
}


# TODO: Repeat everything (or just b) for restricted data 





