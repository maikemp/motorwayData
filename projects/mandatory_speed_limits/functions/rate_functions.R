


get_colnames <- function(df_f){
   y_cols = colnames(df_f)[ifelse(endsWith(colnames(df_f), '_rate'), 
                                  TRUE, FALSE)]
   y_hat_cols = colnames(df_f)[ifelse(startsWith(colnames(df_f), 'Y_hat'), 
                                      TRUE, FALSE)]
   w_hat_cols = colnames(df_f)[ifelse(startsWith(colnames(df_f), 'W_hat'), 
                                      TRUE, FALSE)]
   sl_cols = colnames(df_f)[ifelse(startsWith(colnames(df_f), 'maxspeed'), 
                                   TRUE, FALSE)]
   w_cols = c("maxspeed_100", "maxspeed_120", "maxspeed_130" )
   cl = 'location'
   x_cols = colnames(df_f)[ifelse(!colnames(df_f) %in% c(y_cols, sl_cols, cl, 
                                                         y_hat_cols, w_hat_cols, 
                                                         'overtaking_ht'),
                                  TRUE, FALSE)]
  
   # print(x_cols)
   
   return(list("y_cols"=y_cols, "y_hat_cols"=y_hat_cols, "w_hat_cols"=w_hat_cols,
               "sl_cols"=sl_cols, "w_cols"=w_cols, "cluster"=cl, "x_cols"=x_cols))
}


