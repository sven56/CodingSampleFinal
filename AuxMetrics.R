# ----------------------------------------------------------------------------------------------------
# MSPE
# ----------------------------------------------------------------------------------------------------
MSPE <- function(df) {
  y_hat = df %>% as.data.frame()
  y_real = stat_df[oos_start:oos_end, target_var] %>% as.data.frame()
  
  # NA Error prevention from RW
  y_comb = cbind(df, y_real) %>% 
    drop_na()
  
  MSPE = mean((y_comb[,1] - y_comb[,2])^2)
  
  return(MSPE)
}




# ----------------------------------------------------------------------------------------------------
# MAE
# ----------------------------------------------------------------------------------------------------
MAE <- function(df) {
  y_hat = df %>% as.data.frame()
  y_real = stat_df[oos_start:oos_end, target_var] %>% as.data.frame()
  
  # NA Error prevention from RW
  y_comb = cbind(df, y_real) %>% 
    drop_na()
  
  MAE = mean(abs(y_comb[,1] - y_comb[,2]))
  
  return(MAE)
}





# ----------------------------------------------------------------------------------------------------
# MAE
# ----------------------------------------------------------------------------------------------------
MAD <- function(df) {
  y_hat = df %>% as.data.frame()
  y_real = stat_df[oos_start:oos_end, target_var] %>% as.data.frame()
  
  # NA Error prevention from RW
  y_comb = cbind(df, y_real) %>% 
    drop_na()
  
  residual = y_comb[,2] - y_comb[,1]
  
  MAD = median(abs(residual - median(residual)))
  
  return(MAD)
}





# ----------------------------------------------------------------------------------------------------
# VaR Loss function
# ----------------------------------------------------------------------------------------------------
RLF = function(vector, boundary) {
  sum = 0
  for (i in 1:length(vector)) {
    if (vector[i] < boundary[i]) {
      sum = sum + (vector[i]-boundary[i])^2
    }
  }
  
  return(sum)
}





# ----------------------------------------------------------------------------------------------------
# MSPE - Vintage
# ----------------------------------------------------------------------------------------------------
MSPE_vintage <- function(df, real) {
  y_hat = df %>% as.data.frame()
  y_real = real %>% as.data.frame()
  
  # NA Error prevention from RW
  y_comb = cbind(df, y_real) %>% 
    drop_na()
  
  MSPE = mean((y_comb[,1] - y_comb[,2])^2)
  
  return(MSPE)
}





# ----------------------------------------------------------------------------------------------------
# MAE - Vintage
# ----------------------------------------------------------------------------------------------------
MAE_vintage <- function(df, real) {
  y_hat = df %>% as.data.frame()
  y_real = real %>% as.data.frame()
  
  # NA Error prevention from RW
  y_comb = cbind(df, y_real) %>% 
    drop_na()
  
  MAE = mean(abs(y_comb[,1] - y_comb[,2]))
  
  return(MAE)
}





# ----------------------------------------------------------------------------------------------------
# MAE - Vintage
# ----------------------------------------------------------------------------------------------------
MAD_vintage <- function(df, real) {
  y_hat = df %>% as.data.frame()
  y_real = real %>% as.data.frame()
  
  # NA Error prevention from RW
  y_comb = cbind(df, y_real) %>% 
    drop_na()
  
  residual = y_comb[,2] - y_comb[,1]
  
  MAD = median(abs(residual - median(residual)))
  
  return(MAD)
}

