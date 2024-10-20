# # ----------------------------------------------------------------------------------------------------
# # BENCHMARK: Realized Variance (intra-month)
# # ----------------------------------------------------------------------------------------------------
# library(RSQLite)
# library(dbplyr)
# library(RPostgres)
# wrds <- dbConnect(Postgres(),
#                   host='wrds-pgdata.wharton.upenn.edu',
#                   port=9737,
#                   sslmode='require',
#                   dbname='wrds',
#                   user='svenvanholten')
# ex_r <- dbFetch(dbSendQuery(wrds, "select datadate,exratd from comp_global_daily.g_exrt_dly WHERE fromcurd = 'GBP' AND tocurd = 'USD'"))
# ex_r[,2] = c(NA,diff(ex_r[,2])) # First differencing
# 
# intra_returns <- function(vector) { avg = mean(vector); return(sqrt(sum((vector - avg)^2))) }
# realized_vol = ex_r %>%
#   group_by(month = floor_date(datadate, "month")) %>%
#   summarize(RV_estimate = intra_returns(exratd)) %>%
#   filter(month >= "2014-01-01" & month <= "2022-01-01") %>%
#   as.data.frame()
# 
# ggplot(realized_vol, aes(x = month, y = RV_estimate)) +
#   geom_line(color = "#1b9e77") +
#   labs(x = "Date", y = "Value", title = "Time Series Plot") +
#   theme_classic() + theme(axis.title.y = element_blank(), axis.title.x=element_blank(),plot.title = element_blank())
# 
# # Load first stuff from Main regarding GARCH-family forecasting
# volatility_horizon_m = cbind(arch_mat[,h_1],GARCH_forecast,ARCH_forecast,ARMA_GARCH_forecast[,1], realized_vol[,2])
# cheap_benchmark_plotting(volatility_horizon_m, c("date", "MRF", "GARCH", "ARCH", "ARMA-GARCH", "RV"))
# 




# ----------------------------------------------------------------------------------------------------
# BENCHMARK: RW
# ----------------------------------------------------------------------------------------------------
rw_function <- function(df,horizon) { #df should be real values: stat_df$target_var
  new_vec = c(rep(NA,horizon),head(df,-horizon))
  return(new_vec)
}





# ----------------------------------------------------------------------------------------------------
# BENCHMARK: AR(p)
# ----------------------------------------------------------------------------------------------------
arP_function <- function(df, horizon, p) {
  output_vector = matrix(NA, oos_end-oos_start+1)
  for (oos_i in oos_start:oos_end) {
    data_expanding = head(df, oos_i)
    
    y = data_expanding %>%
      head(-horizon)
    
    ar_setup = arima(y, order=c(p,0,0), include.mean = TRUE)
    ar_p = predict(ar_setup, n.ahead = horizon)
    prediction = ar_p$pred[[horizon]] #First-differenced
    
    output_vector[(oos_i-oos_start+1)] = prediction
  }
  return(output_vector)
}





# ----------------------------------------------------------------------------------------------------
# VaR BENCHMARK: GARCH/ARCH/ARMA-GARCH, returns boundary vector
# ----------------------------------------------------------------------------------------------------
volatility_benchmarks_vector <- function(GARCH, ARMA, horizon) {
  temp_mat =  data.frame(matrix(NA,oos_end-oos_start+1, 8), row.names = stat_df[oos_start:oos_end,1])
  colnames(temp_mat) = c("sigma","mu","omega","alpha","t_df","skew_t_df","skew_t_xi", "pred")
  
  last_successful_iteration <- NULL # Failsave
  for (oos_i in oos_start:oos_end) {
    tryCatch({
      data_expanding <- head(stat_df[,target_var], oos_i) # Expanding
      
      model.fit.GARCH <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = GARCH),
                                    mean.model= list(armaOrder = ARMA, include.mean = TRUE),
                                    distribution.model= "norm")
      fit.GARCH <- ugarchfit(spec=model.fit.GARCH, data=data_expanding, out.sample = horizon, solver ='hybrid') 
      GARCH_forecast <- ugarchforecast(fit.GARCH, n.ahead = horizon, out.sample = horizon)
      
      temp_mat[oos_i - oos_start +1,8] <- as.numeric(GARCH_forecast@forecast$seriesFor[horizon])
      
      temp_mat[oos_i - oos_start + 1,1] <- as.numeric(GARCH_forecast@forecast$sigmaFor[horizon])
      temp_mat[oos_i - oos_start + 1,2] <-  as.numeric(GARCH_forecast@model$pars[1]) # Mu
      temp_mat[oos_i - oos_start + 1,3] <-  as.numeric(GARCH_forecast@model$pars["omega",1]) #omega
      temp_mat[oos_i - oos_start + 1,4] <- as.numeric(GARCH_forecast@model$pars["alpha1",1]) # alpha
      
      # Vector for epsilons
      resid = GARCH_forecast@model$modeldata$residuals
      temp_mat[oos_i - oos_start + 1,5] = stdFit(resid)[[1]][3] # fit t, df
      suppressWarnings(temp_mat[oos_i - oos_start + 1,6] <- sstdFit(resid)[[2]][3]) # fit skewed t, df
      suppressWarnings(temp_mat[oos_i - oos_start + 1,7] <- sstdFit(resid)[[2]][4]) # fit skewed t, xi
      
      last_successful_iteration <- temp_mat[oos_i - oos_start + 1,]  # Update last successful iteration values
    }, error = function(e) {
      temp_mat[oos_i - oos_start + 1,] = last_successful_iteration
    })
  }
  temp_mat = na.locf(temp_mat) #Fwd fill incase of error
  
 return(temp_mat)
}
