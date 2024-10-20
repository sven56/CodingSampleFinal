#
# Code from "MRF-ARCH: A MACHINE LEARNING APPROACH TO FORECAST USD/GBP TAIL RISK IN A HIGH DIMENSIONAL DATASET"
# The code is subdivided by functionality into various R files, e.g. Plotting, Benchmarks, Metrics, and Aux. Functions
# There is one non-connected R file called VintageAnalysis in the folder, which performs the robustness check on the real-time data
# All code is hand-crafted, except for an imported fredMD.R function, as for some reason the package doesn't always load for all versions of R  
# You can comment / un-comment the different "modules" in the MRF function depending on the required purpose
# ~ Sven van Holten Charria
#





# ----------------------------------------------------------------------------------------------------
# Reproducibility of benchmarks and plots - Load in the pre-produced data and skip the 'MRF'-section, go to 'Benchmarks' or 'Plotting'
# ----------------------------------------------------------------------------------------------------
output_mat = readRDS("output_mat_example.rds")

# Reproducibility and basic settings
set.seed(2001)
setwd("/Users/svenvanholten/Documents/Econometrics Erasmus/Year 5/Thesis/Forecasting Seminar/Code Seminar") 





# ----------------------------------------------------------------------------------------------------
# Packages and Sourcing
# ----------------------------------------------------------------------------------------------------
# Download necessary packages - Useful if running on parallel computer setup
if (!requireNamespace("lars", quietly = TRUE)) {
  install.packages("devtools")
  library(devtools)
  install_github("philgoucou/macrorf"); install.packages("lars"); install.packages("vars"); install.packages("tseriers"); install.packages("dplyr");   install.packages("tidyverse"); install.packages("readr")
  install.packages("reshape2")
  }

# Import general libraries
packages <- c("zoo", "tseries", "lars", "MacroRF", "magrittr", "dplyr", 
              "tidyverse", "readr", "vars", "tidyr", "reshape2", 
              "rugarch", "GAS", "forecast","fGarch")

# Function to check if a package is installed, and install it if not
install_if_missing <- function(pkg) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}

# Apply the function to each package in the list
invisible(lapply(packages, install_if_missing))


# Source the auxilliary functions files
source("AuxFunctions.R")
source("AuxMetrics.R")
source("AuxPlotting.R")
source("AuxBenchmarks.R")
source("fredMD.R")





# ----------------------------------------------------------------------------------------------------
# Definitions 
# ----------------------------------------------------------------------------------------------------
# Forecasting horizon 
horizon = 1

### Hyperparameter tuning definitions - Current values are found to be 'optimal'
# Max number of LARS variables
lars_n = 28

# Lags in the Vector Autoregression specification (for ease: just the number of lags of DV, PCAs and MAFs)
lag_k = 8 

# PCA and MAF dimensions
factor_n = 5
MAF_n = 5

# Ridge Lambda and Random Walk Regularization
RL = .5
RWR = .95





# ----------------------------------------------------------------------------------------------------
# Data importing
# ----------------------------------------------------------------------------------------------------
# Obtain MAF DF
MAF_df = MAF_data() #60 months prior / 5Y

# Import monthly FRED data; doResource = TRUE imports the target variable without the standard transformations. 
# After, we only perform a single first-diff as per our methodology.
start_date_import = "1998-01-02"; end_date_import = "2023-05-01"
target_var =  "EXUSUKx"
input_df = getFRED(start_date_import, end_date_import, doResource = TRUE, target_var, wantUK = TRUE)
input_df$EXUSUKx = c(NA, diff(input_df$EXUSUKx))

# Import non-transformed target variable, s.t. we can create 
fred_v_notrans = getFredRaw(start_date_import, end_date_import, target_var)





# ----------------------------------------------------------------------------------------------------
# Data Preperation
# ----------------------------------------------------------------------------------------------------
# Interpolate NAs in middle, and make series that are non-stationary -> stationary. Throw error if a series still stays non-stationary.
stat_df = interpolate_stationary(input_df)

# Out-of-sample definitions
oos_start <- head(which(stat_df$date >= "2014-01-01"),1) 
oos_end <- tail(which(stat_df$date <= "2022-01-01"),1)

# The 'real' vector of differenced data (see "Data Importing") and its 1L. Used as evaluation tool (VaR and metrics)
real_v = stat_df[oos_start:oos_end,target_var]
real_v_shifted = stat_df[(oos_start-1):(oos_end-1),target_var]




# ----------------------------------------------------------------------------------------------------
# MRF - The 'magic'
# ----------------------------------------------------------------------------------------------------
### OUTPUT MATRICES
# Output matrix definition for MRF point forecasting
output_mat = data.frame(matrix(NA,oos_end-oos_start+1, 10), row.names = stat_df[oos_start:oos_end,1])
colnames(output_mat) = c("no_diff","diff","pred_vol_h1", "pred_vol_h2", "pred_vol_h3", "omega", "alpha", "t_df", "sk_t_df", "sk_t_xi")

# Output matrix definitions for volatility forecasting under different horizons
volatility_mat_h1 =  data.frame(matrix(NA,oos_end-oos_start+1, 12), row.names = stat_df[oos_start:oos_end,1]); colnames(volatility_mat_h1) = c(1:12)
volatility_mat_h2 =  data.frame(matrix(NA,oos_end-oos_start+1, 12), row.names = stat_df[oos_start:oos_end,1]); colnames(volatility_mat_h2) = c(1:12)
volatility_mat_h3 =  data.frame(matrix(NA,oos_end-oos_start+1, 12), row.names = stat_df[oos_start:oos_end,1]); colnames(volatility_mat_h3) = c(1:12)



### THE MRF FUNCTION EXPANDING WINDOW
for (oos_i in oos_start:oos_end) { 
  ### STEP 0: Preperation
  data_expanding = head(stat_df, oos_i)
  
  
  
  ### STEP 1: LARS - Selects 'lars_n' variable to include in STEP 2
  var_s = LARS(data_expanding, lars_n)
  
  # Dividing DV and IV sets
  dv_set = data_expanding %>%
    dplyr::select(date, target_var)
  iv_set = data_expanding %>%
    dplyr::select(date, all_of(var_s))
  
  
  
  ### STEP 2.1: PCA
  # Performing PCA
  factors = PCA_function(iv_set)
  factors = factors[,1:factor_n]
  
  ### STEP 2.2: MAF - One year backward looking
  MAF_start = head(data_expanding[,"date"],1) - years(1) 
  MAF_end = tail(data_expanding[,"date"],1)
  MAF_factors = MAF_function(MAF_df, MAF_n, MAF_start, MAF_end)
  
  ### STEP 2.3: Combining the PCAs and MAFs, then lagging
  # Exclude MAF:
  # Y = cbind(dv_set, factors)
  
  # Include MAF:
  Y = cbind(dv_set, factors, MAF_factors)
  
  # Using VAR to construct matrix
  Y_temp = Y[c(1:nrow(Y), rep(nrow(Y),lag_k)), -1]
  step2_out <- vars::VAR(Y_temp, p = lag_k, type = "trend")[["datamat"]] %>%
    as.data.frame() %>%
    dplyr::select(contains(paste(target_var)), contains("PC"),contains(".l1"),contains(".l2"), trend)
  
  
  
  ### STEP 3.1: MRF
  ### Creating train and test matrices
  # X - Full covariate set
  x_all = step2_out %>% 
    dplyr::select(-target_var) 
  
  # Y - Select full Y-vector until t-h, defining the "out-of-sample" forecast
  y_train = step2_out %>% 
    dplyr::select(target_var) %>%
    head(-horizon) %>%
    as.matrix()
  y_test = matrix(NA, nrow = horizon)
  y_all = rbind(y_train, y_test)
  
  # Combining X and Y - This matrix is the 'input' matrix of the MRF
  MRF_in = cbind(y_all,x_all)
  
  # Indices of elements to include in MRF's general model specification: dv ~ dv_L1 + dv_L2 + PCA1_L1 + PCA2_L1
  index_dv1 = which(colnames(MRF_in)==paste(target_var,".l1",sep=""))
  index_dv2 = which(colnames(MRF_in)==paste(target_var,".l2",sep=""))
  index_f1 = which(colnames(MRF_in)=="PC1.l1")
  index_f2 = which(colnames(MRF_in)=="PC2.l1")
  x_position = c(index_dv1, index_dv2, index_f1, index_f2)
  
  # Running the MRF - Incl. hyperparameter values
  mrf_model = MRF(MRF_in, 
                  x.pos = x_position,
                  oos.pos=((nrow(MRF_in)-horizon+1):nrow(MRF_in)), 
                  mtry.frac = 1/3,
                  ridge.lambda = RL,
                  rw.regul = RWR,
                  trend.push = ncol(MRF_in),
                  B = 500,
                  quantile.rate = 0.3,
                  minsize = 15,
                  fast.rw = TRUE) 
  
  
  
  ### STEP 3.2: MRF Result extraction / Forecasting
  # Retrieve last non-transformed value in in-sample (t-h) to do the out-of-sample prediction (of non-diff)
  no_diff_init_index = data_expanding %>% 
    dplyr::select(date) %>% 
    head(oos_i-horizon) %>%
    tail(1) %>%
    as.double()
  notrans_value = fred_v_notrans[fred_v_notrans$date == no_diff_init_index, target_var]
  
  # Extracting predictions of MRF (of size h)
  preds = mrf_model[["pred"]]
  
  # Recording forecasts: Non-transformed & First-differenced
  output_mat[oos_i - oos_start + 1,1] = notrans_value + sum(preds) # NT
  output_mat[oos_i - oos_start + 1,2] = sum(preds) # FD
  
  
  
  ### STEP 4.1: ARCH PARAMETER EXTRACTION
  # Estimating the epsilon matrix. 
  AR_X <- mrf_model$YandX[,-1]
  beta <- mrf_model$betas
  epsilon_m <- matrix(NA, nrow(AR_X))
  
  # For in-sample epsilons, use known betas and values of X. 
  # For out-of-sample epsilons, use estimated betas and last in-sample (known) value of X
  for (i in 1:(nrow(AR_X)-horizon)) { # In-sample
    epsilon_m[i] = mrf_model$YandX[i,1] - c(1,AR_X[i,]) %*% beta[i,]
  } 
  for (i in 1:horizon) { #Out-of-sample
    epsilon_m[nrow(AR_X)-horizon+i] = sum(preds[1:i]) - c(1,AR_X[nrow(AR_X)-horizon,]) %*% beta[nrow(AR_X)-horizon+i,]
  }
  
  # Epsilon-matrix transformations to prepare for AR(p) -> ARCH(p)
  epsilon_sq_m = epsilon_m^2 # Epsilon-squared matrix
  mu_sq = mean(epsilon_sq_m)
  epsilon_sq_hat = epsilon_sq_m - mu_sq # Standardized epsilon-squared matrix
  
  
  
  ### MODULE A: MRF-ARCH(1) to extract time-varying omega / alpha
  ar_1 = arima(epsilon_sq_hat, c(1,0,0), include.mean=FALSE)
  alpha = ar_1$coef[1]
  omega = mu_sq * (1 - alpha)
  output_mat[oos_i - oos_start + 1,6] = omega
  output_mat[oos_i - oos_start + 1,7] = alpha
  
  
  
  ### MODULE B: MRF-ARCH(p) volatility forecasting (given horizon). See Methodology for the forecasting formulas
  ### NOTE: We defined 'p' as 'g' in the following formulas
  # Horizon = 1
  # for (g in 1:12) {
  #   ar_g = arima(epsilon_sq_hat, c(g,0,0), include.mean=FALSE, method = "ML")
  #   sum = 0
  #   for (i in 1:g) {
  #     sum = sum + ar_g$coef[i] * (tail(epsilon_sq_hat,i+1) %>% head(1))
  #   }
  #   omega = mu_sq * (1 - sum(ar_g$coef))
  #   h1_vol = omega + sum
  #   volatility_mat_h1[oos_i - oos_start + 1,g] = sqrt(abs(h1_vol))
  # }
  # 
  # # Horizon = 2
  # if (oos_i >= oos_start + 1) {
  #   for (g in 1:12) {
  #     ar_g = arima(epsilon_sq_hat, c(g,0,0), include.mean=FALSE, method = "ML")
  #     sum = 0
  #     if (g >= 2){
  #       for (i in 2:g) {
  #         sum = sum + ar_g$coef[i] * (tail(epsilon_sq_hat,i+1) %>% head(1))
  #       }
  #     }
  #     omega = mu_sq * (1 - sum(ar_g$coef))
  #     h2_vol = omega + ar_g$coef[1]*volatility_mat_h1[oos_i-oos_start,g]^2 + sum
  #     volatility_mat_h2[oos_i - oos_start + 1,g] = sqrt(abs(h2_vol))
  #   }
  # } 
  # 
  # # Horizon = 3
  # if (oos_i >= oos_start + 2) {
  #   for (g in 1:12) {
  #     ar_g = arima(epsilon_sq_hat, c(g,0,0), include.mean=FALSE, method = "ML")
  #     omega = mu_sq * (1 - sum(ar_g$coef))
  #     if (g == 1) {
  #       h3_vol = omega * (ar_g$coef[1]+1) + ar_g$coef[1]^2*volatility_mat_h1[oos_i-oos_start-1,g]^2
  #       volatility_mat_h3[oos_i - oos_start + 1,g] = sqrt(abs(h3_vol))
  #     } else if (g == 2) {
  #       h3_vol = omega * (ar_g$coef[1]+1) + (ar_g$coef[1]^2+ar_g$coef[2])*volatility_mat_h1[oos_i-oos_start-1,g]^2 + ar_g$coef[2]*(tail(epsilon_sq_hat,4) %>% head(1))
  #       volatility_mat_h3[oos_i - oos_start + 1,g] = sqrt(abs(h3_vol))
  #     } else {
  #       h3_vol_pt1 = omega * (ar_g$coef[1]+1) + (ar_g$coef[1]^2+ar_g$coef[2])*volatility_mat_h1[oos_i-oos_start-1,g]^2
  #       sum = 0
  #       for (i in 1:(g-2)) {
  #         sum = sum + (ar_g$coef[1]*ar_g$coef[i+1] + ar_g$coef[i+2])*(tail(epsilon_sq_hat,3+i) %>% head(1))
  #       }
  #       h3_vol = h3_vol_pt1 + sum + (ar_g$coef[1]*ar_g$coef[g])*(tail(epsilon_sq_hat,g+1) %>% head(1))
  #       volatility_mat_h3[oos_i - oos_start + 1,g] = sqrt(abs(h3_vol))
  #     }
  #   }
  # }
  
  
  
  ### MODULE B.2: Using optimal p = X, forecast volatility 
  g_opt = 4
  
  # Horizon = 1
  ar_g = arima(epsilon_sq_hat, c(g_opt,0,0), include.mean=FALSE, method = "ML")
  sum = 0
  for (i in 1:g_opt) {
    sum = sum + ar_g$coef[i] * (tail(epsilon_sq_hat,i+1) %>% head(1))
  }
  omega = mu_sq * (1 - sum(ar_g$coef))
  h1_vol = omega + sum
  output_mat[oos_i - oos_start + 1,3] = sqrt(abs(h1_vol))
  
  # Horizon = 2
  if (oos_i >= oos_start + 1) {
    ar_g = arima(epsilon_sq_hat, c(g_opt,0,0), include.mean=FALSE, method = "ML")
    sum = 0
    for (i in 2:g_opt) {
      sum = sum + ar_g$coef[i] * (tail(epsilon_sq_hat,i+1) %>% head(1))
    }
    omega = mu_sq * (1 - sum(ar_g$coef))
    h2_vol = omega + ar_g$coef[1]*output_mat[oos_i-oos_start,3]^2 + sum
    output_mat[oos_i - oos_start + 1,4] = sqrt(abs(h2_vol))
  }
  
  
  # Horizon = 3
  if (oos_i >= oos_start + 2) {
    ar_g = arima(epsilon_sq_hat, c(g_opt,0,0), include.mean=FALSE, method = "ML")
    omega = mu_sq * (1 - sum(ar_g$coef))
    h3_vol_pt1 = omega * (ar_g$coef[1]+1) + (ar_g$coef[1]^2+ar_g$coef[2])*output_mat[oos_i-oos_start-1,3]^2
    sum = 0
    for (i in 1:(g_opt-2)) {
      sum = sum + (ar_g$coef[1]*ar_g$coef[i+1] + ar_g$coef[i+2])*(tail(epsilon_sq_hat,3+i) %>% head(1))
    }
    h3_vol = h3_vol_pt1 + sum + (ar_g$coef[1]*ar_g$coef[g_opt])*(tail(epsilon_sq_hat,g_opt+1) %>% head(1))
    output_mat[oos_i - oos_start + 1,5] = sqrt(abs(h3_vol))
  }
  
  
  
  ### MODULE C: Fitted t / skewed-t distribution parameters
  # t-distribution: df
  output_mat[oos_i - oos_start + 1,8] = stdFit(epsilon_m)[[1]][3]
  
  # skewed t-distribution: df & xi (skewing parameter)
  suppressWarnings(output_mat[oos_i - oos_start + 1,9] <- sstdFit(epsilon_m)[[2]][3])
  suppressWarnings(output_mat[oos_i - oos_start + 1,10] <- sstdFit(epsilon_m)[[2]][4])
  
  
  
  ### MODULE Z: Reporting and failsave 
  print(paste("Currently at out-of-sample index:", oos_i-oos_start+1, "Prediction: ", sum(preds)))
  if (oos_i %% 5 == 0) {
    saveRDS(output_mat, "output_mat.rds")
  }
  
}





# ----------------------------------------------------------------------------------------------------
# Benchmarks + Metrics - Forecasting
# ----------------------------------------------------------------------------------------------------
### MRF
MSPE(output_mat[,"diff"])
MAE(output_mat[,"diff"])
MAD(output_mat[,"diff"])



### Random Walk
RW = rw_function(real_v, horizon = 1)
MSPE(RW)
MAE(RW)
MAD(RW)



### AR(p)
# Tuning p (on MSPE)
for (i in 1:8){
  print(paste(i,": ",MSPE(arP_function(stat_df[,target_var], horizon = 1, p = i))))
}

# Metrics on optimal AR(p): AR(3)
AR_p = arP_function(stat_df[,target_var], horizon = 1, p = 3 )
MSPE(AR_p)
MAE(AR_p)
MAD(AR_p)





# ----------------------------------------------------------------------------------------------------
# Benchmark + Metrics - Volatility
# ----------------------------------------------------------------------------------------------------
# GARCH family rolling estimation
GARCH_h = 1 # Horizon
GARCH_forecast = volatility_benchmarks_vector(c(1,1),c(0,0), horizon = GARCH_h)
ARCH_forecast = volatility_benchmarks_vector(c(4,0),c(0,0), horizon = GARCH_h) # Optimal p = 4
ARMA_GARCH_forecast = volatility_benchmarks_vector(c(1,1),c(1,1), horizon = GARCH_h)



# ARCH(p) tuning
for (i in 1:12) {
  # Specify ARCH(p), re-tune for each p (rugarch)
  GARCH_VaR = volatility_benchmarks_vector(c(1,1),c(0,0), horizon = GARCH_h)
  ChrisLoss(pred_vol = GARCH_VaR[,1], alpha = 0.05, horizon = 1, z = i, t_fit = median(GARCH_VaR[,5]), skew_t_fit_nu = median(GARCH_VaR[,6]), skew_t_fit_xi= median(GARCH_VaR[,7]))
}



# MRF Tuning using volatility matrix created - remember to change horizon / volatility matrix
# for (i in 1:12){
#   ChrisLoss(pred_vol = volatility_mat_h1[,i], horizon = 3, alpha = 0.15,  z = i, t_fit = median(output_mat[,5]), skew_t_fit_nu = median(output_mat[,6]), skew_t_fit_xi= median(output_mat[,7]))
# }





# ----------------------------------------------------------------------------------------------------
# Plotting - In order of the paper
# ----------------------------------------------------------------------------------------------------
# Importing created data sets for plotting
omega_alpha_df = read_rds("output_mat_h1.rds")
volatility_mat_h1 = read_rds("volatility_mat_h1.rds")

# Forecast plotting, non-transformed & first-differenced: Figures 1 & 2
cheap_plotting(output_mat, difBool = FALSE)
cheap_plotting(output_mat, difBool = TRUE)

# Volatility plotting for MRF-ARCH + benchmarjs, change GARCH_h and re-run for each horizon: Figure 4 + E.3
volatility_horizon_m = cbind(GARCH_forecast[,1],ARCH_forecast[,1], ARMA_GARCH_forecast[,1], output_mat[,"pred_vol_h1"])
volatility_horizon_m = cbind(GARCH_forecast[,1],ARCH_forecast[,1], ARMA_GARCH_forecast[,1], output_mat[,"pred_vol_h2"])
volatility_horizon_m = cbind(GARCH_forecast[,1],ARCH_forecast[,1], ARMA_GARCH_forecast[,1], output_mat[,"pred_vol_h3"])
cheap_benchmark_plotting(volatility_horizon_m)

# VaR plotting vs benchmarks (V1: ARCH, V2: GARCH/ARMA-GARCH), for different GARCH_h and re-run: Figure 5 + E.5 + E.6 
VaR_benchmark_plotting(q = 0.05, MRF = output_mat, GARCH = GARCH_forecast, ARCH = ARCH_forecast, ARMA_GARCH = ARMA_GARCH_forecast, archBool = TRUE)
VaR_benchmark_plotting(q = 0.05, MRF = output_mat, GARCH = GARCH_forecast, ARCH = ARCH_forecast, ARMA_GARCH = ARMA_GARCH_forecast, archBool = FALSE)

# ARCH parameter plotting: Figure 6
arch_parameter_plotting(MRF_ARCH = omega_alpha_df, ARCH = ARCH_forecast, onlyAlpha = TRUE)
arch_parameter_plotting(MRF_ARCH = omega_alpha_df, ARCH = ARCH_forecast, onlyAlpha = FALSE)

# Forecasting MRF vs AR(3) (benchmark). Change horizon (and p-order) in Benchmark section: Figures E.2
benchmark_plotting(output_mat[,"diff"], stat_df[oos_start:oos_end,target_var], AR_p)


