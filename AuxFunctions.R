# ----------------------------------------------------------------------------------------------------
# FRED Data Collection function
# ----------------------------------------------------------------------------------------------------
getFRED <- function(start_date, end_date, doResource, target_var, wantUK) {
  # FRED US, through website sourcing
  fred_us <- fredmd(file = "https://files.stlouisfed.org/files/htdocs/fred-md/monthly/2024-02.csv",
                    transform = TRUE,
                    date_start = ymd(as.Date("1998-01-01")))
  fred_us_names <- read_csv("https://files.stlouisfed.org/files/htdocs/fred-md/monthly/2024-02.csv")
  colnames(fred_us) <- colnames(fred_us_names)
  
  # Replace the target variable in the dataset with the non-transformed data
  if (doResource) {
    fred_mat = fredmd(file = "https://files.stlouisfed.org/files/htdocs/fred-md/monthly/2024-02.csv",
                      transform = FALSE,
                      date_start = ymd(as.Date("1998-01-01")))
    
    # Collect the non-transformed variable
    sourced_variable = fred_mat[,target_var]
    
    # Replace the transformed by non-transformed
    fred_us[,target_var] = sourced_variable
  }
  
  # Do you want to include the fredMD-UK dataset (only 'no' for robustness check)
  if (wantUK) {
    # FRED UK, through local file. Data taken from Coulombe (2019) 
    fred_uk <- read.csv("balanced_uk_md.csv", header = TRUE, row.names = 1)
    colnames(fred_uk) <- paste("UK_",colnames(fred_uk))
    fred_uk %<>% dplyr::select(-`UK_ GBP_US`)
    
    # Merge FRED US and FRED UK 
    fred_merged <- merge(fred_us, fred_uk, by.x="sasdate", by.y="UK_ Date")
    
    # Filter dataset on dates (+ rename date column for nicer formatting)
    fred_merged = fred_merged %>% 
      rename(date=sasdate) %>%
      filter(date >= start_date & date <= end_date) 
    
    return(fred_merged)
  } else {
    fred_us = fred_us %>% 
      rename(date=sasdate) %>%
      filter(date >= start_date & date <= end_date) 
    
    return(fred_us)
  }
}



getFredRaw <- function(start_date, end_date, target_var) {
  # FRED US, through website sourcing
  fred_us <- fredmd(file = "https://files.stlouisfed.org/files/htdocs/fred-md/monthly/2024-02.csv",
                    transform = FALSE,
                    date_start = ymd(as.Date("1998-01-01")))
  fred_us_names <- read_csv("https://files.stlouisfed.org/files/htdocs/fred-md/monthly/2024-02.csv")
  colnames(fred_us) <- colnames(fred_us_names)
  
  fred_vector <- fred_us %>%
    data.frame() %>%
    rename(date = sasdate) %>%
    dplyr::filter(date >= start_date & date <= end_date) %>%
    dplyr::select(date, target_var)
  
  return(fred_vector)
}


# ----------------------------------------------------------------------------------------------------
# Interpolation for NAs in the middle + Stationarity conversion, using alpha = 0.05
# ----------------------------------------------------------------------------------------------------
# Define function to make stationarity
make_stationary <- function(vector) {
  
  # Check for stationary vector
  if (adf.test(as.matrix(vector[-1]))[4] > 0.05) {
    
    # If yes: First Difference
    diff_vector = c(NA,diff(vector))
    
    #Check: Still stationarity, if yes: error!
    if (adf.test(as.matrix(diff_vector[-c(1,2)]))[4] > 0.05) {
      print("SVEN SAYS: STATIONARITY PROBLEM IN A SERIES!")
    }
    return(diff_vector)
  }
  return(vector)
}

# Tidyverse function, collective to interpolate and create stationary series
interpolate_stationary <- function(df) {
  # Interpolation, then prepare for stationarity test & conversion
  interpolated_df = df %>% 
    mutate(across(-date, ~na.approx(.,na.rm = F))) %>% # Interpolate missing NAs, NO leading/trailing interpolation
    mutate(across(everything(), make_stationary)) %>%
    drop_na() %>%
    suppressWarnings()
  return(interpolated_df)
} 





# ----------------------------------------------------------------------------------------------------
# STEP 1: LARS
# ----------------------------------------------------------------------------------------------------
# Define and run LARS function
LARS <- function(df,n) {
  df <- df %>% dplyr::select(-date)
  
  x <- df %>% dplyr::select(-target_var) %>% as.matrix()
  y <- df %>% dplyr::select(target_var) %>% as.matrix()
  l1 <- lars(x=x, y=y, type="lar")
  
  # have to change this, because its copied from nowcasting paper
  out <- as.data.frame(coef(l1)) %>%
    summarise(across(everything(),~sum(.==0))) %>%
    t() %>%
    as.data.frame() %>%
    arrange(V1) %>%
    filter(V1 <= n) %>%
    rownames()
  
  return(out)
}





# ----------------------------------------------------------------------------------------------------
# STEP 2: PCA
# ----------------------------------------------------------------------------------------------------
#PCA Function
PCA_function <- function(mat) {
  # Select X-matrix
  x_mat <- mat %>%
    dplyr::select(-date)
  
  # Perform PCA
  x_pca <- prcomp(x_mat, scale. = TRUE) 
  
  # Select number of factors with eigenvalue > 1
  pca_select <- data.frame(sd = x_pca$sdev) %>%
    filter(sd^2 > 1)
  
  # Select the factors from set
  x_factors <- x_pca$rotation[, 1:nrow(pca_select)]
  
  # Matrix multiply X by eigenvectors => rotation matrix, which we want as output
  factors <- as.matrix(x_mat) %*% x_factors
  
  #return(x_factors)
  return(factors)
}





# ----------------------------------------------------------------------------------------------------
# STEP 2: MAF
# ----------------------------------------------------------------------------------------------------
# Retrieve the MAF data: the expanding dataset from oos_start to oos_end, but including an extra 60M = 5Y for the MAF implementation
MAF_data <- function() {
  # Retrieve data that is - 5 years from the starting point of the
  fred_us_MAF <- fredmd(file = "https://files.stlouisfed.org/files/htdocs/fred-md/monthly/2024-02.csv",
                        transform = TRUE,
                        date_start = ymd(as.Date("1997-01-01"))) # Note that this is exactly 60 months / 5 years before -> MAF-length = 5
  fred_us_names <- read_csv("https://files.stlouisfed.org/files/htdocs/fred-md/monthly/2024-02.csv")
  colnames(fred_us_MAF) <- colnames(fred_us_names)
  MAF_m = fred_us_MAF %>% cbind() %>%
    rename(date=sasdate) %>%
    filter(date >= "1997-02-01" & date <= "2023-05-01") %>% # Dates are chosen, so its perfectly 5Y before start of regular data
    interpolate_stationary()
  
  return(MAF_m)
}



MAF_function <- function(df, k, start_date, end_date) {
  # Select X-matrix
  x_mat <- df %>%
    filter( date >= start_date & date <= end_date) %>%
    dplyr::select(-date) %>%
    as.matrix() %>%
    # embed(12)
    embed(13)
  
  # Perform PCA on combined mat
  x_pca <- prcomp(x_mat, scale. = TRUE, rank. = k)$x
  #x_pca = trlan.svd(x_mat, neig = 10)$u # Using Lanczos algorithm, truncate PCA at 10
  colnames(x_pca) = paste("MAF_",1:k,sep="")
  return(x_pca)
}



MAF_function_vintage <- function(df, k, start_date, end_date) {
  # Select X-matrix
  x_mat <- df %>%
    filter( date >= start_date & date <= end_date) %>%
    dplyr::select(-date) %>%
    as.matrix() %>%
    embed(12)
  
  # Perform PCA on combined mat
  x_pca <- prcomp(x_mat, scale. = TRUE, rank. = k)$x
  #x_pca = trlan.svd(x_mat, neig = 10)$u # Using Lanczos algorithm, truncate PCA at 10
  colnames(x_pca) = paste("MAF_",1:k,sep="")
  return(x_pca)
}




# ----------------------------------------------------------------------------------------------------
# Christophsen + Loss function
# ----------------------------------------------------------------------------------------------------
ChrisLoss = function(pred_vol, alpha, z, t_fit, skew_t_fit_nu, skew_t_fit_xi, horizon) {
  # Fitting distribution
  quantile_n = qnorm(alpha)
  quantile_t = qstd(alpha, nu = t_fit)
  quantile_sk_t = qsstd(alpha, nu = skew_t_fit_nu, xi = skew_t_fit_xi)
  quantiles = c(quantile_n, quantile_t, quantile_sk_t)
  quantile_names = c("norm", "t", "skewed-t")
  
  for (i in 1:3) {
    boundary = real_v_shifted + quantiles[i] * pred_vol
    
    test_stat = BacktestVaR(data = real_v, VaR = boundary, alpha = alpha)
    uc_stat = test_stat$LRuc[[2]]
    cc_stat = test_stat$LRcc[[2]]
    
    indep_stat = 1-pchisq(test_stat$LRcc[[1]] - test_stat$LRuc[[1]],1)
    
    if (uc_stat > alpha & cc_stat > alpha & indep_stat > alpha) {
      print(paste(z, alpha, quantile_names[i], RLF(real_v[-c(1:horizon)], boundary[-c(1:horizon)])))
    }
  }
}