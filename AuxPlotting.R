# ----------------------------------------------------------------------------------------------------
# MRF Cheap-Plotting of pred / log-dif pred
# ----------------------------------------------------------------------------------------------------
# Plotting the output, adaptive to how many rows are fulfilled, for either real vs transformed data
cheap_plotting <- function(output, difBool) {
  # Acquire starting matrix
  nonEmpty = output %>% 
    as.data.frame %>%
    rownames_to_column() %>%
    mutate(date = as.Date(rowname)) %>%
    dplyr::select("date","no_diff", "diff")
  
  # Get date of comparison matrix through indices
  start_date = nonEmpty[,"date"] %>% head(1)
  len = nrow(nonEmpty) - 1 # How many time periods after start to consider for comparison matrix
  
  # Create divide whether using differenced data or not  
  if (difBool) { #First difference
    comparison_m_start = which(stat_df$date == start_date)
    comparison_m = stat_df[comparison_m_start:(comparison_m_start+len),target_var]
    
    combined_m = cbind(nonEmpty[,c("date","diff")], comparison_m)
    melted = melt(as.data.frame(combined_m), id.vars="date")
    
    # Plotting using melted data
    ggplot(melted, aes(x = date, y = value, color = variable)) +
      geom_line(aes(y = 0), color = "grey", linetype = "longdash") +
      geom_line(show.legend = FALSE) +
      labs(x = "Year", y = "Δ ExR USD/GBP", color = "Variable") +
      scale_color_brewer(palette = "Dark2")+ theme_classic() +
      theme(axis.title = element_text(size = 12), axis.title.x=element_blank())
  } else { #Non-first difference
    
    comparison_m_start = which(fred_v_notrans$date == start_date)
    comparison_m = fred_v_notrans[comparison_m_start:(comparison_m_start+len),target_var]
    
    combined_m = cbind(nonEmpty[,c("date","no_diff")], comparison_m)
    melted = melt(as.data.frame(combined_m), id.vars="date")
    
    # Plotting using melted data
    ggplot(melted, aes(x = date, y = value, color = variable)) +
      geom_line(show.legend = FALSE) +
      labs(x = "Year", y = "ExR USD/GBP", color = "Variable") +
      scale_color_brewer(palette = "Dark2")+ theme_classic() +
      theme(axis.title = element_text(size = 12), axis.title.x=element_blank())
  }
} 



# ----------------------------------------------------------------------------------------------------
# ARCH Parameter plotting
# ----------------------------------------------------------------------------------------------------
# Plotting the output, adaptive to how many rows are fulfilled, for either real vs transformed data
arch_parameter_plotting <- function(MRF_ARCH, ARCH, onlyAlpha) {
  
  if (!onlyAlpha) {
    # Scaling (done manually) for two-axis
    combined = cbind(ARCH[,"omega"]/6+0.000635, MRF_ARCH[,"omega"])
    rownames(combined) = rownames(MRF_ARCH) 
    combined = as.data.frame(combined) %>% rownames_to_column()
    combined$rowname <- as.Date(combined$rowname, format = "%Y-%m-%d")
    plot_df = melt(as.data.frame(combined), id.vars="rowname")
    
    # Plotting
    ggplot(plot_df, aes(x = rowname, y = value, color = variable)) +
      geom_line(show.legend = FALSE, size=0.65) +
      labs(x = "Year", y = "Volatility (%)", color = "Variable") +
      theme_classic() +
      scale_color_manual(values = c("#fabb22", "#1b9e77"), breaks = c("V1", "V2")) +
      theme(axis.title = element_text(size = 12), axis.title.x=element_blank(), axis.title.y=element_blank()) +
      scale_y_continuous(
        name = "x",
        sec.axis = sec_axis(~ (. - 0.000635)*6 , name = "VaR (%)")
      )
  } else {
    # Scaling (done manually) for two-axis
    combined = cbind(ARCH[,"alpha"]/5.5, MRF_ARCH[,"alpha"])
    rownames(combined) = rownames(MRF_ARCH) 
    combined = as.data.frame(combined) %>% rownames_to_column()
    combined$rowname <- as.Date(combined$rowname, format = "%Y-%m-%d")
    plot_df = melt(as.data.frame(combined), id.vars="rowname")
    
    # Plotting
    ggplot(plot_df, aes(x = rowname, y = value, color = variable)) +
      geom_line(show.legend = FALSE, size = 0.65) +
      labs(x = "Year", y = "Volatility (%)", color = "Variable") +
      theme_classic() +
      scale_color_manual(values = c("#fabb22", "#1b9e77"), breaks = c("V1", "V2")) +
      theme(axis.title = element_text(size = 12), axis.title.x=element_blank(), axis.title.y=element_blank()) +
      scale_y_continuous(
        name = "x",
        sec.axis = sec_axis(~ . * 5.5  , name = "VaR (%)")
      )
  }
}



# ----------------------------------------------------------------------------------------------------
# VaR plotting Fixed
# ----------------------------------------------------------------------------------------------------
VaR_plotting_fixed <- function(pred_vol) {
  # Creating the boundaries
  boundary_5 = real_v_shifted + qnorm(0.05) * pred_vol
  boundary_10 = real_v_shifted + qnorm(0.1) * pred_vol
  boundary_15 = real_v_shifted + qnorm(0.15) * pred_vol
  time = stat_df[(oos_start):(oos_start+length(pred_vol)-1),"date"]
  combined_m = data.frame(time, real_v, boundary_5, boundary_10, boundary_15)
  colnames(combined_m) = c("time","Forecasted volatility (h=1)","alpha 0.05", "alpha 0.1", "alpha 0.15")
  melted = melt(as.data.frame(combined_m), id.vars="time")
  
  # Creating the points that indicate boundary crossings
  boundary_5_points <- data.frame(x = time[real_v < boundary_5], y = boundary_5[real_v < boundary_5])
  boundary_10_points <- data.frame(x = time[real_v < boundary_10], y = boundary_10[real_v < boundary_10])
  boundary_15_points <- data.frame(x = time[real_v < boundary_15], y = boundary_15[real_v < boundary_15])
  
  # Plotting
  ggplot(data = melted, aes(x = time, y = value)) +
    geom_line(aes(color = variable), show.legend = FALSE) +
    geom_line(aes(y = 0), color = "grey", linetype = "longdash") +
    labs(x = "Year", y = "Volatility") +
    theme_classic() +
    scale_color_brewer(palette = "Dark2") +
    theme(axis.title = element_text(size = 12), axis.title.x = element_blank()) +
    geom_point(data = boundary_5_points, aes(x = x, y = y), color = "orange", size = 4, alpha=0.6)  +
    geom_point(data = boundary_10_points, aes(x = x, y = y), color = "purple", size = 4, alpha=0.6) +
    geom_point(data = boundary_15_points, aes(x = x, y = y), color = "pink", size = 4, alpha=0.6) 
}



# ----------------------------------------------------------------------------------------------------
# Benchmark comparison plotting
# ----------------------------------------------------------------------------------------------------
# Plotting the output, adaptive to how many rows are fulfilled, for either real vs transformed data
benchmark_plotting <- function(MRF, real, AR) {

  time = stat_df[(oos_start):(oos_start+length(MRF)-1),"date"]
  combined_m = cbind(time, MRF, real, AR)
  colnames(combined_m) = c("date","MRF","Real","AR(3)")
  
  melted = melt(as.data.frame(combined_m), id.vars="date")
  
  
  # Plotting using melted data
  ggplot(melted, aes(x = date, y = value, color = variable)) +
    geom_line(aes(y = 0), color = "grey", linetype = "longdash") +
    geom_line(show.legend = FALSE) +
    labs(x = "Year", y = "Δ ExR USD/GBP", color = "Variable") +
    scale_color_brewer(palette = "Dark2")+ theme_classic() +
    theme(axis.title = element_text(size = 12), axis.title.x=element_blank())
} 



# ----------------------------------------------------------------------------------------------------
# Cheap benchmark comparison plotting (V2)
# ----------------------------------------------------------------------------------------------------
# Plotting the output, adaptive to how many rows are fulfilled, for either real vs transformed data
cheap_benchmark_plotting <- function(df1) {
  # Acquire starting matrix
  date = stat_df[oos_start:oos_end,1]
  combined_m = data.frame(date,df1)
  colnames(combined_m) = c("date", "GARCH","ARCH","ARMA-GARCH", "MRF")
  melted = melt(combined_m, id.vars="date")
  
  # Plotting using melted data
  ggplot(melted, aes(x = date, y = value, color = variable)) +
    geom_line(show.legend = FALSE) +
    labs(x = "Year", y = "Volatility", color = "Variable") +
    theme_classic() +
    scale_color_brewer(palette = "Dark2") +
    theme(axis.title = element_text(size = 12), axis.title.x=element_blank())
} 





# ----------------------------------------------------------------------------------------------------
# VaR of MRF-ARCH vs Benchmarks
# ----------------------------------------------------------------------------------------------------
VaR_benchmark_plotting <- function(q, MRF, GARCH, ARCH, ARMA_GARCH, archBool) {
  MRF_boundary = output_mat[,"pred_vol_h1"] * qt(q, 3) + real_v_shifted
  GARCH_boundary = output_mat[,"pred_vol_h1"] * qt(q, 4) + real_v_shifted
  ARCH_boundary = output_mat[,"pred_vol_h1"] * qnorm(q) + real_v_shifted
  ARMA_GARCH_boundary = output_mat[,"pred_vol_h1"] * qnorm(q) + real_v_shifted
  
  boundary_MRF_points <- data.frame(x = as.Date(rownames(output_mat))[real_v < MRF_boundary], y = MRF_boundary[real_v < MRF_boundary])
  boundary_GARCH_points <- data.frame(x = as.Date(rownames(output_mat))[real_v < GARCH_boundary], y = GARCH_boundary[real_v < GARCH_boundary])
  boundary_ARMA_GARCH_points <- data.frame(x = as.Date(rownames(output_mat))[real_v < ARMA_GARCH_boundary], y = ARMA_GARCH_boundary[real_v < ARMA_GARCH_boundary])
  boundary_ARCH_points <- data.frame(x = as.Date(rownames(output_mat))[real_v < ARCH_boundary], y = ARCH_boundary[real_v < ARCH_boundary])
  
  if (!archBool) {
    combined = data.frame("time" = stat_df[oos_start:oos_end,"date"], real_v, MRF_boundary, GARCH_boundary, ARMA_GARCH_boundary  )
    melted = melt(combined, id.vars="time")
    ggplot(melted, aes(x = time, y = value, color = variable)) +
      geom_line(aes(y = 0), color = "grey", linetype = "longdash") +
      geom_line(show.legend = FALSE) +
      scale_color_brewer(palette = "Dark2") + theme_classic() +
      labs(x = "Year", y = "Volatility", color = "Variable") +
      theme(axis.title = element_text(size = 12), axis.title.x=element_blank()) +
      geom_point(data = boundary_MRF_points, aes(x = x, y = y), color = "orange", size = 4, alpha=0.5) +
      geom_point(data = boundary_GARCH_points, aes(x = x, y = y), color = "purple", size = 4, alpha=0.5) +
      geom_point(data = boundary_ARMA_GARCH_points, aes(x = x, y = y), color = "magenta", size = 4, alpha=0.5)
  } else {
    combined = data.frame("time" = stat_df[oos_start:oos_end,"date"], real_v, MRF_boundary, ARCH_boundary  )
    melted = melt(combined, id.vars="time")
    ggplot(melted, aes(x = time, y = value, color = variable)) +
      geom_line(aes(y = 0), color = "grey", linetype = "longdash") +
      geom_line(show.legend = FALSE) +
      scale_color_brewer(palette = "Dark2")+ theme_classic() +
      labs(x = "Year", y = "Volatility", color = "Variable") +
      theme(axis.title = element_text(size = 12), axis.title.x=element_blank()) +
      geom_point(data = boundary_MRF_points, aes(x = x, y = y), color = "orange", size = 4, alpha=0.5) +
      geom_point(data = boundary_ARCH_points, aes(x = x, y = y), color = "purple", size = 4, alpha=0.5) 
  }
}




# ----------------------------------------------------------------------------------------------------
# Cheap Vintage plotting 
# ----------------------------------------------------------------------------------------------------
# Plotting the output, adaptive to how many rows are fulfilled, for either real vs transformed data
cheap_vintage_plotting <- function(df1, names) {
  # Acquire starting matrix
  date = seq(as.Date("2014-12-01"), as.Date("2022-01-01"), by = "month")
  combined_m = data.frame(date,df1)
  colnames(combined_m) = names
  melted = melt(combined_m, id.vars="date")
  
  # Plotting using melted data
  ggplot(melted, aes(x = date, y = value, color = variable)) +
    geom_line(show.legend = FALSE) +
    #geom_line() +
    labs(x = "Year", y = "Volatility (%)", color = "Variable") +
    theme_classic() +
    scale_color_brewer(palette = "Dark2") +
    labs(x = "Year", y = "Δ CPI (All items)", color = "Variable") +
    theme(axis.title = element_text(size = 12), axis.title.x=element_blank())
} 
