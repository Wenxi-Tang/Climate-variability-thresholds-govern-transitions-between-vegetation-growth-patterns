
setwd('/data/share/tangwenxi')
group_by_n <- 100
RUN_name <- paste0('GroupBy_',group_by_n,'_SlidingWindow_Claude-all')

CLI_names <- c("TMP","PRE","SR","VPD","CO2","SMroot","SMsurf") 
VI_name <- "GPP"
year_start <- 1982
year_end <- 2020
window_size <- 10  # Sliding window size is 10 years

# Parallel settings
max_cores <- 80

# Load packages -----
library(parallel)
library(doParallel)
library(foreach)
library(vars) # VAR
library(tseries) # adf.test
library(urca) # ur.df() & vecm
library(dplyr)
library(tidyr)
library(reshape2)
library(ggplot2)
library(missForest)
library(zoo)
library(purrr)

# Memory and system monitoring functions -----
check_memory <- function() {
  mem_info <- system("free -g", intern = TRUE)
  available_mem <- as.numeric(strsplit(mem_info[2], "\\s+")[[1]][7])
  return(available_mem)
}

clear_memory <- function() {
  gc()
  system("sync && echo 3 | sudo tee /proc/sys/vm/drop_caches > /dev/null 2>&1", ignore.stderr = TRUE)
  cat("Memory cleared\n")
}

get_available_cores <- function(max_cores = 80) {
  total_cores <- detectCores()
  load_avg <- as.numeric(system("cat /proc/loadavg | awk '{print $1}'", intern = TRUE))
  available_cores <- min(max_cores, max(1, floor(total_cores - load_avg)))
  cat("Total cores:", total_cores, "| Load average:", load_avg, "| Using cores:", available_cores, "\n")
  return(available_cores)
}

# Progress bar function -----
print_progress <- function(current, total, prefix = "Progress") {
  percent <- round((current / total) * 100, 1)
  bar_length <- 50
  filled_length <- round((current / total) * bar_length)
  
  bar <- paste(c(rep("█", filled_length), rep("░", bar_length - filled_length)), collapse = "")
  
  cat(sprintf("\r%s: [%s] %s%% (%d/%d)", prefix, bar, percent, current, total))
  
  if (current == total) {
    cat("\n")
  }
  flush.console()
}

# Create directories -----
dir.create('./Out-Table-and-Figure-VPD et al', showWarnings = FALSE)
Datedir_main <- file.path('./Out-Table-and-Figure-VPD et al', Sys.Date())
dir.create(Datedir_main, showWarnings = FALSE)
Datedir <- file.path(Datedir_main, RUN_name)
dir.create(Datedir, showWarnings = FALSE)
VI_dir <- file.path(Datedir, VI_name)
dir.create(VI_dir, showWarnings = FALSE)

# Error log file
error_file <- file.path(VI_dir, "error_log.txt")
cat("Error Log Started at:", as.character(Sys.time()), "\n", file = error_file)

# Error logging function -----
log_error <- function(error_msg, error_file) {
  cat(paste0(rep("=", 80), collapse = ""), "\n", file = error_file, append = TRUE)
  cat("Time:", as.character(Sys.time()), "\n", file = error_file, append = TRUE)
  cat("Error:", error_msg, "\n", file = error_file, append = TRUE)
  cat(paste0(rep("=", 80), collapse = ""), "\n\n", file = error_file, append = TRUE)
}

# Define functions -----
# 1.ADFtest [Stabilization: test<1pct|5pct|10pct]
ADF_test <- function(data){
  output_ADF <- c()
  
  for(typename in c("none","trend","drift")){
    output <- c()
    for(i in 1:ncol(data)){
      urt.i <- ur.df(data[,i],type=typename,selectlags='AIC')
      pct_min <- max(summary(urt.i)@cval[1,])
      DF_test <- summary(urt.i)@teststat[1]
      if(DF_test<pct_min){
        output <- c(output,1)
      }else{
        output <- c(output,0)
      }
    }
    
    table_output = length(table(output))
    if(table_output == 1){
      output_ADF = c(output_ADF, typename)
    }else{
      output_ADF = "Non-stationary"
    }
  }
  
  if(length(output_ADF) > 1 && any(output_ADF %in% c("trend","drift"))){
    output_ADF = "both"
  } else if(length(output_ADF) == 1 && output_ADF == "drift"){
    output_ADF = "const"
  } else if(length(output_ADF) == 0 || output_ADF[1] == "Non-stationary"){
    output_ADF = "both"
  } else {
    output_ADF = output_ADF[1]
  }
  
  return(output_ADF)
}

# Read data -----
data_dir <- './Data'
load("./Data/combined_list-VPD+CO2+SM.RData")
df.GLSP_GS <- read.csv('./Data/LSP_month.csv')

############################################-
#                VAR model                 #
############################################-
n_ahead = 12
CLI_i = 3:9

# Group processing - modified to 100 grids per group ----
lab_groups <- split(1:length(combined_list), 
                    rep(1:ceiling(length(combined_list)/group_by_n), 
                        each = group_by_n))

# Get time and clear memory before starting program
start_time <- Sys.time()
cat("Program started at:", as.character(start_time), "\n")
# clear_memory()

# Define sliding windows ----
window_starts <- seq(year_start, year_end - window_size + 1, by = 1)
windows <- map(window_starts, function(start) {
  list(
    start = start,
    end = start + window_size - 1,
    years = start:(start + window_size - 1)
  )
})

# Define processing function for a single group ----
VI_i = 10
w = 1
group_i = 7
g_i = 51
process_group <- function(group_i, VI_i, current_window, lab_groups, 
                          combined_list, df.GLSP_GS, CLI_names, 
                          n_ahead, year_start, VI_dir, error_file) {
  
  tryCatch({
    # Load necessary packages (in each worker) -----
    library(vars)
    library(tseries)
    library(urca)
    library(dplyr)
    library(tidyr)
    library(reshape2)
    library(ggplot2)
    library(missForest)
    library(zoo)
    
    VI_name <- colnames(combined_list[[1]])[VI_i]
    window_label <- paste0(current_window$start, "_", current_window$end)
    
    lab_group <- lab_groups[[group_i]]
    
    VAR_all_result_list <- list()
    VAR_data_list <- list()
    
    lag_RS_all <- data.frame()
    var_R2.all <- data.frame()
    duration_all <- data.frame()
    causal_all <- data.frame()
    q_all <- data.frame()
    ADF_all <- data.frame()
    johansen_all <- data.frame()
    var_root_all <- data.frame()
    irf_all <- data.frame()
    fevd_all <- data.frame()
    
    # Store results for current window
    window_gpp_response <- data.frame()
    window_gpp_irf <- data.frame()
    
    for (g_i in seq_along(lab_group)) {
      i <- lab_group[g_i]
      tryCatch({
        VAR_all_result_list[[g_i]] <- NA
        result_list <- list()
        x <- i
        
        # (1) Extract data for the model -----
        data <- combined_list[[x]]
        
        data <- data %>%
          separate(Date, into = c("Year", "Month"), sep = "_")
        data$Year <- as.numeric(data$Year)
        data$Month <- as.numeric(data$Month)
        
        # Processing cross-year and non-cross-year data ------
        df.GS_lab <- df.GLSP_GS[which(df.GLSP_GS$lab == unique(data$lab)),]
        SOS <- df.GS_lab$SOS_Month
        EOS <- df.GS_lab$EOS_Month
        target_sequence <- as.numeric(unlist(strsplit(df.GS_lab$GS_Months, ",")))
        seq_length <- length(target_sequence)
        
        # Determine if a complete growing season data exists: Mark target sequence
        data.test <- data %>%
          mutate(
            is_sequence_start = rollapply(
              Month, 
              width = seq_length, 
              function(x) all(x == target_sequence), 
              align = "left", 
              fill = NA
            ),
            Year_lab = NA
          )
        
        # Determine start row positions to mark
        start_pos <- which(data.test$is_sequence_start == TRUE)
        
        # If target sequence is found, mark seq_length cells as "T" from start position
        if (length(start_pos) > 0) {
          for(pos_i in 1:length(start_pos)){
            start_pos_i <- start_pos[pos_i]
            data.test$Year_lab[start_pos_i:(start_pos_i + seq_length - 1)] <- "T"
          }
        }
        
        # Remove temporary columns
        data.test <- data.test %>% dplyr::select(-is_sequence_start)
        
        # Determine if incomplete growing season data exists: Check for rows where year_lab is NA
        if(anyNA(data.test$Year_lab)){
          data.test$Across_years <- T
          data <- data.test %>%
            filter(!is.na(Year_lab) & Year_lab != "")
          season_year <- 38
          season_n <- nrow(data)/38
          
          # Re-add year column
          season_year_a <- 0
          for (season_year_i in 1:season_year) {
            season_year_a <- season_year_a + 1
            season_year_b <- season_year_a + (season_n -1)
            data$year_new[season_year_a:season_year_b] <- year_start + season_year_i
            season_year_a <- season_year_b
          }
        }else{
          data.test$Across_years <- F
          data <- data.test %>%
            filter(!is.na(Year_lab) & Year_lab != "")
          season_year <- 39
          season_n <- nrow(data)/39
          data$year_new <- data$Year
        }
        
        VAR_data_list[[g_i]] <- data
        
        # Filter data for VAR -----
        data <- data %>% subset(year_new %in% current_window$years)
        data <- data[, c(CLI_names, VI_name)]
        data[which(data[,VI_name] == 0), VI_name] <- NA
        
        na_num <- sum(is.na(data))
        season_n_actual <- nrow(data)/length(current_window$years)
        
        # Check if any column has all identical values
        constant_columns <- names(data)[sapply(data, function(col) all(col == col[1]))]
        
        if (length(constant_columns) > 0) {
          VAR_all_result_list[[g_i]] <- paste(
            "Data contains constant column(s):", 
            paste(constant_columns, collapse = ", "), 
            "| Skipping file:", i, 
            "| GS month:", season_n_actual, 
            "| NAs number:", na_num
          )
          next
        }
        
        # Skip if more than 2 NAs
        if(na_num >= 2){
          VAR_all_result_list[[g_i]] <- paste("Data contains more than 2 NAs, skipping file:", i, 
                                              "GS month: ", season_n_actual, "NAs number: ", na_num)
          next
        }
        
        lab <- unique(combined_list[[x]][,2])
        
        # Skip if growing season is shorter than 3 months
        if(season_n_actual < 3){
          VAR_all_result_list[[g_i]] <- paste("The growing season is shorter than 3 months, skipping file:", i,  
                                              "GS month: ", season_n_actual)
          next
        }
        
        # Fill missing values using Random Forest
        if(any(is.na(data))){
          data <- missForest(data)$ximp
        }
        
        result_list[[1]] <- data
        
        # (2) Unit Root Test -----
        lag_type <- ADF_test(data)
        if(lag_type == "Non-stationary"){
          lag_type = "both"
        }
        
        # Determine season parameter
        season_param <- as.integer(season_n_actual)
        
        # Detailed ADF test results
        for (variable_i in 1:ncol(data)) {
          for(typename in c("none","trend","drift")){
            urt.i <- ur.df(data[,variable_i], type=typename, selectlags='AIC')
            
            adf_results <- data.frame(
              VI_name = VI_name,
              lab = lab,
              Variable_name = colnames(data)[variable_i],
              typename = typename,
              Level = c("1%", "5%", "10%"),
              Test_Statistic = round(summary(urt.i)@teststat[1], 3),
              Critical_Value = c(summary(urt.i)@cval[1,]),
              window_start = current_window$start,
              window_end = current_window$end,
              group_i = group_i,
              grid_id = i
            )
            
            adf_results$Stationary <- adf_results$Test_Statistic < adf_results$Critical_Value
            ADF_all <- rbind(ADF_all, adf_results)
          } 
        }
        result_list[[2]] <- ADF_all
        
        # (3) Johansen Cointegration Test -----
        var_model <- VAR(data, p = 1, type = "both", season = as.integer(nrow(data)/length(current_window$years)))
        coeff_array <- as.data.frame(var_model$varresult[[length(var_model$varresult)]]$coefficients)
        if(any(is.na(coeff_array))){
          VAR_all_result_list[[g_i]] <- paste0("VAR model result have NA, skipping file: ", i)
          next
        }
        
        johansen_test <- tryCatch({
          ca.jo(data, type = "trace", ecdet = "trend", K = 2)
        }, error = function(e) {
          return(list(error = TRUE, message = e$message))
        })
        
        # Check if an error occurred
        if(inherits(johansen_test, "list") && !is.null(johansen_test$error) && johansen_test$error) {
          johansen_result_df <- data.frame(
            VI_name = VI_name,
            lab = lab,
            error = TRUE,
            message = johansen_test$message,
            window_start = current_window$start,
            window_end = current_window$end,
            group_i = group_i,
            grid_id = i
          )
        } else {
          johansen_test_summary <- summary(johansen_test)
          test_stats <- johansen_test_summary@teststat
          critical_vals_5pct <- johansen_test_summary@cval[, "5pct"]
          
          johansen_result_df <- data.frame(
            VI_name = VI_name,
            lab = lab,
            typename = "trend",
            r = paste("r <=", (length(critical_vals_5pct)-1):0),
            Test_Statistic = test_stats,
            Critical_Value_5pct = critical_vals_5pct,
            window_start = current_window$start,
            window_end = current_window$end,
            group_i = group_i,
            grid_id = i
          )
          
          johansen_result_df$Cointegration <- johansen_result_df$Test_Statistic > johansen_result_df$Critical_Value_5pct
        }
        
        johansen_all <- rbind(johansen_all, johansen_result_df)
        result_list[[3]] <- johansen_result_df
        
        # (4) Determine lag order -----
        lag_num <- 1
        
        # Store
        lag_RS <- data.frame(VARselect(data, type = "both", lag.max = 10)$criteria)
        lag_RS$VI_name <- VI_name
        lag_RS$lab <- lab
        lag_RS$window_start <- current_window$start
        lag_RS$window_end <- current_window$end
        lag_RS$group_i <- group_i
        lag_RS$grid_id <- i
        lag_RS_all <- rbind(lag_RS_all, lag_RS)
        result_list[[4]] <- lag_RS
        
        # (5) Build VAR model -----
        var_model <- VAR(data, p = 1, type = "both", season = as.integer(nrow(data)/length(current_window$years)))
        coeff_array <- as.data.frame(var_model$varresult[[length(var_model$varresult)]]$coefficients)
        if(any(is.na(coeff_array))){
          VAR_all_result_list[[g_i]] <- paste0("VAR model result have NA, skipping file: ", i)
          next
        }
        var_result <- summary(var_model)
        
        # Overall p-value of VAR
        var_summary_text <- capture.output(summary(var_model))
        var_p_value_line <- grep("p-value:", var_summary_text, value = TRUE)[1]
        var_p_value <- sub(".*p-value: ", " ", var_p_value_line)
        
        # VAR regression R2
        var_R2 <- data.frame(
          VI_name = VI_name,
          lab = lab, 
          R2 = var_result$varresult[[length(var_result$varresult)]]$r.squared,
          AdjR2 = var_result$varresult[[length(var_result$varresult)]]$adj.r.squared,
          p_value = var_p_value,
          window_start = current_window$start,
          window_end = current_window$end,
          group_i = group_i,
          grid_id = i
        )
        var_R2.all <- rbind(var_R2.all, var_R2)
        result_list[[5]] <- var_R2
        
        # (6) Model stability test -----
        var_root_df <- data.frame(
          VI_name = VI_name,
          lab = lab, 
          roots = var_result$roots,
          window_start = current_window$start,
          window_end = current_window$end,
          group_i = group_i,
          grid_id = i
        )
        var_root_df$VARstability <- abs(max(var_root_df$roots)) < 1
        var_root_all <- rbind(var_root_all, var_root_df)
        result_list[[6]] <- var_root_df
        
        # (7) Impulse Response Analysis -----
        irf_result <- irf(var_model, response = VI_name, n.ahead = as.integer(n_ahead))
        
        # IRF estimated values
        irf_est <- as.data.frame(irf_result$irf)
        colnames(irf_est) <- names(irf_result$irf)
        irf_est$Time <- as.factor(c(1:(n_ahead+1)))
        irf_est$VI <- as.factor(VI_name)
        irf_est.m <- melt(irf_est, value.name = 'irf')
        irf_est.m$irf <- round(irf_est.m$irf, 2)
        
        # IRF upper limit
        irf_Upper <- as.data.frame(irf_result$Upper)
        colnames(irf_Upper) <- names(irf_result$Upper)
        irf_Upper$Time <- as.factor(c(1:(n_ahead+1)))
        irf_Upper.m <- melt(irf_Upper,value.name = 'Upper')
        
        # IRF lower limit
        irf_Lower <- as.data.frame(irf_result$Lower)
        colnames(irf_Lower) <- names(irf_result$Lower)
        irf_Lower$Time <- as.factor(c(1:(n_ahead+1)))
        irf_Lower.m <- melt(irf_Lower,value.name = 'Lower')
        
        # Merge IRF results
        irf_est_low_up <- data.frame(
          irf_est.m, 
          Upper = irf_Upper.m$Upper,
          Lower = irf_Lower.m$Lower,
          lab = lab,
          window_start = current_window$start,
          window_end = current_window$end,
          group_i = group_i,
          grid_id = i
        )
        irf_all <- rbind(irf_all, irf_est_low_up)
        result_list[[7]] <- irf_est_low_up
        
        # Extract Impulse Response Curve [IRF] of GPP to itself
        if(VI_name == "GPP" && exists("irf_est_low_up") && nrow(irf_est_low_up) > 0){
          gpp_irf <- irf_est_low_up %>% 
            filter(variable == "GPP") %>%
            mutate(
              window_id = window_label
            )
          
          if(nrow(gpp_irf) > 0){
            window_gpp_irf <- rbind(window_gpp_irf, gpp_irf)
          }
        }
        
        # (8) Plot Impulse Response Analysis -----
        # p <- ggplot(data = irf_est_low_up) + 
        #   geom_line(aes(x = as.numeric(Time), y = irf)) +
        #   geom_line(aes(x = as.numeric(Time), y = Upper), col = 'red') + 
        #   geom_line(aes(x = as.numeric(Time), y = Lower), col = 'red') +
        #   geom_hline(yintercept = 0, linetype = "dashed", size = 0.5) +
        #   labs(x = 'Lag period', y = VI_name, 
        #        title = paste("Grid", lab, "Window", window_label)) + 
        #   theme_bw() +
        #   facet_grid(rows = vars(variable), scales = "free_y")
        # result_list[[8]] <- p
        
        # (9) Lag response order [duration] -----
        for (ncol in 1:(ncol(irf_est)-2)) {
          irf <- round(abs(irf_est[,ncol]),5)
          
          for(lagn in 1:length(irf)){
            if(VI_name %in% c("NDVI", "SIF", "LAI", "GPP")){
              q = max(irf)*0.1
            }
            
            lag_a <- abs(irf[(lagn+1)] - irf[lagn])
            lag_b <- abs(irf[(lagn+2)] - irf[(lagn+1)])
            
            if(is.na(lag_a <= q) || is.na(lag_b <= q)){
              q = max(irf)*0.1
              time_lag <- which(irf < q)[1]
              break()
            }
            
            if(lag_a <= q & lag_b <= q){
              time_lag <- lagn
              break()
            }
          }
          
          df.q <- data.frame(
            VI = VI_name,
            variable = colnames(irf_est)[ncol],
            lab = lab,
            q = q,
            window_start = current_window$start,
            window_end = current_window$end,
            group_i = group_i,
            grid_id = i
          )
          q_all <- rbind(q_all, df.q)
          
          if(is.na(time_lag)){
            time_lag_max <- NA
            time_lag_intensity <- NA
            time_lag_Std <- NA
            time_lag_Stderror <- NA
          }else{
            time_lag_max <- which(abs(irf_est[1:time_lag, ncol]) == 
                                    max(abs(irf_est[1:time_lag, ncol])))
            time_lag_intensity <- round(irf_est[time_lag_max, ncol], 5)
            time_lag_Std <- round(sd(irf_est[1:time_lag, ncol]), 5)
            time_lag_Stderror <- round(time_lag_Std/sqrt(length(irf_est[1:time_lag, ncol])), 5)
          }
          
          time_lag_c <- data.frame(
            VI = VI_name, 
            variable = colnames(irf_est)[ncol], 
            lab = lab,  
            duration = time_lag, 
            intensity = time_lag_intensity, 
            std = time_lag_Std,
            std_error = time_lag_Stderror,
            max_time = time_lag_max,
            window_start = current_window$start,
            window_end = current_window$end,
            group_i = group_i,
            grid_id = i
          )
          
          duration_all <- rbind(duration_all, time_lag_c)
        }
        result_list[[8]] <- duration_all
        
        # (10) Variance Decomposition Analysis -----
        fevd_result <- fevd(var_model, n.ahead = as.integer(n_ahead))
        fevd_VI <- as.data.frame(fevd_result[[length(fevd_result)]])
        fevd_VI$Time <- as.factor(c(1:n_ahead))
        fevd_VI$VI <- as.factor(VI_name)
        fevd_VI.m <- melt(fevd_VI, value.name = "fevd")
        fevd_VI.m$lab <- lab
        fevd_VI.m$window_start <- current_window$start
        fevd_VI.m$window_end <- current_window$end
        fevd_VI.m$group_i <- group_i
        fevd_VI.m$grid_id <- i
        fevd_all <- rbind(fevd_all, fevd_VI.m)
        result_list[[9]] <- fevd_VI.m
        
        # (11) Causality Analysis -----
        causal_result <- causality(var_model)
        causal_df <- data.frame(
          Granger_method = causal_result$Granger$method,
          Granger_p_value = causal_result$Granger$p.value,
          Instant_method = causal_result$Instant$method,
          Instant_p_value = causal_result$Instant$p.value,
          lab = lab,
          window_start = current_window$start,
          window_end = current_window$end,
          group_i = group_i,
          grid_id = i
        )
        causal_all <- rbind(causal_all, causal_df)
        result_list[[10]] <- causal_df
        
        # (12) Save VAR model
        result_list[[11]] <- var_model
        
        VAR_all_result_list[[g_i]] <- result_list
        
      }, error = function(e) {
        error_msg <- paste("Error in grid", i, ":", e$message)
        log_error(error_msg, error_file)
        VAR_all_result_list[[g_i]] <- paste("Error:", e$message)
      })
    }
    
    # Extract GPP response to itself [Duration]
    if(VI_name == "GPP" && exists("duration_all") && nrow(duration_all) > 0){
      gpp_self <- duration_all %>% 
        filter(variable == "GPP") %>%
        mutate(
          window_id = window_label
        )
      
      if(nrow(gpp_self) > 0){
        window_gpp_response <- rbind(window_gpp_response, gpp_self)
      }
    }
    
    # Return results
    return(list(
      VAR_all_result_list = VAR_all_result_list,
      VAR_data_list = VAR_data_list,
      var_R2.all = var_R2.all,
      lag_RS_all = lag_RS_all,
      duration_all = duration_all,
      causal_all = causal_all,
      q_all = q_all,
      ADF_all = ADF_all,
      johansen_all = johansen_all,
      var_root_all = var_root_all,
      irf_all = irf_all,
      fevd_all = fevd_all,
      window_gpp_response = window_gpp_response,
      window_gpp_irf = window_gpp_irf,
      group_i = group_i,
      window_label = window_label,
      success = TRUE
    ))
    
  }, error = function(e) {
    error_msg <- paste("Error in process_group for group", group_i, "window", window_label, ":", e$message)
    log_error(error_msg, error_file)
    return(list(success = FALSE, group_i = group_i, window_label = window_label, error = e$message))
  })
}

# Main loop ----
for(VI_i in c(10)){
  VI_name <- colnames(combined_list[[1]])[VI_i]
  
  cat("\n", paste0(rep("=", 60), collapse = ""), "\n")
  cat("Starting analysis for", VI_name, "\n")
  cat(paste0(rep("=", 60), collapse = ""), "\n")
  
  # Create indicator result directory
  VI_dir <- file.path(Datedir, VI_name)
  dir.create(VI_dir, showWarnings = FALSE)
  
  # Iterate through each sliding window
  for(w in 1:10){ # seq_along(windows)
    # Get available cores
    n_cores <- get_available_cores(max_cores)
    
    # Setup parallel cluster -----
    cl <- makeCluster(n_cores, outfile = file.path(VI_dir, "parallel_log.txt"))
    registerDoParallel(cl)
    
    # Export necessary variables and functions to worker nodes
    clusterExport(cl, c("process_group", "ADF_test", "log_error", "VI_i", 
                        "lab_groups", "combined_list", "df.GLSP_GS", "CLI_names", 
                        "n_ahead", "year_start", "VI_dir", "error_file"))
    
    # Calculate total tasks for progress display
    total_tasks <- length(windows) * length(lab_groups)
    completed_tasks <- 0
    
    current_window <- windows[[w]]
    window_label <- paste0(current_window$start, "_", current_window$end)
    
    cat("\nProcessing window:", window_label, "(", w, "/", length(windows), ")\n")
    
    # Check memory
    available_mem <- check_memory()
    if(available_mem < 10) {
      cat("Low memory detected (", available_mem, "GB). Clearing memory...\n")
      clear_memory()
    }
    
    # Parallel processing all groups
    cat("Starting parallel processing for", length(lab_groups), "groups...\n")
    
    # 1:length(lab_groups)
    parallel_results <- foreach(group_i = 1:length(lab_groups),  
                                .packages = c("vars", "tseries", "urca", "dplyr", 
                                              "tidyr", "reshape2", "ggplot2", 
                                              "missForest", "zoo"),
                                .export = c("process_group", "ADF_test", "log_error", "combined_list", 
                                            "df.GLSP_GS", "CLI_names", "n_ahead", 
                                            "year_start", "VI_dir", "error_file"),
                                .errorhandling = "pass") %dopar% {
                                  process_group(group_i, VI_i, current_window, lab_groups, 
                                                combined_list, df.GLSP_GS, CLI_names, 
                                                n_ahead, year_start, VI_dir, error_file)
                                }
    
    cat("Parallel processing completed for window:", window_label, "\n")
    
    # Process and save results for each group (segmented saving)
    successful_groups <- 0
    failed_groups <- 0
    
    for(result_idx in seq_along(parallel_results)){
      result <- parallel_results[[result_idx]]
      completed_tasks <- completed_tasks + 1
      
      # Update progress bar
      print_progress(completed_tasks, total_tasks, "Overall Progress")
      
      if(class(result) != "try-error" && !is.null(result) && result$success){
        group_i <- result$group_i
        window_label <- result$window_label
        
        # Create group directory
        group_dir <- file.path(VI_dir, paste0("Window_", window_label), paste0("Group_", group_i))
        dir.create(group_dir, recursive = TRUE, showWarnings = FALSE)
        
        # Save results for each group segmentally
        tryCatch({
          # Save list to .RData file
          a1 <- result$VAR_all_result_list
          save(a1, 
               file = file.path(group_dir, paste0(VI_name, "_G", group_i, "_window_", window_label, "_VARresult_list.RData")))
          
          a2 <- result$VAR_data_list
          save(a2, 
               file = file.path(group_dir, paste0(VI_name, "_G", group_i, "_window_", window_label, "_VAR_data_list.RData")))
          
          # Save data.frame to .csv file
          write.csv(result$var_R2.all, 
                    file.path(group_dir, paste0(VI_name, "_G", group_i, "_window_", window_label, "_var_R2.csv")), row.names = FALSE)
          write.csv(result$lag_RS_all,
                    file.path(group_dir, paste0(VI_name, "_G", group_i, "_window_", window_label, "_lag_num_RS.csv")), row.names = FALSE)
          write.csv(result$duration_all, 
                    file.path(group_dir, paste0(VI_name, "_G", group_i, "_window_", window_label, "_duration.csv")), row.names = FALSE)
          write.csv(result$causal_all, 
                    file.path(group_dir, paste0(VI_name, "_G", group_i, "_window_", window_label, "_causal.csv")), row.names = FALSE)
          write.csv(result$q_all, 
                    file.path(group_dir, paste0(VI_name, "_G", group_i, "_window_", window_label, "_q.csv")), row.names = FALSE)
          write.csv(result$ADF_all, 
                    file.path(group_dir, paste0(VI_name, "_G", group_i, "_window_", window_label, "_ADF.csv")), row.names = FALSE)
          write.csv(result$johansen_all, 
                    file.path(group_dir, paste0(VI_name, "_G", group_i, "_window_", window_label, "_johansen.csv")), row.names = FALSE)
          write.csv(result$var_root_all, 
                    file.path(group_dir, paste0(VI_name, "_G", group_i, "_window_", window_label, "_var_root.csv")), row.names = FALSE)
          write.csv(result$irf_all, 
                    file.path(group_dir, paste0(VI_name, "_G", group_i, "_window_", window_label, "_irf.csv")), row.names = FALSE)
          write.csv(result$fevd_all, 
                    file.path(group_dir, paste0(VI_name, "_G", group_i, "_window_", window_label, "_fevd.csv")), row.names = FALSE)
          
          # Save GPP specific results (if exist)
          if(nrow(result$window_gpp_response) > 0){
            write.csv(result$window_gpp_response, 
                      file.path(group_dir, paste0("GPP_self_response_G", group_i, "_window_", window_label, ".csv")), row.names = FALSE)
          }
          if(nrow(result$window_gpp_irf) > 0){
            write.csv(result$window_gpp_irf, 
                      file.path(group_dir, paste0("GPP_irf_response_G", group_i, "_window_", window_label, ".csv")), row.names = FALSE)
          }
          
          successful_groups <- successful_groups + 1
          
        }, error = function(e) {
          error_msg <- paste("Error saving results for group", group_i, "window", window_label, ":", e$message)
          log_error(error_msg, error_file)
          failed_groups <- failed_groups + 1
        })
        
      } else {
        failed_groups <- failed_groups + 1
        if(class(result) == "try-error") {
          error_msg <- paste("Try-error in group", result_idx, "window", window_label, ":", as.character(result))
          log_error(error_msg, error_file)
        } else if(!is.null(result) && !result$success) {
          error_msg <- paste("Processing failed for group", result$group_i, "window", result$window_label, ":", result$error)
          log_error(error_msg, error_file)
        }
      }
    }
    
    cat("Window", window_label, "completed:", successful_groups, "successful,", failed_groups, "failed\n")
    
    # Clear memory
    rm(parallel_results)
    gc()
    
    # Stop parallel cluster
    stopCluster(cl)
  }
  
  cat("\n", VI_name, "processing completed!\n")
}

# Program ended
end_time <- Sys.time()
run_time <- end_time - start_time
cat("\nProgram completed at:", as.character(end_time), "\n")
cat("Total runtime:", run_time, "\n")

# Final memory clearance
clear_memory()

cat("All processing completed successfully!\n")
cat("Results saved in:", Datedir, "\n")
cat("Error log saved in:", error_file, "\n")