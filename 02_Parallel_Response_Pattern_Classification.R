
setwd('/data/share/tangwenxi')
RUN_name <- 'Response-SlidingWindow-10yr-ByGroup-Integrated-Periodic'

# Create directories -----
dir.create('./Out-Table-and-Figure-VPD et al', showWarnings = FALSE)
Datedir_main <- file.path('./Out-Table-and-Figure-VPD et al', Sys.Date())
dir.create(Datedir_main, showWarnings = FALSE)
Datedir <- file.path(Datedir_main, RUN_name)
dir.create(Datedir, showWarnings = FALSE)

# Create figure output directory
figure_date_dir <- file.path('./Figure', Sys.Date())
dir.create('./Figure', showWarnings = FALSE)
dir.create(figure_date_dir, showWarnings = FALSE)

# Parallel settings
max_cores <- 40

# Load packages -----
library(parallel)
library(doParallel)
library(foreach)
library(dplyr)
library(reshape2)
library(tidyr)
library(ggplot2)
library(ggthemes)
library(patchwork)
library(data.table)
library(raster)
library(viridis)
library(RColorBrewer)
library(scales)
library(rnaturalearth)
library(rnaturalearthdata)
library(gganimate)
library(grid)
library(gridExtra)

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

get_available_cores <- function(max_cores = 40) {
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

# Natural sort function -----
natural_sort <- function(x) {
  nums <- as.numeric(gsub("\\D", "", x))
  x[order(nums)]
}

# Get all group info -----
get_group_info <- function(base_path, vi_name = "GPP") {
  window_dirs <- list.dirs(file.path(base_path, vi_name), recursive = FALSE)
  window_dirs <- window_dirs[grepl("Window_", basename(window_dirs))]
  
  if(length(window_dirs) > 0) {
    first_window <- window_dirs[1]
    group_dirs <- list.dirs(first_window, recursive = FALSE)
    group_dirs <- group_dirs[grepl("Group_", basename(group_dirs))]
    
    group_nums <- as.numeric(gsub("Group_", "", basename(group_dirs)))
    group_info <- data.frame(
      group_num = sort(group_nums),
      group_name = paste0("Group_", sort(group_nums))
    )
    return(group_info)
  }
  return(NULL)
}

# Define time period creation function -----
create_time_periods <- function(years, period_length = 5) {
  
  # Create time period segmentation
  # period_length: 5 means 5-year segments, 10 means 10-year segments
  
  min_year <- min(years, na.rm = TRUE)
  max_year <- max(years, na.rm = TRUE)
  
  # Create start years for periods
  start_years <- seq(from = min_year, to = max_year, by = period_length)
  
  # Ensure the last period contains all remaining years
  if(max(start_years) < max_year) {
    start_years <- c(start_years, max_year)
  }
  
  periods <- data.frame()
  for(i in 1:(length(start_years)-1)) {
    periods <- rbind(periods, data.frame(
      period_id = i,
      start_year = start_years[i],
      end_year = start_years[i+1] - 1,
      period_label = paste0(start_years[i], "-", start_years[i+1] - 1)
    ))
  }
  
  # Handle the last period
  if(nrow(periods) > 0 && max(periods$end_year) < max_year) {
    last_period <- data.frame(
      period_id = nrow(periods) + 1,
      start_year = max(periods$end_year) + 1,
      end_year = max_year,
      period_label = paste0(max(periods$end_year) + 1, "-", max_year)
    )
    periods <- rbind(periods, last_period)
  }
  
  return(periods)
}

# Assign year to period -----
assign_year_to_period <- function(year, periods) {
  for(i in 1:nrow(periods)) {
    if(year >= periods$start_year[i] && year <= periods$end_year[i]) {
      return(periods$period_id[i])
    }
  }
  return(NA)
}

# Define equivalence test function -----
perform_equivalence_test <- function(ts1, ts2, Delta) {
  d <- ts2 - ts1
  n <- length(d)
  mean_d <- mean(d)
  sd_d <- sd(d)
  se_d <- sd_d / sqrt(n)
  df <- n - 1
  t_critical <- qt(0.95, df)
  lower_ci <- mean_d - t_critical * se_d
  upper_ci <- mean_d + t_critical * se_d
  max_abs_d <- max(abs(d))
  max_d <- d[which(max_abs_d == abs(d))]
  
  t_lower <- (mean_d - (-Delta)) / se_d
  p1 <- pt(t_lower, df, lower.tail = FALSE)
  t_upper <- (mean_d - Delta) / se_d
  p2 <- pt(t_upper, df, lower.tail = TRUE)
  p_value <- max(p1, p2)
  
  significance <- p_value <= 0.05
  Equ_result <- ifelse(significance, "Equivalence", "Un-Equivalence")
  
  return(list(
    max_abs_d = max_abs_d,
    max_d = max_d,
    mean_d = mean_d,
    lower_ci = lower_ci,
    upper_ci = upper_ci,
    p_value = p_value,
    significance = significance,
    Equ_result = Equ_result
  ))
}

# Modified single grid processing function (added period analysis) -----
process_single_grid_periodic <- function(lab_i, df, vi_name, periods, period_length) {
  tryCatch({
    result_yearly <- data.frame()
    detailed_periods <- data.frame()
    type_change_df <- data.frame()
    
    df.irf <- df %>% filter(lab == lab_i)
    
    if(nrow(df.irf) == 0) {
      return(list(
        result_yearly = data.frame(
          lab = lab_i, Year = NA, max_d = NA, mean_d = NA, 
          p_value = NA, significance = NA, Equ_result = NA,
          type_name = NA, type_code = NA, window_pair = NA, 
          period_id = NA, period_label = NA, note = "no_data"
        ),
        detailed_periods = data.frame(
          lab = lab_i, period_start = NA, period_end = NA, 
          transition_year = NA, response_type = NA, max_d = NA, 
          mean_d = NA, p_value = NA, significance = NA,
          window_pair = NA, period_id = NA, period_label = NA, note = "no_data"
        ),
        type_change_df = data.frame(
          lab = lab_i, change_year = NA, from = NA, to = NA,
          change_id = 0, persistent_type = NA, period_id = NA, 
          period_label = NA, note = "no_data"
        ),
        success = TRUE
      ))
    }
    
    window_pairs <- unique(df.irf[, c("window_start", "window_end")])
    window_pairs <- window_pairs[order(window_pairs$window_start), ]
    years_common <- intersect(window_pairs$window_end, window_pairs$window_start)
    
    if (length(years_common) == 0) {
      result_yearly <- rbind(result_yearly, data.frame(
        lab = lab_i, Year = NA, max_d = NA, mean_d = NA,
        p_value = NA, significance = NA, Equ_result = NA,
        type_name = NA, type_code = NA, window_pair = NA,
        period_id = NA, period_label = NA, note = "no_matched_windows"
      ))
      
      detailed_periods <- rbind(detailed_periods, data.frame(
        lab = lab_i, period_start = NA, period_end = NA, 
        transition_year = NA, response_type = NA, max_d = NA, 
        mean_d = NA, p_value = NA, significance = NA,
        window_pair = NA, period_id = NA, period_label = NA, note = "no_matched_windows"
      ))
      
      type_change_df <- data.frame(
        lab = lab_i, change_year = NA, from = NA, to = NA,
        change_id = 0, persistent_type = NA, period_id = NA,
        period_label = NA, note = "no_matched_windows"
      )
    } else {
      for (y in years_common) {
        w1 <- df.irf %>% filter(window_end == y)
        w2 <- df.irf %>% filter(window_start == y)
        
        if (nrow(w1) > 0 & nrow(w2) > 0) {
          period1_start <- unique(w1$window_start)[1]
          period1_end <- unique(w1$window_end)[1]
          period2_start <- unique(w2$window_start)[1]
          period2_end <- unique(w2$window_end)[1]
          
          # Assign to period
          period_id <- assign_year_to_period(y, periods)
          period_label <- ifelse(is.na(period_id), NA, periods$period_label[periods$period_id == period_id])
          
          Delta <- ifelse(vi_name == "GPP", 0.2,
                          ifelse(vi_name == "SIF", 0.0006,
                                 ifelse(vi_name == "NDVI", 0.0006,
                                        ifelse(vi_name == "LAI", 0.004, 0.3))))
          
          test_result <- perform_equivalence_test(w1$irf, w2$irf, Delta)
          
          if (test_result$Equ_result == "Equivalence") {
            type_name <- "Stable"
          } else if (test_result$max_d > 0) {
            type_name <- "VGC"
          } else {
            type_name <- "VSO"
          }
          type_code <- ifelse(type_name == "VSO", 1,
                              ifelse(type_name == "VGC", 2, 3))
          
          window_pair_str <- paste(unique(w1$window_id), "&", unique(w2$window_id))
          
          result_yearly <- rbind(result_yearly, data.frame(
            lab = lab_i, Year = y, max_d = test_result$max_d,
            mean_d = test_result$mean_d, p_value = test_result$p_value,
            significance = test_result$significance, 
            Equ_result = test_result$Equ_result,
            type_name = type_name, type_code = type_code,
            window_pair = window_pair_str, period_id = period_id,
            period_label = period_label, note = NA
          ))
          
          detailed_periods <- rbind(detailed_periods, data.frame(
            lab = lab_i,
            period_start = paste0(period1_start, "-", period1_end),
            period_end = paste0(period2_start, "-", period2_end),
            transition_year = y,
            response_type = type_name,
            max_d = test_result$max_d,
            mean_d = test_result$mean_d,
            p_value = test_result$p_value,
            significance = test_result$significance,
            window_pair = window_pair_str,
            period_id = period_id,
            period_label = period_label,
            note = NA
          ))
        }
      }
      
      # Analyze changes by period - This is a key modification
      if(nrow(result_yearly) > 1) {
        result_yearly <- result_yearly[order(result_yearly$Year), ]
        
        # Statistics for each period separately
        type_change_df <- data.frame()
        
        for(pid in unique(periods$period_id)) {
          period_data <- result_yearly %>% 
            filter(period_id == pid) %>%
            arrange(Year)
          
          if(nrow(period_data) > 1) {
            # Find type changes within the period
            type_change_idx <- which(period_data$type_name[-1] != period_data$type_name[-nrow(period_data)])
            
            if(length(type_change_idx) > 0) {
              # Changes occurred in this period
              for(idx in type_change_idx) {
                change_record <- data.frame(
                  lab = period_data$lab[idx + 1],
                  change_year = period_data$Year[idx + 1],
                  from = period_data$type_name[idx],
                  to = period_data$type_name[idx + 1],
                  change_id = length(type_change_idx),  # Number of changes in this period
                  persistent_type = NA,
                  period_id = pid,
                  period_label = periods$period_label[periods$period_id == pid],
                  note = NA
                )
                type_change_df <- rbind(type_change_df, change_record)
              }
            } else {
              # No change in this period
              persistent_state <- unique(period_data$type_name)[1]
              no_change_record <- data.frame(
                lab = unique(period_data$lab),
                change_year = NA,
                from = NA,
                to = NA,
                change_id = 0,
                persistent_type = persistent_state,
                period_id = pid,
                period_label = periods$period_label[periods$period_id == pid],
                note = "no_changed_point"
              )
              type_change_df <- rbind(type_change_df, no_change_record)
            }
          } else if(nrow(period_data) == 1) {
            # Only one data point in this period
            persistent_state <- period_data$type_name[1]
            single_data_record <- data.frame(
              lab = lab_i,
              change_year = NA,
              from = NA,
              to = NA,
              change_id = 0,
              persistent_type = persistent_state,
              period_id = pid,
              period_label = periods$period_label[periods$period_id == pid],
              note = "single_data_point"
            )
            type_change_df <- rbind(type_change_df, single_data_record)
          }
        }
        
        # If no data in any period, add default record
        if(nrow(type_change_df) == 0) {
          type_change_df <- data.frame(
            lab = lab_i, change_year = NA, from = NA, to = NA,
            change_id = 0, persistent_type = NA, period_id = NA,
            period_label = NA, note = "insufficient_data"
          )
        }
      } else if(nrow(result_yearly) == 1) {
        persistent_state <- result_yearly$type_name[1]
        pid <- result_yearly$period_id[1]
        type_change_df <- data.frame(
          lab = lab_i, change_year = NA, from = NA, to = NA,
          change_id = 0, persistent_type = persistent_state,
          period_id = pid,
          period_label = ifelse(is.na(pid), NA, periods$period_label[periods$period_id == pid]),
          note = "single_data_point"
        )
      } else {
        type_change_df <- data.frame(
          lab = lab_i, change_year = NA, from = NA, to = NA,
          change_id = 0, persistent_type = NA, period_id = NA,
          period_label = NA, note = "insufficient_data"
        )
      }
    }
    
    return(list(
      result_yearly = result_yearly,
      detailed_periods = detailed_periods,
      type_change_df = type_change_df,
      success = TRUE
    ))
    
  }, error = function(e) {
    return(list(
      result_yearly = data.frame(
        lab = lab_i, Year = NA, max_d = NA, mean_d = NA,
        p_value = NA, significance = NA, Equ_result = NA,
        type_name = NA, type_code = NA, window_pair = NA,
        period_id = NA, period_label = NA, note = paste("error:", e$message)
      ),
      detailed_periods = data.frame(
        lab = lab_i, period_start = NA, period_end = NA, 
        transition_year = NA, response_type = NA, max_d = NA, 
        mean_d = NA, p_value = NA, significance = NA,
        window_pair = NA, period_id = NA, period_label = NA, note = paste("error:", e$message)
      ),
      type_change_df = data.frame(
        lab = lab_i, change_year = NA, from = NA, to = NA,
        change_id = 0, persistent_type = NA, period_id = NA,
        period_label = NA, note = paste("error:", e$message)
      ),
      success = FALSE,
      error = e$message
    ))
  })
}

# Modified group processing function -----
process_group_data_periodic <- function(group_num, group_name, base_path, vi_name, output_dir, periods, period_length) {
  tryCatch({
    cat("\n=== Processing", group_name, "===\n")
    
    cat("Step 1: Combining IRF data for", group_name, "\n")
    
    window_dirs <- list.dirs(file.path(base_path, vi_name), recursive = FALSE)
    window_dirs <- window_dirs[grepl("Window_", basename(window_dirs))]
    window_dirs <- window_dirs[order(basename(window_dirs))]
    
    group_data <- NULL
    
    for(window_dir in window_dirs) {
      window_name <- basename(window_dir)
      group_dir <- file.path(window_dir, group_name)
      
      if(dir.exists(group_dir)) {
        irf_files <- list.files(group_dir, pattern = "*_irf.csv$", full.names = TRUE)
        
        for(irf_file in irf_files) {
          tryCatch({
            irf_data <- fread(irf_file)
            vi_data <- irf_data[irf_data$variable == vi_name, ]
            
            if(nrow(vi_data) > 0) {
              vi_data$window_id <- gsub("Window_", "", window_name)
              group_data <- rbind(group_data, vi_data)
            }
          }, error = function(e) {
            cat("Warning: Error reading", irf_file, ":", e$message, "\n")
          })
        }
      }
    }
    
    if(is.null(group_data) || nrow(group_data) == 0) {
      cat("No data found for", group_name, "\n")
      return(list(
        success = FALSE,
        group_name = group_name,
        note = "no_data_found"
      ))
    }
    
    cat("Step 2: Saving combined data for", group_name, "\n")
    group_output_dir <- file.path(output_dir, group_name)
    dir.create(group_output_dir, showWarnings = FALSE, recursive = TRUE)
    
    combined_file <- file.path(group_output_dir, paste0(group_name, "_", vi_name, "_combined_irf.csv"))
    fwrite(group_data, combined_file)
    cat("Combined data saved:", combined_file, "\n")
    
    cat("Step 3: Analyzing response types for", group_name, "with", period_length, "year periods\n")
    
    lab_list <- unique(group_data$lab)
    cat("Found", length(lab_list), "grids in", group_name, "\n")
    
    if(length(lab_list) == 0) {
      return(list(
        success = FALSE,
        group_name = group_name,
        note = "no_valid_labs"
      ))
    }
    
    group_yearly_results <- data.frame()
    group_detailed_periods <- data.frame()
    group_change_results <- data.frame()
    
    for(lab_i in lab_list) {
      result <- process_single_grid_periodic(lab_i, group_data, vi_name, periods, period_length)
      
      if(result$success) {
        group_yearly_results <- rbind(group_yearly_results, result$result_yearly)
        group_detailed_periods <- rbind(group_detailed_periods, result$detailed_periods)
        group_change_results <- rbind(group_change_results, result$type_change_df)
      }
    }
    
    if(nrow(group_detailed_periods) > 0) {
      group_detailed_periods$group_name <- group_name
      group_detailed_periods$group_num <- group_num
    }
    
    if(nrow(group_yearly_results) > 0) {
      group_yearly_results$group_name <- group_name
      group_yearly_results$group_num <- group_num
    }
    
    if(nrow(group_change_results) > 0) {
      group_change_results$group_name <- group_name
      group_change_results$group_num <- group_num
    }
    
    cat("Step 4: Saving analysis results for", group_name, "\n")
    
    yearly_file <- file.path(group_output_dir, paste0(group_name, "_response_result_yearly_", period_length, "yr.csv"))
    detailed_file <- file.path(group_output_dir, paste0(group_name, "_detailed_periods_", period_length, "yr.csv"))
    change_file <- file.path(group_output_dir, paste0(group_name, "_response_type_change_", period_length, "yr.csv"))
    
    write.csv(group_yearly_results, yearly_file, row.names = FALSE)
    write.csv(group_detailed_periods, detailed_file, row.names = FALSE)
    write.csv(group_change_results, change_file, row.names = FALSE)
    
    cat("Analysis results saved for", group_name, "\n")
    
    return(list(
      success = TRUE,
      group_name = group_name,
      group_num = group_num,
      n_grids = length(lab_list),
      n_yearly_records = nrow(group_yearly_results),
      n_detailed_records = nrow(group_detailed_periods),
      n_change_records = nrow(group_change_results),
      yearly_file = yearly_file,
      detailed_file = detailed_file,
      change_file = change_file,
      combined_file = combined_file,
      yearly_data = group_yearly_results,
      detailed_periods_data = group_detailed_periods,
      change_data = group_change_results
    ))
    
  }, error = function(e) {
    cat("Error processing", group_name, ":", e$message, "\n")
    return(list(
      success = FALSE,
      group_name = group_name,
      error = e$message
    ))
  })
}

# Function to combine results from all groups -----
combine_all_group_results_periodic <- function(processing_results, output_dir, vi_name, period_length) {
  cat("\n=== Combining All Group Results (Yearly & Change) ===\n")
  
  all_yearly_results <- data.frame()
  all_change_results <- data.frame()
  all_detailed_periods <- data.frame()
  successful_groups <- 0
  
  for(result in processing_results) {
    if(class(result) != "try-error" && !is.null(result) && result$success) {
      if(!is.null(result$yearly_data) && nrow(result$yearly_data) > 0) {
        all_yearly_results <- rbind(all_yearly_results, result$yearly_data)
      }
      
      if(!is.null(result$change_data) && nrow(result$change_data) > 0) {
        all_change_results <- rbind(all_change_results, result$change_data)
      }
      
      if(!is.null(result$detailed_periods_data) && nrow(result$detailed_periods_data) > 0) {
        all_detailed_periods <- rbind(all_detailed_periods, result$detailed_periods_data)
      }
      
      successful_groups <- successful_groups + 1
      cat("✓ Combined results from", result$group_name, "\n")
    }
  }
  
  # Save results
  if(nrow(all_yearly_results) > 0) {
    yearly_output <- file.path(output_dir, paste0("Combined_All_Groups_response_result_yearly_", vi_name, "_", period_length, "yr.csv"))
    write.csv(all_yearly_results, yearly_output, row.names = FALSE)
    cat("Yearly results saved to:", yearly_output, "\n")
  }
  
  if(nrow(all_change_results) > 0) {
    change_output <- file.path(output_dir, paste0("Combined_All_Groups_response_type_change_", vi_name, "_", period_length, "yr.csv"))
    write.csv(all_change_results, change_output, row.names = FALSE)
    cat("Change results saved to:", change_output, "\n")
  }
  
  if(nrow(all_detailed_periods) > 0) {
    detailed_output <- file.path(output_dir, paste0("Combined_All_Groups_detailed_periods_", vi_name, "_", period_length, "yr.csv"))
    write.csv(all_detailed_periods, detailed_output, row.names = FALSE)
    cat("Detailed periods saved to:", detailed_output, "\n")
  }
  
  return(list(
    yearly_output_file = if(nrow(all_yearly_results) > 0) yearly_output else NULL,
    change_output_file = if(nrow(all_change_results) > 0) change_output else NULL,
    detailed_output_file = if(nrow(all_detailed_periods) > 0) detailed_output else NULL,
    yearly_data = all_yearly_results,
    change_data = all_change_results,
    detailed_data = all_detailed_periods,
    yearly_records = nrow(all_yearly_results),
    change_records = nrow(all_change_results),
    detailed_records = nrow(all_detailed_periods),
    successful_groups = successful_groups
  ))
}

# ===== Plotting Functions Section =====

# Raster processing function
saveToRsater <- function(labs, values){
  r_test <- raster("/data/share/tangwenxi/Data/resample_example_0.5degree.tif")
  values(r_test) <- NA
  values(r_test)[labs] <- values
  return(r_test)
}

cutByLULC <- function(raster){
  LULC <- raster("/data/share/tangwenxi/Data/LULC_1to15_2000.tif")
  LULC_lab <- which(values(LULC) %in% c(1:15))
  r_LULC <- saveToRsater(labs = LULC_lab, 
                         values = values(raster)[LULC_lab])
  return(r_LULC)
}

# Get world map
world_shp <- ne_coastline()

# Data preprocessing: filter out no_matched_windows data
filter_valid_data <- function(df) {
  df %>% filter(note != "no_matched_windows" | is.na(note))
}

# Persistent Type Visualization
plot_persistent_type_comprehensive_periodic <- function(df_change, period_label, output_dir) {
  persistent_data <- df_change %>%
    filter_valid_data() %>%
    filter(change_id == 0 & 
             (note == "no_changed_point" | note == "single_data_point") &
             !is.na(persistent_type))
  
  if(nrow(persistent_data) == 0) {
    cat("Warning: No persistent type data found for period", period_label, "\n")
    return(NULL)
  }
  
  cat("Found", nrow(persistent_data), "grids with persistent types in period", period_label, "\n")
  
  persistent_data$type_code <- as.numeric(as.factor(persistent_data$persistent_type))
  r_persistent <- saveToRsater(persistent_data$lab, persistent_data$type_code)
  r_persistent <- cutByLULC(r_persistent)
  
  df_map <- as.data.frame(r_persistent, xy = TRUE)
  names(df_map) <- c("x", "y", "type_code")
  df_map <- df_map[!is.na(df_map$type_code), ]
  
  type_labels <- persistent_data %>% 
    dplyr::select(type_code, persistent_type) %>% 
    distinct()
  df_map <- left_join(df_map, type_labels, by = "type_code")
  
  persistent_colors <- c(
    "Stable" = "#2E8B57",
    "VGC" = "#FF6347",
    "VSO" = "#4169E1"
  )
  
  stats_data <- persistent_data %>%
    count(persistent_type) %>%
    mutate(
      percentage = n / sum(n) * 100,
      area_km2 = n * 0.25
    ) %>%
    arrange(desc(n))
  
  count_plot <- ggplot(stats_data, aes(x = reorder(persistent_type, n), y = n, fill = persistent_type)) +
    geom_col(alpha = 0.8, color = "white", linewidth = 0.5) +
    geom_text(aes(label = scales::comma(n)), hjust = -0.1, size = 3, fontface = "bold") +
    scale_fill_manual(values = persistent_colors) +
    coord_flip() +
    labs(title = paste("Grid Count by Persistent Type -", period_label),
         x = "Persistent Type", y = "Number of Grids") +
    theme_minimal() +
    theme(
      legend.position = "none",
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),
      plot.title = element_text(face = "bold", size = 10),
      plot.background = element_rect(fill = "white", color = "black", linewidth = 0.5)
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15)))
  
  area_plot <- ggplot(stats_data, aes(x = reorder(persistent_type, percentage), y = percentage, fill = persistent_type)) +
    geom_col(alpha = 0.8, color = "white", linewidth = 0.5) +
    geom_text(aes(label = sprintf("%.1f%%", percentage)), hjust = -0.1, size = 3, fontface = "bold") +
    scale_fill_manual(values = persistent_colors) +
    coord_flip() +
    labs(title = paste("Area Percentage by Persistent Type -", period_label),
         x = "Persistent Type", y = "Percentage (%)") +
    theme_minimal() +
    theme(
      legend.position = "none",
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),
      plot.title = element_text(face = "bold", size = 10),
      plot.background = element_rect(fill = "white", color = "black", linewidth = 0.5)
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15)))
  
  combined_stats <- count_plot / area_plot
  combined_stats_grob <- ggplotGrob(combined_stats)
  
  map_plot <- ggplot(df_map) +
    geom_hline(yintercept = 0, color = "gray50", linetype = "dashed", linewidth = 0.5) +
    geom_tile(aes(x = x, y = y, fill = persistent_type)) +
    geom_sf(data = world_shp, fill = NA, color = "gray30", linewidth = 0.2) +
    coord_sf(ylim = c(-60, 90), expand = FALSE) +
    scale_fill_manual(name = "Persistent\nType", values = persistent_colors) +
    labs(title = paste("Persistent Vegetation Response Types -", period_label),
         subtitle = "Areas with no significant changes (stable vegetation response patterns)",
         x = "Longitude", y = "Latitude") +
    theme_bw() +
    theme(
      plot.title = element_text(face = "bold", size = 12),
      plot.subtitle = element_text(size = 10),
      legend.position = "bottom",
      panel.background = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank()
    ) +
    annotation_custom(combined_stats_grob, xmin = -180, xmax = -85, ymin = -60, ymax = 20)
  
  # Save figure
  filename <- paste0("persistent_type_comprehensive_", gsub("-", "_", period_label), ".png")
  filepath <- file.path(output_dir, filename)
  ggsave(filepath, map_plot, width = 14, height = 8, dpi = 300)
  cat("Persistent type plot saved:", filepath, "\n")
  
  return(list(
    map = map_plot,
    stats = stats_data,
    count_plot = count_plot,
    area_plot = area_plot,
    filepath = filepath
  ))
}

# Change Frequency Map
plot_change_frequency_map_periodic <- function(df_change, period_label, output_dir) {
  df_filtered <- filter_valid_data(df_change)
  
  all_labs <- unique(df_filtered$lab[!is.na(df_filtered$lab)])
  
  change_freq <- df_filtered %>%
    filter(!is.na(change_id)) %>%
    group_by(lab) %>%
    summarise(max_changes = max(change_id, na.rm = TRUE)) %>%
    ungroup()
  
  complete_freq <- data.frame(lab = all_labs) %>%
    left_join(change_freq, by = "lab") %>%
    mutate(max_changes = ifelse(is.na(max_changes), 0, max_changes))
  
  complete_freq <- complete_freq %>%
    mutate(change_category = case_when(
      max_changes == 0 ~ "No Change",
      max_changes == 1 ~ "1 Change",
      max_changes == 2 ~ "2 Changes", 
      max_changes == 3 ~ "3 Changes",
      max_changes == 4 ~ "4 Changes",
      max_changes == 5 ~ "5 Changes",
      max_changes >= 6 ~ "6+ Changes"
    ))
  
  r_freq <- saveToRsater(complete_freq$lab, as.numeric(as.factor(complete_freq$change_category)))
  r_freq <- cutByLULC(r_freq)
  
  df_map <- as.data.frame(r_freq, xy = TRUE)
  names(df_map) <- c("x", "y", "category_code")
  df_map <- df_map[!is.na(df_map$category_code), ]
  
  category_labels <- complete_freq %>% 
    dplyr::select(change_category) %>% 
    distinct() %>%
    mutate(category_code = as.numeric(as.factor(change_category)))
  
  df_map <- left_join(df_map, category_labels, by = "category_code")
  
  fresh_colors <- c("No Change" = "#ed1c24",        
                    "1 Change" = "#00c7f2",         
                    "2 Changes" = "#6aaed6",        
                    "3 Changes" = "#b870ab",         
                    "4 Changes" = "#fb6ca0",        
                    "5 Changes" = "#ffd700",         
                    "6+ Changes" = "#ffff4d")
  
  area_stats <- df_map %>%
    count(change_category) %>%
    mutate(percentage = n / sum(n) * 100) %>%
    mutate(change_category = factor(change_category, 
                                    levels = c("No Change", "1 Change", "2 Changes", 
                                               "3 Changes", "4 Changes", "5 Changes", "6+ Changes")))
  
  bar_plot <- ggplot(area_stats, aes(x = change_category, y = percentage, fill = change_category)) +
    geom_col(alpha = 0.8) +
    scale_fill_manual(values = fresh_colors) +
    labs(x = NULL, y = "Area %") +
    theme_minimal() +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      axis.text.y = element_text(size = 8),
      axis.title.y = element_text(size = 9),
      panel.grid.minor = element_blank(),
      plot.background = element_rect(fill = "white", color = "black", linewidth = 0.5)
    )
  
  bar_grob <- ggplotGrob(bar_plot)
  
  p1 <- ggplot(df_map) +
    geom_hline(yintercept = 0, color = "gray50", linetype = "dashed", linewidth = 0.5) +
    geom_tile(aes(x = x, y = y, fill = change_category)) +
    geom_sf(data = world_shp, fill = NA, color = "gray30", linewidth = 0.2) +
    coord_sf(ylim = c(-60, 90), expand = FALSE) +
    scale_fill_manual(name = "Change\nFrequency", values = fresh_colors) +
    labs(title = paste("Vegetation Response Change Frequency -", period_label),
         x = "Longitude", y = "Latitude") +
    theme_bw() +
    theme(
      plot.title = element_text(face = "bold", size = 12),
      legend.position = "bottom",
      panel.background = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank()
    ) +
    annotation_custom(bar_grob, xmin = -180, xmax = -85, ymin = -60, ymax = 10)
  
  # Save figure
  filename <- paste0("vegetation_change_frequency_map_", gsub("-", "_", period_label), ".png")
  filepath <- file.path(output_dir, filename)
  ggsave(filepath, p1, width = 12, height = 8, dpi = 300)
  cat("Change frequency plot saved:", filepath, "\n")
  
  return(list(plot = p1, filepath = filepath))
}

# First Change Year Map (Continuous Version)
plot_first_change_year_map1_periodic <- function(df_change, period_label, output_dir) {
  df_filtered <- filter_valid_data(df_change)
  
  first_change <- df_filtered %>%
    filter(!is.na(change_year)) %>%
    group_by(lab) %>%
    summarise(first_year = min(change_year, na.rm = TRUE)) %>%
    ungroup()
  
  if(nrow(first_change) == 0) {
    cat("Warning: No change year data found for period", period_label, "\n")
    return(NULL)
  }
  
  r_year <- saveToRsater(first_change$lab, first_change$first_year)
  r_year <- cutByLULC(r_year)
  
  df_map <- as.data.frame(r_year, xy = TRUE)
  names(df_map) <- c("x", "y", "year")
  df_map <- df_map[!is.na(df_map$year), ]
  
  p2 <- ggplot(df_map) +
    geom_hline(yintercept = 0, color = "gray50", linetype = "dashed", linewidth = 0.5) +
    geom_tile(aes(x = x, y = y, fill = year)) +
    geom_sf(data = world_shp, fill = NA, color = "gray30", linewidth = 0.2) +
    coord_sf(ylim = c(-60, 90), expand = FALSE) +
    scale_fill_gradientn(name = "First Change\nYear", 
                         colors = c("#00c7f2", "#fb6ca0", "#ffff4d"),
                         values = scales::rescale(c(0, 0.5, 1)),
                         breaks = scales::pretty_breaks(n = 5)) +
    labs(title = paste("First Change Year in Vegetation Response -", period_label),
         x = "Longitude", y = "Latitude") +
    theme_bw() +
    theme(
      plot.title = element_text(face = "bold", size = 12),
      legend.position = "bottom",
      panel.background = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank()
    )
  
  # Save figure
  filename <- paste0("first_change_year_map1_", gsub("-", "_", period_label), ".png")
  filepath <- file.path(output_dir, filename)
  ggsave(filepath, p2, width = 12, height = 8, dpi = 300)
  cat("First change year plot (continuous) saved:", filepath, "\n")
  
  return(list(plot = p2, filepath = filepath))
}

# First Change Year Map (Classified Version)
plot_first_change_year_map2_periodic <- function(df_change, period_label, output_dir) {
  df_filtered <- filter_valid_data(df_change)
  
  first_change <- df_filtered %>%
    filter(!is.na(change_year)) %>%
    group_by(lab) %>%
    summarise(first_year = min(change_year, na.rm = TRUE)) %>%
    ungroup()
  
  if(nrow(first_change) == 0) {
    cat("Warning: No change year data found for period", period_label, "\n")
    return(NULL)
  }
  
  year_breaks <- quantile(first_change$first_year, probs = c(0, 0.2, 0.4, 0.6, 0.8, 1), na.rm = TRUE)
  
  first_change <- first_change %>%
    mutate(year_period = case_when(
      first_year <= year_breaks[2] ~ paste0("Early (", round(year_breaks[1]), "-", round(year_breaks[2]), ")"),
      first_year <= year_breaks[3] ~ paste0("Early-Mid (", round(year_breaks[2])+1, "-", round(year_breaks[3]), ")"),
      first_year <= year_breaks[4] ~ paste0("Mid (", round(year_breaks[3])+1, "-", round(year_breaks[4]), ")"),
      first_year <= year_breaks[5] ~ paste0("Mid-Late (", round(year_breaks[4])+1, "-", round(year_breaks[5]), ")"),
      TRUE ~ paste0("Late (", round(year_breaks[5])+1, "-", round(year_breaks[6]), ")")
    ))
  
  first_change$period_code <- as.numeric(as.factor(first_change$year_period))
  r_year <- saveToRsater(first_change$lab, first_change$period_code)
  r_year <- cutByLULC(r_year)
  
  df_map <- as.data.frame(r_year, xy = TRUE)
  names(df_map) <- c("x", "y", "period_code")
  df_map <- df_map[!is.na(df_map$period_code), ]
  
  period_labels <- first_change %>% 
    dplyr::select(period_code, year_period) %>% 
    distinct()
  df_map <- left_join(df_map, period_labels, by = "period_code")
  
  period_colors <- c(
    "#D7191C",  
    "#FDAE61",    
    "#FFFFBF",  
    "#ABD9E9",  
    "#2C7BB6"   
  )
  
  period_order <- sort(unique(df_map$year_period))
  names(period_colors) <- period_order
  
  period_stats <- df_map %>%
    count(year_period) %>%
    mutate(percentage = n / sum(n) * 100) %>%
    arrange(year_period)
  
  bar_plot <- ggplot(period_stats, aes(x = reorder(year_period, percentage), y = percentage)) +
    geom_col(aes(fill = year_period), alpha = 0.8) +
    scale_fill_manual(values = period_colors) +
    labs(x = NULL, y = "Area %") +
    coord_flip() +
    theme_minimal() +
    theme(
      legend.position = "none",
      axis.text.x = element_text(size = 7),
      axis.text.y = element_text(size = 6),
      axis.title.x = element_text(size = 8),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),
      plot.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
      plot.margin = margin(5, 5, 5, 5)
    )
  
  bar_grob <- ggplotGrob(bar_plot)
  
  p2 <- ggplot(df_map) +
    geom_hline(yintercept = 0, color = "gray50", linetype = "dashed", linewidth = 0.5) +
    geom_tile(aes(x = x, y = y, fill = year_period)) +
    geom_sf(data = world_shp, fill = NA, color = "gray30", linewidth = 0.2) +
    coord_sf(ylim = c(-60, 90), expand = FALSE) +
    scale_fill_manual(name = "First Change\nPeriod", 
                      values = period_colors,
                      guide = guide_legend(
                        ncol = 1,
                        byrow = TRUE
                      )) +
    labs(title = paste("First Change Year in Vegetation Response -", period_label),
         x = "Longitude", y = "Latitude") +
    theme_bw() +
    theme(
      plot.title = element_text(face = "bold", size = 12),
      legend.position = "bottom",
      panel.background = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank()
    ) +
    annotation_custom(bar_grob, xmin = -180, xmax = -85, ymin = -60, ymax = 10)
  
  # Save figure
  filename <- paste0("first_change_year_map2_", gsub("-", "_", period_label), ".png")
  filepath <- file.path(output_dir, filename)
  ggsave(filepath, p2, width = 12, height = 8, dpi = 300)
  cat("First change year plot (classified) saved:", filepath, "\n")
  
  return(list(plot = p2, filepath = filepath))
}

# Change Type Spatial Distribution Map
plot_change_type_map_periodic <- function(df_change, period_label, output_dir) {
  df_filtered <- filter_valid_data(df_change)
  
  change_type <- df_filtered %>%
    filter(!is.na(from) & !is.na(to) & from != to) %>%
    mutate(change_type = paste(from, "→", to)) %>%
    group_by(lab) %>%
    slice(1) %>%
    ungroup()
  
  if(nrow(change_type) == 0) {
    cat("Warning: No change type data found for period", period_label, "\n")
    return(NULL)
  }
  
  change_type$type_code <- as.numeric(as.factor(change_type$change_type))
  r_type <- saveToRsater(change_type$lab, change_type$type_code)
  r_type <- cutByLULC(r_type)
  
  df_map <- as.data.frame(r_type, xy = TRUE)
  names(df_map) <- c("x", "y", "type_code")
  df_map <- df_map[!is.na(df_map$type_code), ]
  
  type_labels <- change_type %>% 
    dplyr::select(type_code, change_type) %>% 
    distinct()
  df_map <- left_join(df_map, type_labels, by = "type_code")
  
  area_stats <- df_map %>%
    filter(!is.na(change_type)) %>%
    count(change_type) %>%
    mutate(percentage = n / sum(n) * 100) %>%
    arrange(desc(percentage))
  
  all_change_types <- c("Stable → VGC", "Stable → VSO", "VGC → Stable", 
                        "VGC → VSO", "VSO → Stable", "VSO → VGC")
  
  area_stats <- area_stats %>%
    right_join(data.frame(change_type = all_change_types), by = "change_type") %>%
    mutate(
      n = ifelse(is.na(n), 0, n),
      percentage = ifelse(is.na(percentage), 0, percentage)
    ) %>%
    filter(percentage > 0) %>%
    arrange(desc(percentage))
  
  bar_plot <- ggplot(area_stats, aes(x = reorder(change_type, percentage), y = percentage, fill = change_type)) +
    geom_col(alpha = 0.8) +
    scale_fill_manual(
      values = c(
        "Stable → VGC" = "#2E8B57",   
        "Stable → VSO" = "#48D1CC",   
        "VGC → Stable" = "#DAA520",   
        "VGC → VSO" = "#CD853F",        
        "VSO → Stable" = "#4682B4",   
        "VSO → VGC" = "#9370DB")
    ) +
    labs(x = NULL, y = "Area %") +
    coord_flip() +
    theme_minimal() +
    theme(
      legend.position = "none",
      axis.text.x = element_text(size = 7),
      axis.text.y = element_text(size = 7),
      axis.title.x = element_text(size = 8),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),
      plot.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
      plot.margin = margin(5, 5, 5, 5)
    )
  
  bar_grob <- ggplotGrob(bar_plot)
  
  p3 <- ggplot(df_map) +
    geom_hline(yintercept = 0, color = "gray50", linetype = "dashed", linewidth = 0.5) +
    geom_tile(aes(x = x, y = y, fill = change_type)) +
    geom_sf(data = world_shp, fill = NA, color = "gray30", linewidth = 0.2) +
    coord_sf(ylim = c(-60, 90), expand = FALSE) +
    scale_fill_manual(
      name = "Change Type", 
      values = c(
        "Stable → VGC" = "#2E8B57",   
        "Stable → VSO" = "#48D1CC",   
        "VGC → Stable" = "#DAA520",   
        "VGC → VSO" = "#CD853F",        
        "VSO → Stable" = "#4682B4",   
        "VSO → VGC" = "#9370DB"        
      ),
      na.value = "transparent",
      guide = guide_legend(
        ncol = 3,
        byrow = TRUE
      )
    ) +
    labs(title = paste("Vegetation Response Change Types -", period_label),
         x = "Longitude", y = "Latitude") +
    theme_bw() +
    theme(
      plot.title = element_text(face = "bold", size = 12),
      legend.position = "bottom",
      panel.background = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank()
    ) +
    annotation_custom(bar_grob, xmin = -180, xmax = -100, ymin = -60, ymax = -15)
  
  # Save figure
  filename <- paste0("change_type_distribution_map_", gsub("-", "_", period_label), ".png")
  filepath <- file.path(output_dir, filename)
  ggsave(filepath, p3, width = 12, height = 8, dpi = 300)
  cat("Change type distribution plot saved:", filepath, "\n")
  
  return(list(plot = p3, filepath = filepath))
}

# Change Year Time Series Bar Chart
plot_change_timeline_periodic <- function(df_change, period_label, output_dir) {
  df_filtered <- filter_valid_data(df_change)
  
  yearly_changes <- df_filtered %>%
    filter(!is.na(change_year)) %>%
    count(change_year, name = "count")
  
  if(nrow(yearly_changes) == 0) {
    cat("Warning: No change year data found for period", period_label, "\n")
    return(NULL)
  }
  
  p4 <- ggplot(yearly_changes, aes(x = change_year, y = count)) +
    geom_col(fill = "steelblue", alpha = 0.7) +
    geom_smooth(method = "loess", se = TRUE, color = "red", linewidth = 1) +
    labs(title = paste("Temporal Trend of Vegetation Response Changes -", period_label),
         x = "Year", y = "Number of Changes",
         subtitle = "Red line shows smoothed trend; no_matched_windows data excluded") +
    theme_bw() +
    theme(
      plot.title = element_text(face = "bold", size = 12),
      panel.grid.minor = element_blank()
    )
  
  # Save figure
  filename <- paste0("change_timeline_", gsub("-", "_", period_label), ".png")
  filepath <- file.path(output_dir, filename)
  ggsave(filepath, p4, width = 10, height = 6, dpi = 300)
  cat("Change timeline plot saved:", filepath, "\n")
  
  return(list(plot = p4, filepath = filepath))
}

# State Transition Matrix Heatmap
plot_transition_matrix_periodic <- function(df_change, period_label, output_dir) {
  df_filtered <- filter_valid_data(df_change)
  
  valid_states <- c("Stable", "VGC", "VSO")
  
  transition_data <- df_filtered %>%
    filter(!is.na(from) & !is.na(to) & 
             from %in% valid_states & to %in% valid_states) %>%
    count(from, to, name = "count")
  
  if(nrow(transition_data) == 0) {
    cat("Warning: No transition data found for period", period_label, "\n")
    return(NULL)
  }
  
  complete_matrix <- expand.grid(from = valid_states, to = valid_states) %>%
    left_join(transition_data, by = c("from", "to")) %>%
    mutate(count = ifelse(is.na(count), 0, count))
  
  total_changes <- sum(complete_matrix$count)
  complete_matrix$percentage <- complete_matrix$count / total_changes * 100
  
  p5 <- ggplot(complete_matrix, aes(x = from, y = to, fill = percentage)) +
    geom_tile(color = "white", linewidth = 1) +
    geom_text(aes(label = sprintf("%.1f%%", percentage)), 
              color = "white", fontface = "bold", size = 4) +
    scale_fill_gradient(name = "Percentage", 
                        low = "#90EE90", high = "#FF69B4",
                        breaks = scales::pretty_breaks(n = 5)) +
    labs(title = paste("Vegetation Response State Transition Matrix -", period_label),
         x = "From State", y = "To State",
         subtitle = "Percentage of total transitions (Stable, VGC, VSO only); no_matched_windows excluded") +
    theme_bw() +
    theme(
      plot.title = element_text(face = "bold", size = 12),
      plot.subtitle = element_text(size = 10),
      axis.text.x = element_text(angle = 0, hjust = 0.5),
      axis.text.y = element_text(angle = 0),
      panel.grid = element_blank(),
      legend.position = "right"
    )
  
  # Save figure
  filename <- paste0("transition_matrix_", gsub("-", "_", period_label), ".png")
  filepath <- file.path(output_dir, filename)
  ggsave(filepath, p5, width = 8, height = 6, dpi = 300)
  cat("Transition matrix plot saved:", filepath, "\n")
  
  return(list(plot = p5, filepath = filepath))
}

# Change Type Bar Chart
plot_change_type_bar_periodic <- function(df_change, period_label, output_dir) {
  df_filtered <- filter_valid_data(df_change)
  
  valid_states <- c("Stable", "VGC", "VSO")
  
  change_summary <- df_filtered %>%
    filter(!is.na(from) & !is.na(to) & from != to &
             from %in% valid_states & to %in% valid_states) %>%
    mutate(change_type = paste(from, "→", to)) %>%
    count(change_type, name = "count") %>%
    mutate(percentage = count / sum(count) * 100) %>%
    arrange(desc(count))
  
  if(nrow(change_summary) == 0) {
    cat("Warning: No change type data found for period", period_label, "\n")
    return(NULL)
  }
  
  all_change_types <- c("Stable → VGC", "Stable → VSO", "VGC → Stable", 
                        "VGC → VSO", "VSO → Stable", "VSO → VGC")
  
  change_summary <- change_summary %>%
    right_join(data.frame(change_type = all_change_types), by = "change_type") %>%
    mutate(
      count = ifelse(is.na(count), 0, count),
      percentage = ifelse(is.na(percentage), 0, percentage)
    ) %>%
    arrange(desc(percentage))
  
  p6 <- ggplot(change_summary, aes(x = reorder(change_type, percentage), y = percentage, fill = change_type)) +
    geom_col(alpha = 0.8, color = "white", linewidth = 0.5) +
    scale_fill_manual(
      name = "Change Type", 
      values = c(
        "Stable → VGC" = "#7FCDCD",
        "Stable → VSO" = "#98D8C8",
        "VGC → Stable" = "#F7DC6F", 
        "VGC → VSO" = "#F8C471", 
        "VSO → Stable" = "#85C1E9", 
        "VSO → VGC" = "#D7BDE2")
    ) +
    geom_text(aes(label = sprintf("%.1f%%", percentage)), 
              hjust = -0.1, color = "black", fontface = "bold", size = 3.5) +
    coord_flip() +
    labs(title = paste("Global Distribution of Change Types -", period_label),
         subtitle = "Proportion of vegetation response transitions (no_matched_windows excluded)",
         x = "Change Type", y = "Percentage (%)") +
    theme_bw() +
    theme(
      plot.title = element_text(face = "bold", size = 12),
      plot.subtitle = element_text(size = 10),
      legend.position = "none",
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),
      axis.text.y = element_text(size = 10),
      axis.text.x = element_text(size = 10),
      axis.title = element_text(size = 11)
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1)))
  
  # Save figure
  filename <- paste0("change_type_bar_", gsub("-", "_", period_label), ".png")
  filepath <- file.path(output_dir, filename)
  ggsave(filepath, p6, width = 8, height = 6, dpi = 300)
  cat("Change type bar plot saved:", filepath, "\n")
  
  return(list(plot = p6, filepath = filepath))
}

# Main Plotting Function
create_all_plots_for_period <- function(combined_data, period_label, output_dir) {
  cat(paste("\n=== Creating plots for period", period_label, "===\n"))
  
  plot_results <- list()
  
  # First output data filtering statistics
  original_count <- nrow(combined_data)
  filtered_count <- nrow(filter_valid_data(combined_data))
  filtered_out <- original_count - filtered_count
  cat("Period", period_label, "- Original data:", original_count, "rows\n")
  cat("Period", period_label, "- Filtered data:", filtered_count, "rows\n")
  cat("Period", period_label, "- Filtered out no_matched_windows:", filtered_out, "rows\n\n")
  
  # Create all charts
  plot_results$persistent <- plot_persistent_type_comprehensive_periodic(combined_data, period_label, output_dir)
  plot_results$frequency <- plot_change_frequency_map_periodic(combined_data, period_label, output_dir)
  plot_results$first_year_continuous <- plot_first_change_year_map1_periodic(combined_data, period_label, output_dir)
  plot_results$first_year_classified <- plot_first_change_year_map2_periodic(combined_data, period_label, output_dir)
  plot_results$change_type <- plot_change_type_map_periodic(combined_data, period_label, output_dir)
  plot_results$timeline <- plot_change_timeline_periodic(combined_data, period_label, output_dir)
  plot_results$transition <- plot_transition_matrix_periodic(combined_data, period_label, output_dir)
  plot_results$bar_chart <- plot_change_type_bar_periodic(combined_data, period_label, output_dir)
  
  return(plot_results)
}

# Main Program - Integrated Version -----
main_analysis_periodic <- function(period_length = 5) {  # Default 5-year period
  start_time <- Sys.time()
  cat("=== Integrated Vegetation Response Analysis (Periodic Version) ===\n")
  cat("Program started at:", as.character(start_time), "\n")
  cat("Period length:", period_length, "years\n\n")
  
  # Set data path
  base_path <- './Out-Table-and-Figure-VPD et al/2025-07-17/GroupBy_100_SlidingWindow_Claude-all'
  vi_name <- "GPP"
  
  # Get all group info
  cat("\n=== Getting Group Information ===\n")
  group_info <- get_group_info(base_path, vi_name)
  
  if(is.null(group_info)) {
    stop("No group information found!")
  }
  
  cat("Found", nrow(group_info), "groups to process\n")
  print(head(group_info, 10))
  
  # Create output folders
  results_dir <- file.path(Datedir, paste0("Group_Results_", period_length, "yr"))
  dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Create time periods
  # Assuming analysis years range, adjust according to your actual data
  analysis_years <- 1991:2020  # Adjust according to your actual data range
  periods <- create_time_periods(analysis_years, period_length)
  
  cat("\n=== Time Periods Created ===\n")
  print(periods)
  
  # Set parallel processing
  n_cores <- get_available_cores(max_cores)
  
  cat("\n=== Starting Parallel Group Processing ===\n")
  cat("Using", n_cores, "cores for parallel processing\n")
  
  # Process all groups in parallel
  cl <- makeCluster(n_cores, outfile = file.path(Datedir, paste0("group_processing_log_", period_length, "yr.txt")))
  registerDoParallel(cl)
  
  processing_results <- foreach(i = 1:nrow(group_info),
                                .packages = c("data.table", "dplyr"),
                                .export = c("process_group_data_periodic", "process_single_grid_periodic", 
                                            "perform_equivalence_test", "assign_year_to_period",
                                            "base_path", "vi_name", "results_dir", "group_info", 
                                            "periods", "period_length"),
                                .errorhandling = "pass") %dopar% {
                                  group_num <- group_info$group_num[i]
                                  group_name <- group_info$group_name[i]
                                  
                                  process_group_data_periodic(group_num, group_name, base_path, 
                                                              vi_name, results_dir, periods, period_length)
                                }
  
  stopCluster(cl)
  
  # Summarize processing results
  cat("\n=== Processing Summary ===\n")
  successful_groups <- 0
  failed_groups <- 0
  
  for(i in seq_along(processing_results)) {
    result <- processing_results[[i]]
    group_name <- group_info$group_name[i]
    
    if(class(result) != "try-error" && !is.null(result) && result$success) {
      successful_groups <- successful_groups + 1
      cat("✓", group_name, "- Grids:", result$n_grids, 
          "| Yearly:", result$n_yearly_records, 
          "| Detailed:", result$n_detailed_records, 
          "| Changes:", result$n_change_records, "\n")
    } else {
      failed_groups <- failed_groups + 1
      error_msg <- if(class(result) == "try-error") as.character(result) else result$error
      cat("✗", group_name, "- Error:", error_msg, "\n")
    }
    
    print_progress(i, nrow(group_info), "Group Processing")
  }
  
  cat("\nGroup processing completed:", successful_groups, "successful,", failed_groups, "failed\n")
  
  # Combine results from all groups
  cat("\n=== Combining All Group Results ===\n")
  combined_results <- combine_all_group_results_periodic(processing_results, results_dir, vi_name, period_length)
  
  # Create plots for each period
  if(!is.null(combined_results$change_data) && nrow(combined_results$change_data) > 0) {
    cat("\n=== Creating Plots for Each Period ===\n")
    
    all_plot_results <- list()
    
    # Group data by period and plot
    unique_periods <- unique(combined_results$change_data$period_label)
    unique_periods <- unique_periods[!is.na(unique_periods)]
    
    if(length(unique_periods) > 0) {
      for(period_label in unique_periods) {
        period_data <- combined_results$change_data %>% 
          filter(period_label == !!period_label)
        
        if(nrow(period_data) > 0) {
          period_plot_results <- create_all_plots_for_period(period_data, period_label, figure_date_dir)
          all_plot_results[[period_label]] <- period_plot_results
        }
      }
    }
    
    # Also create plots for overall data
    overall_plot_results <- create_all_plots_for_period(combined_results$change_data, 
                                                        paste0("Overall_", period_length, "yr"), 
                                                        figure_date_dir)
    all_plot_results[["Overall"]] <- overall_plot_results
  }
  
  # Program ended
  end_time <- Sys.time()
  run_time <- end_time - start_time
  cat("\n=== Analysis Completed Successfully ===\n")
  cat("End time:", as.character(end_time), "\n")
  cat("Total runtime:", run_time, "\n")
  
  cat("\n=== Final Statistics ===\n")
  cat("Successful groups:", successful_groups, "\n")
  cat("Failed groups:", failed_groups, "\n")
  cat("Period length:", period_length, "years\n")
  cat("Number of periods:", nrow(periods), "\n")
  if(!is.null(combined_results)) {
    cat("Total yearly records:", combined_results$yearly_records, "\n")
    cat("Total detailed records:", combined_results$detailed_records, "\n")
    cat("Total change records:", combined_results$change_records, "\n")
  }
  
  cat("\n=== Output Files Structure ===\n")
  cat("├── Results (", results_dir, ")\n", sep = "")
  cat("│   ├── Group_1/ ... Group_N/\n")
  cat("│   ├── Combined_All_Groups_response_result_yearly_", vi_name, "_", period_length, "yr.csv\n", sep = "")
  cat("│   ├── Combined_All_Groups_response_type_change_", vi_name, "_", period_length, "yr.csv\n", sep = "")
  cat("│   └── Combined_All_Groups_detailed_periods_", vi_name, "_", period_length, "yr.csv\n", sep = "")
  cat("├── Plots (", figure_date_dir, ")\n", sep = "")
  cat("│   ├── persistent_type_comprehensive_[period].png\n")
  cat("│   ├── vegetation_change_frequency_map_[period].png\n")
  cat("│   ├── first_change_year_map1_[period].png\n")
  cat("│   ├── first_change_year_map2_[period].png\n")
  cat("│   ├── change_type_distribution_map_[period].png\n")
  cat("│   ├── change_timeline_[period].png\n")
  cat("│   ├── transition_matrix_[period].png\n")
  cat("│   └── change_type_bar_[period].png\n")
  cat("\n=== All Functions Completed Successfully ===\n")
  
  return(list(
    combined_results = combined_results,
    plot_results = if(exists("all_plot_results")) all_plot_results else NULL,
    periods = periods,
    run_time = run_time,
    successful_groups = successful_groups,
    failed_groups = failed_groups
  ))
}

# Run Main Program
cat("=== Integrated Vegetation Response Analysis System (Periodic Version) ===\n")
cat("This script combines statistical analysis and visualization:\n")
cat("1. Group-wise IRF data processing by time periods\n")
cat("2. Detailed period analysis with transition years\n")
cat("3. Response type classification (VGC/VSO/Stable) by periods\n")
cat("4. Persistent type tracking within each period\n")
cat("5. Global result combination and statistics\n")
cat("6. Comprehensive visualization for each period\n")
cat("\nYou can choose period length (5 or 10 years):\n")

# Run Analysis - 5-year period version
cat("\n=== Running 5-year period analysis ===\n")
results_5yr <- main_analysis_periodic(period_length = 5)

# Run Analysis - 10-year period version  
cat("\n=== Running 10-year period analysis ===\n")
results_10yr <- main_analysis_periodic(period_length = 10)

cat("\n=== Both analyses completed! ===\n")
cat("Check the output folders for results and plots organized by period length.\n")
