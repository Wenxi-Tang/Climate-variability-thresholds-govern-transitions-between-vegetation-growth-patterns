
setwd("E:/GlobalVeg")
RUN_name <- 'Threshold_analysis_segment'

variable_name <- "VPD+CO2+SM-CO2TRUE"
RunCodeData_path <- paste0('1.RunCodeData-',variable_name)
VAR_path <- paste0('2.VARmodelResult-',variable_name,'/ALL_1982_2020')

# Create directories -----
# Folder for output tables and figures
dir.create('./Out-Table-and-Figure-VPD et al', showWarnings = T)

# Folder for the current date
Datedir_main <- file.path('./Out-Table-and-Figure-VPD et al', Sys.Date())
dir.create(Datedir_main, showWarnings = T)

# Folder for the RUN_name
Datedir <- file.path(Datedir_main, RUN_name)
dir.create(Datedir, showWarnings = FALSE)

# Load packages -----
library(ggplot2)
library(dplyr)
library(reshape2)
library(tidyr)
library(segmented)
library(raster)

# Define functions -----
saveToRsater <- function(labs, values){
  r_test <- raster("./Data/resample_example_0.5degree.tif")
  values(r_test) <- NA
  values(r_test)[labs] <- values
  return(r_test)
}

cutByLULC <- function(raster){
  LULC <- raster("./Data/GLC_2000/LULC_1to15_2000.tif")
  LULC_lab <- which(values(LULC) %in% c(1:15))
  r_LULC <- saveToRsater(labs = LULC_lab, 
                         values = values(raster)[LULC_lab])
  return(r_LULC)
}

# sd(veg) ~ sd(cli) -----
cli_dirs <- list.dirs(file.path('./Out-Table-and-Figure-VPD et al',RunCodeData_path,
                                'RunThreshold/CLI_VI_data_group'),
                      recursive = F)
vi_paths <- list.files(file.path('./Out-Table-and-Figure-VPD et al',VAR_path),
                       '.RData$', full.names = T)


main_cli_name = "CO2"
path_i = 1
i = 1
j = 1
# "TMP", "PRE", 
for(main_cli_name in c("TMP", "PRE", "SR", "VPD", "SMroot", "SMsurf", "CO2")){
  for (path_i in 1) { # :length(vi_paths)
    cli_paths <- list.files(cli_dirs[path_i], '.Rdata', 
                            full.names = T)
    
    vi_path <- vi_paths[path_i]
    load(vi_path)
    vi_labs <- names(VARdata_combined_list)
    
    file_names <- basename(cli_paths)
    tmp_files <- cli_paths[grepl(main_cli_name, file_names)]
    
    VI_name <- basename(cli_dirs[path_i])
    CLI_name <- main_cli_name
    
    VI_CLI_dir <- file.path(Datedir, paste0(VI_name, "_",CLI_name))
    dir.create(VI_CLI_dir, showWarnings = FALSE)
    
    df.result_all <- data.frame()
    for (i in 1:length(tmp_files)) {
      load(tmp_files[i])
      
      NAMEs <- unlist(strsplit(basename(tmp_files[i]), ".Rdata"))[1]
      Group_name <- unlist(strsplit(NAMEs, "_"))[3]
      
      df.result_group <- data.frame()
      veg_cli_SDdata_list <- list()
      segments_plot_list <- list()
      
      for (j in 1:length(veg_cli_list)) { 
        lab <- names(veg_cli_list)[j]
        
        df.all_cli <- veg_cli_list[[j]]
        df.all_vi <- VARdata_combined_list[[which(vi_labs == lab)]]
        cat(lab, j, which(vi_labs == lab), vi_labs[which(vi_labs == lab)] == lab,"\n")
        
        df.sd_vi <- df.all_vi %>% 
          group_by(year_new) %>% 
          summarise(mean_GS = mean(GPP),
                    sd_GS = sd(GPP),
                    max_GS = max(GPP))
        
        df.sd_cli <- df.all_cli %>% 
          group_by(year_new) %>% 
          summarise(mean_cli = mean(LCE_data),
                    sd_cli = sd(LCE_data))
        dftest <- df.all_cli %>% subset(year_new == 1982)
        
        # if(CLI_name == "PRE"){
        #   df.sd_cli <- df.all_cli %>% 
        #     group_by(year_new) %>% 
        #     summarise(mean_cli = sum(LCE_data),
        #               sd_cli = sd(LCE_data))
        # }else{
        #   df.sd_cli <- df.all_cli %>% 
        #     group_by(year_new) %>% 
        #     summarise(mean_cli = mean(LCE_data),
        #               sd_cli = sd(LCE_data))
        # }
        
        df.sd <- merge(df.sd_vi, df.sd_cli, by = "year_new")
        df.sd <- na.omit(df.sd)
        
        df <- data.frame(x = df.sd$sd_cli, 
                         y = df.sd$sd_GS)
        
        # Segmented regression model 
        fit <- lm(y ~ x, data = df)
        
        # Try the segmented function
        result <- tryCatch({
          # Fit segmented regression model
          segmented.fit <- segmented(fit, seg.Z = ~x)
          summary(segmented.fit)
          
          # If successful, store the result
          list(
            success = TRUE,
            model = segmented.fit
          )
        }, error = function(e) {
          # Catch error and return NA
          list(success = FALSE, message = e$message)
        })
        
        if(result$success == T){
          segmented.fit <- segmented(fit, seg.Z = ~x)
          
          if(is.null(segmented.fit$psi[2])){
            signif = "F"
            df.result <- data.frame(lab, 
                                    Group_name,
                                    j,
                                    signif,
                                    t(data.frame(rep(NA, (ncol(df.result_group) - 4)))))
            colnames(df.result) <- colnames(df.result_group)
            row.names(df.result) <- 1
            
            veg_cli_SDdata_list[[j]] <- df.sd
            segments_plot_list[[j]] <- NA
            
            df.result_group <- rbind(df.result_group, df.result)
            df.result_all <- rbind(df.result_all, df.result)
            cat(NAMEs, lab, j, signif, "no breakponit! \n")
            
            next()
          }else{
            seg_result <- summary(segmented.fit)
            df$fit <- broken.line(segmented.fit)$fit
            df.sd$fit <- broken.line(segmented.fit)$fit
            veg_cli_SDdata_list[[j]] <- df.sd
            
            # Plotting
            segments_plot_list[[j]] <- ggplot(df, aes(x = x, y = y)) +
              geom_point(color = "steelblue", size = 3) +
              geom_line(aes(x = x, y = fit), 
                        color = "red") +
              labs(title = paste("Segmented Regression: ", lab),
                   x = CLI_name, y = VI_name)
            
            # Extract goodness of fit
            R2 <- round(seg_result$r.squared, 4)
            
            # Extract coefficients/slopes/p-values for two segments
            # Intercept
            intercept <- intercept(segmented.fit)$x
            intercept_df <- as.data.frame(t(intercept))
            # Slope
            slope <- slope(segmented.fit)$x
            colnames(slope)[2:5] <- c("St.Err", "t.value", "CI.l", "CI.u")
            slope.m <- melt(slope)
            slope_df <- data.frame(t(slope.m$value))
            colnames(slope_df) <- c(as.character(slope.m$Var1[1:2]), 
                                    paste0(slope.m$Var1[3:nrow(slope.m)], 
                                           "_",
                                           slope.m$Var2[3:nrow(slope.m)]))
            # Significance
            coef_summary <- round(seg_result$coefficients, 2)
            pvalue1 <- coef_summary["x", "Pr(>|t|)"]
            pvalue2 <- coef_summary["U1.x", "Pr(>|t|)"]
            
            fit1 <- if(slope_df$slope1 >= 0){
              paste0("y = ", intercept_df$intercept1, " + ", slope_df$slope1, " * x")
            }else{
              paste0("y = ", intercept_df$intercept1, " - ", abs(slope_df$slope1), " * x")
            }
            
            fit2 <- if((slope_df$slope2) >= 0){
              paste0("y = ", intercept_df$intercept2, " + ", slope_df$slope2, " * x")
            }else{
              paste0("y = ", intercept_df$intercept2, " - ", abs(slope_df$slope2), " * x")
            }
            
            signif <- if(pvalue1 < 0.1 || pvalue2 < 0.1){
              "T"
            }else{
              "F"
            }
            
            # Breakpoint / Significance of breakpoint / Confidence interval
            # Breakpoint
            breakpoint_x <- round(segmented.fit$psi[2], 4)
            breakpoint_y <- intercept_df$intercept1 + slope_df$slope1 * breakpoint_x
            # Significance
            davies_test <- davies.test(fit, k = 10)
            breakpoint_p <- davies_test$p.value
            # Confidence interval
            breakpoint_CI <- as.data.frame(confint(segmented.fit))
            colnames(breakpoint_CI)[2:3] <- c("psi1.x_CI.l", "psi1.x_CI.u")
            
            df.result <- data.frame(lab, 
                                    Group_name,
                                    j,
                                    signif, 
                                    R2,
                                    breakpoint_x, 
                                    breakpoint_y,
                                    breakpoint_p,
                                    breakpoint_CI[,2:3],
                                    intercept_df, 
                                    slope_df,
                                    pvalue1, 
                                    pvalue2, 
                                    seg1 = paste0("x <= ", breakpoint_x), 
                                    fit1, 
                                    seg2 = paste0("x > ", breakpoint_x),
                                    fit2)
            
            df.result_group <- rbind(df.result_group, df.result)
            df.result_all <- rbind(df.result_all, df.result)
            cat(NAMEs, lab, j, signif, "\n")
          }
        }else{
          signif = "F"
          df.result <- data.frame(lab, 
                                  Group_name,
                                  j,
                                  signif,
                                  t(data.frame(rep(NA, (ncol(df.result_group) - 4)))))
          colnames(df.result) <- colnames(df.result_group)
          row.names(df.result) <- 1
          
          veg_cli_SDdata_list[[j]] <- df.sd
          segments_plot_list[[j]] <- NA
          
          df.result_group <- rbind(df.result_group, df.result)
          df.result_all <- rbind(df.result_all, df.result)
          cat(NAMEs, lab, j, signif, "can't fit segment model! \n")
          
          next()
        }
      }
      
      names(segments_plot_list) <- names(veg_cli_list)
      save(segments_plot_list, 
           file = file.path(VI_CLI_dir, paste0(NAMEs, '_plot.RData')))
      
      names(veg_cli_SDdata_list) <- names(veg_cli_list)
      save(veg_cli_SDdata_list,
           file = file.path(VI_CLI_dir, paste0(NAMEs, '_SDdata.RData')))
      
      row.names(df.result_group) <- 1:nrow(df.result_group)
      write.csv(df.result_group, 
                file.path(VI_CLI_dir, paste0(NAMEs, '.csv')),
                row.names = F)
      cat(NAMEs, "OK!", "\n")
    }
    
    row.names(df.result_all) <- 1:nrow(df.result_all)
    write.csv(df.result_all, 
              file.path(VI_CLI_dir, paste0(VI_name,"_", CLI_name, '_all.csv')),
              row.names = F)
    cat(VI_name, CLI_name, "OK!", "\n")
  }
}

# Merge SD data -----------------------
# Set file directory
# SD_dirs <- list.dirs('./Out-Table-and-Figure-VPD et al/4.ThresholdResult-VPD+CO2+SM-CO2TRUE/GPP',
#                      recursive = FALSE)
names <- basename(SD_dirs)

# Output path
Datedir <- './Merged_SDdata'
if (!dir.exists(Datedir)) dir.create(Datedir)

# Iterate through each climate factor
for (i in seq_along(SD_dirs)) {
  dir_name <- names[i]
  SD_group_paths <- list.files(SD_dirs[i], full.names = TRUE, pattern = '_SDdata\\.RData$')
  
  message("Processing climate factor: ", dir_name)
  
  # Create an empty list for all SD data of the current climate factor
  SD_list <- list()
  
  # Iterate through each _SDdata.RData file
  for (j in seq_along(SD_group_paths)) {
    load(SD_group_paths[j])  # Load file, assuming it contains veg_cli_SDdata_list
    print(head(veg_cli_SDdata_list[[1]]))
    
    if (exists("veg_cli_SDdata_list")) {
      # Flatten merge: Append each element directly to SD_list
      SD_list <- c(SD_list, veg_cli_SDdata_list)
    } else {
      warning("Object veg_cli_SDdata_list not found in file: ", SD_group_paths[j])
    }
    
    rm(veg_cli_SDdata_list)  # Clean up environment variables
  }
  
  # Save merged climate factor SD data
  # save(SD_list, file = file.path(Datedir, paste0(dir_name, "_SDdata.RData")))
  
  message("✅ ", dir_name, " Merge completed, total ", length(SD_list), " elements.")
}

message("🎉 Flat merge of SD data for all climate factors completed!")

# Merge grouped breakpoint data + Visualization + Save as TIFF -----
df_dirs <- list.dirs(Datedir, recursive = F)
# dir_i <- 1
# i <- 1
for (dir_i in 1:length(df_dirs)) {
  df_paths <- list.files(df_dirs[dir_i], '.csv$', full.names = T)
  df_paths_filtered <- subset(df_paths, !grepl("all.csv", df_paths)) # Filter out paths containing all.csv (all.csv does not contain lab column)
  VI_CLI_name <- basename(df_dirs[dir_i])
  
  
  df_all <- data.frame()
  for (i in 1:length(df_paths_filtered)) {
    df <- read.csv(df_paths_filtered[i])
    df_all <- rbind(df_all, df)
    print(VI_CLI_name)
    print(head(df))
  }
  write.csv(df_all, file.path(Datedir, paste0(VI_CLI_name, '_breakpoint_all.csv')))
  
  r_mean_break.x <- cutByLULC(saveToRsater(as.numeric(df_all$lab), df_all$breakpoint_x))
  plot(r_mean_break.x, main = paste0(VI_CLI_name,"_breakpoint_x"))
  writeRaster(r_mean_break.x, 
              filename = file.path(Datedir, paste0(VI_CLI_name,"_breakpoint_x")), 
              format = "GTiff")
  
  r_mean_break.y <- cutByLULC(saveToRsater(as.numeric(df_all$lab), df_all$breakpoint_y))
  plot(r_mean_break.y, main = paste0(VI_CLI_name,"_breakpoint_y"))
  writeRaster(r_mean_break.y, 
              filename = file.path(Datedir, paste0(VI_CLI_name,"_breakpoint_y")), 
              format = "GTiff")
}