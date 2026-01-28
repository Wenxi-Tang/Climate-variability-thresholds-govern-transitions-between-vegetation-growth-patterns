
# ==============================================================================
# Comprehensive Vegetation Response Visualization Script - Direct Execution Version
# Scientific visualization combining threshold crossing analysis and time series analysis
# ==============================================================================

# Set working directory
setwd('E:/GlobalVeg')

# Create output folders -----
RUN_name <- "Comprehensive-Visualization-5yr"
Datedir_main <- file.path('./Out-Table-and-Figure-VPD et al', Sys.Date())
dir.create(Datedir_main, showWarnings = FALSE, recursive = TRUE)
Datedir <- file.path(Datedir_main, RUN_name)
dir.create(Datedir, showWarnings = FALSE, recursive = TRUE)

# Create figure output folder
figure_date_dir <- file.path('./Figure', Sys.Date())
dir.create('./Figure', showWarnings = FALSE, recursive = TRUE)
dir.create(figure_date_dir, showWarnings = FALSE, recursive = TRUE)

cat("=== Comprehensive Vegetation Response Visualization Analysis ===\n")
cat("Output folder:", Datedir, "\n")
cat("Figure folder:", figure_date_dir, "\n\n")

# Load required packages -----
library(dplyr)
library(tidyr)
library(raster)
library(ggplot2)
library(pROC)
library(caret)
library(gridExtra)
library(grid)
library(patchwork)
library(viridis)
library(RColorBrewer)
library(scales)
library(rnaturalearth)
library(rnaturalearthdata)
library(sf)
library(cowplot)
library(data.table)

cat("All packages loaded!\n\n")

# Define helper functions -----
saveToRaster <- function(labs, values){
  r_test <- raster("./Data/resample_example_0.5degree.tif")
  values(r_test) <- NA
  values(r_test)[labs] <- values
  return(r_test)
}

cutByLULC <- function(raster){
  LULC <- raster("./Data/GLC_2000/LULC_1to15_2000.tif")
  LULC_lab <- which(values(LULC) %in% c(1:15))
  r_LULC <- saveToRaster(labs = LULC_lab, 
                         values = values(raster)[LULC_lab])
  return(r_LULC)
}

# Get world map
world_shp <- ne_coastline()

# Data filtering function
filter_valid_data <- function(df) {
  df %>% filter(note != "no_matched_windows" | is.na(note))
}

# ==============================================================================
# Step 0: Read all required data
# ==============================================================================

cat("Step 0: Reading data...\n")

# 1. Read threshold data
Threshold_paths <- list.files('./Out-Table-and-Figure-VPD et al/4.ThresholdResult-VPD+CO2+SM-CO2TRUE/GPP',
                              full.names = TRUE, 
                              pattern = '_breakpoint_x.tif$')

# 2. Read SD data
SD_paths <- list.files('./Out-Table-and-Figure-VPD et al/4.ThresholdResult-VPD+CO2+SM-CO2TRUE/Merged_SDdata',
                       full.names = TRUE, 
                       pattern = '_SDdata.RData$')
load(SD_paths[1])
SD_list[1]
SD_list[[1]]

# 3. Read vegetation response type classification data
r_type <- raster('./Out-Table-and-Figure-VPD et al/3.VegResponseType-VPD+CO2+SM-CO2TRUE/GPP_VPD+CO2+SM-CO2TRUE_3type.tif')
df_type <- data.frame(
  lab = 1:ncell(r_type),
  type = values(r_type)
) %>%
  mutate(type_name = case_when(
    type == 1 ~ "VSO",
    type == 2 ~ "VGC",
    type == 3 ~ "Stable",
    TRUE ~ NA_character_
  )) %>%
  na.omit()

# 4. Read sliding window data
df_transitions <- read.csv("E:/Global-SlidingWindow/Out-Figure-and-Table/5yr统计/Combined_All_Groups_response_type_change_GPP_5yr.csv")

cat("Sliding window data dimensions:", nrow(df_transitions), "x", ncol(df_transitions), "\n")
cat("Sliding window data structure:\n")
print(head(df_transitions))

# 5. Read SHAP threshold data
df_trippingtype <- read.csv('./Out-Table-and-Figure-VPD et al/5.SHAP-VPD+CO2+SM-CO2TRUE/SHAP_ClimateThreshold_TrippingPoint.csv')

# Define climate factor groups 
energy_vars <- c("CO2", "TMP", "SR")            # Energy-limited -> Stable↔VGC
water_vars <- c("VPD", "PRE", "SMroot", "SMsurf")  # Water-limited -> VGC↔VSO
climate_vars <- c(energy_vars, water_vars)

single_threshold_vars <- c("CO2", "VPD", "PRE", "SMroot", "SMsurf")
double_threshold_vars <- c("TMP", "SR")

cat("Data reading completed!\n")
cat("Energy-limited factors:", paste(energy_vars, collapse=", "), "\n")
cat("Water-limited factors:", paste(water_vars, collapse=", "), "\n\n")

# Check data quality
cat("Data quality check:\n")
cat("- Valid sliding window data:", nrow(filter_valid_data(df_transitions)), "\n")
cat("- Filtered out no_matched_windows:", sum(df_transitions$note == "no_matched_windows", na.rm = TRUE), "\n")
cat("- Number of unique periods:", length(unique(df_transitions$period_label[!is.na(df_transitions$period_label)])), "\n\n")

# == Step 1: Identify Climate Threshold Crossing Events (Group Coding + Trend Judgment) ====================================
cat("Step 1: Identifying climate threshold crossing and trend events (Group Coding)...\n")

# Function: Calculate single threshold crossing and trend events (Group Coding)
calculate_single_threshold_cross <- function(sd_list, var_name, tripping_data, energy_vars) {
  cat(paste0("  Processing single threshold factor: ", var_name, "\n"))
  
  # Determine factor type
  is_energy <- var_name %in% energy_vars
  
  # Extract threshold for this factor from df_trippingtype
  threshold_info <- tripping_data %>%
    filter(Climate_Variable == var_name)
  
  if(nrow(threshold_info) == 0) {
    cat(paste0("    Warning: ", var_name, " not found in SHAP threshold table\n"))
    return(NULL)
  }
  
  shap_threshold <- threshold_info$X_Position[1]
  left_type <- threshold_info$Left_contribution[1]
  right_type <- threshold_info$Right_contribution[1]
  
  cat(paste0("    SHAP Threshold = ", round(shap_threshold, 4), 
             ", Left:", left_type, " -> Right:", right_type, "\n"))
  cat(paste0("    Factor Type: ", ifelse(is_energy, "Energy", "Water"), "\n"))
  
  cross_events_list <- list()
  
  pixel_names <- names(sd_list)
  if(is.null(pixel_names)) {
    pixel_names <- as.character(1:length(sd_list))
  }
  
  for(i in 1:length(sd_list)) {
    pixel_data <- sd_list[[i]]
    pixel_lab <- pixel_names[i]
    
    if(is.null(pixel_data) || nrow(pixel_data) < 2) next
    
    if(!"sd_cli" %in% names(pixel_data)) {
      cat(paste0("    Warning: Pixel ", pixel_lab, " missing sd_cli column\n"))
      next
    }
    
    cli_values <- pixel_data$sd_cli
    years <- pixel_data$year_new
    
    position <- ifelse(cli_values < shap_threshold, "left", "right")
    
    # Core modification: Consider crossing + non-crossing but with trend cases
    cross_dir <- rep(0, length(years))
    event_type <- rep("none", length(years))  # New: Event type marker
    
    for(j in 2:length(years)) {
      cli_change <- cli_values[j] - cli_values[j-1]  # Climate value change
      
      # Determine if threshold is crossed
      is_cross <- position[j-1] != position[j]
      
      if(is_cross) {
        # Case 1: Crossing threshold
        if(position[j-1] == "left" && position[j] == "right") {
          # Crossing from left to right
          cross_dir[j] <- ifelse(is_energy, 1, 3)  # Energy: Stable→VGC, Water: VSO→VGC
          event_type[j] <- "cross_L2R"
        } else if(position[j-1] == "right" && position[j] == "left") {
          # Crossing from right to left
          cross_dir[j] <- ifelse(is_energy, 2, 3)  # Energy: VGC→Stable, Water: VSO→VGC
          event_type[j] <- "cross_R2L"
        }
      } else {
        # Case 2: Not crossing threshold, but has change trend
        if(abs(cli_change) > 0) {  # Ensure there is change
          if(position[j] == "left") {
            # On the left side of threshold
            if(cli_change > 0) {
              # Moving right (approaching threshold)
              cross_dir[j] <- ifelse(is_energy, 2, 3)  # Energy: VGC→Stable, Water: VSO→VGC
              event_type[j] <- "left_increasing"
            } else {
              # Moving left (away from threshold)
              cross_dir[j] <- ifelse(is_energy, 1, 4)  # Energy: Stable→VGC, Water: VGC→VSO
              event_type[j] <- "left_decreasing"
            }
          } else {
            # On the right side of threshold
            if(cli_change > 0) {
              # Moving right (away from threshold)
              cross_dir[j] <- ifelse(is_energy, 1, 4)  # Energy: Stable→VGC, Water: VGC→VSO
              event_type[j] <- "right_increasing"
            } else {
              # Moving left (approaching threshold)
              cross_dir[j] <- ifelse(is_energy, 2, 3)  # Energy: VGC→Stable, Water: VSO→VGC
              event_type[j] <- "right_decreasing"
            }
          }
        }
      }
    }
    
    cross_events_list[[i]] <- data.frame(
      lab = pixel_lab,
      year = years,
      sd_cli = cli_values,
      threshold = shap_threshold,
      position = position,
      cross_dir = cross_dir,
      event_type = event_type,
      left_type = left_type,
      right_type = right_type,
      var = var_name,
      var_group = ifelse(is_energy, "Energy", "Water"),
      stringsAsFactors = FALSE
    )
  }
  
  return(dplyr::bind_rows(cross_events_list))
}

# Function: Calculate double threshold crossing and trend events (Group Coding)
calculate_double_threshold_cross <- function(sd_list, var_name, tripping_data, energy_vars) {
  cat(paste0("  Processing double threshold factor: ", var_name, "\n"))
  
  # Determine factor type (TMP and SR are both energy type)
  is_energy <- var_name %in% energy_vars
  
  threshold_info <- tripping_data %>%
    filter(Climate_Variable == var_name) %>%
    arrange(X_Position)
  
  if(nrow(threshold_info) != 2) {
    cat(paste0("    Warning: ", var_name, " should have 2 thresholds, actually found ", nrow(threshold_info), "\n"))
    return(NULL)
  }
  
  threshold1 <- threshold_info$X_Position[1]
  threshold2 <- threshold_info$X_Position[2]
  
  cat(paste0("    SHAP Threshold 1 = ", round(threshold1, 4), 
             " (", threshold_info$Left_contribution[1], "->", threshold_info$Right_contribution[1], ")\n"))
  cat(paste0("    SHAP Threshold 2 = ", round(threshold2, 4),
             " (", threshold_info$Left_contribution[2], "->", threshold_info$Right_contribution[2], ")\n"))
  
  cross_events_list <- list()
  
  pixel_names <- names(sd_list)
  if(is.null(pixel_names)) {
    pixel_names <- as.character(1:length(sd_list))
  }
  
  for(i in 1:length(sd_list)) {
    pixel_data <- sd_list[[i]]
    pixel_lab <- pixel_names[i]
    
    if(is.null(pixel_data) || nrow(pixel_data) < 2) next
    
    if(!"sd_cli" %in% names(pixel_data)) {
      cat(paste0("    Warning: Pixel ", pixel_lab, " missing sd_cli column\n"))
      next
    }
    
    cli_values <- pixel_data$sd_cli
    years <- pixel_data$year_new
    
    # Define three zones
    zone <- rep(NA, length(cli_values))
    zone[cli_values < threshold1] <- "Zone1"       # Low state -> Stable
    zone[cli_values >= threshold1 & cli_values < threshold2] <- "Zone2"  # Mid state -> VGC
    zone[cli_values >= threshold2] <- "Zone3"      # High state -> Stable
    
    # Core modification: Consider crossing + non-crossing but with trend cases
    cross_dir <- rep(0, length(years))
    event_type <- rep("none", length(years))
    
    for(j in 2:length(years)) {
      if(is.na(zone[j-1]) || is.na(zone[j])) next
      
      cli_change <- cli_values[j] - cli_values[j-1]
      
      # Determine if threshold is crossed
      is_cross <- zone[j-1] != zone[j]
      
      if(is_cross) {
        # Case 1: Crossing threshold
        if(zone[j-1] == "Zone1" && zone[j] == "Zone2") {
          cross_dir[j] <- 1  # Stable → VGC (crossing threshold1 from left to right)
          event_type[j] <- "cross_T1_L2R"
        } else if(zone[j-1] == "Zone2" && zone[j] == "Zone1") {
          cross_dir[j] <- 2  # VGC → Stable (crossing threshold1 from right to left)
          event_type[j] <- "cross_T1_R2L"
        } else if(zone[j-1] == "Zone2" && zone[j] == "Zone3") {
          cross_dir[j] <- 2  # VGC → Stable (crossing threshold2 from left to right)
          event_type[j] <- "cross_T2_L2R"
        } else if(zone[j-1] == "Zone3" && zone[j] == "Zone2") {
          cross_dir[j] <- 1  # Stable → VGC (crossing threshold2 from right to left)
          event_type[j] <- "cross_T2_R2L"
        }
      } else {
        # Case 2: Not crossing threshold, but has change trend
        if(abs(cli_change) > 0) {
          if(zone[j] == "Zone1") {
            # In Zone1 (left of threshold1)
            if(cli_change > 0) {
              # Moving right (approaching threshold1)
              cross_dir[j] <- 2  # VGC → Stable
              event_type[j] <- "Zone1_increasing"
            } else {
              # Moving left (away from threshold1)
              cross_dir[j] <- 1  # Stable → VGC
              event_type[j] <- "Zone1_decreasing"
            }
          } else if(zone[j] == "Zone2") {
            # In Zone2 (between threshold1 and threshold2)
            if(cli_change > 0) {
              # Moving right (approaching threshold2)
              cross_dir[j] <- 1  # Stable → VGC
              event_type[j] <- "Zone2_increasing"
            } else {
              # Moving left (approaching threshold1)
              cross_dir[j] <- 2  # VGC → Stable
              event_type[j] <- "Zone2_decreasing"
            }
          } else if(zone[j] == "Zone3") {
            # In Zone3 (right of threshold2)
            if(cli_change > 0) {
              # Moving right (away from threshold2)
              cross_dir[j] <- 2  # VGC → Stable
              event_type[j] <- "Zone3_increasing"
            } else {
              # Moving left (approaching threshold2)
              cross_dir[j] <- 1  # Stable → VGC
              event_type[j] <- "Zone3_decreasing"
            }
          }
        }
      }
    }
    
    cross_events_list[[i]] <- data.frame(
      lab = pixel_lab,
      year = years,
      sd_cli = cli_values,
      threshold1 = threshold1,
      threshold2 = threshold2,
      zone = zone,
      cross_dir = cross_dir,
      event_type = event_type,
      var = var_name,
      var_group = "Energy",
      stringsAsFactors = FALSE
    )
  }
  
  return(dplyr::bind_rows(cross_events_list))
}

# Process all climate factors
all_cross_events <- list()

for(var in climate_vars) {
  sd_path <- SD_paths[grepl(paste0("GPP_", var, "_"), SD_paths)]
  
  if(length(sd_path) == 0) {
    cat(paste0("  Warning: SD data for ", var, " not found\n"))
    next
  }
  
  load(sd_path)
  
  if(var %in% single_threshold_vars) {
    cross_events <- calculate_single_threshold_cross(SD_list, var, df_trippingtype, energy_vars)
  } else if(var %in% double_threshold_vars) {
    cross_events <- calculate_double_threshold_cross(SD_list, var, df_trippingtype, energy_vars)
  }
  
  if(!is.null(cross_events)) {
    all_cross_events[[var]] <- cross_events
  }
}

df_cross_all <- dplyr::bind_rows(all_cross_events)

head(df_cross_all)
if(nrow(df_cross_all) > 0) {
  write.csv(df_cross_all, file.path(Datedir, "Step1_CrossEvent_All.csv"), row.names = FALSE)
  cat(paste0("  Identified ", nrow(df_cross_all), " event records\n"))
  cat(paste0("  Involving ", length(unique(df_cross_all$lab)), " pixels\n"))
  
  # Statistics of cross/trend events by direction
  cat("\nEvent direction statistics:\n")
  print(table(df_cross_all$cross_dir, df_cross_all$var_group))
  
  cat("\nEvent type statistics:\n")
  event_summary <- df_cross_all %>%
    filter(cross_dir != 0) %>%
    group_by(var, event_type) %>%
    summarise(count = n(), .groups = 'drop') %>%
    arrange(var, event_type)
  print(event_summary)
  
  # Statistics of Crossing Events vs Trend Events
  cross_count <- sum(grepl("cross_", df_cross_all$event_type))
  trend_count <- sum(!grepl("cross_|none", df_cross_all$event_type))
  cat(paste0("\nNumber of crossing events: ", cross_count, " (", round(cross_count/nrow(df_cross_all)*100, 2), "%)\n"))
  cat(paste0("Number of trend events: ", trend_count, " (", round(trend_count/nrow(df_cross_all)*100, 2), "%)\n"))
} else {
  cat("  Warning: No events identified\n")
}

cat("Step 1 Completed!\n\n")




# ==============================================================================
# Data Preprocessing: Group data by period
# ==============================================================================

cat("Data Preprocessing: Grouping data by period...\n")

# Get all unique period labels
unique_periods <- sort(unique(df_transitions$period_label[!is.na(df_transitions$period_label)]))
cat("Found", length(unique_periods), "periods:", paste(unique_periods, collapse=", "), "\n")

# Create data subset for each period
df_by_period <- list()
for(period in unique_periods) {
  df_by_period[[period]] <- df_transitions %>%
    filter(period_label == period)
  cat("Period", period, ":", nrow(df_by_period[[period]]), "records\n")
}

# Create period label mapping (for panel annotation)
period_mapping <- data.frame(
  period_label = unique_periods,
  panel_letter = letters[1:length(unique_periods)]
)

cat("Period label mapping:\n")
print(period_mapping)
cat("\n")

# ==============================================================================
# Visualization 1: Change Frequency Map (5-year stats, column group plot)
# ==============================================================================

cat("=== Visualization 1: Creating Change Frequency Column Group Plot ===\n")
df_area_stats <- data.frame()
plot_list_freq <- list()
for(i in 1:4) {
  period_label <- unique_periods[i]
  panel_letter <- period_mapping$panel_letter[i]
  period_data <- df_by_period[[period_label]]
  
  cat("Processing period", period_label, "...\n")
  
  df_filtered <- filter_valid_data(period_data)
  
  # Calculate change frequency
  all_labs <- unique(df_filtered$lab[!is.na(df_filtered$lab)])
  
  change_freq <- df_filtered %>%
    filter(!is.na(change_id)) %>%
    group_by(lab) %>%
    summarise(max_changes = max(change_id, na.rm = TRUE)) %>%
    ungroup()
  
  complete_freq <- data.frame(lab = all_labs) %>%
    left_join(change_freq, by = "lab") %>%
    mutate(max_changes = ifelse(is.na(max_changes), 0, max_changes)) %>%
    mutate(change_category = case_when(
      max_changes == 0 ~ "No Change",
      max_changes == 1 ~ "1 Change",
      max_changes == 2 ~ "2 Changes", 
      max_changes == 3 ~ "3 Changes",
      max_changes == 4 ~ "4 Changes",
      max_changes == 5 ~ "5 Changes",
      max_changes >= 6 ~ "6+ Changes"
    ))
  
  cat("  - Change frequency stats:\n")
  print(table(complete_freq$change_category))
  
  # Create raster
  r_freq <- saveToRaster(complete_freq$lab, as.numeric(as.factor(complete_freq$change_category)))
  r_freq <- cutByLULC(r_freq)
  
  df_map <- as.data.frame(r_freq, xy = TRUE)
  names(df_map) <- c("x", "y", "category_code")
  df_map <- df_map[!is.na(df_map$category_code), ]
  
  category_labels <- complete_freq %>% 
    dplyr::select(change_category) %>% 
    distinct() %>%
    mutate(category_code = as.numeric(as.factor(change_category)))
  
  df_map <- left_join(df_map, category_labels, by = "category_code")
  
  # Define colors
  fresh_colors <- c("No Change" = "#ed1c24",        
                    "1 Change" = "#00c7f2",         
                    "2 Changes" = "#6aaed6",        
                    "3 Changes" = "#b870ab",         
                    "4 Changes" = "#fb6ca0",        
                    "5 Changes" = "#ffd700",         
                    "6+ Changes" = "#ffff4d")
  
  # Add bar plot
  area_stats <- df_map %>%
    count(change_category) %>%
    mutate(percentage = n / sum(n) * 100) %>%
    mutate(change_category = factor(change_category, 
                                    levels = c("No Change", "1 Change", "2 Changes", 
                                               "3 Changes", "4 Changes", "5 Changes", "6+ Changes")))
  area_stats$period <- period_label
  df_area_stats <- rbind(df_area_stats, area_stats)
  
  bar_plot <- ggplot(area_stats, 
                     aes(x = change_category, 
                         y = percentage, 
                         fill = change_category)) +
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
  
  
  # Create map
  p <- ggplot(df_map) +
    geom_hline(yintercept = 0, color = "gray50", linetype = "dashed", linewidth = 0.3) +
    geom_tile(aes(x = x, y = y, fill = change_category)) +
    geom_sf(data = world_shp, fill = NA, color = "gray30", linewidth = 0.2) +
    coord_sf(ylim = c(-60, 90), expand = FALSE) +
    scale_fill_manual(name = "Transition\nFrequency", values = fresh_colors) +
    labs(title = paste0("(", panel_letter, ") ", period_label)) +
    theme_bw() +
    theme(
      plot.title = element_text(face = "bold", size = 20, 
                                hjust = 0),
      legend.position = "none",
      panel.background = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.title = element_blank(),
      axis.text = element_text(size = 15)
      # plot.margin = margin(5, 5, 5, 5)
    ) +
    annotation_custom(bar_grob, 
                      xmin = -180, 
                      xmax = -85, 
                      ymin = -60, 
                      ymax = 10)
  
  plot_list_freq[[i]] <- p
}

# Create unified legend
legend_plot_freq <- ggplot(df_map) +
  geom_tile(aes(x = x, y = y, fill = change_category)) +
  scale_fill_manual(name = "Transition Frequency",
                    values = fresh_colors) +
  guides(fill = guide_legend(nrow = 2)) +  # Key: Force 2 rows
  theme_void() +
  theme(legend.position = "bottom",
        legend.title = element_text(size = 18, face = "bold"),
        legend.text = element_text(size = 18))


g <- ggplotGrob(legend_plot_freq)
legend_freq <- g$grobs[[which(sapply(g$grobs, function(x) x$name) == "guide-box")]]

# Combine plots
combined_freq <- plot_grid(plotlist = plot_list_freq[1:4], 
                           ncol = 2, align = "v")
final_freq <- plot_grid(combined_freq, 
                        legend_freq, 
                        ncol = 1, 
                        rel_heights = c(0.92, 0.08))

# Save image
ggsave(file.path(figure_date_dir, "Fig1_change_frequency_column_5yr.png"), 
       final_freq, width = 12, height = 8, dpi = 600, bg = "white")
ggsave(file.path(figure_date_dir, "Fig1_change_frequency_column_5yr.pdf"), 
       final_freq, width = 12, height = 8, device = "pdf")
write.csv(df_area_stats, file.path(figure_date_dir, "df_area_stat.csv"))

cat("Change frequency column group plot saved!\n\n")

# ==============================================================================
# Visualization 2: First Change Year Map (with frequency bar chart)
# ==============================================================================

cat("=== Visualization 2: Creating Comprehensive First Change Year Map ===\n")

df_filtered_all <- filter_valid_data(df_transitions)

first_change <- df_filtered_all %>%
  filter(!is.na(change_year)) %>%
  group_by(lab) %>%
  summarise(first_year = min(change_year, na.rm = TRUE)) %>%
  ungroup()

cat("First change statistics:\n")
cat("- Number of pixels with change:", nrow(first_change), "\n")
cat("- First change year range:", min(first_change$first_year), "-", max(first_change$first_year), "\n")

# Create year classification
year_breaks <- quantile(first_change$first_year, probs = c(0, 0.2, 0.4, 0.6, 0.8, 1), na.rm = TRUE)
cat("Year quantiles:", round(year_breaks), "\n")

first_change <- first_change %>%
  mutate(year_period = case_when(
    first_year <= year_breaks[2] ~ paste0("Early (", round(year_breaks[1]), "-", round(year_breaks[2]), ")"),
    first_year <= year_breaks[3] ~ paste0("Early-Mid (", round(year_breaks[2])+1, "-", round(year_breaks[3]), ")"),
    first_year <= year_breaks[4] ~ paste0("Mid (", round(year_breaks[3])+1, "-", round(year_breaks[4]), ")"),
    first_year <= year_breaks[5] ~ paste0("Mid-Late (", round(year_breaks[4])+1, "-", round(year_breaks[5]), ")"),
    TRUE ~ paste0("Late (", round(year_breaks[5])+1, "-", round(year_breaks[6]), ")")
  ))

first_change$period_code <- as.numeric(as.factor(first_change$year_period))
r_year <- saveToRaster(first_change$lab, first_change$period_code)
r_year <- cutByLULC(r_year)

df_map_year <- as.data.frame(r_year, xy = TRUE)
names(df_map_year) <- c("x", "y", "period_code")
df_map_year <- df_map_year[!is.na(df_map_year$period_code), ]

period_labels_year <- first_change %>% 
  dplyr::select(period_code, year_period) %>% 
  distinct()
df_map_year <- left_join(df_map_year, period_labels_year, by = "period_code")

# Define colors
period_colors <- c("#D7191C", "#FDAE61", "#FFFFBF", "#ABD9E9", "#2C7BB6")
period_order <- sort(unique(df_map_year$year_period))
names(period_colors) <- period_order

# Statistics
period_stats <- df_map_year %>%
  count(year_period) %>%
  mutate(percentage = n / sum(n) * 100) %>%
  arrange(year_period)

cat("Year period distribution:\n")
print(period_stats)

# Create bar chart
period_stats$year_period <- factor(period_stats$year_period,
                                   levels = rev(c("Early (1992-1992)",
                                                  "Early-Mid (1993-1994)",
                                                  "Mid (1995-1997)",
                                                  "Mid-Late (1998-2000)",
                                                  "Late (2001-2010)")))

bar_plot_year <- ggplot(period_stats, 
                        aes(x = year_period, 
                            y = percentage)) +
  geom_col(aes(fill = year_period),
           alpha = 0.8, width = 0.7) +
  scale_fill_manual(values = period_colors) +
  labs(x = NULL, y = "Area (%)") +
  coord_flip() +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 12),
    axis.text.y = element_blank(),
    axis.title.x = element_text(size = 12),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    plot.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
    plot.margin = margin(5, 5, 5, 5)
  )

bar_grob_year <- ggplotGrob(bar_plot_year)

# Main map
df_map_year$year_period <- factor(df_map_year$year_period,
                                  levels = c("Early (1992-1992)",
                                             "Early-Mid (1993-1994)",
                                             "Mid (1995-1997)",
                                             "Mid-Late (1998-2000)",
                                             "Late (2001-2010)"))

main_plot_year <- ggplot(df_map_year) +
  geom_hline(yintercept = 0, color = "gray50",
             linetype = "dashed", linewidth = 0.5) +
  geom_tile(aes(x = x, y = y, fill = year_period)) +
  geom_sf(data = world_shp, fill = NA,
          color = "gray30", linewidth = 0.2) +
  coord_sf(ylim = c(-60, 90), expand = FALSE) +
  scale_fill_manual(name = "First Change Period",
                    values = period_colors) +
  guides(fill = guide_legend(nrow = 2)) +  # Add this line
  labs(title = "First Change Year in Vegetation Response (5-year periods)",
       x = "Longitude", y = "Latitude") +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold", size = 22),
    legend.position = "bottom",
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 18),
    legend.title = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 18)
  ) +
  annotation_custom(bar_grob_year, xmin = -180, xmax = -85, ymin = -60, ymax = 10)

# Save image
ggsave(file.path(figure_date_dir, "Fig2_first_change_year_comprehensive_5yr.png"), 
       main_plot_year, width = 11, height = 6, dpi = 600, bg = "white")
ggsave(file.path(figure_date_dir, "Fig2_first_change_year_comprehensive_5yr.pdf"), 
       main_plot_year, width = 11, height = 6, device = "pdf")

cat("First change year map saved!\n\n")

# ==============================================================================
# Visualization 3: Persistent Type Visualization (Column group plot)
# ==============================================================================

cat("=== Visualization 3: Creating Persistent Type Column Group Plot ===\n")

plot_list_persistent <- list()

for(i in 1:4) {
  period_label <- unique_periods[i]
  panel_letter <- period_mapping$panel_letter[i]
  period_data <- df_by_period[[period_label]]
  
  cat("Processing period", period_label, "...\n")
  
  persistent_data <- period_data %>%
    filter_valid_data() %>%
    filter(change_id == 0 & 
             (note == "no_changed_point" | note == "single_data_point") &
             !is.na(persistent_type))
  persistent_data$persistent_type <- factor(x = persistent_data$persistent_type, 
                                            levels = c("VGC","Stable","VSO"),
                                            labels = c("δVGC","Stable","δVSO"))
  
  
  if(nrow(persistent_data) == 0) {
    cat("  Warning: No persistent type data for this period\n")
    next
  }
  
  cat("  - Persistent pixels:", nrow(persistent_data), "\n")
  cat("  - Type distribution:\n")
  print(table(persistent_data$persistent_type))
  
  persistent_data$type_code <- as.numeric(as.factor(persistent_data$persistent_type))
  r_persistent <- saveToRaster(persistent_data$lab, persistent_data$type_code)
  r_persistent <- cutByLULC(r_persistent)
  
  df_map_pers <- as.data.frame(r_persistent, xy = TRUE)
  names(df_map_pers) <- c("x", "y", "type_code")
  df_map_pers <- df_map_pers[!is.na(df_map_pers$type_code), ]
  
  type_labels <- persistent_data %>% 
    dplyr::select(type_code, persistent_type) %>% 
    distinct()
  df_map_pers <- left_join(df_map_pers, type_labels, by = "type_code")
  df_map_pers$persistent_type <- factor(x = df_map_pers$persistent_type, 
                                        levels = c("δVGC","Stable","δVSO"),
                                        labels = c("δVGC","Stable","δVSO"))
  
  persistent_colors <- c(
    "Stable" = "#56a0d3",
    "δVGC" = "#DA5CA3",
    "δVSO" = "#4DAF4A"
  )
  
  stats_data <- persistent_data %>%
    count(persistent_type) %>%
    mutate(
      percentage = n / sum(n) * 100
    ) %>%
    arrange(desc(n))
  
  count_plot <- ggplot(stats_data, 
                       aes(x = reorder(persistent_type, n), 
                           y = n, fill = persistent_type)) +
    geom_col(alpha = 0.8, color = "white", linewidth = 0.5) +
    geom_text(aes(label = scales::comma(n)),
              hjust = -0.1, size = 3, fontface = "bold") +
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
  
  area_plot <- ggplot(stats_data, 
                      aes(x = reorder(persistent_type, percentage), 
                          y = percentage, fill = persistent_type)) +
    geom_col(alpha = 0.8, color = "white", linewidth = 0.5) +
    geom_text(aes(label = sprintf("%.1f%%", percentage)), 
              hjust = -0.05, size = 4.5, fontface = "bold") +
    scale_fill_manual(values = persistent_colors) +
    coord_flip() +
    labs(title = NULL,
         x = NULL, 
         y = NULL) +
    theme_bw() +
    theme(
      legend.position = "none",
      panel.background = element_blank(),  
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.title = element_text(size = 9),
      axis.ticks.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.x = element_blank(),
      axis.text = element_text(size = 12)
      # panel.border = element_blank()
      # axis.line.x = element_line(size = 0.5, linetype = "solid", colour = "black"),
      # axis.line.y = element_line(size = 0.5, linetype = "solid", colour = "black")
    ) +
    scale_y_continuous(limits = c(0, 100),
                       expand = expansion(mult = c(0, 0.15)))
  
  combined_stats_grob <- ggplotGrob(area_plot)
  
  
  p <- ggplot(df_map_pers) +
    geom_hline(yintercept = 0, color = "gray50", linetype = "dashed", linewidth = 0.3) +
    geom_tile(aes(x = x, y = y, fill = persistent_type)) +
    geom_sf(data = world_shp, fill = NA, color = "gray30", linewidth = 0.2) +
    coord_sf(ylim = c(-60, 90), expand = FALSE) +
    scale_fill_manual(name = "Persistent\nType",
                      values = persistent_colors) +
    labs(title = paste0("(", panel_letter, ") ", period_label)) +
    theme_bw() +
    theme(
      plot.title = element_text(face = "bold",
                                size = 20, hjust = 0),
      legend.position = "none",
      panel.background = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.title = element_blank(),
      axis.text = element_text(size = 15)
    ) +
    annotation_custom(combined_stats_grob, 
                      # xmin = -180, 
                      # xmax = -82, 
                      # ymin = -60, 
                      # ymax = 8)
                      xmin = -180,
                      xmax = -82,
                      ymin = -60,
                      ymax = 8)
  p
  plot_list_persistent[[i]] <- p
}

# Create unified legend
legend_plot_pers <- ggplot(df_map_pers) +
  geom_tile(aes(x = x, y = y, fill = persistent_type)) +
  scale_fill_manual(name = "Persistent Type",
                    values = persistent_colors) +
  guides(fill = guide_legend(nrow = 1)) +  # Key: Force 1 row
  theme_void() +
  theme(legend.position = "bottom",
        legend.title = element_text(size = 18, face = "bold"),
        legend.text = element_text(size = 18))


g <- ggplotGrob(legend_plot_pers)
legend_pers <- g$grobs[[which(sapply(g$grobs, function(x) x$name) == "guide-box")]]

# Combine plots
combined_pers <- plot_grid(plotlist = plot_list_persistent, ncol = 1, align = "v")
final_pers <- plot_grid(combined_pers, legend_pers, 
                        ncol = 1, rel_heights = c(0.92, 0.08))

# Save image
ggsave(file.path(figure_date_dir, "Fig3_persistent_type_column_5yr.png"), 
       final_pers, width = 6, height = 14, dpi = 600, bg = "white")
ggsave(file.path(figure_date_dir, "Fig3_persistent_type_column_5yr.pdf"), 
       final_pers, width = 6, height = 14, device = "pdf")

cat("Persistent Type column group plot saved!\n\n")

# ==============================================================================
# Visualization 4: Change Type Spatial Distribution Map (Column group plot)
# ==============================================================================

cat("=== Visualization 4: Creating Change Type Column Group Plot ===\n")

df_area_change_stat <- data.frame()
plot_list_change_type <- list()
panel_letters <- letters[5:8]
for(i in 1:4) {
  period_label <- unique_periods[i]
  panel_letter <- panel_letters[i]
  period_data <- df_by_period[[period_label]]
  
  cat("Processing period", period_label, "...\n")
  
  df_filtered <- filter_valid_data(period_data)
  
  change_type_data <- df_filtered %>%
    filter(!is.na(from) & !is.na(to) & from != to) %>%
    mutate(change_type = paste0(from, "->", to)) %>%
    group_by(lab) %>%
    slice(1) %>%
    ungroup()
  
  if(nrow(change_type_data) == 0) {
    cat("  Warning: No change type data for this period\n")
    next
  }
  
  cat("  - Changed pixels:", nrow(change_type_data), "\n")
  cat("  - Change type distribution:\n")
  print(table(change_type_data$change_type))
  
  # Calculate area percentage
  area_stats <- change_type_data %>%
    group_by(change_type) %>%
    summarise(count = n()) %>%
    mutate(
      percentage = count / sum(count) * 100,
      change_type = factor(change_type,
                           levels = c("VGC->VSO", "VGC->Stable",
                                      "VSO->Stable", "VSO->VGC",
                                      "Stable->VGC", "Stable->VSO"),
                           labels = c("δVGC->δVSO", "δVGC->Stable",
                                      "δVSO->Stable", "δVSO->δVGC",
                                      "Stable->δVGC", "Stable->δVSO"))
    ) %>%
    arrange(change_type)
  df_area_change_stat <- rbind(df_area_change_stat, area_stats)
  
  change_type_data$type_code <- as.numeric(as.factor(change_type_data$change_type))
  r_change <- saveToRaster(change_type_data$lab, 
                           change_type_data$type_code)
  r_change <- cutByLULC(r_change)
  
  df_map_change <- as.data.frame(r_change, xy = TRUE)
  names(df_map_change) <- c("x", "y", "type_code")
  df_map_change <- df_map_change[!is.na(df_map_change$type_code), ]
  
  type_labels_change <- change_type_data %>% 
    dplyr::select(type_code, change_type) %>% 
    distinct()
  df_map_change <- left_join(df_map_change,
                             type_labels_change, by = "type_code")
  
  df_map_change$change_type <-factor(df_map_change$change_type,
                                     levels = c("VGC->VSO", "VGC->Stable",
                                                "VSO->Stable", "VSO->VGC",
                                                "Stable->VGC", "Stable->VSO"),
                                     labels = c("δVGC->δVSO", "δVGC->Stable",
                                                "δVSO->Stable", "δVSO->δVGC",
                                                "Stable->δVGC", "Stable->δVSO"))
  
  change_colors <- c(
    "Stable->δVGC" = "#2E8B57",
    "Stable->δVSO" = "#48D1CC",
    "δVGC->Stable" = "#DAA520",
    "δVGC->δVSO" = "#CD853F",
    "δVSO->Stable" = "#4682B4",
    "δVSO->δVGC" = "#9370DB"
  )
  
  # Create area percentage bar chart (similar style to area_plot)
  area_plot <- ggplot(area_stats, 
                      aes(x = reorder(change_type, percentage), 
                          y = percentage, 
                          fill = change_type)) +
    geom_col(alpha = 0.8, color = "white", linewidth = 0.5) +
    geom_text(aes(label = sprintf("%.1f%%", percentage)),
              hjust = -0.05, size = 4.5, fontface = "bold") +
    scale_fill_manual(values = change_colors) +
    coord_flip() +
    labs(title = NULL,
         x = NULL, 
         y = NULL) +
    theme_bw() +
    theme(
      legend.position = "none",
      panel.background = element_blank(),  
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.title = element_text(size = 9),
      axis.ticks.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.x = element_blank(),
      axis.text = element_text(size = 8)
      # panel.border = element_blank()
      # axis.line.x = element_line(size = 0.5, linetype = "solid", colour = "black"),
      # axis.line.y = element_line(size = 0.5, linetype = "solid", colour = "black")
    ) +
    scale_y_continuous(limits = c(0, 100),
                       expand = expansion(mult = c(0, 0.15)))
  
  # Convert to grob object
  combined_stats_grob <- ggplotGrob(area_plot)
  
  p <- ggplot(df_map_change) +
    geom_hline(yintercept = 0, 
               color = "gray50", linetype = "dashed", 
               linewidth = 0.3) +
    geom_tile(aes(x = x, y = y, fill = change_type)) +
    geom_sf(data = world_shp, fill = NA, color = "gray30",
            linewidth = 0.2) +
    coord_sf(ylim = c(-60, 90), expand = FALSE) +
    scale_fill_manual(name = "Transition Type", 
                      values = change_colors, 
                      na.value = "transparent") +
    labs(title = paste0("(", panel_letter, ") ",
                        period_label)) +
    theme_bw() +
    theme(
      plot.title = element_text(face = "bold",
                                size = 20, hjust = 0),
      legend.position = "none",
      panel.background = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.title = element_blank(),
      axis.text = element_text(size = 15)
    ) +
    annotation_custom(combined_stats_grob, 
                      # xmin = -180, 
                      # xmax = -108, 
                      # ymin = -60, 
                      # ymax = 22)
                      
                      xmin = -180,
                      xmax = -82,
                      ymin = -60,
                      ymax = 8)
  p
  plot_list_change_type[[i]] <- p
}

# Create unified legend
legend_plot_change <- ggplot(df_map_change) +
  geom_tile(aes(x = x, y = y, fill = change_type)) +
  scale_fill_manual(name = "Transition Type",
                    values = change_colors) +
  guides(fill = guide_legend(nrow = 3)) +  # Key: Force 3 rows
  theme_void() +
  theme(legend.position = "bottom",
        legend.title = element_text(size = 18, face = "bold"),
        legend.text = element_text(size = 18))


g <- ggplotGrob(legend_plot_change)
legend_change <- g$grobs[[which(sapply(g$grobs, function(x) x$name) == "guide-box")]]

# Combine plots
combined_change <- plot_grid(plotlist = plot_list_change_type[1:4], 
                             ncol = 1, align = "v")
final_change <- plot_grid(combined_change, 
                          legend_change,
                          ncol = 1, 
                          rel_heights = c(0.92, 0.08))

# Save image
ggsave(file.path(figure_date_dir, "Fig4_change_type_column_5yr.png"), 
       final_change, width = 6, height = 14, dpi = 600, bg = "white")
ggsave(file.path(figure_date_dir, "Fig4_change_type_column_5yr.pdf"), 
       final_change, width = 6, height = 14, device = "pdf")
write.csv(df_area_change_stat, 
          file.path(figure_date_dir, "Fig4_change_type_stat.csv"))

cat("Change type column group plot saved!\n\n")

# ==============================================================================
# Visualization 5: Significance Test Analysis and Heatmap by Year Period
# ==============================================================================

cat("=== Visualization 5: Significance Test by Year Period ===\n")

# Prepare event matrix data
cat("Preparing transition event data...\n")

library(dplyr)
library(tidyr)

# Define climate factors
vars <- c("CO2", "TMP", "SR", "VPD", "PRE", "SMroot", "SMsurf")

# Convert crossing events to wide format
df_cross_wide <- df_cross_all %>%
  dplyr::select(lab, year, var, cross_dir) %>%
  tidyr::pivot_wider(
    names_from = var,
    values_from = cross_dir,
    names_prefix = "Cross_"
  )

# Store test results for all periods
test_results_by_period <- list()
i = 4
# Loop through each period
for(i in 1:4) {
  period_label <- unique_periods[i]
  period_data <- df_by_period[[period_label]]
  
  cat("Processing period", period_label, "...\n")
  
  # Process change events and persistent types simultaneously
  df_combined_events <- period_data %>%
    mutate(
      # Assign change_dir for transition events
      change_dir = case_when(
        !is.na(change_year) & from == "Stable" & to == "VGC" ~ 1,    
        !is.na(change_year) & from == "VGC" & to == "Stable" ~ 2,    
        !is.na(change_year) & from == "VSO" & to == "VGC" ~ 3,       
        !is.na(change_year) & from == "VGC" & to == "VSO" ~ 4,
        # Assign corresponding change_dir for persistent types
        !is.na(persistent_type) & persistent_type == "VGC" ~ 13,  # VGC Persistent - matches both 1 and 3
        !is.na(persistent_type) & persistent_type == "VSO" ~ 4,   # VSO Persistent - matches 4
        !is.na(persistent_type) & persistent_type == "Stable" ~ 2, # Stable Persistent - matches 2
        TRUE ~ 0
      ),
      # Assign groups
      change_group = case_when(
        change_dir %in% c(1, 2) ~ "Energy",
        change_dir %in% c(3, 4) ~ "Water",
        change_dir == 13 ~ "Both",  # VGC persistent type belongs to both groups
        TRUE ~ NA_character_
      ),
      # Create event type label
      event_type = case_when(
        !is.na(change_year) ~ paste0(from, "->", to, " (Change)"),
        !is.na(persistent_type) ~ paste0(persistent_type, " (Persistent)"),
        TRUE ~ NA_character_
      )
    ) %>%
    filter(change_dir != 0)
  
  cat("  Combined event stats:\n")
  cat("  Total events:", nrow(df_combined_events), "\n")
  if(nrow(df_combined_events) > 0) {
    print(table(df_combined_events$event_type))
    cat("  Stats by group:\n")
    print(table(df_combined_events$change_group))
  }
  
  df_cross_wide$lab <- as.character(df_cross_wide$lab)
  df_combined_events$lab <- as.character(df_combined_events$lab)
  print(nrow(df_cross_wide))
  print(nrow(df_combined_events))
  
  # Find intersection
  common_labs <- intersect(unique(df_cross_wide$lab), unique(df_combined_events$lab))
  cat(paste0("  Intersecting pixels: ", length(common_labs), "\n"))
  
  df_cross_wide <- df_cross_wide %>% filter(lab %in% common_labs)
  df_combined_events <- df_combined_events %>% filter(lab %in% common_labs)
  
  lags <- 0
  event_matrix_list <- list()
  
  for(lag in lags) {
    cat(paste0("  Processing lag: ", lag, " years\n"))
    
    matched_events <- df_combined_events %>%
      left_join(
        df_cross_wide %>% mutate(match_year = year + lag),
        by = c("lab" = "lab", "change_year" = "match_year")
      ) %>%
      mutate(lag = lag)
    
    event_matrix_list[[lag + 1]] <- matched_events
  }
  
  event_matrix <- bind_rows(event_matrix_list)
  print(nrow(event_matrix))
  
  # Core modification: Group calculation of directional consistency (1234 coding, direct matching)
  # Assuming event_matrix already has:
  # change_dir, change_group, Cross_CO2, Cross_TMP, Cross_SR,
  # Cross_VPD, Cross_PRE, Cross_SMroot, Cross_SMsurf
  
  event_matrix <- event_matrix %>%
    mutate(
      match_CO2      = Cross_CO2      == change_dir,
      match_TMP      = Cross_TMP      == change_dir,
      match_SR       = Cross_SR       == change_dir,
      match_VPD      = Cross_VPD      == change_dir,
      match_PRE      = Cross_PRE      == change_dir,
      match_SMroot   = Cross_SMroot   == change_dir,
      match_SMsurf   = Cross_SMsurf   == change_dir
    )
  
  event_matrix <- event_matrix %>%
    mutate(
      match_energy = if_any(c(match_CO2, match_TMP, match_SR), ~ .x == TRUE),
      match_water  = if_any(c(match_VPD, match_PRE, match_SMroot, match_SMsurf), ~ .x == TRUE)
    )
  
  # By climate factor
  test_results <- data.frame()
  unique(event_matrix$event_type)
  
  for (v in vars) {
    col_cross <- paste0("Cross_", v)
    col_match <- paste0("match_", v)
    
    if(v %in% c("CO2", "TMP", "SR")){
      df <- event_matrix %>% 
        subset(event_type %in% c("VGC->Stable (Change)",
                                 "Stable->VGC (Change)"))
    }else{
      df <- event_matrix %>% 
        subset(event_type %in% c("VGC->VSO (Change)",
                                 "VSO->VGC (Change)"))
    }
    
    df <- df %>%
      filter(!!sym(col_cross) != 0)
    
    n <- nrow(df)
    k <- sum(df[[col_match]], na.rm = TRUE)
    
    test <- prop.test(k, n, p = 0.25)
    
    test_results <- rbind(test_results, data.frame(
      factor = v,
      period = period_label,
      n = n,
      match_rate = test$estimate,
      CI_low = test$conf.int[1],
      CI_high = test$conf.int[2],
      p_value = test$p.value
    ))
  }
  
  # ---- match_energy ----
  df_energy <- event_matrix %>%
    subset(event_type %in% c("Stable (Persistent)",
                             "VGC->Stable (Change)",
                             "Stable->VGC (Change)"))
  
  n_energy <- nrow(df_energy)
  k_energy <- sum(df_energy$match_energy, na.rm = TRUE)
  
  test_energy <- prop.test(k_energy, n_energy, p = 0.25)
  
  test_results <- rbind(test_results, data.frame(
    factor = "Energy_factors",
    period = unique(df_energy$period_label),
    n = n_energy,
    match_rate = test_energy$estimate,
    CI_low = test_energy$conf.int[1],
    CI_high = test_energy$conf.int[2],
    p_value = test_energy$p.value
  ))
  
  # ---- match_water ----
  df_water <- event_matrix %>%
    subset(event_type %in% c("VGC->VSO (Change)",
                             "VSO->VGC (Change)"))
  
  n_water <- nrow(df_water)
  k_water <- sum(df_water$match_water, na.rm = TRUE)
  
  test_water <- prop.test(k_water, n_water, p = 0.25)
  
  test_results <- rbind(test_results, data.frame(
    factor = "Water_factors",
    period = unique(df_water$period_label),
    n = n_water,
    match_rate = test_water$estimate,
    CI_low = test_water$conf.int[1],
    CI_high = test_water$conf.int[2],
    p_value = test_water$p.value
  ))
  
  print(test_results)
  # write.csv(test_results, file.path(Datedir, "Step4_SignificanceTest_all.csv"), row.names = FALSE)
  
  test_results_by_period[[i]] <- test_results
  
  
  head(event_matrix)
}

# Merge results from all periods
df_test_by_period <- bind_rows(test_results_by_period)
df111 <- df_test_by_period %>% subset(factor == "Water_factors") %>%
  dplyr::select(match_rate)
mean(df111$match_rate)

# Save significance test results
write.csv(df_test_by_period, 
          file.path(Datedir, "Significance_Test_by_Period.csv"), 
          row.names = FALSE)

cat("Significance test results saved to:", 
    file.path(Datedir, "Significance_Test_by_Period.csv"), "\n\n")


# Create significance test heatmap
cat("\nCreating significance test heatmap...\n")

df_plot_sig <- df_test_by_period %>%
  mutate(
    # Add significance judgment (based on p-value)
    significant = p_value < 0.05,
    
    # Create period label (if not exists)
    period_label = period,
    
    # Create factor group
    factor_group = case_when(
      factor %in% c("Energy_factors", "CO2", "TMP", "SR") ~ "Energy",
      factor %in% c("Water_factors", "VPD", "PRE", "SMroot", "SMsurf") ~ "Water",
      TRUE ~ "Other"
    ),
    
    # Create significance label
    sig_label = case_when(
      is.na(match_rate) ~ "NA",
      p_value < 0.001 ~ paste0(round(match_rate * 100, 1), "%***"),
      p_value < 0.01 ~ paste0(round(match_rate * 100, 1), "%**"),
      p_value < 0.05 ~ paste0(round(match_rate * 100, 1), "%*"),
      TRUE ~ paste0(round(match_rate * 100, 1), "%")
    ),
    
    # Factor sorting
    factor = factor(factor, 
                    levels = c("Energy_factors", "TMP", "SR", "CO2", 
                               "Water_factors", "VPD", "PRE", "SMroot", "SMsurf")),
    
    # Period sorting (chronological reverse for y-axis display)
    period_label = factor(period_label, 
                          levels = rev(c("1991-1995", "1996-2000", "2001-2005", "2006-2010"))),
    
    # Add panel labels for facets
    factor_group = factor(factor_group, 
                          levels = c("Energy", "Water"),
                          labels = c("(a)", "(b)"))
  )

# Create Heatmap
df_plot_sig$period <- factor(df_plot_sig$period,
                             levels = c("1991-1995", "1996-2000", "2001-2005", "2006-2010"))

df_plot_sig$factor <- factor(df_plot_sig$factor,
                             levels = c("Energy_factors", "TMP", "SR", "CO2", 
                                        "Water_factors", "VPD", "PRE", "SMroot", "SMsurf"))


p_heatmap <- ggplot(df_plot_sig, 
                    aes(x = factor, 
                        y = period_label, 
                        fill = match_rate)) +
  geom_tile(color = "white", linewidth = 0.8) +
  geom_text(aes(label = sig_label), 
            size = 4.5,
            color = "black",
            fontface = "bold") +
  scale_fill_gradientn(
    colors = c("#FAD689", "#f3eed9","#F07881"),  # Blue-White-Red gradient
    name = "Match rate",
    labels = scales::percent_format(accuracy = 1)
    # limits = c(0.45, 0.55)
  ) +
  facet_wrap(~factor_group, scales = "free_x", nrow = 1) +
  labs(
    x = NULL,
    y = NULL
  ) +
  theme_bw() +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = 18, face = "bold", hjust = 0),
    axis.text.x = element_text(size = 15, color = "black", angle = 0),
    axis.text.y = element_text(size = 15, color = "black"),
    axis.title = element_text(size = 15, face = "bold"),
    legend.position = "right",
    legend.title = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 12),
    legend.key.height = unit(1, "cm"),
    legend.key.width = unit(0.5, "cm"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 0.5)
  )
p_heatmap

# Save significance heatmap
ggsave(file.path(figure_date_dir, "Fig5_significance_heatmap_by_period.png"), 
       p_heatmap, width = 10, height = 5, dpi = 600, bg = "white")
ggsave(file.path(figure_date_dir, "Fig5_significance_heatmap_by_period.pdf"), 
       p_heatmap, width = 10, height = 5, device = "pdf")

cat("Heatmap saved!\n")

# ==============================================================
# Create Line Plot
# ==============================================================

cat("\nCreating line plot...\n")

# Prepare data
df_plot_line <- df_test_by_period %>%
  mutate(
    # Add significance judgment
    significant = p_value < 0.05,
    
    # Create factor group
    factor_group = case_when(
      factor %in% c("Energy_factors", "Water_factors") ~ "ALL",
      factor %in% c("CO2", "TMP", "SR") ~ "Energy",
      factor %in% c("VPD", "PRE", "SMroot", "SMsurf") ~ "Water",
      TRUE ~ "Other"
    ),
    
    # Factor sorting
    factor = factor(factor, 
                    levels = c("Energy_factors", 
                               "Water_factors", 
                               "TMP", "SR", "CO2", 
                               "VPD", "PRE", "SMroot", "SMsurf"),
                    labels = c("Energy-related factors", 
                               "Moisture-related factors", 
                               "TMP", "SR", "CO2", 
                               "VPD", "PRE", "SMroot", "SMsurf")),
    
    # Period sorting (chronological)
    period = factor(period, 
                    levels = c("1991-1995", "1996-2000", "2001-2005", "2006-2010")),
    
    # Add panel labels for facets
    lab_title = factor(factor_group, 
                       levels = c("ALL", "Energy", "Water"),
                       labels = c("(i) Combination of energy and\nwater limiting factors", 
                                  "(j) Energy-related factors", 
                                  "(k) Moisture-related factors"))
  )


# Add panel labels for facets
df_hline <- data.frame(
  lab_title = factor(c("(j)", "(k)"), 
                     levels = c("(j)", 
                                "(k)"),
                     labels = c("(j) Energy-related factors", 
                                "(k) Moisture-related factors")),
  yintercept = 50
)

df_hline_2 <- data.frame(
  lab_title = factor(c("(i)"),
                     levels = c("(i)"),
                     labels = c("(i) Combination of energy and\nwater limiting factors")
  ),
  yintercept = 87.5
)


p_lineplot_all <- df_plot_line %>%
  ggplot(aes(x = period, 
             y = match_rate * 100,  # Convert to percentage
             group = factor,
             color = factor,
             shape = factor)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 4, fill = "white", stroke = 1.5) +
  geom_hline(data = df_hline, 
             aes(yintercept = yintercept),
             linetype = "dashed", 
             color = "gray50", 
             linewidth = 0.8) +
  geom_hline(data = df_hline_2, 
             aes(yintercept = yintercept),
             linetype = "dashed", 
             color = "gray50", 
             linewidth = 0.8) +
  
  # Color scheme
  scale_color_manual(
    values = c("Energy-related factors" = "#FAD689",
               "Moisture-related factors" = "#F07881",
               "TMP" = "#E41A1C", "SR" = "#377EB8", "CO2" = "#4DAF4A",
               "VPD" = "#984EA3", "PRE" = "#FF7F00",
               "SMroot" = "#A65628", "SMsurf" = "#F781BF"),
    name = "Climate factor"
  ) +
  # Point shapes
  scale_shape_manual(
    values = c("Energy-related factors" = 16, 
               "Moisture-related factors" = 17,
               "TMP" = 16, "SR" = 17, "CO2" = 15,
               "VPD" = 16, "PRE" = 17, "SMroot" = 15, "SMsurf" = 18),
    name = "Climate factor"
  ) +
  
  # Axis settings
  scale_y_continuous(
    name = "Match rate (%)"
  ) +
  facet_wrap(~lab_title, scales = "free", nrow = 1) +
  
  labs(x = NULL) +
  
  theme_bw() +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = 22, face = "bold", hjust = 0),
    axis.text.x = element_text(size = 18, color = "black", 
                               angle = 45, hjust = 1),
    axis.text.y = element_text(size = 18, color = "black"),
    axis.title.y = element_text(size = 18, face = "bold"),
    legend.position = "bottom",
    legend.title = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 15),
    legend.key.height = unit(1, "cm"),
    legend.key.width = unit(1.5, "cm"),
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

p_lineplot_all

# Save plot
ggsave(file.path(figure_date_dir, "Fig5_final5.png"), 
       p_lineplot_all, width = 20, height = 5, dpi = 600, bg = "white")
ggsave(file.path(figure_date_dir, "Fig5_final5.pdf"), 
       p_lineplot_all, width = 20, height = 5, device = "pdf")


cat("Line plot saved!\n")


# == SD Raw Data and Classification Data ==============================================

# Extract climate variable name
extract_climate_var <- function(path) {
  basename(path) %>%
    gsub("GPP_", "", .) %>%
    gsub("_SDdata.RData", "", .)
}

# Classification function
classify_by_threshold <- function(sd_cli,
                                  climate_var, 
                                  df_threshold) {
  # Get threshold info for this climate variable
  threshold_info <- df_threshold %>%
    filter(Climate_Variable == climate_var)
  
  if(nrow(threshold_info) == 0) {
    return(NA)
  }
  
  if(climate_var %in% single_threshold_vars) {
    # Single threshold: 2 classes
    threshold <- threshold_info$X_Position[1]
    if(threshold_info$Transition_type[1] == "VGC-VSO") {
      return(ifelse(sd_cli < threshold, 1, 2))
    } else if(threshold_info$Transition_type[1] == "VGC-Stable") {
      return(ifelse(sd_cli < threshold, 1, 2)) 
    }
  } else if(climate_var %in% double_threshold_vars) {
    # Double threshold: 3 classes
    if(nrow(threshold_info) >= 2) {
      thresholds <- sort(threshold_info$X_Position)
      if(sd_cli < thresholds[1]) {
        return(1) 
      } else if(sd_cli <= thresholds[2]) {
        return(2)
      } else {
        return(3)
      }
    }
  }
  
  return(NA)
}

# Store all results
all_results <- list()

# == Step 3: Process each climate factor individually =========================================
sd_path <- SD_paths[1]
for(sd_path in SD_paths) {
  
  climate_var <- extract_climate_var(sd_path)
  cat("\nProcessing climate factor:", climate_var, "\n")
  cat("File path:", sd_path, "\n")
  
  # Load SD data
  load(sd_path)
  
  # Get all pixel IDs
  pixel_ids <- as.numeric(names(SD_list))
  
  # Get year range
  years <- SD_list[[1]]$year_new
  cat("Year range:", min(years), "-", max(years), "\n")
  cat("Number of pixels:", length(pixel_ids), "\n")
  
  # Create classified rasters for each year
  yearly_rasters_original <- list()
  yearly_rasters_classified <- list()
  
  for(year in years) {
    # Extract data for this year
    year_data <- lapply(names(SD_list), function(id) {
      df <- SD_list[[id]]
      year_row <- df[df$year_new == year, ]
      if(nrow(year_row) > 0) {
        data.frame(
          pixel_id = as.numeric(id),
          sd_cli = year_row$sd_cli,
          sd_GS = year_row$sd_GS,
          stringsAsFactors = FALSE
        )
      } else {
        NULL
      }
    }) %>% bind_rows()
    
    # Classification
    year_data$class <- sapply(year_data$sd_cli, function(x) {
      classify_by_threshold(x, climate_var, df_trippingtype)
    })
    
    # Create original value raster
    r_original <- saveToRaster(
      labs = year_data$pixel_id,
      values = year_data$sd_cli
    )
    plot(r_original)
    
    # Create classified raster
    r_classified <- saveToRaster(
      labs = year_data$pixel_id,
      values = year_data$class
    )
    
    yearly_rasters_original[[as.character(year)]] <- r_original
    yearly_rasters_classified[[as.character(year)]] <- r_classified
    
    # Save raster files
    writeRaster(r_original, 
                file.path(Datedir, paste0(climate_var, "_", year, "_original.tif")),
                overwrite = TRUE)
    writeRaster(r_classified, 
                file.path(Datedir, paste0(climate_var, "_", year, "_classified.tif")),
                overwrite = TRUE)
  }
  
  cat("Raster saving completed!\n")
  
  # Store results
  all_results[[climate_var]] <- list(
    original = yearly_rasters_original,
    classified = yearly_rasters_classified,
    years = years
  )
}
# load(SD_paths[1])
# SD_list[1]

raster('./Out-Table-and-Figure-VPD et al/2025-11-04/Comprehensive-Visualization-5yr/CO2_1982_classified.tif')
raster('./Out-Table-and-Figure-VPD et al/2025-11-04/Comprehensive-Visualization-5yr/CO2_1982_original.tif')
# == Step 4: Create Visualization Charts ======================================================
cat("\n\nStep 4: Creating visualization charts...\n")

# Define color schemes
colors_3class <- c("1" = "#2166ac", "2" = "#d6604d", "3" = "#92c5de")
labels_3class <- c("1" = "< T1", "2" = "T1 ≥ & < T2", "3" = "> T2")

colors_2class_vgc_vso <- c("1" = "#2166ac", "2" = "#d6604d")
labels_2class_vgc_vso <- c("1" = "< T", "2" = "> T")

colors_2class_vgc_stable <- c("1" = "#2166ac", "3" = "#92c5de")
labels_2class_vgc_stable <- c("1" = "< T", "2" = "> T")

# Create single year map
create_year_map <- function(raster_obj, year, type = "original", climate_var = NULL, 
                            colors = NULL, labels = NULL) {
  
  # Convert to data frame
  df_plot <- as.data.frame(rasterToPoints(raster_obj))
  colnames(df_plot) <- c("x", "y", "value")
  
  # Create base map
  p <- ggplot() +
    geom_raster(data = df_plot, aes(x = x, y = y, fill = value)) +
    geom_sf(data = st_as_sf(world_shp), fill = NA, color = "gray30", size = 0.2) +
    coord_sf(expand = FALSE) +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      axis.title = element_blank(),
      axis.text = element_text(size = 6),
      plot.title = element_text(size = 8, face = "bold", hjust = 0.5),
      legend.position = "none",
      plot.margin = margin(2, 2, 2, 2)
    ) +
    labs(title = year)
  
  # Set color by type
  if(type == "original") {
    p <- p + scale_fill_viridis_c(option = "D", na.value = "transparent")
  } else {
    # Classified map
    df_plot$value <- factor(df_plot$value)
    p <- ggplot() +
      geom_raster(data = df_plot, aes(x = x, y = y, fill = value)) +
      geom_sf(data = st_as_sf(world_shp), fill = NA, color = "gray30", size = 0.2) +
      coord_sf(expand = FALSE) +
      scale_fill_manual(values = colors, labels = labels, na.value = "transparent") +
      theme_minimal() +
      theme(
        panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_text(size = 6),
        plot.title = element_text(size = 8, face = "bold", hjust = 0.5),
        legend.position = "none",
        plot.margin = margin(2, 2, 2, 2)
      ) +
      labs(title = year)
  }
  
  return(p)
}

# Create combined charts for each climate factor
for(climate_var in names(all_results)) {
  
  cat("Plotting chart for", climate_var, "...\n")
  
  result <- all_results[[climate_var]]
  years <- result$years
  
  # Determine color scheme
  threshold_info <- df_trippingtype %>%
    filter(Climate_Variable == climate_var)
  
  if(climate_var %in% double_threshold_vars) {
    colors <- colors_3class
    labels <- labels_3class
  } else {
    if(threshold_info$Transition_type[1] == "VGC-VSO") {
      colors <- colors_2class_vgc_vso
      labels <- labels_2class_vgc_vso
    } else {
      colors <- colors_2class_vgc_stable
      labels <- labels_2class_vgc_stable
    }
  }
  
  # Create original value maps (5 rows 8 columns = 40 plots)
  plots_original <- lapply(years, function(year) {
    create_year_map(result$original[[as.character(year)]], 
                    year, 
                    type = "original")
  })
  
  # Create classified maps
  plots_classified <- lapply(years, function(year) {
    create_year_map(result$classified[[as.character(year)]], 
                    year, 
                    type = "classified", 
                    climate_var = climate_var, 
                    colors = colors, 
                    labels = labels)
  })
  
  # Combine charts - Original values
  cat("  - Combining original value charts...\n")
  
  # Create legend
  legend_original <- ggplot() +
    geom_blank() +
    scale_fill_viridis_c(option = "D", name = paste0(climate_var, "\n(Original)")) +
    theme_minimal() +
    theme(legend.position = "right")
  
  legend_grob_original <- cowplot::get_legend(legend_original)
  
  # Combine all years
  combined_original <- wrap_plots(plots_original, ncol = 8, nrow = 5) +
    plot_annotation(
      title = paste0(climate_var, " - Original Climate Values (1982-2020)"),
      theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))
    )
  
  # Save
  ggsave(
    filename = file.path(figure_date_dir, paste0(climate_var, "_Original_Values.png")),
    plot = combined_original,
    width = 24,
    height = 15,
    dpi = 300,
    bg = "white"
  )
  
  # Combine charts - Classified results
  cat("  - Combining classified charts...\n")
  
  # Create legend
  legend_df <- data.frame(
    x = 1,
    y = 1:length(colors),
    class = factor(names(colors), levels = names(colors))
  )
  
  legend_classified <- ggplot(legend_df, aes(x = x, y = y, fill = class)) +
    geom_tile() +
    scale_fill_manual(
      values = colors,
      labels = labels,
      name = "Class"
    ) +
    theme_minimal() +
    theme(legend.position = "right")
  
  legend_grob_classified <- cowplot::get_legend(legend_classified)
  
  # Combine all years
  combined_classified <- wrap_plots(plots_classified, ncol = 8, nrow = 5) +
    plot_annotation(
      title = paste0(climate_var, " - Threshold Classification Results (1982-2020)"),
      theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))
    )
  
  # Add legend
  final_plot <- plot_grid(
    combined_classified,
    legend_grob_classified,
    ncol = 2,
    rel_widths = c(1, 0.1)
  )
  
  # Save
  ggsave(
    filename = file.path(figure_date_dir, paste0(climate_var, "_Classified.png")),
    plot = final_plot,
    width = 26,
    height = 15,
    dpi = 300,
    bg = "white"
  )
  
  cat("  - Done!\n")
}

cat("\n\n=== All processing completed! ===\n")
cat("Processed", length(all_results), "climate factors\n")
cat("Generated 2 charts for each factor: original values and classified results\n")
cat("Total generated", length(all_results) * 2, "charts\n")
cat("\nCharts saved at:", figure_date_dir, "\n")
cat("Rasters saved at:", Datedir, "\n")

# ==============================================================================
# Step 5: Generate Statistical Summary
# ==============================================================================

cat("\n\nStep 5: Generating statistical summary...\n")

# Count pixels for each category for each climate factor
summary_stats <- list()

for(climate_var in names(all_results)) {
  
  result <- all_results[[climate_var]]
  years <- result$years
  
  yearly_stats <- lapply(years, function(year) {
    r <- result$classified[[as.character(year)]]
    freq_table <- freq(r, useNA = "no")
    
    if(!is.null(freq_table) && nrow(freq_table) > 0) {
      data.frame(
        climate_var = climate_var,
        year = year,
        class = freq_table[, 1],
        count = freq_table[, 2],
        stringsAsFactors = FALSE
      )
    } else {
      NULL
    }
  }) %>% bind_rows()
  
  summary_stats[[climate_var]] <- yearly_stats
}

# Merge all stats
all_stats <- bind_rows(summary_stats)

# Add class label
all_stats <- all_stats %>%
  mutate(
    class_label = case_when(
      class == 1 ~ "VGC",
      class == 2 ~ "VSO",
      class == 3 ~ "Stable",
      TRUE ~ NA_character_
    )
  )

# Save statistical results
write.csv(all_stats, 
          file.path(Datedir, "Classification_Statistics.csv"),
          row.names = FALSE)

cat("Statistical summary saved!\n")

# Print summary
cat("\nClassification statistics summary:\n")
print(all_stats %>%
        group_by(climate_var, class_label) %>%
        summarise(
          mean_count = mean(count, na.rm = TRUE),
          total_count = sum(count, na.rm = TRUE),
          .groups = "drop"
        ))

cat("\n\n=== Script execution completed! ===\n")