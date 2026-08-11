

setwd("G:/03-Program/GlobalVeg")

RUN_name <- "Fig4-new"

Datedir <- file.path('./Out-Table-and-Figure-VPD et al', RUN_name)
dir.create(Datedir, showWarnings = FALSE, recursive = TRUE)
figure_date_dir <- Datedir

# --- KEY INPUT PATHS (edit to match your machine) ---
TRANSITIONS_CSV <- "./Nas/NCC-ResultData/Result_4/Combined_All_Groups_response_type_change_GPP_5yr.csv"
TRIPPING_CSV    <- "./Nas/NCC-ResultData/Result_3/SHAP_ClimateThreshold_TrippingPoint.csv"
SD_DIR          <- "./Nas/NCC-ResultData/Result_2/Threshold_Pipeline_Server/GPP_breakpoint"

LULC_GRID_TIF   <- "./Nas/resample_example_0.5degree.tif"
LULC_CLASS_TIF  <- "./Nas/GLC_FCS30D_stable_0p5deg_1985_2020_8class.tif"
LULC_KEEP_CODES <- 1:8

TARGET_PERIODS <- c("1991-1995", "1996-2000", "2001-2005", "2006-2010")

cat("=== 05 Map + ring-inset + lines redesign ===\n")
cat("Output folder:", figure_date_dir, "\n\n")


# ============================================================================ #
# 1. PACKAGES
# ============================================================================ #
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(cowplot)
  library(scales)
  library(data.table)
  library(raster)
  library(sf)
  library(rnaturalearth)
})

world_shp <- ne_coastline()


# ============================================================================ #
# 2. ISOLATED RASTER HELPERS  -- UNCHANGED
# ============================================================================ #
saveToRaster <- function(labs, values) {
  r_test <- raster(LULC_GRID_TIF)
  values(r_test) <- NA
  values(r_test)[labs] <- values
  r_test
}

cutByLULC <- function(r) {
  LULC <- raster(LULC_CLASS_TIF)
  keep_lab <- which(values(LULC) %in% LULC_KEEP_CODES)
  saveToRaster(labs = keep_lab, values = values(r)[keep_lab])
}

rasterToMapDF <- function(lab_code_df, code_col = "code") {
  r <- saveToRaster(lab_code_df$lab, lab_code_df[[code_col]])
  r <- cutByLULC(r)
  dfm <- as.data.frame(r, xy = TRUE)
  names(dfm) <- c("x", "y", code_col)
  dfm[!is.na(dfm[[code_col]]), ]
}

filter_valid_data <- function(df) {
  df %>% filter(note != "no_matched_windows" | is.na(note))
}


# ============================================================================ #
# 3. READ DATA  -- UNCHANGED
# ============================================================================ #
cat("Step 0: reading data...\n")

df_transitions <- read.csv(TRANSITIONS_CSV)
cat("Transitions:", nrow(df_transitions), "rows x", ncol(df_transitions), "cols\n")

df_trippingtype <- read.csv(TRIPPING_CSV)

SD_paths <- list.files(SD_DIR, full.names = TRUE, pattern = "_SDdata\\.RData$")

energy_vars <- c("CO2", "TMP", "SR")
water_vars  <- c("VPD", "PRE", "SMroot", "SMsurf")
climate_vars <- c(energy_vars, water_vars)
single_threshold_vars <- c("CO2", "VPD", "PRE", "SMroot", "SMsurf")
double_threshold_vars <- c("TMP", "SR")

cat("Valid transitions:", nrow(filter_valid_data(df_transitions)), "\n\n")


# ============================================================================ #
# 4. STEP 1 -- THRESHOLD-CROSSING EVENTS  -- UNCHANGED
# ============================================================================ #
calculate_single_threshold_cross <- function(sd_list, var_name, tripping_data, energy_vars) {
  is_energy <- var_name %in% energy_vars
  ti <- tripping_data %>% filter(Climate_Variable == var_name)
  if (nrow(ti) == 0) return(NULL)
  shap_threshold <- ti$X_Position[1]
  left_type  <- ti$Left_contribution[1]
  right_type <- ti$Right_contribution[1]
  
  pixel_names <- names(sd_list); if (is.null(pixel_names)) pixel_names <- as.character(seq_along(sd_list))
  out <- list()
  for (i in seq_along(sd_list)) {
    pd <- sd_list[[i]]; lb <- pixel_names[i]
    if (is.null(pd) || nrow(pd) < 2 || !("sd_cli" %in% names(pd))) next
    cli <- pd$sd_cli; yrs <- pd$year_new
    position <- ifelse(cli < shap_threshold, "left", "right")
    cross_dir <- rep(0, length(yrs)); event_type <- rep("none", length(yrs))
    for (j in 2:length(yrs)) {
      dch <- cli[j] - cli[j - 1]
      is_cross <- position[j - 1] != position[j]
      if (is_cross) {
        if (position[j - 1] == "left" && position[j] == "right") {
          cross_dir[j] <- ifelse(is_energy, 1, 3); event_type[j] <- "cross_L2R"
        } else {
          cross_dir[j] <- ifelse(is_energy, 2, 3); event_type[j] <- "cross_R2L"
        }
      } else if (abs(dch) > 0) {
        if (position[j] == "left") {
          if (dch > 0) { cross_dir[j] <- ifelse(is_energy, 2, 3); event_type[j] <- "left_increasing" }
          else         { cross_dir[j] <- ifelse(is_energy, 1, 4); event_type[j] <- "left_decreasing" }
        } else {
          if (dch > 0) { cross_dir[j] <- ifelse(is_energy, 1, 4); event_type[j] <- "right_increasing" }
          else         { cross_dir[j] <- ifelse(is_energy, 2, 3); event_type[j] <- "right_decreasing" }
        }
      }
    }
    out[[i]] <- data.frame(lab = lb, year = yrs, sd_cli = cli, threshold = shap_threshold,
                           position = position, cross_dir = cross_dir, event_type = event_type,
                           left_type = left_type, right_type = right_type, var = var_name,
                           var_group = ifelse(is_energy, "Energy", "Water"),
                           stringsAsFactors = FALSE)
  }
  dplyr::bind_rows(out)
}

calculate_double_threshold_cross <- function(sd_list, var_name, tripping_data, energy_vars) {
  ti <- tripping_data %>% filter(Climate_Variable == var_name) %>% arrange(X_Position)
  if (nrow(ti) != 2) return(NULL)
  t1 <- ti$X_Position[1]; t2 <- ti$X_Position[2]
  pixel_names <- names(sd_list); if (is.null(pixel_names)) pixel_names <- as.character(seq_along(sd_list))
  out <- list()
  for (i in seq_along(sd_list)) {
    pd <- sd_list[[i]]; lb <- pixel_names[i]
    if (is.null(pd) || nrow(pd) < 2 || !("sd_cli" %in% names(pd))) next
    cli <- pd$sd_cli; yrs <- pd$year_new
    zone <- rep(NA, length(cli))
    zone[cli < t1] <- "Zone1"; zone[cli >= t1 & cli < t2] <- "Zone2"; zone[cli >= t2] <- "Zone3"
    cross_dir <- rep(0, length(yrs)); event_type <- rep("none", length(yrs))
    for (j in 2:length(yrs)) {
      if (is.na(zone[j - 1]) || is.na(zone[j])) next
      dch <- cli[j] - cli[j - 1]; is_cross <- zone[j - 1] != zone[j]
      if (is_cross) {
        if (zone[j-1]=="Zone1" && zone[j]=="Zone2") { cross_dir[j] <- 1; event_type[j] <- "cross_T1_L2R" }
        else if (zone[j-1]=="Zone2" && zone[j]=="Zone1") { cross_dir[j] <- 2; event_type[j] <- "cross_T1_R2L" }
        else if (zone[j-1]=="Zone2" && zone[j]=="Zone3") { cross_dir[j] <- 2; event_type[j] <- "cross_T2_L2R" }
        else if (zone[j-1]=="Zone3" && zone[j]=="Zone2") { cross_dir[j] <- 1; event_type[j] <- "cross_T2_R2L" }
      } else if (abs(dch) > 0) {
        if (zone[j]=="Zone1") { if (dch>0){cross_dir[j]<-2;event_type[j]<-"Zone1_increasing"} else {cross_dir[j]<-1;event_type[j]<-"Zone1_decreasing"} }
        else if (zone[j]=="Zone2") { if (dch>0){cross_dir[j]<-1;event_type[j]<-"Zone2_increasing"} else {cross_dir[j]<-2;event_type[j]<-"Zone2_decreasing"} }
        else if (zone[j]=="Zone3") { if (dch>0){cross_dir[j]<-2;event_type[j]<-"Zone3_increasing"} else {cross_dir[j]<-1;event_type[j]<-"Zone3_decreasing"} }
      }
    }
    out[[i]] <- data.frame(lab = lb, year = yrs, sd_cli = cli, threshold1 = t1, threshold2 = t2,
                           zone = zone, cross_dir = cross_dir, event_type = event_type,
                           var = var_name, var_group = "Energy", stringsAsFactors = FALSE)
  }
  dplyr::bind_rows(out)
}

cat("Step 1: identifying threshold-crossing events...\n")
all_cross_events <- list()
for (v in climate_vars) {
  sd_path <- SD_paths[grepl(paste0("GPP_", v, "_"), SD_paths)]
  if (length(sd_path) == 0) { cat("  missing SD for", v, "\n"); next }
  load(sd_path[1])  # loads sd_all
  ev <- if (v %in% single_threshold_vars) calculate_single_threshold_cross(sd_all, v, df_trippingtype, energy_vars)
  else calculate_double_threshold_cross(sd_all, v, df_trippingtype, energy_vars)
  if (!is.null(ev)) all_cross_events[[v]] <- ev
}
df_cross_all <- dplyr::bind_rows(all_cross_events)
if (nrow(df_cross_all) > 0)
  write.csv(df_cross_all, file.path(Datedir, "Step1_CrossEvent_All.csv"), row.names = FALSE)
cat("  events:", nrow(df_cross_all), "| pixels:", length(unique(df_cross_all$lab)), "\n\n")


# ============================================================================ #
# 5. GROUP TRANSITIONS BY PERIOD  -- UNCHANGED
# ============================================================================ #
unique_periods <- sort(unique(df_transitions$period_label[!is.na(df_transitions$period_label)]))
df_by_period <- list()
for (p in unique_periods) df_by_period[[p]] <- df_transitions %>% filter(period_label == p)
period_mapping <- data.frame(period_label = unique_periods,
                             panel_letter = letters[seq_along(unique_periods)])
cat("Periods found:", paste(unique_periods, collapse = ", "), "\n\n")


# ============================================================================ #
# 6. MATCH RATE  -- UNCHANGED
# ============================================================================ #
cat("=== Match-rate computation ===\n")
vars <- c("CO2", "TMP", "SR", "VPD", "PRE", "SMroot", "SMsurf")

df_cross_wide <- df_cross_all %>%
  dplyr::select(lab, year, var, cross_dir) %>%
  tidyr::pivot_wider(names_from = var, values_from = cross_dir, names_prefix = "Cross_")
df_cross_wide$lab <- as.character(df_cross_wide$lab)

test_results_by_period <- list()
for (i in 1:4) {
  period_label <- TARGET_PERIODS[i]
  period_data <- df_by_period[[period_label]]
  if (is.null(period_data)) { cat("  (no data for", period_label, ")\n"); next }
  
  df_combined_events <- period_data %>%
    mutate(
      change_dir = case_when(
        !is.na(change_year) & from == "Stable" & to == "VGC" ~ 1,
        !is.na(change_year) & from == "VGC" & to == "Stable" ~ 2,
        !is.na(change_year) & from == "VSO" & to == "VGC" ~ 3,
        !is.na(change_year) & from == "VGC" & to == "VSO" ~ 4,
        !is.na(persistent_type) & persistent_type == "VGC" ~ 13,
        !is.na(persistent_type) & persistent_type == "VSO" ~ 4,
        !is.na(persistent_type) & persistent_type == "Stable" ~ 2,
        TRUE ~ 0),
      change_group = case_when(
        change_dir %in% c(1, 2) ~ "Energy",
        change_dir %in% c(3, 4) ~ "Water",
        change_dir == 13 ~ "Both",
        TRUE ~ NA_character_),
      event_type = case_when(
        !is.na(change_year) ~ paste0(from, "->", to, " (Change)"),
        !is.na(persistent_type) ~ paste0(persistent_type, " (Persistent)"),
        TRUE ~ NA_character_)
    ) %>%
    filter(change_dir != 0)
  df_combined_events$lab <- as.character(df_combined_events$lab)
  
  common_labs <- intersect(unique(df_cross_wide$lab), unique(df_combined_events$lab))
  dcw <- df_cross_wide %>% filter(lab %in% common_labs)
  dce <- df_combined_events %>% filter(lab %in% common_labs)
  
  event_matrix <- dce %>%
    left_join(dcw %>% mutate(match_year = year + 0),
              by = c("lab" = "lab", "change_year" = "match_year")) %>%
    mutate(lag = 0)
  
  event_matrix <- event_matrix %>%
    mutate(
      match_CO2    = Cross_CO2    == change_dir,
      match_TMP    = Cross_TMP    == change_dir,
      match_SR     = Cross_SR     == change_dir,
      match_VPD    = Cross_VPD    == change_dir,
      match_PRE    = Cross_PRE    == change_dir,
      match_SMroot = Cross_SMroot == change_dir,
      match_SMsurf = Cross_SMsurf == change_dir
    ) %>%
    mutate(
      match_energy = if_any(c(match_CO2, match_TMP, match_SR), ~ .x == TRUE),
      match_water  = if_any(c(match_VPD, match_PRE, match_SMroot, match_SMsurf), ~ .x == TRUE)
    )
  
  test_results <- data.frame()
  for (v in vars) {
    col_cross <- paste0("Cross_", v); col_match <- paste0("match_", v)
    if (v %in% c("CO2", "TMP", "SR")) {
      df <- event_matrix %>% subset(event_type %in% c("VGC->Stable (Change)", "Stable->VGC (Change)"))
    } else {
      df <- event_matrix %>% subset(event_type %in% c("VGC->VSO (Change)", "VSO->VGC (Change)"))
    }
    df <- df %>% filter(!!sym(col_cross) != 0)
    n <- nrow(df); k <- sum(df[[col_match]], na.rm = TRUE)
    test <- prop.test(k, n, p = 0.25)
    test_results <- rbind(test_results, data.frame(
      factor = v, period = period_label, n = n,
      match_rate = test$estimate, CI_low = test$conf.int[1],
      CI_high = test$conf.int[2], p_value = test$p.value))
  }
  
  df_energy <- event_matrix %>% subset(event_type %in% c("Stable (Persistent)", "VGC->Stable (Change)", "Stable->VGC (Change)"))
  te <- prop.test(sum(df_energy$match_energy, na.rm = TRUE), nrow(df_energy), p = 0.25)
  test_results <- rbind(test_results, data.frame(
    factor = "Energy_factors", period = period_label, n = nrow(df_energy),
    match_rate = te$estimate, CI_low = te$conf.int[1], CI_high = te$conf.int[2], p_value = te$p.value))
  
  df_water <- event_matrix %>% subset(event_type %in% c("VGC->VSO (Change)", "VSO->VGC (Change)"))
  tw <- prop.test(sum(df_water$match_water, na.rm = TRUE), nrow(df_water), p = 0.25)
  test_results <- rbind(test_results, data.frame(
    factor = "Water_factors", period = period_label, n = nrow(df_water),
    match_rate = tw$estimate, CI_low = tw$conf.int[1], CI_high = tw$conf.int[2], p_value = tw$p.value))
  
  test_results_by_period[[i]] <- test_results
  cat("  period", period_label, "done\n")
}

df_test_by_period <- bind_rows(test_results_by_period)
write.csv(df_test_by_period, file.path(Datedir, "Significance_Test_by_Period.csv"), row.names = FALSE)
cat("Saved match-rate table.\n\n")


# ============================================================================ #
# 7. SHARED STYLE CONSTANTS
# ============================================================================ #
persistent_colors <- c("stable" = "#56a0d3", "δVGC" = "#DA5CA3", "δVSO" = "#4DAF4A")
change_colors <- c(
  "stable->δVGC" = "#2E8B57", "stable->δVSO" = "#48D1CC",
  "δVGC->stable" = "#DAA520", "δVGC->δVSO" = "#CD853F",
  "δVSO->stable" = "#4682B4", "δVSO->δVGC" = "#9370DB")
ch_levels <- c("VGC->VSO", "VGC->Stable", "VSO->Stable", "VSO->VGC", "Stable->VGC", "Stable->VSO")
ch_labels <- c("δVGC->δVSO", "δVGC->stable", "δVSO->stable", "δVSO->δVGC", "stable->δVGC", "stable->δVSO")

combined_colors <- c(persistent_colors, change_colors)
combined_levels <- names(combined_colors)   # 9 categories

# --- Sunburst category order (also decides what starts at 12 o'clock) -------
SUNBURST_PERS_ORDER <- c("δVSO", "stable", "δVGC")
SUNBURST_CHG_ORDER  <- c("δVGC->δVSO", "δVSO->δVGC", "δVGC->stable",
                         "stable->δVGC", "δVSO->stable", "stable->δVSO")
sunburst_levels <- c(SUNBURST_PERS_ORDER, SUNBURST_CHG_ORDER)

# Inner-ring (aggregate) colours: neutral, so they read as a summary layer.
INNER_RING_COLORS <- c("Persistent" = "#4A90E2", "Transition"  = "#D9534F")
INNER_RING_COLORS <- c("Persistent" = "#3B7EA1", "Transition"  = "#E05A47")
INNER_RING_COLORS <- c("Persistent" = "#248c2a", "Transition" = "#7555b8")
sunburst_fill_values <- c(combined_colors, INNER_RING_COLORS)

# geom_text() `size` is in mm; theme() text sizes are in pt. This converts.
PT_TO_MM <- 25.4 / 72.27

# --- Line-plot colors / shapes ---
line_color_map <- c(
  "Energy-related factors" = "#FAD689", "Moisture-related factors" = "#F07881",
  "TMP" = "#E41A1C", "SR" = "#377EB8", "CO2" = "#4DAF4A",
  "VPD" = "#984EA3", "PRE" = "#FF7F00", "SMroot" = "#A65628", "SMsurf" = "#F781BF")
line_shape_map <- c(
  "Energy-related factors" = 16, "Moisture-related factors" = 17,
  "TMP" = 16, "SR" = 17, "CO2" = 15, "VPD" = 16, "PRE" = 17, "SMroot" = 15, "SMsurf" = 18)

df_plot_line <- df_test_by_period %>%
  mutate(
    factor_group = case_when(
      factor %in% c("Energy_factors", "Water_factors") ~ "ALL",
      factor %in% c("CO2", "TMP", "SR") ~ "Energy",
      factor %in% c("VPD", "PRE", "SMroot", "SMsurf") ~ "Water",
      TRUE ~ "Other"),
    factor = factor(factor,
                    levels = c("Energy_factors", "Water_factors", "TMP", "SR", "CO2",
                               "VPD", "PRE", "SMroot", "SMsurf"),
                    labels = c("Energy-related factors", "Moisture-related factors",
                               "TMP", "SR", "CO2", "VPD", "PRE", "SMroot", "SMsurf")),
    period = factor(period, levels = TARGET_PERIODS),
    lab_title = factor(factor_group,
                       levels = c("ALL", "Energy", "Water"),
                       labels = c("Combination of energy and\nwater limiting factors",
                                  "Energy-related factors",
                                  "Moisture-related factors"))
  )
write.csv(df_plot_line, file.path(Datedir, "MatchPercentage.csv"), row.names = FALSE)


# ============================================================================ #
# 8. COMBINED RASTER PREP: persistent + transition per period
# ----------------------------------------------------------------------------
# [BUGFIX] recode persistent "Stable" -> "stable" (lower-case) so it matches
# the colour/level vector; otherwise those pixels became NA and vanished.
# ============================================================================ #
prep_combined <- function(period_label) {
  period_data <- df_by_period[[period_label]]
  if (is.null(period_data)) return(NULL)
  
  pers <- period_data %>%
    filter_valid_data() %>%
    filter(change_id == 0 &
             (note == "no_changed_point" | note == "single_data_point") &
             !is.na(persistent_type)) %>%
    mutate(rtype = dplyr::recode(persistent_type,
                                 "VGC" = "δVGC", "VSO" = "δVSO", "Stable" = "stable"),  # [BUGFIX]
           .src = 0L) %>%
    dplyr::select(lab, rtype, .src)
  
  chg <- filter_valid_data(period_data) %>%
    filter(!is.na(from) & !is.na(to) & from != to) %>%
    mutate(change_type = paste0(from, "->", to)) %>%
    group_by(lab) %>% slice(1) %>% ungroup() %>%
    mutate(rtype = as.character(factor(change_type, levels = ch_levels, labels = ch_labels)),
           .src = 1L) %>%
    dplyr::select(lab, rtype, .src)
  
  combined <- bind_rows(pers, chg)
  if (nrow(combined) == 0) return(NULL)
  combined <- combined %>% arrange(lab, desc(.src)) %>%
    distinct(lab, .keep_all = TRUE) %>% dplyr::select(lab, rtype)
  
  combined$rtype <- factor(combined$rtype, levels = combined_levels)
  if (any(is.na(combined$rtype)))
    warning("prep_combined(", period_label, "): ",
            sum(is.na(combined$rtype)), " pixels fell outside the 9 categories.")
  combined$type_code <- as.numeric(combined$rtype)
  
  df_map <- rasterToMapDF(data.frame(lab = combined$lab, code = combined$type_code), "code")
  lab_map <- combined %>% dplyr::select(type_code, rtype) %>% distinct()
  df_map <- left_join(df_map, lab_map, by = c("code" = "type_code"))
  
  bar_stats <- combined %>% count(rtype, .drop = FALSE) %>%
    mutate(percentage = n / sum(n) * 100, period = period_label)
  
  pers_stats <- pers %>% count(rtype) %>% mutate(percentage = n / sum(n) * 100,
                                                 period = period_label) %>%
    rename(persistent_type = rtype)
  chg_stats  <- chg  %>% count(rtype) %>% mutate(percentage = n / sum(n) * 100,
                                                 period = period_label) %>%
    rename(change_type = rtype)
  
  list(map_df = df_map, bar_stats = bar_stats,
       pers_stats = pers_stats, chg_stats = chg_stats)
}

cat("Building combined (persistent + transition) data per period...\n")
combined_prep <- list()
df_pers_stats_all <- data.frame(); df_change_stats_all <- data.frame()
for (p in TARGET_PERIODS) {
  combined_prep[[p]] <- prep_combined(p)
  if (!is.null(combined_prep[[p]])) {
    df_pers_stats_all   <- rbind(df_pers_stats_all,   combined_prep[[p]]$pers_stats)
    df_change_stats_all <- rbind(df_change_stats_all, combined_prep[[p]]$chg_stats)
  }
}
write.csv(df_pers_stats_all,   file.path(Datedir, "Fig4_persistent_stat.csv"),  row.names = FALSE)
write.csv(df_change_stats_all, file.path(Datedir, "Fig4_change_type_stat.csv"), row.names = FALSE)

bar_stats_all <- bind_rows(lapply(TARGET_PERIODS, function(p) combined_prep[[p]]$bar_stats))
write.csv(bar_stats_all, file.path(Datedir, "Fig4_area_share_all_periods.csv"), row.names = FALSE)
cat("  done.\n\n")


# ============================================================================ #
# 9. MAP PANEL RENDERER + TYPE-LEGEND BUILDER  (functions only)
# ============================================================================ #
render_map_panel <- function(prep, period_label, panel_letter = NULL,
                             title_size = 26, axis_text_size = 24,
                             show_title_period = TRUE) {
  if (is.null(prep) || nrow(prep$map_df) == 0) return(ggplot() + theme_void())
  ttl <- if (show_title_period) period_label else ""
  if (!is.null(panel_letter)) ttl <- paste0("(", panel_letter, ") ", ttl)
  ggplot(prep$map_df) +
    geom_hline(yintercept = 0, color = "gray50", linetype = "dashed", linewidth = 0.3) +
    geom_tile(aes(x = x, y = y, fill = rtype)) +
    geom_sf(data = world_shp, fill = NA, color = "gray30", linewidth = 0.2) +
    coord_sf(ylim = c(-60, 90), expand = FALSE) +
    scale_fill_manual(values = combined_colors, na.value = "transparent", drop = FALSE) +
    labs(title = ttl) +
    theme_bw() +
    theme(plot.title = element_text(face = "bold", size = title_size, hjust = 0),
          legend.position = "none", panel.background = element_blank(),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.title = element_blank(), axis.text = element_text(size = axis_text_size),
          plot.margin = margin(4, 6, 4, 6))
}

build_type_legend <- function(color_vec, legend_name, nrow = 1,
                              title_size = 26, text_size = 24) {
  cowplot::get_legend(
    ggplot(data.frame(x = 1, y = factor(names(color_vec), levels = names(color_vec)))) +
      geom_tile(aes(x = x, y = y, fill = y)) +
      scale_fill_manual(name = legend_name, values = color_vec) +
      guides(fill = guide_legend(nrow = nrow)) + theme_void() +
      theme(legend.position = "bottom",
            legend.title = element_text(size = title_size, face = "bold"),
            legend.text = element_text(size = text_size)))
}


# ============================================================================ #
# 10. TWO-RING NESTED DONUT (RING/SUNBURST)  + inset wrapper
# ----------------------------------------------------------------------------
#  inner ring : Persistent vs Transition total    (colour only, via legend)
#  outer ring : the 9 sub-categories               (colour only, via legends)
#  show_pct = FALSE  ->  NO numbers printed anywhere (used for the map inset).
# ============================================================================ #
R_HOLE_OUT <- 1.00
R_IN_OUT   <- 2.40
R_OUT_IN   <- 2.40
R_OUT_OUT  <- 3.70
R_ELBOW    <- 3.88
R_LABEL    <- 4.06
R_XMAX     <- 4.30
RING_BORDER_COLOR <- NA
RING_BORDER_WIDTH <- 0.15

LABEL_MIN_INSIDE <- 4.0
LABEL_GAP        <- 4.2
INNER_LABEL_MODE <- "pct"

fmt_pct <- function(p) {
  ifelse(p <= 0,     "",
         ifelse(p >= 1,     sprintf("%.1f%%", p),
                ifelse(p >= 0.1,   sprintf("%.2f%%", p),
                       sprintf("%.3f%%", p))))
}

spread_y <- function(y, gap, lo, hi) {
  n <- length(y)
  if (n <= 1) return(y)
  o  <- order(y)
  ys <- y[o]
  for (i in 2:n) ys[i] <- max(ys[i], ys[i - 1] + gap)
  if (ys[n] > hi) {
    ys[n] <- hi
    for (i in (n - 1):1) ys[i] <- min(ys[i], ys[i + 1] - gap)
  }
  ys[1] <- max(ys[1], lo)
  out <- numeric(n); out[o] <- ys
  out
}

build_sunburst_data <- function(bar_stats) {
  d <- bar_stats %>%
    mutate(rtype = as.character(rtype)) %>%
    filter(!is.na(rtype))
  
  miss <- setdiff(sunburst_levels, d$rtype)
  if (length(miss) > 0)
    d <- bind_rows(d, data.frame(rtype = miss, n = 0L, percentage = 0,
                                 period = d$period[1], stringsAsFactors = FALSE))
  
  d <- d %>%
    mutate(group = ifelse(rtype %in% SUNBURST_PERS_ORDER, "Persistent", "Transition"),
           rtype = factor(rtype, levels = sunburst_levels)) %>%
    arrange(rtype)
  
  d$ymax <- cumsum(d$percentage)
  d$ymin <- d$ymax - d$percentage
  
  outer <- d %>%
    mutate(ring = "outer", xmin = R_OUT_IN, xmax = R_OUT_OUT,
           fill_key = as.character(rtype), lab_name = as.character(rtype))
  
  inner <- d %>%
    group_by(group) %>%
    summarise(percentage = sum(percentage), n = sum(n),
              ymin = min(ymin), ymax = max(ymax), .groups = "drop") %>%
    mutate(ring = "inner", xmin = R_HOLE_OUT, xmax = R_IN_OUT,
           fill_key = group, lab_name = group) %>%
    arrange(ymin)
  
  list(outer = outer, inner = inner, total_n = sum(d$n),
       period = as.character(d$period[1]))
}

build_sunburst_panel <- function(bar_stats, letter = NULL,
                                 title_text     = NULL,
                                 title_size     = 26,
                                 inner_lab_size = 24,
                                 outer_lab_size = 24,
                                 center_size    = 24,
                                 show_center    = TRUE,
                                 show_pct       = TRUE,          # <- NEW
                                 inner_label_mode = INNER_LABEL_MODE) {
  
  sb <- build_sunburst_data(bar_stats)
  o  <- sb$outer
  i  <- sb$inner
  
  o$ymid <- (o$ymin + o$ymax) / 2
  i$ymid <- (i$ymin + i$ymax) / 2
  
  i$lab_txt <- switch(inner_label_mode,
                      "pct"  = fmt_pct(i$percentage),
                      "name" = i$lab_name,
                      "both" = paste0(i$lab_name, "\n", fmt_pct(i$percentage)),
                      fmt_pct(i$percentage))
  
  inside  <- o %>% filter(percentage >= LABEL_MIN_INSIDE)
  outside <- o %>% filter(percentage > 0 & percentage < LABEL_MIN_INSIDE)
  
  if (show_pct && nrow(outside) > 0) {
    outside$side  <- ifelse(outside$ymid < 50, "R", "L")
    outside$y_adj <- outside$ymid
    for (s in unique(outside$side)) {
      k  <- outside$side == s
      lo <- if (s == "R") 1.5  else 51.5
      hi <- if (s == "R") 48.5 else 98.5
      outside$y_adj[k] <- spread_y(outside$ymid[k], LABEL_GAP, lo, hi)
    }
    outside$hj <- ifelse(outside$side == "R", 0, 1)
  }
  
  ttl <- ""
  if (!is.null(title_text)) ttl <- title_text
  if (!is.null(letter))     ttl <- paste0("(", letter, ") ", ttl)
  
  ring_layer <- function(dat) {
    if (is.na(RING_BORDER_COLOR)) {
      geom_rect(data = dat,
                aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
                    fill = fill_key, colour = fill_key),
                linewidth = RING_BORDER_WIDTH)
    } else {
      geom_rect(data = dat,
                aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
                    fill = fill_key),
                colour = RING_BORDER_COLOR, linewidth = RING_BORDER_WIDTH)
    }
  }
  
  p <- ggplot() + ring_layer(i) + ring_layer(o)
  
  # ---- all numeric labels are gated behind show_pct ----
  if (show_pct && nrow(outside) > 0) {
    p <- p +
      geom_segment(data = outside,
                   aes(x = R_OUT_OUT, xend = R_ELBOW, y = ymid, yend = ymid),
                   color = "gray25", linewidth = 0.4) +
      geom_segment(data = outside,
                   aes(x = R_ELBOW, xend = R_LABEL - 0.05, y = ymid, yend = y_adj),
                   color = "gray25", linewidth = 0.4) +
      geom_text(data = outside,
                aes(x = R_LABEL, y = y_adj, label = fmt_pct(percentage), hjust = hj),
                size = outer_lab_size * PT_TO_MM, color = "gray10")
  }
  
  if (show_pct && nrow(inside) > 0)
    p <- p + geom_text(data = inside,
                       aes(x = (R_OUT_IN + R_OUT_OUT) / 2, y = ymid,
                           label = fmt_pct(percentage)),
                       size = outer_lab_size * PT_TO_MM,
                       fontface = "bold", color = "white")
  
  if (show_pct)
    p <- p + geom_text(data = i,
                       aes(x = (R_HOLE_OUT + R_IN_OUT) / 2, y = ymid, label = lab_txt),
                       size = inner_lab_size * PT_TO_MM,
                       fontface = "bold", color = "white", lineheight = 0.95)
  
  if (show_center)
    p <- p + annotate("text", x = 0, y = 0,
                      label = paste0("n = ", format(sb$total_n, big.mark = ",")),
                      size = center_size * PT_TO_MM, fontface = "bold",
                      color = "gray15")
  
  p +
    coord_polar(theta = "y", start = 0, direction = 1, clip = "off") +
    scale_fill_manual(values = sunburst_fill_values, guide = "none") +
    scale_colour_manual(values = sunburst_fill_values, guide = "none") +
    scale_x_continuous(limits = c(0, R_XMAX), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) +
    labs(title = ttl, x = NULL, y = NULL) +
    theme_void() +
    theme(plot.title = element_text(face = "bold", size = title_size, hjust = 0),
          plot.margin = margin(6, 10, 6, 10))
}

# --- the ring chart as it appears INSIDE a map corner: no numbers, no title,
#     translucent white backing so it reads over the ocean/coastlines. --------
sunburst_inset <- function(bar_stats) {
  build_sunburst_panel(bar_stats, letter = NULL, title_text = NULL,
                       show_pct = FALSE, show_center = FALSE) +
    theme(plot.background = element_rect(fill = NA, colour = NA),
          plot.margin = margin(0, 0, 0, 0))
}


# ============================================================================ #
# 11. MATCH-RATE LINE-PLOT BUILDER + LINE LEGEND BUILDER  (functions only)
# ============================================================================ #
LINE_TITLE_SIZE <- 26
LINE_AXIS_SIZE  <- 24

build_lineplot_panel <- function(df_line, group_name, letter, y_line, title_text,
                                 title_size = LINE_TITLE_SIZE,
                                 axis_text_size = LINE_AXIS_SIZE,
                                 axis_title_size = LINE_AXIS_SIZE) {
  d <- df_line %>% filter(factor_group == group_name)
  full_title <- paste0("(", letter, ") ", title_text)
  
  ggplot(d, aes(x = period, y = match_rate * 100,
                group = factor, color = factor, shape = factor)) +
    geom_hline(yintercept = y_line, linetype = "dashed", color = "gray50", linewidth = 0.8) +
    geom_line(linewidth = 1.2) +
    geom_point(size = 4, fill = "white", stroke = 1.5) +
    scale_color_manual(values = line_color_map, name = "Climate factor", drop = FALSE) +
    scale_shape_manual(values = line_shape_map, name = "Climate factor", drop = FALSE) +
    scale_y_continuous(name = "Match rate (%)") +
    labs(title = full_title, x = NULL) +
    theme_bw() +
    theme(
      plot.title = element_text(size = title_size, face = "bold", hjust = 0),
      axis.text.x = element_text(size = axis_text_size, color = "black", angle = 45, hjust = 1),
      axis.text.y = element_text(size = axis_text_size, color = "black"),
      axis.title.y = element_text(size = axis_title_size, face = "bold"),
      legend.position = "none",
      panel.background = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.margin = margin(6, 10, 6, 6)
    )
}

build_line_legend <- function(title_size = 26, text_size = 24, nrow = 3) {
  # Build from the FULL data so every factor level has real line+point rows and
  # its legend key renders a glyph. (Filtering to one group left the individual
  # factors with no data, so their keys drew as bare text with no line.)
  p_leg <- ggplot(df_plot_line,
                  aes(x = period, y = match_rate * 100,
                      group = factor, color = factor, shape = factor)) +
    geom_line(linewidth = 1.2) +
    geom_point(size = 4, fill = "white", stroke = 1.5) +
    scale_color_manual(values = line_color_map, name = "Climate factor", drop = FALSE) +
    scale_shape_manual(values = line_shape_map, name = "Climate factor", drop = FALSE) +
    guides(color = guide_legend(nrow = nrow), shape = guide_legend(nrow = nrow)) +
    theme(legend.position = "bottom",
          legend.title = element_text(size = title_size, face = "bold"),
          legend.text = element_text(size = text_size),
          legend.key.width = unit(1.4, "cm"),
          legend.key = element_rect(fill = "transparent", color = NA),
          legend.background = element_rect(fill = "transparent", color = NA))
  cowplot::get_legend(p_leg)
}


# ============================================================================ #
# 12. THE FIGURE:  LEFT = 4 maps (ring inset in corner)  |  RIGHT = 3 lines
# ============================================================================ #
cat("Assembling the combined maps + lines figure...\n")

MAP_TITLE_SIZE <- 26
MAP_AXIS_SIZE  <- 24

INSET_L <- -0.545      # left edge  -> hug the SW corner (~180 W)
INSET_B <- -0.12   # bottom edge -> map floor (~60 S)
INSET_R <- 0.8   # right edge  -> out toward ~90 W (left quarter)
INSET_T <- 0.85  # top edge    -> up to ~40 N

COL_W_RATIO <- c(1.5, 0.5)     # left(maps) : right(lines) width
FIG_W  <- 14               # inches
BODY_H <- 18               # body (panels) relative height
LEG_H  <- 0                # legend block height (now 3 stacked map rows)
FIG_H  <- BODY_H + LEG_H

# ---- legends (colour only; the figure prints no percentages) ----
legend_pers   <- build_type_legend(persistent_colors, "Persistent type", nrow = 1)
legend_change <- build_type_legend(change_colors,     "Transition type", nrow = 1)
legend_inner  <- build_type_legend(INNER_RING_COLORS, "Response class",  nrow = 1)
legend_line   <- build_line_legend()

# ---- one map with its ring chart inset into the bottom-left corner ----
build_map_inset <- function(prep, period_label, letter) {
  m <- render_map_panel(prep, period_label, letter,
                        title_size = MAP_TITLE_SIZE, axis_text_size = MAP_AXIS_SIZE)
  if (is.null(prep) || is.null(prep$bar_stats)) return(m)
  s <- sunburst_inset(prep$bar_stats)
  m + inset_element(s, left = INSET_L, bottom = INSET_B,
                    right = INSET_R, top = INSET_T,
                    align_to = "panel", on_top = TRUE,
                    clip = FALSE, ignore_tag = TRUE)
}

map_letters <- c("a", "b", "c", "d")
map_panels  <- lapply(seq_along(TARGET_PERIODS), function(k)
  build_map_inset(combined_prep[[TARGET_PERIODS[k]]], TARGET_PERIODS[k], map_letters[k]))
left_col <- Reduce(`/`, map_panels)                    # 4 maps stacked

# ---- the 3 line plots (right column) ----
line_a <- build_lineplot_panel(df_plot_line, "ALL",    "e", 87.5,
                               "Combination of energy\n  and water limiting factors")
line_b <- build_lineplot_panel(df_plot_line, "Energy", "f", 50,
                               "Energy-related factors")
line_c <- build_lineplot_panel(df_plot_line, "Water",  "g", 50,
                               "Moisture-related factors")
right_col <- line_a / line_b / line_c

# ---- body (2 columns) + shared legend row underneath ----
body <- (left_col | right_col) + plot_layout(widths = COL_W_RATIO)

# 3 map legends stacked (one category per row) under the map column; the
# climate-factor legend sits under the line column. Widths match the body 3:1.
map_legends_col <- wrap_elements(full = legend_pers)   /
  wrap_elements(full = legend_change) /
  wrap_elements(full = legend_inner)
legends_row <- (map_legends_col | wrap_elements(full = legend_line)) +
  plot_layout(widths = COL_W_RATIO)

final_fig <- (body / legends_row) + plot_layout(heights = c(BODY_H, LEG_H))
final_fig <- body

ggsave(file.path(figure_date_dir, "Fig4_maps_lines_combined3.png"), final_fig,
       width = FIG_W, height = FIG_H, dpi = 300, bg = "white", limitsize = FALSE)
ggsave(file.path(figure_date_dir, "Fig4_maps_lines_combined3.pdf"), final_fig,
       width = FIG_W, height = FIG_H, device = "pdf", limitsize = FALSE)
cat("  saved: Fig4_maps_lines_combined.png/pdf\n\n")


ggsave(file.path(figure_date_dir, "Fig4_legend.png"), legends_row,
       width = 30, height = 3,dpi = 300, limitsize = FALSE)

# ============================================================================ #
# 13. DIAGNOSTIC STAT-ONLY BAR CHARTS (for checking the numbers only)
# ============================================================================ #
cat("Saving diagnostic check bar charts...\n")

p_check_pers <- ggplot(df_pers_stats_all,
                       aes(x = factor(period, levels = TARGET_PERIODS),
                           y = percentage, fill = persistent_type)) +
  geom_col(position = position_dodge(0.8), width = 0.75) +
  geom_text(aes(label = sprintf("%.1f", percentage)),
            position = position_dodge(0.8), vjust = -0.3, size = 3) +
  scale_fill_manual(values = persistent_colors, name = "Persistent type") +
  labs(x = NULL, y = "Area %", title = "Persistent type composition by period (check)") +
  theme_bw()
ggsave(file.path(figure_date_dir, "check_persistent_bars.png"),
       p_check_pers, width = 9, height = 5, dpi = 200, bg = "white")

p_check_change <- ggplot(df_change_stats_all,
                         aes(x = factor(period, levels = TARGET_PERIODS),
                             y = percentage, fill = change_type)) +
  geom_col(position = position_dodge(0.9), width = 0.85) +
  scale_fill_manual(values = change_colors, name = "Transition type") +
  labs(x = NULL, y = "Area %", title = "Transition type composition by period (check)") +
  theme_bw()
ggsave(file.path(figure_date_dir, "check_transition_bars.png"),
       p_check_change, width = 10, height = 5, dpi = 200, bg = "white")

cat("\n=== Done. Main figure + checks written to:\n    ", figure_date_dir, "\n")
