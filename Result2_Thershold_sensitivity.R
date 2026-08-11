

# ============================================================================ #
#  0. USER CONFIG  -- edit these paths/parameters, then `Rscript` the file
# ============================================================================ #
setwd("/data/share/tangwenxi")              

## ---- 0.1 Indicator / climate factors / period -----------------------------
VI_NAME    <- "GPP"                          
CLI_NAMES  <- c("TMP", "PRE", "SR", "VPD", "CO2", "SMroot", "SMsurf")
YEAR_START <- 1982
YEAR_END   <- 2020
DUR_MAX    <- 12                              

## ---- 0.2 INPUT 1: code #1 (VAR) per-group outputs --------------------------
VAR_GROUPS_DIR <- "./Out-Table-and-Figure-VPD et al/2026-06-26/History_GroupBy_100_SlidingWindow_39yr_1982_2020/Hist_GPP/Window_1982_2020"   ### CHECK ###
URATION_COMBINED_CSV <- NA   

## ---- 0.3 INPUT 2: raw monthly climate tables -------------------------------
CLIMATE_CSV_DIR <- "./Data/Monthly_df_1982_2020"   
CLI_FILE_SUFFIX <- c(TMP = "TMP", PRE = "PRE", SR = "MeanSR", VPD = "VPD",
                     CO2 = "CO2", SMroot = "SMroot", SMsurf = "SMsurf")  

## ---- 0.4 OUTPUT ------------------------------------------------------------
OUT_ROOT <- file.path("./Out-Table-and-Figure-VPD et al",
                      as.character(Sys.Date()),
                      "Threshold_Pipeline_Server")
CLI_WIDE_DIR <- file.path(OUT_ROOT, "CLI_wide")    
RESULT_DIR   <- file.path(OUT_ROOT, paste0(VI_NAME, "_breakpoint"))
SD_DIR       <- file.path(OUT_ROOT, "Threshold_sd") 

## ---- 0.5 Parallel / resource knobs ----------------------------------------
CORE_FRACTION <- 0.80      
CORE_HARD_CAP <- 80        
CORE_MIN      <- 2         
LOW_MEM_GB    <- 10       

for (d in c(OUT_ROOT, CLI_WIDE_DIR, RESULT_DIR, SD_DIR))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
ERROR_FILE <- file.path(OUT_ROOT, "error_log.txt")
cat("Pipeline started:", as.character(Sys.time()), "\n", file = ERROR_FILE)


# ============================================================================ #
#  1. PACKAGES
# ============================================================================ #
suppressPackageStartupMessages({
  library(parallel)
  library(doParallel)
  library(foreach)
  library(dplyr)
  library(tidyr)
  library(reshape2)
  library(segmented)
  library(data.table)
})


# ============================================================================ #
#  2. HELPER FUNCTIONS
# ============================================================================ #

## ---- 2.1 cores: 80% of idle cores -----------------------------------------
get_use_cores <- function(fraction = CORE_FRACTION,
                          hard_cap = CORE_HARD_CAP,
                          min_cores = CORE_MIN) {
  total <- parallel::detectCores()
  load1 <- tryCatch(
    as.numeric(system("cat /proc/loadavg | awk '{print $1}'", intern = TRUE)),
    error = function(e) 0)
  if (is.na(load1)) load1 <- 0
  idle  <- max(0, total - load1)                 
  use   <- floor(idle * fraction)
  use   <- max(min_cores, min(hard_cap, use))
  cat(sprintf("[cores] total=%d | 1min-load=%.2f | idle~=%.1f | using %d (%.0f%% of idle)\n",
              total, load1, idle, use, fraction * 100))
  use
}

## ---- 2.2 free memory (GB) -------------------------------------------------
free_mem_gb <- function() {
  out <- tryCatch(system("free -g", intern = TRUE), error = function(e) NULL)
  if (is.null(out)) return(Inf)
  as.numeric(strsplit(trimws(out[2]), "\\s+")[[1]][7])
}
maybe_gc <- function() {
  if (free_mem_gb() < LOW_MEM_GB) { gc(); cat("[mem] low memory -> gc()\n") }
}

## ---- 2.3 text progress bar ------------------------------------------------
progress_bar <- function(cur, total, prefix = "Progress", width = 40) {
  frac   <- if (total > 0) cur / total else 1
  filled <- round(frac * width)
  bar    <- paste0(strrep("#", filled), strrep("-", width - filled))
  cat(sprintf("\r%s [%s] %5.1f%% (%d/%d)", prefix, bar, frac * 100, cur, total))
  if (cur >= total) cat("\n")
  flush.console()
}

log_err <- function(msg) {
  cat(strrep("=", 70), "\n", as.character(Sys.time()), "\n", msg, "\n",
      strrep("=", 70), "\n\n", file = ERROR_FILE, append = TRUE)
}

## ---- 2.4 lagged-climate-exposure (LCE) for ONE grid -----------------------
extract_LCE_grid <- function(veg_df, cli_df, duration) {
  veg_df$Month <- as.numeric(as.character(veg_df$Month))
  years_all <- unique(veg_df$year_new)
  out <- vector("list", length(years_all))
  
  for (ri in seq_along(years_all)) {
    year     <- years_all[ri]
    veg_row  <- veg_df[veg_df$year_new == year, ]
    tgt      <- veg_row$Month                      
    prev     <- duration
    start_m  <- as.numeric(tgt) - prev
    end_m    <- as.numeric(tgt) - 1
    start_a  <- start_m[1]
    end_a    <- end_m[length(end_m)]
    
    # ---- end month / end year of the lag window ----
    if (tgt[length(tgt)] == 1 && end_a == 0) {
      end_month <- 12 + end_a; end_year <- year - 1
    } else {
      end_month <- end_a;      end_year <- year
    }
    
    # ---- is the GROWING SEASON itself cross-year? ----
    mid_year <- NA; months_range_mid <- NA
    if (12 %in% tgt[1:(length(tgt) - 1)]) {
      gs_cross <- TRUE
      if (start_a <= 0) {
        start_month <- 12 + start_a; start_year <- year - 2
        mid_year <- year - 1;        months_range_mid <- 1:12
      } else {
        start_month <- start_a;      start_year <- year - 1
      }
    } else {
      gs_cross <- FALSE
      if (start_a <= 0) {
        start_month <- 12 + start_a; start_year <- year - 1
      } else {
        start_month <- start_a;      start_year <- year
      }
    }
    
    # ---- pull the climate values over the lag window ----
    if (start_year != end_year) {                   
      if (is.na(months_range_mid[1]) && is.na(mid_year)) {
        mr1 <- start_month:12; mr2 <- 1:end_month
        months_range <- c(mr1, mr2)
        years_range  <- c(rep(start_year, length(mr1)), rep(end_year, length(mr2)))
        if (start_year < YEAR_START && end_year < YEAR_START) {
          cli_vals <- rep(NA_real_, length(months_range))
        } else if (start_year < YEAR_START && end_year >= YEAR_START) {
          cli_vals <- c(rep(NA_real_, length(mr1)),
                        as.numeric(cli_df[cli_df$Year == end_year,   mr2 + 1]))
        } else {
          cli_vals <- c(as.numeric(cli_df[cli_df$Year == start_year, mr1 + 1]),
                        as.numeric(cli_df[cli_df$Year == end_year,   mr2 + 1]))
        }
      } else {                                       
        mr1 <- start_month:12; mr2 <- 1:end_month
        months_range <- c(mr1, months_range_mid, mr2)
        years_range  <- c(rep(start_year, length(mr1)),
                          rep(mid_year,   length(months_range_mid)),
                          rep(end_year,   length(mr2)))
        if (start_year < YEAR_START && end_year < YEAR_START) {
          cli_vals <- rep(NA_real_, length(months_range))
        } else if (start_year < YEAR_START && end_year >= YEAR_START) {
          cli_vals <- c(rep(NA_real_, length(mr1)),
                        as.numeric(cli_df[cli_df$Year == mid_year,  months_range_mid + 1]),
                        as.numeric(cli_df[cli_df$Year == end_year,  mr2 + 1]))
        } else {
          cli_vals <- c(as.numeric(cli_df[cli_df$Year == start_year, mr1 + 1]),
                        as.numeric(cli_df[cli_df$Year == mid_year,    months_range_mid + 1]),
                        as.numeric(cli_df[cli_df$Year == end_year,    mr2 + 1]))
        }
      }
    } else {                                         
      months_range <- start_month:end_month
      years_range  <- rep(start_year, length(months_range))
      if (start_year < YEAR_START) {
        cli_vals <- rep(NA_real_, length(months_range))
      } else {
        cli_vals <- as.numeric(cli_df[cli_df$Year == start_year, months_range + 1])
      }
    }
    
    out[[ri]] <- data.frame(year_new  = year,
                            LCE_year  = years_range,
                            LCE_month = months_range,
                            LCE_data  = cli_vals)
  }
  data.table::rbindlist(out)
}

## ---- 2.5 segmented regression SD(VI) ~ SD(CLI) for ONE grid ---------------
fit_segment_grid <- function(df_sd, lab, group_name) {
  base <- data.frame(lab = lab, Group_name = group_name,
                     stringsAsFactors = FALSE)
  attr_na <- function(sig) {
    cbind(base, data.frame(signif = sig, R2 = NA_real_,
                           breakpoint_x = NA_real_, breakpoint_y = NA_real_,
                           breakpoint_p = NA_real_,
                           psi1.x_CI.l = NA_real_, psi1.x_CI.u = NA_real_,
                           intercept1 = NA_real_, intercept2 = NA_real_,
                           slope1 = NA_real_, slope2 = NA_real_,
                           slope1_SE = NA_real_, slope2_SE = NA_real_,
                           pvalue1 = NA_real_, pvalue2 = NA_real_,
                           stringsAsFactors = FALSE))
  }
  # need enough distinct x to estimate a breakpoint
  if (nrow(df_sd) < 4 || length(unique(round(df_sd$sd_cli, 8))) < 3)
    return(list(row = attr_na(NA_character_), sd = df_sd))
  
  df  <- data.frame(x = df_sd$sd_cli, y = df_sd$sd_GS)
  fit <- lm(y ~ x, data = df)
  seg <- tryCatch(segmented::segmented(fit, seg.Z = ~ x),
                  error = function(e) NULL)
  if (is.null(seg) || is.null(seg$psi) || nrow(seg$psi) < 1)
    return(list(row = attr_na("F"), sd = df_sd))
  
  ss <- summary(seg)
  R2 <- round(ss$r.squared, 4)
  
  ic  <- segmented::intercept(seg)$x          
  i1  <- as.numeric(ic[1, 1]); i2 <- as.numeric(ic[2, 1])
  sl  <- segmented::slope(seg)$x             
  s1  <- as.numeric(sl[1, 1]); s2  <- as.numeric(sl[2, 1])
  s1e <- as.numeric(sl[1, 2]); s2e <- as.numeric(sl[2, 2])
  
  cf  <- ss$coefficients
  pv1 <- if ("x"    %in% rownames(cf)) cf["x",    "Pr(>|t|)"] else NA_real_
  pv2 <- if ("U1.x" %in% rownames(cf)) cf["U1.x", "Pr(>|t|)"] else NA_real_
  
  bx  <- round(as.numeric(seg$psi[1, "Est."]), 4)
  by  <- i1 + s1 * bx
  dv  <- tryCatch(segmented::davies.test(fit, k = 10), error = function(e) NULL)
  bp  <- if (is.null(dv)) NA_real_ else dv$p.value
  ci  <- tryCatch(as.data.frame(confint(seg)), error = function(e) NULL)
  cil <- if (is.null(ci)) NA_real_ else as.numeric(ci[1, 2])
  ciu <- if (is.null(ci)) NA_real_ else as.numeric(ci[1, 3])
  
  sig <- if (!is.na(pv1) && pv1 < 0.1) "T" else if (!is.na(pv2) && pv2 < 0.1) "T" else "F"
  
  row <- cbind(base, data.frame(
    signif = sig, R2 = R2,
    breakpoint_x = bx, breakpoint_y = by, breakpoint_p = bp,
    psi1.x_CI.l = cil, psi1.x_CI.u = ciu,
    intercept1 = round(i1, 4), intercept2 = round(i2, 4),
    slope1 = round(s1, 4), slope2 = round(s2, 4),
    slope1_SE = round(s1e, 4), slope2_SE = round(s2e, 4),
    pvalue1 = pv1, pvalue2 = pv2,
    stringsAsFactors = FALSE))
  
  df_sd$fit <- tryCatch(segmented::broken.line(seg)$fit, error = function(e) NA_real_)
  list(row = row, sd = df_sd)
}


# ============================================================================ #
#  STEP A.  Merge code #1 per-group *_VAR_data_list.RData  ->  veg_min
#           (only year_new / Month / VI columns are kept -- all that the SD step needs)
# ============================================================================ #
cat("\n================ STEP A: merge VAR_data_list ================\n")

var_group_dirs <- list.dirs(VAR_GROUPS_DIR, recursive = FALSE)
var_group_dirs <- var_group_dirs[grepl("Group_", basename(var_group_dirs))]
var_group_dirs <- var_group_dirs[order(as.integer(sub("^Group_", "", basename(var_group_dirs))))]

head(basename(var_group_dirs))
length(var_group_dirs)

if (length(var_group_dirs) == 0)
  stop("No 'Group_*' folders found under VAR_GROUPS_DIR -- check the path.")
cat("Found", length(var_group_dirs), "group folders.\n")

veg_min <- list()        
n_grp   <- length(var_group_dirs)
for (gi in seq_len(n_grp)) {
  f <- list.files(var_group_dirs[gi], "_VAR_data_list\\.RData$", full.names = TRUE)
  for (ff in f) {
    e <- new.env()
    load(ff, envir = e)
    obj_name <- ls(e)[1]                 
    vdl <- get(obj_name, envir = e)
    for (d in vdl) {
      if (is.null(d) || !is.data.frame(d)) next
      if (!all(c("year_new", "Month", VI_NAME, "lab") %in% names(d))) next
      lab <- as.character(unique(d$lab))[1]
      veg_min[[lab]] <- d[, c("year_new", "Month", VI_NAME)]
    }
  }
  progress_bar(gi, n_grp, "  VAR_data merge")
}
cat("VARdata grids loaded:", length(veg_min), "\n")
maybe_gc()

## ---- STEP A.1  overall SD(VI) per grid (CLI-independent) -> GPP_sd.csv ------
cat("Writing overall VI SD per grid ...\n")
vi_sd_df <- data.frame(
  lab = as.numeric(names(veg_min)),
  vi  = VI_NAME,
  sd  = vapply(veg_min, function(d) sd(d[[VI_NAME]], na.rm = TRUE), numeric(1)),
  stringsAsFactors = FALSE)
vi_sd_df$sd2 <- 2 * vi_sd_df$sd
vi_sd_df$sd3 <- 3 * vi_sd_df$sd
data.table::fwrite(vi_sd_df, file.path(SD_DIR, paste0(VI_NAME, "_sd.csv")))
rm(vi_sd_df)


# ============================================================================ #
#  STEP B.  Merge code #1 per-group *_duration.csv  ->  valid grids + per-CLI duration
# ============================================================================ #
cat("\n================ STEP B: merge duration ================\n")

if (!is.na(DURATION_COMBINED_CSV) && file.exists(DURATION_COMBINED_CSV)) {
  dur_all <- as.data.frame(data.table::fread(DURATION_COMBINED_CSV))
  cat("Loaded combined duration:", DURATION_COMBINED_CSV, "\n")
} else {
  dur_files <- list.files(var_group_dirs, "_duration\\.csv$", full.names = TRUE)
  cat("Reading", length(dur_files), "per-group duration files...\n")
  dur_all <- data.table::rbindlist(lapply(dur_files, data.table::fread),
                                   use.names = TRUE, fill = TRUE)
  dur_all <- as.data.frame(dur_all)
}
# `variable` holds the impulse name (TMP/PRE/.../GPP); `duration` is the lag length.
labs_big  <- unique(dur_all$lab[which(dur_all$duration > DUR_MAX)])
labs_na   <- unique(dur_all$lab[is.na(dur_all$duration)])
labs_excl <- union(as.character(labs_big), as.character(labs_na))
cat("Grids excluded (duration > ", DUR_MAX, " or NA): ", length(labs_excl), "\n", sep = "")


# ============================================================================ #
#  STEP C.  Raw monthly climate CSV  ->  wide (Year x Month) per CLI, cached to disk
#           (re-implements code #7 part 2; vectorized + parallel over grids)
# ============================================================================ #
cat("\n================ STEP C: climate -> wide per CLI ================\n")

build_cli_wide <- function(df_cli) {
  # df_cli: col1 = lab, remaining cols = <prefix>_YYYY_MM
  id   <- as.character(df_cli[[1]])
  vcol <- names(df_cli)[-1]
  yr   <- as.integer(sub(".*_(\\d{4})_\\d{2}$", "\\1", vcol))
  mo   <- as.integer(sub(".*_(\\d{4})_(\\d{2})$", "\\2", vcol))
  ord  <- order(yr, mo)                              # year-major ordering
  years   <- sort(unique(yr)); months <- sort(unique(mo))
  M    <- as.matrix(df_cli[, -1, drop = FALSE][, ord, drop = FALSE])
  ny   <- length(years); nm <- length(months)
  setNames(lapply(seq_len(nrow(M)), function(i) {
    w <- as.data.frame(matrix(as.numeric(M[i, ]), nrow = ny, ncol = nm, byrow = TRUE))
    names(w) <- as.character(months)
    cbind(Year = years, w)
  }), id)
}

n_cores_C <- get_use_cores()
cl <- makeCluster(n_cores_C, outfile = file.path(OUT_ROOT, "parallel_log.txt"))
registerDoParallel(cl)

for (ci in seq_along(CLI_NAMES)) {
  cli  <- CLI_NAMES[ci]
  cache <- file.path(CLI_WIDE_DIR, paste0(cli, "_wide.RData"))
  if (file.exists(cache)) { cat("  [", cli, "] cached, skip.\n", sep = ""); next }
  
  suffix <- CLI_FILE_SUFFIX[[cli]]
  fpath  <- list.files(CLIMATE_CSV_DIR, paste0("_", suffix, "\\.csv$"), full.names = TRUE)
  if (length(fpath) == 0) {
    log_err(paste0("Climate file for ", cli, " (suffix _", suffix, ".csv) not found."))
    cat("  [", cli, "] climate file NOT FOUND -- skipped.\n", sep = ""); next
  }
  df_cli <- as.data.frame(data.table::fread(fpath[1]))
  df_cli <- df_cli[complete.cases(df_cli), ]        
  
  # parallel over grid-chunks
  idx    <- seq_len(nrow(df_cli))
  chunks <- split(idx, cut(seq_along(idx), n_cores_C, labels = FALSE))
  pieces <- foreach(rows = chunks, .export = "build_cli_wide",
                    .packages = character(0)) %dopar% {
                      build_cli_wide(df_cli[rows, , drop = FALSE])
                    }
  CLI_wide <- do.call(c, pieces)
  save(CLI_wide, file = cache)
  cat("  [", cli, "] wide grids: ", length(CLI_wide), " -> ", cache, "\n", sep = "")
  rm(df_cli, pieces, CLI_wide); maybe_gc()
}
stopCluster(cl)


# ============================================================================ #
#  STEP D.  LCE + SD + segmented regression  (fused #8 + #9 part 1, parallel)
#           For each CLI, for each grid:
#             1) extract lagged climate over `duration` months -> LCE
#             2) SD(VI) per year   &   SD(CLI) per year
#             3) segmented regression  SD(VI) ~ SD(CLI)  -> breakpoint
# ============================================================================ #
cat("\n================ STEP D: LCE + SD + segmented regression ================\n")

n_cli <- length(CLI_NAMES)
for (ci in seq_along(CLI_NAMES)) {
  cli <- CLI_NAMES[ci]
  cat(sprintf("\n---- [%d/%d] %s_%s ----\n", ci, n_cli, VI_NAME, cli))
  
  cache <- file.path(CLI_WIDE_DIR, paste0(cli, "_wide.RData"))
  if (!file.exists(cache)) { cat("  no wide climate cache, skip.\n"); next }
  load(cache)                                        # -> CLI_wide
  
  # per-grid duration for THIS climate factor (excluding bad grids)
  dcli <- dur_all[dur_all$variable == cli & !(as.character(dur_all$lab) %in% labs_excl),
                  c("lab", "duration")]
  dcli <- dcli[!duplicated(dcli$lab), ]
  dur_vec <- setNames(dcli$duration, as.character(dcli$lab))
  
  # grids that have veg data AND climate data AND a valid duration
  valid <- Reduce(intersect, list(names(veg_min), names(CLI_wide), names(dur_vec)))
  cat("  valid grids:", length(valid), "\n")
  if (length(valid) == 0) next
  
  # build per-grid input slices (only what the worker needs) and chunk them
  n_cores_D <- get_use_cores()
  chunks <- split(valid, cut(seq_along(valid), n_cores_D, labels = FALSE))
  veg_ch <- lapply(chunks, function(ls) veg_min[ls])
  cli_ch <- lapply(chunks, function(ls) CLI_wide[ls])
  dur_ch <- lapply(chunks, function(ls) dur_vec[ls])
  lab_ch <- chunks
  rm(CLI_wide); gc()
  
  cl <- makeCluster(n_cores_D, outfile = file.path(OUT_ROOT, "parallel_log.txt"))
  registerDoParallel(cl)
  
  res <- foreach(vc = veg_ch, cc = cli_ch, dc = dur_ch, lc = lab_ch,
                 .packages = c("segmented", "data.table"),
                 .export   = c("extract_LCE_grid", "fit_segment_grid",
                               "VI_NAME", "YEAR_START", "cli"),
                 .errorhandling = "pass") %dopar% {
                   
                   rows <- vector("list", length(lc))
                   sds  <- vector("list", length(lc))
                   osd  <- vector("list", length(lc))  
                   grp  <- "G1"
                   for (k in seq_along(lc)) {
                     lab <- lc[k]
                     veg <- vc[[lab]]; cliw <- cc[[lab]]; dur <- dc[[lab]]
                     if (is.null(veg) || is.null(cliw) || is.na(dur)) next
                     
                     # --- SD(VI) per growing-season year ---
                     vi_vals  <- as.numeric(veg[[VI_NAME]])
                     grp_year <- veg$year_new
                     sd_vi <- data.frame(
                       year_new = as.numeric(names(tapply(vi_vals, grp_year, mean))),
                       mean_GS  = as.numeric(tapply(vi_vals, grp_year, mean)),
                       sd_GS    = as.numeric(tapply(vi_vals, grp_year, sd)),
                       max_GS   = as.numeric(tapply(vi_vals, grp_year, max)))
                     
                     # --- lagged climate (LCE) then SD(CLI) per year ---
                     lce <- tryCatch(extract_LCE_grid(veg, cliw, dur), error = function(e) NULL)
                     if (is.null(lce) || nrow(lce) == 0) next
                     
                     # overall SD(CLI) over ALL lagged months/years (single value per grid):
                     # this is exactly what Fig3.Threshold_plot uses on the x-axis of b/c/d.
                     osd[[k]] <- data.frame(lab = as.numeric(lab), cli = cli, vi = VI_NAME,
                                            sd = sd(lce$LCE_data, na.rm = TRUE),
                                            stringsAsFactors = FALSE)
                     
                     sd_cli <- data.frame(
                       year_new = as.numeric(names(tapply(lce$LCE_data, lce$year_new,
                                                          function(z) mean(z, na.rm = TRUE)))),
                       mean_cli = as.numeric(tapply(lce$LCE_data, lce$year_new,
                                                    function(z) mean(z, na.rm = TRUE))),
                       sd_cli   = as.numeric(tapply(lce$LCE_data, lce$year_new,
                                                    function(z) sd(z, na.rm = TRUE))))
                     
                     df_sd <- merge(sd_vi, sd_cli, by = "year_new")
                     df_sd <- df_sd[complete.cases(df_sd[, c("sd_GS", "sd_cli")]), ]
                     if (nrow(df_sd) == 0) next
                     
                     fitres <- fit_segment_grid(df_sd, lab, grp)
                     rows[[k]] <- fitres$row
                     sds[[k]]  <- fitres$sd
                   }
                   list(rows = data.table::rbindlist(Filter(Negate(is.null), rows),
                                                     use.names = TRUE, fill = TRUE),
                        sds  = setNames(sds, lc),
                        osd  = data.table::rbindlist(Filter(Negate(is.null), osd),
                                                     use.names = TRUE, fill = TRUE))
                 }
  stopCluster(cl)
  
  # collect results (drop any worker that errored)
  ok       <- Filter(function(x) is.list(x) && !is.null(x$rows), res)
  rows_all <- data.table::rbindlist(lapply(ok, `[[`, "rows"),
                                    use.names = TRUE, fill = TRUE)
  sd_all   <- do.call(c, lapply(ok, function(x) Filter(Negate(is.null), x$sds)))
  osd_all  <- data.table::rbindlist(lapply(ok, `[[`, "osd"),
                                    use.names = TRUE, fill = TRUE)
  
  out_csv <- file.path(RESULT_DIR, paste0(VI_NAME, "_", cli, "_breakpoint_all.csv"))
  # `j` is a passthrough integer index expected by Fig3.Threshold_plot's e-g block
  # (must be a real integer, not NA, or that block's na.omit() drops every row).
  if (nrow(rows_all) > 0) rows_all$j <- seq_len(nrow(rows_all))
  data.table::fwrite(rows_all, out_csv)
  save(sd_all, file = file.path(RESULT_DIR, paste0(VI_NAME, "_", cli, "_SDdata.RData")))
  # overall per-grid climate SD -> matches Fig3.Threshold_plot's <VI>_<CLI>_sd.csv
  osd_all$sd2 <- 2 * osd_all$sd; osd_all$sd3 <- 3 * osd_all$sd
  data.table::fwrite(osd_all, file.path(SD_DIR, paste0(VI_NAME, "_", cli, "_sd.csv")))
  cat(sprintf("  -> %s  (rows: %d, with breakpoint: %d)\n",
              basename(out_csv), nrow(rows_all),
              sum(!is.na(rows_all$breakpoint_x))))
  progress_bar(ci, n_cli, "  overall (CLI)")
  rm(res, ok, rows_all, sd_all, osd_all, veg_ch, cli_ch); gc()
}


# ============================================================================ #
#  DONE
# ============================================================================ #
cat("\n========================================================\n")
cat("All done.\n")
cat("Outputs to copy to your LOCAL machine for the figure:\n")
cat("  (a) per-grid breakpoint   -> ", file.path(RESULT_DIR, paste0(VI_NAME, "_<CLI>_breakpoint_all.csv")), "\n")
cat("  (b-d) overall per-grid SD -> ", file.path(SD_DIR,     paste0(VI_NAME, "_<CLI>_sd.csv")), "\n")
cat("                              ", file.path(SD_DIR,     paste0(VI_NAME, "_sd.csv")), "\n")
cat("Then locally: rasterize breakpoint_x (Local_Rasterize_and_Cut.R) and run Fig3.Threshold_plot.\n")
cat("Error log:", ERROR_FILE, "\n")
cat("Finished:", as.character(Sys.time()), "\n")
