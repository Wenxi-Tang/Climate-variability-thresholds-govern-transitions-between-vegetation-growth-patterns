

library(data.table)
library(terra)
library(ggplot2)
library(patchwork)
library(sf)
library(rnaturalearth)
library(dplyr)

set.seed(20240101)

# ---- Section 0: Configuration -----------------------------------------------

IRF_ROOT     <- "H:/2026-06-18"
TEMPLATE_TIF <- "G:/03-Program/GlobalVeg/NCC-一审/Nas/resample_example_0.5degree.tif"
LULC_TIF     <- "G:/03-Program/GlobalVeg/NCC-一审/Nas/GLC_FCS30D_stable_0p5deg_1985_2020_8class.tif"
OUT_DIR      <- "G:/03-Program/GlobalVeg/Out-Table-and-Figure-VPD et al/FigS5-3"

PERIOD_EARLY <- "1982_2000"
PERIOD_LATE  <- "2001_2020"
YEAR_TAG     <- "2000"
VI_LIST      <- c("GPP", "LAI", "SIF", "NDVI")

ANCHOR_VI    <- "GPP"   # index whose crossing point defines the percentile
ROUND_ANCHOR <- 1       # decimals for the anchor threshold
SIGNIF_OTHER <- 3       # significant digits for the other indices
N_CURVE_PT   <- 22      # markers per sensitivity curve

CLS_KEY <- c("VSO", "Stable", "VGC")
CLS_COL <- c(VSO = "#4DAF4A", Stable = "#56a0d3", VGC = "#DA5CA3")
DISP    <- c(VSO = "\u03b4VSO", Stable = "Stable", VGC = "\u03b4VGC")
CODE_OF <- c(VSO = 1L, VGC = 2L, Stable = 3L)

COL_THR <- "#B22222"    # threshold reference line
TAGS    <- c("(a)", "(b)", "(c)", "(d)")

PS      <- 26           # minimum font size for base graphics
RES     <- 200
FIG_DPI <- 400

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

r_tpl  <- terra::rast(TEMPLATE_TIF)
r_lulc <- terra::rast(LULC_TIF)
stopifnot(terra::ncell(r_tpl) == terra::ncell(r_lulc))
v_lulc <- terra::values(r_lulc, mat = FALSE)


# ---- Section 0b: Core computation -------------------------------------------

read_irf <- function(vi, period) {
  f <- list.files(IRF_ROOT, pattern = paste0(vi, "_Window_", period, "_irf\\.csv$"),
                  recursive = TRUE, full.names = TRUE)
  if (!length(f))
    f <- list.files(IRF_ROOT, pattern = paste0(vi, "_G.*_window_", period, "_irf\\.csv$"),
                    recursive = TRUE, full.names = TRUE)
  if (!length(f)) stop("No IRF files for ", vi, " ", period)
  dt <- rbindlist(lapply(f, fread), use.names = TRUE, fill = TRUE)
  if ("V1" %in% names(dt)) dt[, V1 := NULL]
  dt[variable == vi][, Time := as.integer(as.character(Time))][]
}

# Per-cell statistics. Delta_min is the smallest equivalence limit at which a
# cell passes the TOST, i.e. |mean_d| + t(0.95, n-1) * se_d.
grid_stats <- function(vi) {
  e <- read_irf(vi, PERIOD_EARLY)[, .(irf_e = irf), by = .(lab, Time)]
  l <- read_irf(vi, PERIOD_LATE )[, .(irf_l = irf), by = .(lab, Time)]
  m <- merge(e, l, by = c("lab", "Time"))
  m[, lab := as.integer(lab)][, d_val := irf_l - irf_e]
  setorder(m, lab, Time)
  
  gs <- m[, {
    d <- d_val; i <- which.max(abs(d))
    .(n_time = .N, mean_d = mean(d), sd_d = sd(d), max_d = d[i])
  }, by = lab]
  gs <- gs[n_time >= 2 & is.finite(mean_d) & is.finite(sd_d)]
  gs[, Delta_min := abs(mean_d) + qt(0.95, n_time - 1) * sd_d / sqrt(n_time)]
  gs[, dir_code  := fifelse(max_d > 0, CODE_OF[["VGC"]], CODE_OF[["VSO"]])]
  gs[!is.na(v_lulc[lab])][]                    # land-cover mask
}

classify <- function(gs, Delta)
  fifelse(gs$Delta_min <= Delta, CODE_OF[["Stable"]], gs$dir_code)

proportions_at <- function(gs, Delta) {
  cl <- classify(gs, Delta)
  c(VSO    = 100 * mean(cl == CODE_OF[["VSO"]]),
    Stable = 100 * mean(cl == CODE_OF[["Stable"]]),
    VGC    = 100 * mean(cl == CODE_OF[["VGC"]]))
}

sens_curve <- function(gs, n_pt = N_CURVE_PT) {
  qs <- unname(quantile(gs$Delta_min, seq(0.01, 0.99, length.out = n_pt)))
  rbindlist(lapply(qs, function(D) {
    p <- proportions_at(gs, D)
    data.table(Delta = D, VSO = p[["VSO"]], Stable = p[["Stable"]], VGC = p[["VGC"]])
  }))
}

# Limit where the Stable and VSO curves intersect; both are monotone in Delta.
crossing_point <- function(gs) {
  f <- function(D) mean(gs$Delta_min <= D) -
    mean(gs$Delta_min > D & gs$dir_code == CODE_OF[["VSO"]])
  lo <- min(gs$Delta_min) * 0.5; hi <- max(gs$Delta_min)
  if (f(lo) > 0 || f(hi) < 0) return(NA_real_)
  uniroot(f, c(lo, hi), tol = .Machine$double.eps^0.6)$root
}

message("Computing per-cell statistics ...")
GS <- lapply(VI_LIST, grid_stats); names(GS) <- VI_LIST
for (v in VI_LIST) message("  ", v, ": ", nrow(GS[[v]]), " cells")


# ---- Section 1: Anchor threshold and its percentile -------------------------

gs_a      <- GS[[ANCHOR_VI]]
cross_raw <- crossing_point(gs_a)
if (is.na(cross_raw)) stop("Stable and VSO curves do not intersect for ", ANCHOR_VI)

Delta_anchor <- round(cross_raw, ROUND_ANCHOR)
PCT_ANCHOR   <- mean(abs(gs_a$mean_d) <= Delta_anchor)
p_anchor     <- proportions_at(gs_a, Delta_anchor)

cat(sprintf("\n[Section 1] %s crossing = %.4f -> rounded %.*f\n",
            ANCHOR_VI, cross_raw, ROUND_ANCHOR, Delta_anchor))
cat(sprintf("            percentile in |mean_d|: P%.1f\n", 100 * PCT_ANCHOR))
cat(sprintf("            proportions: VSO %.1f%% | Stable %.1f%% | VGC %.1f%%\n",
            p_anchor[["VSO"]], p_anchor[["Stable"]], p_anchor[["VGC"]]))


# ---- Section 2: Thresholds for all indices and sensitivity panels ------------
DELTA <- sapply(VI_LIST, function(v) {
  if (v == ANCHOR_VI) return(Delta_anchor)
  ad <- abs(GS[[v]]$mean_d); ad <- ad[ad > 0]      # drop exact zeros
  signif(unname(quantile(ad, PCT_ANCHOR)), SIGNIF_OTHER)
})

fmt_d <- function(v) if (v == ANCHOR_VI) sprintf("%.*f", ROUND_ANCHOR, DELTA[[v]]) else
  formatC(DELTA[[v]], format = "g", digits = SIGNIF_OTHER)

thr_tbl <- data.table(
  VI           = VI_LIST,
  n_cells      = sapply(VI_LIST, function(v) nrow(GS[[v]])),
  median_abs_d = sapply(VI_LIST, function(v) median(abs(GS[[v]]$mean_d))),
  Delta        = unname(DELTA),
  percentile   = sapply(VI_LIST, function(v) 100 * mean(abs(GS[[v]]$mean_d) <= DELTA[[v]])))
print(thr_tbl, digits = 4)
fwrite(thr_tbl, file.path(OUT_DIR, "thresholds.csv"))

open_png <- function(fn, w = 4200, h = 3400)
  png(file.path(OUT_DIR, fn), width = w, height = h, res = RES, pointsize = PS)

# par(mfrow) shrinks cex to 0.83; reset it so PS is the true minimum size.
setup_par <- function(mar) par(mfrow = c(2, 2), mar = mar, mgp = c(3.0, 0.9, 0),
                               cex = 1, cex.axis = 1, cex.lab = 1.05,
                               cex.main = 1.05, lend = 1, tcl = -0.4)
add_tag <- function(i) mtext(TAGS[i], side = 3, adj = 0, line = 0.6, font = 2, cex = 1.1)

# Single reference line at the adopted threshold, annotated with its percentile.
thr_line <- function(v) abline(v = DELTA[[v]], col = COL_THR, lwd = 3, lty = 2)
thr_text <- function(v) sprintf("%s = %s (P%.0f)", "\u0394", fmt_d(v), 100 * PCT_ANCHOR)

# Panel 1: stable fraction versus the equivalence limit
open_png("panel_fig1_stable_vs_Delta.png")
setup_par(c(5.2, 6.0, 3.2, 2.0))
for (i in seq_along(VI_LIST)) {
  v <- VI_LIST[i]; gs <- GS[[v]]
  dm <- sort(gs$Delta_min); xm <- unname(quantile(dm, 0.99))
  plot(dm, 100 * seq_along(dm) / length(dm), type = "l", lwd = 3,
       xlim = c(0, xm), ylim = c(0, 100), xaxs = "i", yaxs = "i",
       xlab = expression(Delta), ylab = "Stable grid cells (%)", main = v)
  thr_line(v)
  legend("bottomright", bty = "n", seg.len = 1.4, col = COL_THR, lwd = 3, lty = 2,
         legend = thr_text(v))
  add_tag(i)
}
dev.off()

# Panel 2: class composition versus the equivalence limit
open_png("panel_fig2_composition.png")
setup_par(c(5.2, 6.0, 3.2, 2.0))
for (i in seq_along(VI_LIST)) {
  v <- VI_LIST[i]; gs <- GS[[v]]
  s  <- sens_curve(gs); xm <- unname(quantile(gs$Delta_min, 0.99))
  s  <- s[Delta <= xm]
  plot(s$Delta, s$Stable, type = "b", pch = 16, lwd = 3, cex = 0.9,
       col = CLS_COL[["Stable"]], xlim = c(0, xm), ylim = c(0, 100), yaxs = "i",
       xlab = expression(Delta), ylab = "Proportion (%)", main = v)
  lines(s$Delta, s$VGC, type = "b", pch = 17, lwd = 3, cex = 0.9, col = CLS_COL[["VGC"]])
  lines(s$Delta, s$VSO, type = "b", pch = 15, lwd = 3, cex = 0.9, col = CLS_COL[["VSO"]])
  thr_line(v)
  legend("topright", bty = "n", seg.len = 1.4, lwd = 3, pch = c(16, 17, 15),
         col = c(CLS_COL[["Stable"]], CLS_COL[["VGC"]], CLS_COL[["VSO"]]),
         legend = c("Stable", expression(delta*"VGC"), expression(delta*"VSO")))
  mtext(thr_text(v), side = 3, adj = 1, line = 0.6, cex = 0.85, col = COL_THR)
  add_tag(i)
}
dev.off()

# Panel 3: distribution of |mean_d| relative to the limit
open_png("panel_fig3_meand_hist.png")
setup_par(c(5.2, 6.6, 3.2, 2.0))
for (i in seq_along(VI_LIST)) {
  v <- VI_LIST[i]; x <- abs(GS[[v]]$mean_d); x <- x[is.finite(x)]
  hist(x[x <= quantile(x, 0.99)], breaks = 50, col = "grey86", border = NA,
       xlab = expression(paste("|", bar(d), "|")), ylab = "Grid cells", main = v)
  thr_line(v)
  legend("topright", bty = "n", seg.len = 1.4, col = COL_THR, lwd = 3, lty = 2,
         legend = thr_text(v))
  add_tag(i)
}
dev.off()


# ---- Section 3: Classification maps (Fig. S5) -------------------------------

world_shp <- ne_coastline(scale = "medium", returnclass = "sf")

to_raster <- function(labs, vals) {
  r <- terra::rast(r_tpl); terra::values(r) <- NA; r[labs] <- vals; r
}

pct_all <- data.table(); p_map <- list()

for (i in seq_along(VI_LIST)) {
  v  <- VI_LIST[i]; gs <- GS[[v]]
  cl <- classify(gs, DELTA[[v]])
  
  r <- to_raster(gs$lab, cl)
  terra::writeRaster(r, file.path(OUT_DIR, paste0(v, "-", YEAR_TAG, "-type_code.tif")),
                     overwrite = TRUE, datatype = "INT2S")
  terra::writeRaster(to_raster(gs$lab, gs$mean_d),
                     file.path(OUT_DIR, paste0(v, "-", YEAR_TAG, "-mean_d.tif")),
                     overwrite = TRUE, datatype = "FLT4S")
  
  df <- as.data.frame(r, xy = TRUE, na.rm = TRUE); names(df)[3] <- "value"
  # Legend order VSO, Stable, VGC corresponds to codes 1, 3, 2
  df$group <- factor(df$value, levels = c(1, 3, 2), labels = unname(DISP[CLS_KEY]))
  
  p <- proportions_at(gs, DELTA[[v]])
  pct_all <- rbind(pct_all, data.table(VI = v, Delta = DELTA[[v]],
                                       VSO = p[["VSO"]], Stable = p[["Stable"]],
                                       VGC = p[["VGC"]]))
  sub <- sprintf("%s %.1f%%   %s %.1f%%   %s %.1f%%",
                 DISP[["VSO"]], p[["VSO"]], DISP[["Stable"]], p[["Stable"]],
                 DISP[["VGC"]], p[["VGC"]])
  
  p_map[[i]] <- ggplot(df) +
    geom_hline(yintercept = 0, colour = "gray50", linetype = "dashed", linewidth = 0.5) +
    geom_tile(aes(x, y, fill = group)) +
    geom_sf(data = world_shp, fill = NA, colour = "gray30", linewidth = 0.2,
            inherit.aes = FALSE) +
    # unname(): a named vector would be matched to factor levels by name
    scale_fill_manual(name = NULL, values = unname(CLS_COL[CLS_KEY]),
                      na.value = "transparent", drop = FALSE) +
    coord_sf(ylim = c(-60, 90), expand = FALSE) +
    labs(title = paste(TAGS[i], v), subtitle = sub, x = NULL, y = NULL) +
    theme_bw(base_size = 12) +
    theme(plot.title    = element_text(face = "bold", size = 12),
          plot.subtitle = element_text(size = 10),
          panel.grid    = element_blank(),
          legend.position = "bottom")
}

fwrite(pct_all, file.path(OUT_DIR, "FigS5-class_percentages.csv"))
print(pct_all, digits = 4)

map_plot <- wrap_plots(p_map, ncol = 2, nrow = 2) +
  plot_layout(guides = "collect") & theme(legend.position = "bottom")
ggsave(file.path(OUT_DIR, "FigS5-map.png"), map_plot, width = 24, height = 15,
       units = "cm", dpi = FIG_DPI, bg = "white")
ggsave(file.path(OUT_DIR, "FigS5-map.pdf"), map_plot, width = 24, height = 15,
       units = "cm", device = cairo_pdf, bg = "white")


# ---- Section 4: Association and effect sizes (Table S1) ---------------------

cls_wide <- Reduce(function(a, b) merge(a, b, by = "lab"),
                   lapply(VI_LIST, function(v)
                     data.table(lab = GS[[v]]$lab, cls = classify(GS[[v]], DELTA[[v]]))[
                       , setnames(.SD, "cls", v)]))
cat("\n[Section 4] cells shared by all indices:", nrow(cls_wide), "\n")

as_cls  <- function(x) factor(x, levels = c(1, 2, 3), labels = c("VSO", "VGC", "Stable"))
targets <- setdiff(VI_LIST, "GPP")

X2_of <- function(tab) as.numeric(suppressWarnings(chisq.test(tab, correct = FALSE))$statistic)

# Cramer's V with Bergsma (2013) bias correction
cramers_v <- function(tab) {
  n <- sum(tab); r <- nrow(tab); k <- ncol(tab)
  phi2 <- max(0, X2_of(tab) / n - (k - 1) * (r - 1) / (n - 1))
  r <- r - (r - 1)^2 / (n - 1); k <- k - (k - 1)^2 / (n - 1)
  sqrt(phi2 / min(k - 1, r - 1))
}
cohens_w    <- function(tab) sqrt(X2_of(tab) / sum(tab))
kappa_stats <- function(tab) {
  n <- sum(tab); po <- sum(diag(tab)) / n
  pe <- sum(rowSums(tab) * colSums(tab)) / n^2
  c(po = po, kappa = (po - pe) / (1 - pe))
}
v_magnitude <- function(v, df_star = 2)
  as.character(cut(v, breaks = c(-Inf, c(0.10, 0.30, 0.50) / sqrt(df_star), Inf),
                   right = FALSE, labels = c("negligible", "small", "medium", "large")))
es_row <- function(tab, cmp) {
  ct <- suppressWarnings(chisq.test(tab, correct = FALSE)); kp <- kappa_stats(tab)
  data.frame(Comparison = cmp, n = sum(tab), X2 = as.numeric(ct$statistic),
             df = as.numeric(ct$parameter), p = ct$p.value,
             CramersV = cramers_v(tab), CohensW = cohens_w(tab),
             Agreement = kp[["po"]], Kappa = kp[["kappa"]], row.names = NULL)
}

# Full map: V and kappa are independent of sample size.
full_res <- bind_rows(lapply(targets, function(v)
  es_row(table(as_cls(cls_wide$GPP), as_cls(cls_wide[[v]])), paste0("GPP vs. ", v))))
full_res$V_magnitude <- v_magnitude(full_res$CramersV)
write.csv(full_res, file.path(OUT_DIR, "TableS1-effectsize_fullmap.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

# Subsampling: chi-squared and its p-value scale with n, so report them on
# repeated draws of fixed size and give percentile intervals for effect sizes.
N_DRAW <- 100; N_SUB <- 5000
boot <- bind_rows(lapply(seq_len(N_DRAW), function(it) {
  s <- cls_wide[sample(.N, N_SUB)]
  bind_rows(lapply(targets, function(v)
    cbind(iter = it, es_row(table(as_cls(s$GPP), as_cls(s[[v]])), paste0("GPP vs. ", v)))))
}))
write.csv(boot, file.path(OUT_DIR, "TableS1-boot_raw.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

tableS1 <- boot %>%
  group_by(Comparison) %>%
  summarise(X2_mean = mean(X2), p_max = max(p), p_median = median(p),
            V_mean = mean(CramersV), V_lo = quantile(CramersV, .025),
            V_hi = quantile(CramersV, .975), W_mean = mean(CohensW),
            K_mean = mean(Kappa), K_lo = quantile(Kappa, .025),
            K_hi = quantile(Kappa, .975), A_mean = mean(Agreement), .groups = "drop") %>%
  mutate(mag = v_magnitude(V_mean)) %>%
  transmute(Comparison,
            `Chi-squared (X2)`       = sprintf("%.2f", X2_mean),
            `p-value`                = ifelse(p_max < 0.05, "p < 0.05", sprintf("%.3f", p_median)),
            `Cramer's V (95% CI)`    = sprintf("%.3f (%.3f-%.3f)", V_mean, V_lo, V_hi),
            `Cohen's w`              = sprintf("%.3f", W_mean),
            `Cohen's kappa (95% CI)` = sprintf("%.3f (%.3f-%.3f)", K_mean, K_lo, K_hi),
            `Overall agreement (%)`  = sprintf("%.1f", A_mean * 100),
            Interpretation = sprintf(
              "Cramer's V = %.2f indicates a %s effect; the two classifications agree on %.1f%% of vegetated cells (kappa = %.2f).",
              V_mean, mag, A_mean * 100, K_mean))
write.csv(tableS1, file.path(OUT_DIR, "TableS1-revised.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")
print(as.data.frame(tableS1))

# Cell-wise standardised residuals: which class pairs drive the association.
# VSO / VGC are relabelled with the delta prefix for the exported table only.
cell_effect <- bind_rows(lapply(targets, function(v) {
  tab <- table(as_cls(cls_wide$GPP), as_cls(cls_wide[[v]]))
  ct  <- suppressWarnings(chisq.test(tab, correct = FALSE))
  obs <- as.data.frame(tab / sum(tab));                  names(obs) <- c("GPP_class", "Other_class", "obs_prop")
  exp <- as.data.frame(as.table(ct$expected / sum(tab))); names(exp) <- c("GPP_class", "Other_class", "exp_prop")
  res <- as.data.frame(as.table(ct$stdres));             names(res) <- c("GPP_class", "Other_class", "adj_std_resid")
  obs %>% left_join(exp, by = c("GPP_class", "Other_class")) %>%
    left_join(res, by = c("GPP_class", "Other_class")) %>%
    mutate(diff_pp = (obs_prop - exp_prop) * 100, Comparison = paste0("GPP vs. ", v))
}))
to_disp <- function(x) DISP[as.character(x)]
cell_effect$GPP_class   <- to_disp(cell_effect$GPP_class)
cell_effect$Other_class <- to_disp(cell_effect$Other_class)
write.csv(cell_effect, file.path(OUT_DIR, "TableS2-cellwise_effect.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

message("Done -> ", OUT_DIR)