# -*- coding: utf-8 -*-
"""
============================================================================
 05_1_Prepare_Threshold_SHAP_data.py
----------------------------------------------------------------------------
 STAGE 1 of 2  --  DATA PREPARATION ONLY (no model, no SHAP here)

 Purpose
 -------
 Assemble the per-grid threshold values produced by the R pipeline
 (03_1_Thershold_sensitivity.R) into a single wide table that can be fed
 directly into the SHAP analysis script (05_2_RandomForest_SHAP_Thresholds.py).

 What the R pipeline gives us
 ----------------------------
 In  Z:\\tangwenxi\\NCC-ResultData\\Result_2\\Threshold_Pipeline_Server\\GPP_breakpoint
 there is one CSV per climate factor:
     GPP_TMP_breakpoint_all.csv
     GPP_PRE_breakpoint_all.csv
     GPP_SR_breakpoint_all.csv
     GPP_VPD_breakpoint_all.csv
     GPP_CO2_breakpoint_all.csv
     GPP_SMroot_breakpoint_all.csv
     GPP_SMsurf_breakpoint_all.csv

 Each file has one row per grid cell with (among others) the columns:
     lab           : grid cell index (the 0.5 deg raster cell id)
     breakpoint_x  : the threshold on SD(climate)  <-- THIS is the value we want
     breakpoint_y, slope1, slope2, pvalue1, pvalue2, R2, ...

 We keep ONLY 'breakpoint_x' per grid and rename it to '<CLI>_Threshold',
 then join the 7 factors on 'lab' into one wide table.

 The Stable / VGC / VSO label for each grid (ResponseType) is NOT produced by
 the threshold pipeline, so it is read from a separate classification CSV that
 you provide (CLASS_CSV below: must contain columns 'lab' and 'ResponseType').

 Output
 ------
 SHAP_data_Threshold.csv  with columns:
     lab,
     TMP_Threshold, PRE_Threshold, SR_Threshold, VPD_Threshold,
     CO2_Threshold, SMroot_Threshold, SMsurf_Threshold,
     ResponseType
 This is the ONLY file the SHAP stage needs. The SHAP stage does not touch
 the R outputs again.
============================================================================
"""

import os
import pandas as pd

# ============================================================================
# USER CONFIG  -- edit these paths only
# ============================================================================

# Folder that holds the 7 GPP_<CLI>_breakpoint_all.csv files from the R pipeline.
BREAKPOINT_DIR = r"Z:\tangwenxi\NCC-ResultData\Result_2\Threshold_Pipeline_Server\GPP_breakpoint"

# Vegetation index prefix used in the R filenames (GPP_<CLI>_breakpoint_all.csv).
VI_NAME = "GPP"

# The 7 climate factors, in the order we want the columns to appear.
# (Same factor set as the original SHAP code.)
CLI_NAMES = ["TMP", "PRE", "SR", "VPD", "CO2", "SMroot", "SMsurf"]

# Column in each breakpoint CSV that holds the threshold value we keep.
THRESHOLD_COL = "breakpoint_x"

# Grid-cell id column shared by all files.
LAB_COL = "lab"

# Separate classification table: must contain 'lab' and 'ResponseType'
# (ResponseType values expected: 'Stable', 'VGC', 'VSO').
CLASS_CSV = "Z:/tangwenxi/NCC-ResultData/Result_1/Two_Period_Compare_2000_GPP/GPP_self_response_classification.csv"

# Where to write the assembled SHAP-ready table.
OUT_CSV = r"Z:\tangwenxi\NCC-ResultData\Result_3\SHAP_data_Threshold.csv"

# If a grid cell is missing the threshold for ANY of the 7 factors, drop it.
# (SHAP needs a complete feature row.) Set to False to keep partial rows (NaN).
DROP_INCOMPLETE = True


# ============================================================================
# STEP 1: read each breakpoint file, keep lab + breakpoint_x, rename column
# ============================================================================
def load_one_threshold(cli):
    """Read GPP_<cli>_breakpoint_all.csv and return a 2-column frame:
    lab + <cli>_Threshold."""
    fname = f"{VI_NAME}_{cli}_breakpoint_all.csv"
    fpath = os.path.join(BREAKPOINT_DIR, fname)
    if not os.path.exists(fpath):
        raise FileNotFoundError(f"Missing breakpoint file: {fpath}")

    df = pd.read_csv(fpath)

    # Sanity checks on required columns.
    for col in (LAB_COL, THRESHOLD_COL):
        if col not in df.columns:
            raise KeyError(
                f"Column '{col}' not found in {fname}. "
                f"Available columns: {list(df.columns)}"
            )

    out = df[[LAB_COL, THRESHOLD_COL]].copy()
    out = out.rename(columns={THRESHOLD_COL: f"{cli}_Threshold"})

    # A grid can appear only once per factor; guard against accidental dupes.
    out = out.drop_duplicates(subset=[LAB_COL], keep="first")
    print(f"  [{cli}] rows: {len(out)}  (non-NA thresholds: "
          f"{out[f'{cli}_Threshold'].notna().sum()})")
    return out


print("=== STEP 1: load per-factor threshold files ===")
merged = None
for cli in CLI_NAMES:
    part = load_one_threshold(cli)
    merged = part if merged is None else merged.merge(part, on=LAB_COL, how="outer")

print(f"Merged grid table: {merged.shape[0]} grids x {merged.shape[1]} cols")


# ============================================================================
# STEP 2: attach ResponseType (Stable / VGC / VSO) from the classification CSV
# ============================================================================
print("\n=== STEP 2: attach ResponseType labels ===")
if not os.path.exists(CLASS_CSV):
    raise FileNotFoundError(f"Classification CSV not found: {CLASS_CSV}")

cls = pd.read_csv(CLASS_CSV)
for col in (LAB_COL, "type_name"):
    if col not in cls.columns:
        raise KeyError(
            f"Column '{col}' not found in {os.path.basename(CLASS_CSV)}. "
            f"Available columns: {list(cls.columns)}"
        )
cls = cls[[LAB_COL, "type_name"]].drop_duplicates(subset=[LAB_COL], keep="first")

# Inner join: only grids that have BOTH thresholds and a class label survive.
df_out = merged.merge(cls, on=LAB_COL, how="inner")
print(f"After joining labels: {df_out.shape[0]} grids")
print("ResponseType distribution:")
print(df_out["type_name"].value_counts())


# ============================================================================
# STEP 3: order columns, optionally drop incomplete rows, save
# ============================================================================
print("\n=== STEP 3: finalize and save ===")
threshold_cols = [f"{c}_Threshold" for c in CLI_NAMES]
ordered_cols = [LAB_COL] + threshold_cols + ["type_name"]
df_out = df_out[ordered_cols]

if DROP_INCOMPLETE:
    before = len(df_out)
    df_out = df_out.dropna(subset=threshold_cols)
    print(f"Dropped {before - len(df_out)} grids with incomplete thresholds; "
          f"{len(df_out)} remain.")

os.makedirs(os.path.dirname(OUT_CSV), exist_ok=True)
# index=False -> SHAP stage controls its own indexing.
df_out.to_csv(OUT_CSV, index=False)
print(f"\nSaved SHAP-ready table -> {OUT_CSV}")
print(df_out.head())
