# -*- coding: utf-8 -*-
"""
============================================================================
 05_2_RandomForest_SHAP_Thresholds.py
----------------------------------------------------------------------------
 STAGE 2 of 2  --  SHAP ANALYSIS ONLY (runs independently of data prep)

 Prerequisite
 ------------
 Run 05_1_Prepare_Threshold_SHAP_data.py first. It writes:
     SHAP_data_Threshold.csv
 with columns:
     lab,
     TMP_Threshold, PRE_Threshold, SR_Threshold, VPD_Threshold,
     CO2_Threshold, SMroot_Threshold, SMsurf_Threshold,
     ResponseType            (Stable / VGC / VSO)

 This script:
   1) loads ONLY that table (it never re-reads the raw R outputs),
   2) trains a RandomForestClassifier on the 7 Threshold features,
   3) computes SHAP values (TreeExplainer),
   4) saves the per-class SHAP CSVs and the inverse-standardized X_test
      that the R plotting script (Fig4) consumes,
   5) draws the summary / swarm (group-peak) plots and the dependence plots.

 Because stages are separated, you can re-run the plotting block at the
 bottom on a previously trained model without re-doing the data assembly.

 NOTE on "interaction" dependence plots
 --------------------------------------
 Reviewer comment (3) asked us to inspect interaction effects via dependence
 plots. SHAP's dependence_plot can color each point by a second feature
 (interaction_index). Below we produce BOTH:
     - main-effect dependence plots (interaction_index=None), and
     - interaction dependence plots (interaction_index='auto', i.e. SHAP picks
       the strongest interacting feature; plus a full NxN grid).
============================================================================
"""

import os
from datetime import datetime

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import seaborn as sns

from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import (
    confusion_matrix, accuracy_score, f1_score,
    precision_score, recall_score, roc_auc_score,
)
import shap


# ============================================================================
# USER CONFIG
# ============================================================================
name = "Threshold-1"

# The SHAP-ready table written by stage 1.
INPUT_CSV = r"Z:\tangwenxi\NCC-ResultData\Result_3\SHAP_data_Threshold_LULCcut.csv"

# Parent output directory. A dated sub-folder is created underneath.
SHAP_dir =  r"Z:\tangwenxi\NCC-ResultData\Result_3"
print(f"Output directory: {SHAP_dir}")

# The 7 threshold features, in a fixed order (controls feature_names order).
FEATURE_COLS = [
    "TMP_Threshold", "PRE_Threshold", "SR_Threshold", "VPD_Threshold",
    "CO2_Threshold", "SMroot_Threshold", "SMsurf_Threshold",
]
RESPONSE_COL = "ResponseType"
LAB_COL = "lab"


# ============================================================================
# ============================  MODEL TRAINING  ==============================
# ============================================================================
# Folders
output_folder = os.path.join(SHAP_dir, name)
os.makedirs(output_folder, exist_ok=True)
data_output_folder = os.path.join(SHAP_dir, name, "DATA")
os.makedirs(data_output_folder, exist_ok=True)

# --- Load the assembled table (the ONLY input) ---
df = pd.read_csv(INPUT_CSV)
df = df.dropna(subset=FEATURE_COLS + [RESPONSE_COL])
print(df.info())
print(df[RESPONSE_COL].value_counts())

X = df[FEATURE_COLS].copy()
y = df[RESPONSE_COL].copy()
labs = df[LAB_COL].copy()

# --- Train / test split (stratified, identical settings to original) ---
X_train, X_test, y_train, y_test, labs_train, labs_test = train_test_split(
    X, y, labs, test_size=0.3, random_state=42, stratify=y
)

# Save raw splits
X_train.to_csv(os.path.join(data_output_folder, f"X_train_{name}.csv"), index=True)
X_test.to_csv(os.path.join(data_output_folder, f"X_test_{name}.csv"), index=True)
y_train.to_csv(os.path.join(data_output_folder, f"y_train_{name}.csv"), index=True)
y_test.to_csv(os.path.join(data_output_folder, f"y_test_{name}.csv"), index=True)
labs_train.to_csv(os.path.join(data_output_folder, f"labs_train_{name}.csv"), index=True)
labs_test.to_csv(os.path.join(data_output_folder, f"labs_test_{name}.csv"), index=True)

# --- Standardize ---
scaler = StandardScaler()
X_train_scaled = pd.DataFrame(scaler.fit_transform(X_train), columns=X_train.columns)
X_test_scaled = pd.DataFrame(scaler.transform(X_test), columns=X_test.columns)
X_train_scaled.to_csv(os.path.join(data_output_folder, f"X_train_scaled_{name}.csv"), index=True)
X_test_scaled.to_csv(os.path.join(data_output_folder, f"X_test_scaled_{name}.csv"), index=True)

# --- Random Forest ---
model = RandomForestClassifier(n_estimators=100, random_state=42)
model.fit(X_train_scaled, y_train)

# Quick feature-importance bar
importances = pd.Series(model.feature_importances_, index=X.columns)
importances.sort_values(ascending=False).plot(kind="bar", title="Feature Importances")
plt.tight_layout()
plt.savefig(os.path.join(output_folder, f"feature_importances_{name}.png"))
plt.close()

# --- Predict + metrics ---
y_pred = model.predict(X_test_scaled)
pd.DataFrame({"lab": labs_test, "True Values": y_test, "Predictions": y_pred}).to_csv(
    os.path.join(output_folder, f"predictions_{name}.csv"), index=False
)

eval_df = pd.DataFrame([{
    "Accuracy": accuracy_score(y_test, y_pred),
    "F1 Score (weighted)": f1_score(y_test, y_pred, average="weighted"),
    "Precision (weighted)": precision_score(y_test, y_pred, average="weighted"),
    "Recall (weighted)": recall_score(y_test, y_pred, average="weighted"),
    "AUC": roc_auc_score(y_test, model.predict_proba(X_test_scaled),
                         multi_class="ovr", average="weighted"),
}])
eval_df.to_csv(os.path.join(output_folder, f"model_evaluation_metrics_{name}.csv"), index=False)
print(eval_df)

# Confusion matrix
cm = confusion_matrix(y_test, y_pred, labels=model.classes_)
plt.figure(figsize=(8, 6))
sns.heatmap(cm, annot=True, fmt="d", cmap="Blues",
            xticklabels=model.classes_, yticklabels=model.classes_)
plt.xlabel("Predicted"); plt.ylabel("True"); plt.title("Confusion Matrix")
plt.tight_layout()
plt.savefig(os.path.join(output_folder, f"confusion_matrix_{name}.png"))
plt.close()


# ============================================================================
# ===========================  SHAP COMPUTATION  =============================
# ============================================================================
explainer = shap.TreeExplainer(model)
shap_values = explainer.shap_values(X_test_scaled)   # shape: (n, n_features, n_classes)

# Save per-class SHAP values (these feed the R Fig4 script).
for class_idx, class_name in enumerate(model.classes_):
    shap_values_class = shap_values[:, :, class_idx]
    sdf = pd.DataFrame(shap_values_class, columns=X.columns)
    sdf.index = X_test.index
    sdf["lab"] = labs_test
    sdf.to_csv(os.path.join(output_folder, f"shap_values_{class_name}_{name}.csv"), index=True)

# Inverse-standardized X_test (R script merges SHAP with the ORIGINAL x-values).
X_test_original = pd.DataFrame(scaler.inverse_transform(X_test_scaled), columns=X.columns)
X_test_original.index = X_test.index
X_test_original["lab"] = labs_test
X_test_original.to_csv(os.path.join(data_output_folder, f"X_test_original_{name}.csv"), index=True)

print("SHAP values shape:", shap_values.shape)


# ============================================================================
# =====================  PLOTTING (re-runnable block)  =======================
# ----------------------------------------------------------------------------
# Everything below only needs `model`, `shap_values`, `X_test`, `X_test_scaled`.
# If you saved/loaded those, you can re-run from here without retraining.
# ============================================================================

# Pretty axis labels for the 7 threshold features, in FEATURE_COLS order.
new_columns = [
    r"$T_{\Delta \mathrm{TMP}}$",
    r"$T_{\Delta \mathrm{PRE}}$",
    r"$T_{\Delta \mathrm{SR}}$",
    r"$T_{\Delta \mathrm{VPD}}$",
    r"$T_{\Delta \mathrm{CO_2}}$",
    r"$T_{\Delta \mathrm{SMroot}}$",
    r"$T_{\Delta \mathrm{SMsurf}}$",
]

custom_cmap = mcolors.LinearSegmentedColormap.from_list(
    "my_cmap", ["#4DAF4A", "#ffffff", "#d62728"]      # green -> white -> red
)


# ----------------------------------------------------------------------------
# 1. SUMMARY / SWARM PLOTS  (one "group-peak" beeswarm per class)
# ----------------------------------------------------------------------------
labels = ["(c) ", "(b) ", "(a) "]
for class_idx in range(len(model.classes_)):
    feature_names = [
        f"{a}: {b}"
        for a, b in zip(new_columns,
                        np.abs(shap_values[:, :, class_idx]).mean(0).round(3))
    ]
    class_name = model.classes_[class_idx]
    if class_name == "VSO":
        title_name = labels[class_idx] + "\u03b4VSO"
    elif class_name == "VGC":
        title_name = labels[class_idx] + "\u03b4VGC"
    else:
        title_name = labels[class_idx] + class_name

    print(f"Summary (beeswarm) plot for class {class_name} (idx {class_idx})")
    shap.summary_plot(
        shap_values[:, :, class_idx], X_test,
        feature_names=feature_names,
        max_display=X_test.shape[1],
        plot_type="dot", cmap=custom_cmap,
        plot_size=[5, 5], show=False,
    )
    plt.tight_layout()
    plt.text(-0.15, 1.05, title_name, transform=plt.gca().transAxes,
             fontsize=16, va="top", ha="left")
    plt.xlabel("SHAP value", fontsize=16)
    plt.savefig(os.path.join(output_folder, f"{class_idx}_{title_name}_shap_summary_plot.png"),
                dpi=300, bbox_inches="tight")
    plt.close()

# Bar plot across all classes
shap.summary_plot(shap_values, X_test, plot_type="bar", show=False,
                  class_names=model.classes_)
plt.savefig(os.path.join(output_folder, "shap_bar_plot.png"), dpi=300, bbox_inches="tight")
plt.close()


# ----------------------------------------------------------------------------
# 2. MAIN-EFFECT DEPENDENCE PLOTS  (interaction_index=None)
#    One PNG per feature x class, plus a combined NfeaturesxNclasses grid.
# ----------------------------------------------------------------------------
dependence_output_folder = os.path.join(SHAP_dir, name, "Dependence")
os.makedirs(dependence_output_folder, exist_ok=True)

for feature in X.columns:
    for class_idx, class_name in enumerate(model.classes_):
        shap.dependence_plot(
            feature, shap_values[:, :, class_idx], X_test,
            interaction_index=None, show=False,
        )
        plt.title(f"{feature} - {class_name}")
        plt.savefig(os.path.join(dependence_output_folder,
                                 f"{feature}_{class_name}_shap_dependence.png"),
                    dpi=300, bbox_inches="tight")
        plt.close()

n_features = len(X.columns)
n_classes = len(model.classes_)
fig, axes = plt.subplots(n_features, n_classes, figsize=(n_classes * 4, n_features * 3))
for i, feature in enumerate(X.columns):
    for class_idx, class_name in enumerate(model.classes_):
        shap.dependence_plot(
            feature, shap_values[:, :, class_idx], X_test,
            interaction_index=None, show=False, ax=axes[i, class_idx],
        )
        axes[i, class_idx].set_title(f"{feature} - {class_name}")
plt.tight_layout()
plt.savefig(os.path.join(dependence_output_folder, "all_shap_dependence.png"),
            dpi=300, bbox_inches="tight")
plt.close()


# ----------------------------------------------------------------------------
# 3. INTERACTION DEPENDENCE PLOTS  (reviewer comment 3)
#    3a) auto-interaction: SHAP picks the strongest interacting feature and
#        colors the points by it -> reveals 2-way interactions at a glance.
#    3b) full interaction grid per class: for each (feature, class), color by
#        every other feature in turn (N x N panel) so every pairwise
#        interaction can be inspected.
# ----------------------------------------------------------------------------
interaction_output_folder = os.path.join(SHAP_dir, name, "Dependence_interaction")
os.makedirs(interaction_output_folder, exist_ok=True)

# 3a) Auto-selected interacting feature (one PNG per feature x class).
for feature in X.columns:
    for class_idx, class_name in enumerate(model.classes_):
        shap.dependence_plot(
            feature, shap_values[:, :, class_idx], X_test,
            interaction_index="auto",          # SHAP chooses strongest interaction
            show=False,
        )
        plt.title(f"{feature} - {class_name} (auto interaction)")
        plt.savefig(os.path.join(interaction_output_folder,
                                 f"{feature}_{class_name}_auto_interaction.png"),
                    dpi=300, bbox_inches="tight")
        plt.close()

# 3b) Full pairwise interaction grid, one figure PER class.
#     Rows = main feature, Columns = the feature used for coloring.
for class_idx, class_name in enumerate(model.classes_):
    sv_class = shap_values[:, :, class_idx]
    fig, axes = plt.subplots(n_features, n_features,
                             figsize=(n_features * 3.2, n_features * 3.0))
    for i, feat_main in enumerate(X.columns):
        for k, feat_color in enumerate(X.columns):
            ax = axes[i, k]
            shap.dependence_plot(
                feat_main, sv_class, X_test,
                interaction_index=feat_color,   # explicit pairwise interaction
                show=False, ax=ax,
            )
            ax.set_title(f"{feat_main} x {feat_color}", fontsize=8)
            ax.tick_params(labelsize=7)
            ax.set_xlabel(feat_main, fontsize=8)
            if k == 0:
                ax.set_ylabel(f"SHAP ({feat_main})", fontsize=8)
            else:
                ax.set_ylabel("")
    fig.suptitle(f"Pairwise SHAP interaction grid - class {class_name}", fontsize=14)
    plt.tight_layout(rect=[0, 0, 1, 0.98])
    plt.savefig(os.path.join(interaction_output_folder,
                             f"interaction_grid_{class_name}.png"),
                dpi=300, bbox_inches="tight")
    plt.close()


# ----------------------------------------------------------------------------
# 3c) SHAP interaction VALUES heatmap (mean |interaction|) per class.
#     This quantifies which feature PAIRS interact most strongly. It uses the
#     full interaction tensor (can be slow / memory heavy on large test sets;
#     subsample X_test if needed).
# ----------------------------------------------------------------------------
try:
    # Optional subsample to keep the interaction tensor affordable.
    max_rows_for_interaction = 1000
    if X_test_scaled.shape[0] > max_rows_for_interaction:
        idx_sub = np.random.RandomState(42).choice(
            X_test_scaled.shape[0], max_rows_for_interaction, replace=False)
        X_inter = X_test_scaled.iloc[idx_sub]
    else:
        X_inter = X_test_scaled

    # shape: (n, n_features, n_features, n_classes)
    inter_values = explainer.shap_interaction_values(X_inter)

    for class_idx, class_name in enumerate(model.classes_):
        iv = inter_values[..., class_idx] if inter_values.ndim == 4 else inter_values
        mean_abs = np.abs(iv).mean(axis=0)        # n_features x n_features
        plt.figure(figsize=(8, 6.5))
        sns.heatmap(mean_abs, annot=True, fmt=".3f", cmap="PRGn",
                    xticklabels=X.columns, yticklabels=X.columns)
        plt.title(f"Mean |SHAP interaction| - class {class_name}")
        plt.tight_layout()
        plt.savefig(os.path.join(interaction_output_folder,
                                 f"mean_abs_interaction_{class_name}.png"),
                    dpi=300, bbox_inches="tight")
        plt.close()

        # Also export the matrix so the R script / paper table can use it.
        pd.DataFrame(mean_abs, index=X.columns, columns=X.columns).to_csv(
            os.path.join(interaction_output_folder,
                         f"mean_abs_interaction_{class_name}.csv"))
except Exception as e:
    # Interaction values are expensive; never let them break the main run.
    print(f"[warn] SHAP interaction values skipped: {e}")


print("All SHAP plots written under:", os.path.join(SHAP_dir, name))
