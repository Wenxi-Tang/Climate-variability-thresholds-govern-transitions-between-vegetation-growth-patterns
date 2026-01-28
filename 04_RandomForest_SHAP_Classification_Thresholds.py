# -*- coding: utf-8 -*-
"""
Created on Thu Jan 29 04:49:00 2026

@author: ADMIN
"""


import pandas as pd
import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import mean_squared_error, r2_score, confusion_matrix
from sklearn.metrics import accuracy_score, f1_score, precision_score, recall_score, classification_report, roc_auc_score
import shap
import matplotlib.pyplot as plt
import matplotlib
import seaborn as sns
import os
import matplotlib.colors as mcolors
from datetime import datetime


# === Set folder name ===
name = 'Threshold'

# === Create directories ===
# Set parent directory (replace with your actual path)
file_dir = "E:/GlobalVeg/Out-Table-and-Figure-VPD et al/"
file_dir2 = file_dir+"5.SHAP-VPD+CO2+SM-CO2TRUE/"
# Get today's date
today = datetime.today().strftime('%Y-%m-%d')
# Construct full path
SHAP_dir = os.path.join(file_dir2, today)
# Create directory (if it doesn't exist)
os.makedirs(SHAP_dir, exist_ok=True)
print(f"✅ Directory created: {SHAP_dir}")


# === Adjust variables ===
columns_to_select = ['TMP', 'PRE', 'SR', 'VPD', 'CO2', 'SMroot', 'SMsurf']
columns_to_select = ['GPP', 'SIF', 'LAI', 'NDVI']
columns_to_select = ['TMP', 'PRE', 'SR', 'VPD', 'CO2', 'SMroot', 'SMsurf', 'GPP', 'SIF', 'LAI', 'NDVI']
columns_to_select = ['TMP_Threshold', 'PRE_Threshold', 'SR_Threshold', 'VPD_Threshold', 
                     'CO2_Threshold', 'SMroot_Threshold', 'SMsurf_Threshold', 
                     'GPP_sd_1982_2020', 'SIF_sd_1982_2020', 'LAI_sd_1982_2020', 'NDVI_sd_1982_2020']


year_to_selected = ['Threshold']
year_to_selected = ['Lagged']
year_to_selected = ['1982_2020']
year_to_selected = ['1982_2020', 'Threshold']


# === Adjust variables ===
columns_to_select = ['TMP', 'PRE', 'SR', 'VPD', 'CO2', 'SMroot', 'SMsurf']
year_to_selected = ['Threshold']

# === Step 1: Load Data ===
# Create folder if it does not exist
output_folder = os.path.join(SHAP_dir, name)  # Create target folder under SHAP_dir
if not os.path.exists(output_folder):
    os.makedirs(output_folder)


# Assuming data format is as follows:
# | temp | precip | co2 | ... | response_type |
df = pd.read_csv(file_dir2+"/SHAP-Data/SHAP_data.csv")  # Replace with actual file path
df = df.dropna()
df = df.drop(['Unnamed: 0'], axis=1)  # Drop a column, axis=1 means column

# Filter data
df_selected = df[df.columns[df.columns.str.contains('|'.join(columns_to_select))]]
df_selected = df_selected[df_selected.columns[df_selected.columns.str.contains('|'.join(year_to_selected))]]
response_type = 'ResponseType'


# Check for missing values and basic info
print(df_selected.info())
print(df[response_type].value_counts())  # Check response type distribution


# Historical Training ---------------------------------------------------------------------------------------------
# === Step 2: Split Feature and Target Columns ===
X = df_selected
y = df[response_type]
labs = df['lab']  # Save 'lab' column

# === Step 3: Split Train and Test Sets ===
X_train, X_test, y_train, y_test, labs_train, labs_test = train_test_split(
    X, y, labs, test_size=0.3, random_state=42, stratify=y
)

# Set output folder path
data_output_folder = os.path.join(SHAP_dir, name, "DATA")  # Create target folder under SHAP_dir
if not os.path.exists(data_output_folder):
    os.makedirs(data_output_folder)

# Save original split data
X_train.to_csv(data_output_folder + "/X_train_" + name + ".csv", index=True)
X_test.to_csv(data_output_folder + "/X_test_" + name + ".csv", index=True)

y_train.to_csv(data_output_folder + "/y_train_" + name + ".csv", index=True)
y_test.to_csv(data_output_folder + "/y_test_" + name + ".csv", index=True)

labs_train.to_csv(data_output_folder + "/labs_train_" + name + ".csv", index=True)
labs_test.to_csv(data_output_folder + "/labs_test_" + name + ".csv", index=True)

# === Step 4: Data Standardization ===
scaler = StandardScaler()
X_train_scaled = pd.DataFrame(scaler.fit_transform(X_train), columns=X_train.columns)
X_test_scaled = pd.DataFrame(scaler.transform(X_test), columns=X_test.columns)

# === Step 5: Save Standardized Data ===
X_train_scaled.to_csv(data_output_folder + "/X_train_scaled_" + name + ".csv", index=True)
X_test_scaled.to_csv(data_output_folder + "/X_test_scaled_" + name + ".csv", index=True)


# === Step 6: Random Forest Training ===
model = RandomForestClassifier(n_estimators=100, random_state=42)
model.fit(X_train_scaled, y_train)

# Briefly check feature importance
importances = pd.Series(model.feature_importances_, index=X.columns)
importances.sort_values(ascending=False).plot(kind='bar', title='Feature Importances')
plt.tight_layout()
plt.savefig(output_folder + "/feature_importances_"+name+".png")
plt.close()


# === Step 7: Model Prediction ===
y_pred = model.predict(X_test_scaled)

predictions_df = pd.DataFrame({'lab': labs_test, 'True Values': y_test, 'Predictions': y_pred})
predictions_df.to_csv(output_folder + "/predictions_"+name+".csv", index=False)


# === Step 8: Save Evaluation Results to DataFrame ===
# Calculate evaluation metrics
accuracy = accuracy_score(y_test, y_pred)
f1 = f1_score(y_test, y_pred, average='weighted')  # Add average='weighted' for multi-class classification
precision = precision_score(y_test, y_pred, average='weighted')  # Add average='weighted' for multi-class classification
recall = recall_score(y_test, y_pred, average='weighted')  # Add average='weighted' for multi-class classification
# For AUC, calculate micro and macro AUC (suitable for multi-class problems)
auc = roc_auc_score(y_test, model.predict_proba(X_test_scaled), multi_class='ovr', average='weighted')
 
evaluation_metrics = {
    'Accuracy': accuracy,
    'F1 Score (weighted)': f1,
    'Precision (weighted)': precision,
    'Recall (weighted)': recall,
    'AUC': auc
}

# Convert evaluation results to DataFrame
eval_df = pd.DataFrame([evaluation_metrics])

# Save evaluation results to CSV file
eval_df.to_csv(output_folder + "/model_evaluation_metrics_"+name+".csv", index=False)

# Print evaluation results
print(eval_df)


# Confusion Matrix
cm = confusion_matrix(y_test, y_pred, labels=model.classes_)
plt.figure(figsize=(8, 6))
sns.heatmap(cm, annot=True, fmt='d', cmap='Blues', xticklabels=model.classes_, yticklabels=model.classes_)
plt.xlabel('Predicted')
plt.ylabel('True')
plt.title('Confusion Matrix')
plt.tight_layout()
plt.savefig(output_folder + "/confusion_matrix_"+name+".png")
plt.close()


# SHAP Interpretability Analysis -----------------------------------------------------------------------------------
# === Step 1: Initialize SHAP Explainer ===
explainer = shap.TreeExplainer(model)
shap_values = explainer.shap_values(X_test_scaled)

# Save SHAP results
for class_idx, class_name in enumerate(model.classes_):
    # Extract SHAP values for each class
    shap_values_class = shap_values[:, :, class_idx]
    
    # Convert SHAP values to DataFrame
    shap_values_class_df = pd.DataFrame(shap_values_class, columns=X.columns)
    
    # Set row names to X_test index
    shap_values_class_df.index = X_test.index
    
    # Add 'lab' column
    shap_values_class_df['lab'] = labs_test
    
    # Save to CSV file
    shap_values_class_df.to_csv(
        os.path.join(output_folder, f"shap_values_{class_name}_{name}.csv"),
        index=True
    )


# === Step 2: Plot all mainstream SHAP analysis charts ===
# Inverse transform feature data 
X_test_original = pd.DataFrame(scaler.inverse_transform(X_test_scaled), columns=X.columns)

print(shap_values.shape)



# 1. Summary Plot
# new_columns = ['CO$_2$', 'PRE', 'SM$_{root}$', 'SM$_{surf}$', 'SR', 'TMP', 'VPD']
new_columns = [
    r'$T_{\Delta \mathrm{CO_2}}$',
    r'$T_{\Delta \mathrm{PRE}}$',
    r'$T_{\Delta \mathrm{SMroot}}$',
    r'$T_{\Delta \mathrm{SMsurf}}$',
    r'$T_{\Delta \mathrm{SR}}$',
    r'$T_{\Delta \mathrm{TMP}}$',
    r'$T_{\Delta \mathrm{VPD}}$'
]


# Custom color scheme
custom_cmap = mcolors.LinearSegmentedColormap.from_list(
    "my_cmap",
    ["#4DAF4A", "#ffffff", "#d62728"]   # Blue to White to Red
)

labels = ["(c) ", "(b) ", "(a) "]
for class_idx in [0, 1, 2]:
    feature_names = [
        a + ": " + str(b) for a,b in zip(new_columns, np.abs(shap_values[:, :, class_idx]).mean(0).round(3))
    ]

    class_name = model.classes_[class_idx]
    if class_name == 'VSO':
        title_name = labels[class_idx] + 'δVSO'
    elif class_name == 'VGC':
        title_name = labels[class_idx] + 'δVGC'
    else:
        title_name = labels[class_idx] + class_name  # Keep original by default
    
    print(f"Generating SHAP plot for Class {class_name} (Index {class_idx})")

    # Plot SHAP summary plot, reduce dot size
    shap.summary_plot(shap_values[:, :, class_idx], 
                      X_test, 
                      feature_names=feature_names, 
                      max_display=X_test.shape[1], 
                      plot_type="dot", 
                      cmap=custom_cmap, 
                      plot_size=[5,5],
                      show=False)

    plt.tight_layout()

    # Add top-left subtitle (consistent with loop order)
    plt.text(-0.15, 1.05, title_name, transform=plt.gca().transAxes,
             fontsize=16, va='top', ha='left')
    plt.xlabel("SHAP value", fontsize=16)

    # Ensure save path exists
    os.makedirs(output_folder, exist_ok=True)

    # Save
    plt.savefig(os.path.join(output_folder, f"{class_idx}_{title_name}_shap_summary_plot.png"),
                dpi=300, bbox_inches="tight")
    plt.close()



# 2. Bar Plot
print(model.classes_)
shap.summary_plot(shap_values, X_test, plot_type="bar", show=False, class_names=model.classes_)
plt.savefig(output_folder +"/"+ "shap_bar_plot.png", dpi=300, bbox_inches="tight")
plt.close()

shap.summary_plot(shap_values[:, :, 2],  X_test, plot_type="bar", show=True)


# 3. Heat Plot
# Create storage path
heatplot_output_folder = output_folder+"/shap_heatmaps"
if not os.path.exists(heatplot_output_folder):
    os.makedirs(heatplot_output_folder)

# Iterate through SHAP values for each class (0 to 2)
shap_values_sampled = explainer(X_test_scaled)
shap_values_sampled.values = shap_values_sampled.values.round(3)

for class_idx in range(3):  # Corresponds to classes 0, 1, 2
    # Get SHAP values for current class
    shap_values_class = shap_values_sampled[:, :, class_idx]
    
    # Create SHAP Explanation object
    expl = shap.Explanation(values=shap_values_class, 
                             data=X_test_scaled, 
                             feature_names=X_test_scaled.columns)

    # Generate heatmap and add title
    class_name = model.classes_[class_idx]
    shap.plots.heatmap(expl, 
                       max_display=12,  # Added max_display or it might error if not enough features
                       plot_width=16, cmap='PRGn', show=False) # Note: Original code had syntax error here, fixed arguments
    plt.title(f"SHAP Heatmap for {class_name}")
    
    # Save heatmap
    heatmap_filename = os.path.join(heatplot_output_folder, f"shap_heatmap_{class_name}.png")
    plt.savefig(heatmap_filename, dpi=300, bbox_inches="tight")
    plt.close()  # Close current figure to avoid overlap




# === Step 3: Plot dependence plots for three ResponseTypes of the same climate factor in one figure ===
dependence_output_folder = os.path.join(SHAP_dir, name, 'Dependence') 
if not os.path.exists(dependence_output_folder):
    os.makedirs(dependence_output_folder)

plt.figure(figsize=(15, 15))


# Single plot
for i, feature in enumerate(X.columns):
    for class_idx, class_name in enumerate(model.classes_):
        shap_values_class = shap_values[:, :, class_idx] 
        shap.dependence_plot(
            feature,
            shap_values_class,
            X_test,
            interaction_index=None,
            show=False
        )
        plt.title(f"{feature} - {class_name}")
        plt.savefig(dependence_output_folder +"/" + f"{feature}_{class_name}"+ "_shap_dependence"+".png", 
                    dpi=300, bbox_inches="tight")
        plt.close()
        

# Combined plot
# Set plot dimensions
n_features = len(X.columns)  # Number of features
n_classes = len(model.classes_)  # Number of classes
fig, axes = plt.subplots(n_features, n_classes, figsize=(n_classes * 4, n_features * 3))

# Iterate through features and classes to generate each dependence plot
for i, feature in enumerate(X.columns):
    for class_idx, class_name in enumerate(model.classes_):
        # Get SHAP values for each class
        shap_values_class = shap_values[:, :, class_idx]

        # Create SHAP dependence plot and embed it into the corresponding subplot
        shap.dependence_plot(
            feature,                    # Feature name to plot
            shap_values_class,          # SHAP values
            X_test,                     # Feature data
            show=False,                 # Do not show plot directly
            interaction_index=None,
            ax=axes[i, class_idx]       # Plot on corresponding subplot
        )
        
        # Add title to each subplot
        axes[i, class_idx].set_title(f"{feature} - {class_name}")

# Adjust subplot layout to prevent overlap
plt.tight_layout()

# Save final combined image
plt.savefig(dependence_output_folder + "/all_shap_dependence.png", dpi=300, bbox_inches="tight")
plt.close()


# Save inverse standardized results
X_test_original.index = X_test.index
X_test_original['lab'] = labs_test
X_test_original.to_csv(data_output_folder + "/X_test_original_"+name+".csv", index=True)






import seaborn as sns
import statsmodels.api as sm

# Ensure output folder exists
os.makedirs(dependence_output_folder, exist_ok=True)

# Set up canvas
n_features = len(X.columns)  # Number of features
n_classes = len(model.classes_)  # Number of classes
fig, axs = plt.subplots(n_features, n_classes, figsize=(n_classes * 4, n_features * 3))

lowess = sm.nonparametric.lowess

# Iterate through features and classes
for i, feature in enumerate(X.columns):
    for j, class_idx in enumerate([2, 0, 1]):  # Order 2, 0, 1
        class_name = model.classes_[class_idx]

        ax = axs[i, j]

        # Extract X and corresponding SHAP values
        x = X_test[feature].values
        shap_value = shap_values[:, :, class_idx]
        y = shap_value[:, i]

        # Plot scatter
        ax.scatter(x=x, y=y, s=0.1, color='slategrey')

        # Plot density
        sns.kdeplot(x=x, y=y, cmap="Blues", ax=ax, fill=True, levels=6, alpha=1)

        # Plot lowess fit line
        x_low, x_high = np.percentile(x, [0, 100])
        mask = (x >= x_low) & (x <= x_high)
        x_filtered = x[mask]
        y_filtered = y[mask]
        z = lowess(y_filtered, x_filtered, frac=0.3)
        ax.plot(z[:, 0], z[:, 1], color='white', linestyle='--', linewidth=3)

        # Y-axis range is not fixed here! No ax.set_ylim()

        # Axis settings
        ax.set_xlabel(feature, fontsize=14)
        if j == 0:
            ax.set_ylabel("SHAP value", fontsize=14)
        else:
            ax.set_ylabel("")

        # Small annotations
        ax.text(0.05, max(y)*0.9, 'Positive Contribution (PC)', ha='left', va='bottom', fontsize=10)
        ax.text(0.05, min(y)*1.1, 'Negative Contribution (NC)', ha='left', va='top', fontsize=10)

        # Subplot title
        ax.set_title(f"{feature} - {class_name}", fontsize=14)

plt.tight_layout()

# Save
plt.savefig(os.path.join(dependence_output_folder, "all_shap_dependence_auto_ylim.png"), dpi=300, bbox_inches="tight")
plt.close()


# === Step 4: Extract "Threshold" for each feature ===
threshold_dict = {}

# Iterate through each class
for class_idx, class_name in enumerate(model.classes_):  # Iterate through classes
    for feature in X.columns:
        # Get corresponding SHAP values and original values for each feature
        values = X_test[feature].values
        shap_vals = shap_values[:, X.columns.get_loc(feature), class_idx]  # Get SHAP values corresponding to the class

        # Estimate inflection point using segmented average and slope change
        sorted_idx = np.argsort(values)
        sorted_x = values[sorted_idx]
        sorted_shap = shap_vals[sorted_idx]

        # Simple method: Point of maximum SHAP rate of change as "threshold"
        diff = np.gradient(sorted_shap, sorted_x)
        max_idx = np.argmax(np.abs(diff))

        # Check index range to ensure no out-of-bounds
        if max_idx < len(sorted_x):
            threshold_value = sorted_x[max_idx]
        else:
            threshold_value = sorted_x[-1]  # Prevent index out of bounds

        # Save threshold for each class to dictionary
        threshold_dict[f"{feature}_{class_name}"] = threshold_value

# Print thresholds for all variables
print("Extracted response thresholds:")
for k, v in threshold_dict.items():
    print(f"{k}: {v:.3f}")



# === Step 5: Save Threshold Data ===
# Extract climate indicators, response types, and thresholds
climate = []
response_type = []
threshold = []

for key, value in threshold_dict.items():
    parts = key.split('_')
    # If key has multiple parts, take first part as climate_type
    climate_type = parts[0]
    response = '_'.join(parts[2:])  # Join remaining parts as response
    
    climate.append(climate_type)
    response_type.append(response)
    threshold.append(value)

# Create DataFrame
df_threshold = pd.DataFrame({
    'Climate': climate,
    'ResponseType': response_type,
    'Threshold': threshold
})

# Display DataFrame
print(df_threshold)
df_threshold.to_csv(output_folder + "/climate_response_thresholds.csv")