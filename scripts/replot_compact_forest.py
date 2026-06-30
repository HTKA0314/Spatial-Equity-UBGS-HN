import os
import pandas as pd
import matplotlib.pyplot as plt

OUTPUT_DIR = r"d:\KLTN\outputs\step13b_stratified_regression"
PLOTS_DIR = r"d:\KLTN\outputs\plots"

coef_path = os.path.join(OUTPUT_DIR, "stratified_coefficients_all.csv")

if not os.path.exists(coef_path):
    print("Files not found.")
    exit(1)

df = pd.read_csv(coef_path)

# Filter to Model-1A Physical Observed
df_plot = df[df["model_folder"].str.startswith("model_1a_physical_observed")].copy()

# Filter by the preferred model for reporting (this guarantees we use the SLM/SEM direct coefficients as reported in Table 5)
df_plot = df_plot[df_plot["is_preferred_for_reporting"] == True]

# Extract stratum from folder name
def get_stratum(folder):
    if "global_interaction" in folder: return "Global_Interaction"
    if "high_density" in folder: return "High_Density"
    if "medium_density" in folder: return "Medium_Density"
    if "low_density" in folder: return "Low_Density"
    return "Unknown"

df_plot["stratum"] = df_plot["model_folder"].apply(get_stratum)

plot_vars = ["tree_frac", "grass_frac", "crop_frac", "water_frac", "bare_frac"]
df_plot = df_plot[df_plot["variable"].isin(plot_vars)]

variables = [v for v in plot_vars if v in df_plot["variable"].unique()]

fig, ax = plt.subplots(figsize=(8, 6.5))

palette = {
    "Global_Interaction": "#1a1a1a", # Dark grey
    "Low_Density": "#4575b4",      # Blue
    "Medium_Density": "#fdae61",   # Orange
    "High_Density": "#d73027",     # Red
}

var_labels = {
    "tree_frac": "Tree cover",
    "grass_frac": "Grass cover",
    "crop_frac": "Cropland",
    "water_frac": "Water",
    "bare_frac": "Bare land",
}

y_ticks = []
y_labels = []

categories = ["Global_Interaction", "High_Density", "Medium_Density", "Low_Density"]

for i, var in enumerate(reversed(variables)):
    df_var = df_plot[df_plot["variable"] == var].copy()
    
    df_var["stratum"] = pd.Categorical(df_var["stratum"], categories=categories, ordered=True)
    df_var = df_var.sort_values("stratum")
    
    y_center = i * 2.5
    y_ticks.append(y_center)
    y_labels.append(var_labels.get(var, var))
    
    for idx, row in enumerate(df_var.itertuples()):
        offset_map = {
            "Global_Interaction": 0.6,
            "High_Density": 0.2, 
            "Medium_Density": -0.2, 
            "Low_Density": -0.6
        }
        offset = offset_map.get(row.stratum, 0)
        y_pos = y_center + offset
        
        lbl = row.stratum.replace("_Interaction", "").replace("_", " ") if i == len(variables) - 1 else ""
        if lbl == "Global":
            lbl = "Global (All)"
            
        marker = "D" if row.stratum == "Global_Interaction" else "o"
        
        ax.errorbar(
            row.coef,
            y_pos,
            xerr=[[row.coef - row.ci_lower], [row.ci_upper - row.coef]] if pd.notna(row.ci_lower) else None,
            fmt=marker,
            color=palette.get(row.stratum, "black"),
            markersize=7 if marker=="D" else 8,
            capsize=4,
            linewidth=2,
            label=lbl
        )

ax.axvline(0, color="grey", linestyle="--", linewidth=1)
ax.set_yticks(y_ticks)
ax.set_yticklabels(y_labels, fontweight="bold")
ax.set_xlabel("SLM coefficient (°C)", fontweight="bold")
ax.grid(axis="x", linestyle="--", alpha=0.7)

handles, labels = ax.get_legend_handles_labels()
ax.legend(handles, labels, title="Model", loc="lower left", framealpha=0.9)

plt.tight_layout()
out_path = os.path.join(PLOTS_DIR, "fig8_compact_forest_plot.png")
plt.savefig(out_path, dpi=300, bbox_inches="tight")
print(f"Saved compact plot to {out_path}")
