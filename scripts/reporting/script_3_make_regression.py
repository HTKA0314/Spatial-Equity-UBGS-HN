import os
import sys
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))
from scripts.reporting import constants as const
from scripts.reporting.script_2a_make_main_plots import set_style, save_fig

def _sig_color(p):
    if pd.isna(p): return "#AAAAAA"
    if p < 0.01: return "#C0392B"
    if p < 0.05: return "#E67E22"
    return "#AAAAAA"

def supp3_density_dependent_lst_associations():
    # Attempt to load data
    path = os.path.join(const.OUTPUT_DIR, "step13b_density_interactions", "marginal_effects_by_density.csv")
    sig_path = os.path.join(const.OUTPUT_DIR, "step13b_density_interactions", "interaction_significance.csv")
    
    if not os.path.exists(path):
        print("Skipping Supp3: marginal_effects_by_density.csv not found.")
        return
        
    df = pd.read_csv(path)
    df = df[(df["model_label"] == "Model-1A-Physical-Observed") & (df["model_name"] == "SLM")].copy()
    if df.empty: return
    
    sig = pd.read_csv(sig_path) if os.path.exists(sig_path) else pd.DataFrame()
    var_order = ["tree_frac", "grass_frac", "crop_frac", "water_frac"]
    var_labels = {"tree_frac": "Tree cover", "grass_frac": "Grass cover", "crop_frac": "Cropland", "water_frac": "Water cover"}
    strata = ["Low_Density", "Medium_Density", "High_Density"]
    
    fig, axes = plt.subplots(len(var_order), 1, figsize=(8.4, 2.35 * len(var_order)), sharex=True)
    if len(var_order) == 1: axes = [axes]
    
    for ax, var in zip(axes, var_order):
        d = df[df["variable"] == var].copy()
        d["stratum"] = pd.Categorical(d["stratum"], categories=strata, ordered=True)
        d = d.sort_values("stratum")
        labels = []
        for i, row in enumerate(d.itertuples()):
            label = const.DENSITY_CODE_TO_LABEL.get(row.stratum, row.stratum)
            if row.stratum != "Low_Density" and not sig.empty:
                match = sig[(sig["model_label"] == "Model-1A-Physical-Observed") & (sig["variable"] == var)]
                if not match.empty:
                    p_col = "Medium_p_value" if row.stratum == "Medium_Density" else "High_p_value"
                    p = match.iloc[0].get(p_col, np.nan)
                    if pd.notna(p) and p < 0.05: label += " *"
            labels.append(label)
            color = const.DENSITY_COLORS.get(const.DENSITY_CODE_TO_LABEL.get(row.stratum, ""), "#555555")
            ax.errorbar(row.coef, i, xerr=[[row.coef - row.ci_lower], [row.ci_upper - row.coef]], fmt="o", color=color, ecolor=color, markersize=7, capsize=4, linewidth=1.5)
            
        ax.axvline(0, color="#666666", linestyle="--", linewidth=1)
        ax.set_yticks(range(len(d)))
        ax.set_yticklabels(labels)
        ax.set_title(var_labels[var], fontsize=10.8, fontweight="bold")
        ax.grid(axis="x", linestyle="--", alpha=0.25)
        
    fig.supxlabel("SLM coefficient (°C per unit land-cover fraction; not decomposed direct impact)", fontsize=10)
    fig.suptitle("Fig S3. Density-Dependent Land-Cover/LST Associations", fontsize=13, fontweight="bold", y=0.99)
    fig.tight_layout(rect=[0, 0.04, 1, 0.95])
    save_fig(fig, "supp3_density_dependent_lst_associations.png")

def supp4_preferred_spatial_regression_coefficients():
    files = glob.glob(os.path.join(const.OUTPUT_DIR, "step13_regression_results", "**", "coefficients.csv"), recursive=True)
    if not files:
        print("Skipping Supp4: no coefficients.csv found.")
        return
        
    records = []
    for f in files:
        try:
            d = pd.read_csv(f)
            if {"variable", "coef", "ci_lower", "ci_upper"}.issubset(d.columns):
                records.append(d)
        except Exception: pass
        
    if not records: return
    df = pd.concat(records, ignore_index=True)
    
    if "is_preferred_for_reporting" in df.columns:
        pref_mask = df["is_preferred_for_reporting"].astype(str).str.strip().str.lower().isin(["true", "1", "yes", "t"])
        if pref_mask.any(): df = df[pref_mask].copy()
            
    wanted = ["tree_frac", "grass_frac", "crop_frac", "water_frac", "bare_frac", "log_dist_jrc_water"]
    labels = {"tree_frac": "Tree cover", "grass_frac": "Grass cover", "crop_frac": "Cropland", "water_frac": "Water cover", "bare_frac": "Bare soil", "log_dist_jrc_water": "Log distance to water"}
    
    df = df[df["variable"].isin(wanted)].copy().drop_duplicates(subset=["variable"], keep="first")
    if df.empty: return
    if "p_value" not in df.columns: df["p_value"] = np.nan
    df["_color"] = df["p_value"].apply(_sig_color)
    
    y_order = [v for v in wanted if v in set(df["variable"])]
    y_pos = {v: i for i, v in enumerate(reversed(y_order))}
    
    fig, ax = plt.subplots(figsize=(8, 5.4))
    ax.axvline(0, color="#666666", linestyle="--", linewidth=1)
    
    for _, row in df.iterrows():
        v = row["variable"]
        if v not in y_pos: continue
        ax.errorbar(row["coef"], y_pos[v], xerr=[[row["coef"] - row["ci_lower"]], [row["ci_upper"] - row["coef"]]], fmt="o", color=row["_color"], ecolor=row["_color"], markersize=7, capsize=4, linewidth=1.5)
        
    ax.set_yticks([y_pos[v] for v in reversed(y_order)])
    ax.set_yticklabels([labels[v] for v in reversed(y_order)])
    ax.set_xlabel("Estimated coefficient (°C per unit change in predictor)")
    ax.set_title("Fig S4. Preferred Spatial Regression Coefficients", fontsize=12, fontweight="bold")
    ax.grid(axis="x", linestyle="--", alpha=0.25)
    
    handles = [mpatches.Patch(color="#C0392B", label="p < 0.01"), mpatches.Patch(color="#E67E22", label="p < 0.05"), mpatches.Patch(color="#AAAAAA", label="Not significant")]
    ax.legend(handles=handles, loc="lower right", title="Significance", fontsize=8.5)
    fig.tight_layout(rect=[0, 0.04, 1, 1])
    save_fig(fig, "supp4_preferred_spatial_regression_coefficients.png")

def make_regression_figures():
    set_style()
    print("Generating regression figures...")
    supp3_density_dependent_lst_associations()
    supp4_preferred_spatial_regression_coefficients()
    print("Regression figures complete.")

if __name__ == "__main__":
    make_regression_figures()
