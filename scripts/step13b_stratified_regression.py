"""
Step 13B: Stratified spatial regression (Heterogeneity / Robustness check).

Input:
  config.ACCESSIBILITY_GPKG
    hanoi_grid_250m_walking_access.gpkg from Step 12

Output:
  ../outputs/step13b_stratified_regression/
    stratified_model_summary.csv
    stratified_coefficients_all.csv
    <stratum>/<model>/coefficients.csv, ols_summary.txt, etc.

Purpose:
  This step repeats the main Step 13 models separately for three population
  density strata (Low / Medium / High-density wards). Its primary role is as
  a **heterogeneity check and robustness analysis**, NOT a replacement for the
  city-wide Step 13 models.

  Key methodological caveat (for paper reporting):
    Density strata are spatially fragmented (non-contiguous), which may cause
    longer KNN links within each stratum than for the full city-wide sample.
    Interpret stratum-level results with appropriate caution and report KNN
    distance diagnostics alongside coefficients.

Strata:
  Low_Density    : pop density < 5,000 ppl/km²
  Medium_Density : 5,000 – 15,000 ppl/km²
  High_Density   : > 15,000 ppl/km²
"""


import os
import sys
import glob
import warnings
import shutil

import numpy as np
import pandas as pd
import libpysal as ps
import matplotlib.pyplot as plt

import config

from step13_spatial_regression import (
    _read_input_grid,
    _prepare_regression_variables,
    _find_first_existing_column,
    run_regression_submodel,
    MIN_CELL_AREA_M2,
    EXCLUDE_EDGE_CELLS,
    KNN_K,
)

warnings.filterwarnings("ignore", category=RuntimeWarning)
sys.stdout.reconfigure(encoding="utf-8")


OUTPUT_DIR = os.path.join(config.OUTPUT_DIR, "step13b_stratified_regression")

LOW_DENSITY_UPPER = 5000
MEDIUM_DENSITY_UPPER = 15000

# Minimum observations per stratum to attempt regression
MIN_STRATUM_N = getattr(config, "STRATIFIED_REGRESSION_MIN_N", 200)

# KNN distance sanity-check thresholds (metres)
KNN_MEAN_DISTANCE_WARNING_M = getattr(
    config,
    "STRATIFIED_KNN_MEAN_DISTANCE_WARNING_M",
    1000,
)

KNN_MAX_DISTANCE_WARNING_M = getattr(
    config,
    "STRATIFIED_KNN_MAX_DISTANCE_WARNING_M",
    3000,
)

def _prepare_base_sample(gdf):
    base = gdf.copy()

    if "area_m2" in base.columns:
        base = base[
            pd.to_numeric(base["area_m2"], errors="coerce") >= MIN_CELL_AREA_M2
        ].copy()

    if "edge_cell" in base.columns and EXCLUDE_EDGE_CELLS:
        base = base[base["edge_cell"] == 0].copy()

    if "lst_valid" in base.columns:
        base = base[base["lst_valid"] == 1].copy()

    return base


def _add_density_variables(gdf):
    gdf = gdf.copy()

    # Prefer pre-computed column, then derive from population / area.
    if "pop_density_km2" in gdf.columns:
        density = pd.to_numeric(gdf["pop_density_km2"], errors="coerce")
    elif "density" in gdf.columns:
        density = pd.to_numeric(gdf["density"], errors="coerce")
    else:
        if "population" not in gdf.columns or "area_m2" not in gdf.columns:
            raise ValueError(
                "Need pop_density_km2, density, OR both population and area_m2 "
                "to calculate population density."
            )
        pop = pd.to_numeric(gdf["population"], errors="coerce").fillna(0)
        area_km2 = pd.to_numeric(gdf["area_m2"], errors="coerce") / 1_000_000
        density = pop / area_km2

    gdf["pop_density_km2"] = density.replace([np.inf, -np.inf], np.nan)

    conditions = [
        (gdf["pop_density_km2"] < LOW_DENSITY_UPPER),

        (gdf["pop_density_km2"] >= LOW_DENSITY_UPPER)
        & (gdf["pop_density_km2"] <= MEDIUM_DENSITY_UPPER),

        gdf["pop_density_km2"] > MEDIUM_DENSITY_UPPER,
    ]

    choices = ["Low_Density", "Medium_Density", "High_Density"]

    gdf["density_stratum"] = np.select(conditions, choices, default="Unknown")

    return gdf


def _collect_coefficients(output_dir):
    coef_files = glob.glob(os.path.join(output_dir, "**", "coefficients.csv"), recursive=True)

    records = []

    for path in coef_files:
        try:
            df = pd.read_csv(path)
            df["model_folder"] = os.path.basename(os.path.dirname(path))
            records.append(df)
        except Exception as e:
            print(f"WARNING: Could not read coefficients file {path}: {e}")

    if not records:
        return pd.DataFrame()

    return pd.concat(records, ignore_index=True)


def _prepare_interaction_terms(df, x_vars):
    df = df.copy()
    df["dummy_Medium_Density"] = (df["density_stratum"] == "Medium_Density").astype(float)
    df["dummy_High_Density"] = (df["density_stratum"] == "High_Density").astype(float)
    
    interact_vars = list(x_vars) + ["dummy_Medium_Density", "dummy_High_Density"]
    
    for v in x_vars:
        med_col = f"{v}_x_Medium_Density"
        high_col = f"{v}_x_High_Density"
        df[med_col] = df[v] * df["dummy_Medium_Density"]
        df[high_col] = df[v] * df["dummy_High_Density"]
        interact_vars.extend([med_col, high_col])
        
    return df, interact_vars


def _calculate_marginal_effects(model, model_name, x_vars, model_label):
    betas = np.asarray(model.betas).flatten()
    
    if hasattr(model, 'name_x') and len(model.name_x) > 0:
        names = list(model.name_x)
        if len(betas) > len(names):
            extra_n = len(betas) - len(names)
            for i in range(extra_n):
                if model_name == "SLM" and i == extra_n - 1:
                    names.append("W_lag_dependent")
                elif model_name == "SEM" and i == extra_n - 1:
                    names.append("lambda")
                else:
                    names.append(f"extra_coef_{i + 1}")
    else:
        names = ["CONSTANT"] + list(x_vars)

    name_to_idx = {name: idx for idx, name in enumerate(names)}
    
    vm = None
    if hasattr(model, 'vm'):
        vm = np.asarray(model.vm)
        
    records = []
    
    for v in x_vars:
        if v not in name_to_idx:
            continue
        
        i_base = name_to_idx[v]
        
        # Low Density (Reference)
        coef_low = betas[i_base]
        se_low = np.nan
        if vm is not None:
            se_low = np.sqrt(vm[i_base, i_base])
            
        records.append({
            "model_label": model_label,
            "model_name": model_name,
            "variable": v,
            "stratum": "Low_Density",
            "coef": coef_low,
            "std_err": se_low,
            "ci_lower": coef_low - 1.96 * se_low if pd.notna(se_low) else np.nan,
            "ci_upper": coef_low + 1.96 * se_low if pd.notna(se_low) else np.nan,
        })
        
        # Medium Density
        med_interact = f"{v}_x_Medium_Density"
        if med_interact in name_to_idx:
            i_med = name_to_idx[med_interact]
            coef_med = coef_low + betas[i_med]
            se_med = np.nan
            if vm is not None:
                se_med = np.sqrt(vm[i_base, i_base] + vm[i_med, i_med] + 2 * vm[i_base, i_med])
            
            records.append({
                "model_label": model_label,
                "model_name": model_name,
                "variable": v,
                "stratum": "Medium_Density",
                "coef": coef_med,
                "std_err": se_med,
                "ci_lower": coef_med - 1.96 * se_med if pd.notna(se_med) else np.nan,
                "ci_upper": coef_med + 1.96 * se_med if pd.notna(se_med) else np.nan,
            })
            
        # High Density
        high_interact = f"{v}_x_High_Density"
        if high_interact in name_to_idx:
            i_high = name_to_idx[high_interact]
            coef_high = coef_low + betas[i_high]
            se_high = np.nan
            if vm is not None:
                se_high = np.sqrt(vm[i_base, i_base] + vm[i_high, i_high] + 2 * vm[i_base, i_high])
            
            records.append({
                "model_label": model_label,
                "model_name": model_name,
                "variable": v,
                "stratum": "High_Density",
                "coef": coef_high,
                "std_err": se_high,
                "ci_lower": coef_high - 1.96 * se_high if pd.notna(se_high) else np.nan,
                "ci_upper": coef_high + 1.96 * se_high if pd.notna(se_high) else np.nan,
            })
            
    return pd.DataFrame(records)


def _extract_interaction_significance(model, model_name, x_vars, model_label):
    betas = np.asarray(model.betas).flatten()
    
    if hasattr(model, 'name_x') and len(model.name_x) > 0:
        names = list(model.name_x)
        if len(betas) > len(names):
            extra_n = len(betas) - len(names)
            for i in range(extra_n):
                if model_name == "SLM" and i == extra_n - 1:
                    names.append("W_lag_dependent")
                elif model_name == "SEM" and i == extra_n - 1:
                    names.append("lambda")
                else:
                    names.append(f"extra_coef_{i + 1}")
    else:
        names = ["CONSTANT"] + list(x_vars)

    try:
        if hasattr(model, 'z_stat'):
            p_values = [z[1] for z in model.z_stat]
        elif hasattr(model, 't_stat'):
            p_values = [z[1] for z in model.t_stat]
        else:
            p_values = [np.nan] * len(betas)
    except Exception:
        p_values = [np.nan] * len(betas)

    name_to_p = {name: p for name, p in zip(names, p_values)}
    
    records = []
    for v in x_vars:
        med_p = name_to_p.get(f"{v}_x_Medium_Density", np.nan)
        high_p = name_to_p.get(f"{v}_x_High_Density", np.nan)
        
        records.append({
            "model_label": model_label,
            "model_name": model_name,
            "variable": v,
            "Medium_p_value": med_p,
            "High_p_value": high_p,
            "Medium_sig": "Yes" if pd.notna(med_p) and med_p < 0.05 else "No",
            "High_sig": "Yes" if pd.notna(high_p) and high_p < 0.05 else "No"
        })
        
    return pd.DataFrame(records)


def _extract_model_diagnostics(result_dict, model_obj, label):
    preferred_model = result_dict["preferred_model"]
    
    ak_stat = np.nan
    ak_p = np.nan
    if preferred_model in ["SLM", "SEM"] and hasattr(model_obj, "ak_test"):
        try:
            ak_stat = model_obj.ak_test[0]
            ak_p = model_obj.ak_test[1]
        except Exception:
            pass
            
    return {
        "model_label": label,
        "preferred_model": preferred_model,
        "observations": result_dict["observations"],
        "ols_r2": result_dict["ols_r2"],
        "ols_moran_i": result_dict["moran_i"],
        "ols_moran_p": result_dict["moran_p"],
        "robust_lm_lag_p": result_dict["robust_lm_lag_p"],
        "robust_lm_error_p": result_dict["robust_lm_error_p"],
        "preferred_residual_moran_i": result_dict.get("preferred_moran_i", np.nan),
        "preferred_residual_moran_p": result_dict.get("preferred_moran_p", np.nan),
        "anselin_kelejian_stat": ak_stat,
        "anselin_kelejian_p": ak_p,
    }


def _plot_marginal_effects(df_marginal, output_path, sig_df=None):
    df_plot = df_marginal.copy()
    
    # Filter to only the core land-cover variables (exclude log_dist_jrc_water)
    plot_vars = ["tree_frac", "grass_frac", "crop_frac", "water_frac", "bare_frac"]
    df_plot = df_plot[df_plot["variable"].isin(plot_vars)]
    
    # If both SLM and OLS exist, prefer SLM
    if "SLM" in df_plot["model_name"].values:
        df_plot = df_plot[df_plot["model_name"] == "SLM"]
    elif "SEM" in df_plot["model_name"].values:
        df_plot = df_plot[df_plot["model_name"] == "SEM"]
    else:
        df_plot = df_plot[df_plot["model_name"] == "OLS"]
        
    variables = df_plot["variable"].unique()
    if len(variables) == 0:
        return
        
    fig, axes = plt.subplots(len(variables), 1, figsize=(8, 3 * len(variables)), sharex=False)
    if len(variables) == 1:
        axes = [axes]
        
    # sns.set_theme(style="whitegrid")
    
    palette = {
        "Low_Density": "#4575b4",      # Blue
        "Medium_Density": "#fdae61",   # Orange
        "High_Density": "#d73027",     # Red
    }
    
    for ax, var in zip(axes, variables):
        df_var = df_plot[df_plot["variable"] == var].copy()
        
        df_var["stratum"] = pd.Categorical(df_var["stratum"], categories=["Low_Density", "Medium_Density", "High_Density"], ordered=True)
        df_var = df_var.sort_values("stratum")
        
        yticklabels = []
        for idx, row in enumerate(df_var.itertuples()):
            y_pos = idx
            ax.errorbar(
                row.coef,
                y_pos,
                xerr=[[row.coef - row.ci_lower], [row.ci_upper - row.coef]] if pd.notna(row.ci_lower) else None,
                fmt="o",
                color=palette.get(row.stratum, "black"),
                markersize=8,
                capsize=5,
                linewidth=2,
            )
            
            # Check interaction significance
            is_sig = False
            if sig_df is not None and row.stratum != "Low_Density":
                match = sig_df[(sig_df["model_label"] == row.model_label) & (sig_df["variable"] == row.variable)]
                if not match.empty:
                    p_val = match.iloc[0]["Medium_p_value"] if row.stratum == "Medium_Density" else match.iloc[0]["High_p_value"]
                    if pd.notna(p_val) and p_val < 0.05:
                        is_sig = True
            
            label = row.stratum.replace("_", " ")
            if is_sig:
                label += "*"
            yticklabels.append(label)
            
        ax.axvline(0, color="grey", linestyle="--", linewidth=1)
        ax.set_yticks(range(len(df_var)))
        ax.set_yticklabels(yticklabels)
        
        var_labels = {
            "tree_frac": "Tree cover",
            "grass_frac": "Grass cover",
            "crop_frac": "Cropland",
            "water_frac": "Water",
            "bare_frac": "Bare land",
        }
        ax.set_title(var_labels.get(var, var), fontsize=11, fontweight="bold")
        ax.set_xlabel("Estimated association with LST (°C)")
        ax.grid(axis="x", linestyle="--", alpha=0.7)
        
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Forest plot saved: {output_path}")


def run_stratified_regression():
    print("--- STEP 13B: STRATIFIED SPATIAL REGRESSION ---")

    if os.path.exists(OUTPUT_DIR):
        print("Cleaning up old stratified regression output directory...")
        # NOTE: This deletes the entire directory. If run is interrupted, partial results will be lost.
        shutil.rmtree(OUTPUT_DIR)

    os.makedirs(OUTPUT_DIR, exist_ok=True)

    print("Loading Step 12 grid...")
    gdf = _read_input_grid()
    print(f"Rows loaded: {len(gdf)}")

    print("Preparing regression variables...")
    gdf = _prepare_regression_variables(gdf)

    dep_var_obs = _find_first_existing_column(
        gdf,
        ["lst_p75_observed", "lst_p75_mean", "LST_p75_C_mean", "LST_p75"],
        "observed LST variable",
    )

    print(f"Observed dependent variable: {dep_var_obs}")

    # Calculate Total Vegetation for RobustVeg model
    if all(col in gdf.columns for col in ["tree_frac", "grass_frac", "crop_frac"]):
        gdf["total_vegetation_frac"] = gdf["tree_frac"] + gdf["grass_frac"] + gdf["crop_frac"]

    # Density classification
    gdf = _add_density_variables(gdf)

    # Apply the same base filters as Step 13 observed models
    base = _prepare_base_sample(gdf).reset_index(drop=True)
    print(f"Total valid observations for stratified regression: {len(base)}")

    print("\nDensity strata counts:")
    print(base["density_stratum"].value_counts(dropna=False).to_string())

    # Build W globally first
    print("\nBuilding Global KNN weight matrix...")
    coords_global = np.column_stack([
        base.geometry.representative_point().x,
        base.geometry.representative_point().y,
    ])
    W_global_raw = ps.weights.KNN.from_array(coords_global, k=KNN_K)
    
    # Re-key global W with grid_ids so subsetting and alignment match grid_id
    id_map = dict(enumerate(base["grid_id"].tolist()))
    new_neighbors = {id_map[i]: [id_map[j] for j in neighbors] for i, neighbors in W_global_raw.neighbors.items()}
    new_weights = {id_map[i]: list(weights) for i, weights in W_global_raw.weights.items()}
    W_global = ps.weights.W(new_neighbors, new_weights, id_order=base["grid_id"].tolist())
    W_global.transform = "r"

    physical_vars = [
        "tree_frac",
        "grass_frac",
        "crop_frac",
        "water_frac",
        "bare_frac",
        "log_dist_jrc_water",
    ]

    robust_veg_vars = [
        "total_vegetation_frac",
        "water_frac",
        "bare_frac",
        "log_dist_jrc_water",
    ]

    models_to_run = [
        ("Model-1A-Physical-Observed", dep_var_obs, physical_vars),
        ("Model-1C-Physical-RobustVeg", dep_var_obs, robust_veg_vars),
    ]

    summary_records = []
    all_marginal_effects = []
    all_interaction_sigs = []
    all_diagnostics = []

    # ==========================================================================
    # 1. RUN PRIMARY GLOBAL INTERACTION MODELS (OPTION D)
    # ==========================================================================
    print("\n=======================================================")
    print("RUNNING PRIMARY GLOBAL INTERACTION MODELS (OPTION D)")
    print("=======================================================")

    for model_base_name, y_var, x_vars in models_to_run:
        label = f"{model_base_name}_Global_Interaction"
        print(f"\nPreparing interaction terms for {label}...")
        
        base_interact, interact_vars = _prepare_interaction_terms(base, x_vars)
        
        try:
            result = run_regression_submodel(
                base_interact,
                y_var,
                interact_vars,
                label,
                OUTPUT_DIR,
                w=W_global,
            )
            result["stratum"] = "Global_Interaction"
            result["status"] = "ok"
            summary_records.append(result)
            
            # Retrieve model object directly from results (no re-estimation)
            model_obj = result["model_obj"]
            preferred_model = result["preferred_model"]
            sub_dir = result["output_dir"]
            
            # Extract marginal effects
            marginal_df = _calculate_marginal_effects(model_obj, preferred_model, x_vars, model_base_name)
            if marginal_df["coef"].isna().any():
                print(f"  WARNING: Some marginal effects are NaN in {label} — interaction terms may have been dropped by VIF.")
            all_marginal_effects.append(marginal_df)
            
            # Extract significance tables
            sig_df = _extract_interaction_significance(model_obj, preferred_model, x_vars, model_base_name)
            all_interaction_sigs.append(sig_df)
            
            # Extract diagnostics
            diag_dict = _extract_model_diagnostics(result, model_obj, label)
            all_diagnostics.append(diag_dict)
            
            # Plot forest plot with significance stars
            _plot_marginal_effects(marginal_df, os.path.join(sub_dir, f"{model_base_name}_marginal_effects.png"), sig_df=sig_df)
            
        except Exception as e:
            print(f"ERROR in {label}: {e}")
            summary_records.append({
                "stratum": "Global_Interaction",
                "label": label,
                "status": "error",
                "error": str(e),
            })

    # Save outputs of Global Interaction models
    if all_marginal_effects:
        marginal_all_df = pd.concat(all_marginal_effects, ignore_index=True)
        marginal_all_path = os.path.join(OUTPUT_DIR, "marginal_effects_by_density.csv")
        marginal_all_df.to_csv(marginal_all_path, index=False, encoding="utf-8-sig")
        print(f"All marginal effects saved: {marginal_all_path}")
        
    if all_interaction_sigs:
        sig_all_df = pd.concat(all_interaction_sigs, ignore_index=True)
        sig_all_path = os.path.join(OUTPUT_DIR, "interaction_significance.csv")
        sig_all_df.to_csv(sig_all_path, index=False, encoding="utf-8-sig")
        print(f"Interaction significance table saved: {sig_all_path}")
        
    # Plot overall forest plot with stars
    if all_marginal_effects:
        model1a_marginal = marginal_all_df[marginal_all_df["model_label"] == "Model-1A-Physical-Observed"]
        if not model1a_marginal.empty:
            sig_for_plot = sig_all_df if all_interaction_sigs else None
            _plot_marginal_effects(model1a_marginal, os.path.join(OUTPUT_DIR, "supp6_regression_forest_plot.png"), sig_df=sig_for_plot)

    # ==========================================================================
    # 2. RUN SEPARATE STRATIFIED MODELS (OPTION B - ROBUSTNESS CHECKS)
    # ==========================================================================
    print("\n=======================================================")
    print("RUNNING SEPARATE STRATIFIED MODELS (OPTION B)")
    print("=======================================================")

    strata = {
        "Low_Density": base[base["density_stratum"] == "Low_Density"].copy(),
        "Medium_Density": base[base["density_stratum"] == "Medium_Density"].copy(),
        "High_Density": base[base["density_stratum"] == "High_Density"].copy(),
    }

    island_records = []

    for stratum_name, stratum_df in strata.items():
        print(f"\nSTRATUM: {stratum_name} (n={len(stratum_df)})")
        
        if len(stratum_df) < MIN_STRATUM_N:
            print(f"Skipping {stratum_name}; n < {MIN_STRATUM_N}")
            continue

        # Option B: Subset the global weight matrix to stratum grid_ids
        stratum_ids = stratum_df["grid_id"].tolist()
        w_sub = ps.weights.w_subset(W_global, stratum_ids)
        
        # Remove islands
        islands = [k for k, neighbors in w_sub.neighbors.items() if len(neighbors) == 0]
        island_count = len(islands)
        
        if island_count > 0:
            print(f"Removing {island_count} island cells with no neighbors in stratum...")
            valid_ids = [idx for idx in stratum_ids if idx not in islands]
            w_sub = ps.weights.w_subset(W_global, valid_ids)
            stratum_df_filtered = stratum_df[stratum_df["grid_id"].isin(valid_ids)].copy()
        else:
            stratum_df_filtered = stratum_df.copy()

        w_sub.transform = "r"

        island_records.append({
            "stratum": stratum_name,
            "original_N": len(stratum_df),
            "islands_removed": island_count,
            "final_N": len(stratum_df_filtered),
            "pct_removed": (island_count / len(stratum_df)) * 100
        })

        for model_base_name, y_var, x_vars in models_to_run:
            label = f"{model_base_name}_{stratum_name}"
            
            missing_x = [v for v in x_vars if v not in stratum_df_filtered.columns]
            if missing_x:
                print(f"Skipping {label}; missing variables: {missing_x}")
                continue

            try:
                result = run_regression_submodel(
                    stratum_df_filtered,
                    y_var,
                    x_vars,
                    label,
                    OUTPUT_DIR,
                    w=w_sub,
                )
                
                knn_mean_dist = result.get("knn_mean_distance_m", np.nan)
                knn_max_dist = result.get("knn_max_distance_m", np.nan)
                
                flags = []
                if pd.notna(knn_mean_dist) and knn_mean_dist > KNN_MEAN_DISTANCE_WARNING_M:
                    flags.append(f"KNN mean distance {knn_mean_dist:.0f} m > {KNN_MEAN_DISTANCE_WARNING_M} m")
                if pd.notna(knn_max_dist) and knn_max_dist > KNN_MAX_DISTANCE_WARNING_M:
                    flags.append(f"KNN max distance {knn_max_dist:.0f} m > {KNN_MAX_DISTANCE_WARNING_M} m")

                for flag in flags:
                    print(f"WARNING (Subset W): {flag}")

                result["stratum"] = stratum_name
                result["status"] = "ok"
                result["knn_flags"] = " | ".join(flags)
                summary_records.append(result)

                # Extract diagnostics for the stratified robustness submodel
                diag_dict = _extract_model_diagnostics(result, result["model_obj"], label)
                all_diagnostics.append(diag_dict)

            except Exception as e:
                print(f"ERROR in {label}: {e}")
                summary_records.append({
                    "stratum": stratum_name,
                    "label": label,
                    "status": "error",
                    "error": str(e),
                })

    # Save outputs
    summary_df = pd.DataFrame(summary_records)
    summary_path = os.path.join(OUTPUT_DIR, "stratified_model_summary.csv")
    summary_df.to_csv(summary_path, index=False, encoding="utf-8-sig")

    coef_df = _collect_coefficients(OUTPUT_DIR)
    if not coef_df.empty:
        coef_path = os.path.join(OUTPUT_DIR, "stratified_coefficients_all.csv")
        coef_df.to_csv(coef_path, index=False, encoding="utf-8-sig")
        print(f"Stratified coefficients saved: {coef_path}")

    # Save diagnostics for all models
    if all_diagnostics:
        diag_df = pd.DataFrame(all_diagnostics)
        diag_path = os.path.join(OUTPUT_DIR, "global_model_diagnostics.csv")
        diag_df.to_csv(diag_path, index=False, encoding="utf-8-sig")
        print(f"All model diagnostics saved to: {diag_path}")

    # Save island statistics
    if island_records:
        island_df = pd.DataFrame(island_records)
        island_path = os.path.join(OUTPUT_DIR, "island_statistics.csv")
        island_df.to_csv(island_path, index=False, encoding="utf-8-sig")
        print(f"Island statistics table saved: {island_path}")

    print("\nDONE STEP 13B.")
    print(f"Summary saved: {summary_path}")


if __name__ == "__main__":
    run_stratified_regression()
