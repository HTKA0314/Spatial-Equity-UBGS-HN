"""
Step 13: Spatial regression for LST and green-blue cooling/accessibility analysis.

Input:
  config.ACCESSIBILITY_GPKG
    hanoi_grid_250m_walking_access.gpkg from Step 12

Output:
  ../outputs/step13_spatial_regression/

Main dependent variables:
  - lst_p75_observed : observed LST only
  - lst_p75_filled   : IDW-filled LST sensitivity

Model groups:
  1A Physical Observed
  1B Physical Filled
  1C Physical RobustVeg

Modeling logic:
  - Land-cover classes are converted to fractions 0-1.
  - In physical models, built-up is the reference class.
    Therefore built_frac is omitted.
  - log_dist_jrc_water is used instead of legacy log_dist_water.
  - Accessible green-blue variables come from Step 10.
  - Network accessibility variables come from Step 12.
  - OLS is always run first with spatial diagnostics.
  - SLM/SEM are selected using Anselin LM / robust LM decision rules.
  - coefficients.csv is saved for each model so plotting scripts do not parse text summaries.

VIF policy (two-threshold):
  - VIF > 5  : warn in console and flag in vif_results.csv
  - VIF > 10 : drop variable and refit (hard threshold)
  Primary model k=8 (config.REGRESSION_KNN_K). Sensitivity runs k=4.
"""

import os
import sys
import warnings

import numpy as np
import pandas as pd
import geopandas as gpd

import libpysal as ps
from spreg import OLS, GM_Lag, GM_Error
from statsmodels.stats.outliers_influence import variance_inflation_factor

import config

sys.stdout.reconfigure(encoding="utf-8")
warnings.filterwarnings("ignore", category=RuntimeWarning)


# ---------------------------------------------------------------------
# CONFIG
# ---------------------------------------------------------------------

INPUT_GRID = getattr(
    config,
    "ACCESSIBILITY_GPKG",
    os.path.join(config.DATA_DIR, "hanoi_grid_250m_walking_access.gpkg"),
)

OUTPUT_DIR = os.path.join(config.OUTPUT_DIR, "step13_spatial_regression")

# FIX #1: default k=8 (primary model), sensitivity runs k=4.
# Previously defaulted to 6, inconsistent with agreed methodology.
KNN_K = getattr(config, "REGRESSION_KNN_K", 8)

MIN_CELL_AREA_M2 = getattr(config, "REGRESSION_MIN_CELL_AREA_M2", 31250)
EXCLUDE_EDGE_CELLS = getattr(config, "REGRESSION_EXCLUDE_EDGE_CELLS", True)

# FIX #2: Two-threshold VIF policy (warn at 5, drop at 10).
# Previously only had a single warn threshold at 10.
VIF_WARN_THRESHOLD  = getattr(config, "VIF_WARN_THRESHOLD",  5.0)
VIF_DROP_THRESHOLD  = getattr(config, "VIF_DROP_THRESHOLD",  10.0)

# Keep backward-compat alias so any external code using the old name still works.
VIF_WARNING_THRESHOLD = VIF_DROP_THRESHOLD


# ---------------------------------------------------------------------
# BASIC UTILITIES
# ---------------------------------------------------------------------

def _read_input_grid() -> gpd.GeoDataFrame:
    if not os.path.exists(INPUT_GRID):
        raise FileNotFoundError(
            f"Input grid not found: {INPUT_GRID}. Run Step 12 first."
        )

    try:
        gdf = gpd.read_file(INPUT_GRID, layer="grid_walking_access")
    except Exception:
        gdf = gpd.read_file(INPUT_GRID)

    if gdf.crs is None:
        raise ValueError("Input grid has no CRS.")

    if str(gdf.crs).upper() != str(config.PROJ_CRS).upper():
        gdf = gdf.to_crs(config.PROJ_CRS)

    if "grid_id" not in gdf.columns:
        raise ValueError("grid_id column missing from input grid.")

    gdf["grid_id"] = gdf["grid_id"].astype(int)

    return gdf


def _find_first_existing_column(gdf: gpd.GeoDataFrame, candidates: list[str], label: str) -> str:
    col = next((c for c in candidates if c in gdf.columns), None)

    if col is None:
        raise ValueError(
            f"Cannot find {label}. Tried: {candidates}\n"
            f"Available columns: {list(gdf.columns)}"
        )

    return col


def _to_fraction(series: pd.Series) -> pd.Series:
    """
    Convert percentage 0-100 to fraction 0-1 if needed.
    If already 0-1, keep as is.
    """
    s = pd.to_numeric(series, errors="coerce").astype(float)
    max_val = s.max(skipna=True)

    if pd.notna(max_val) and max_val > 1.5:
        s = s / 100.0

    return s.clip(lower=0.0, upper=1.0)


def _safe_log1p(series: pd.Series) -> pd.Series:
    s = pd.to_numeric(series, errors="coerce").astype(float)
    s = s.clip(lower=0)
    return np.log1p(s)


def _model_dir(output_dir: str, label: str) -> str:
    safe = (
        label.lower()
        .replace(" ", "_")
        .replace("-", "_")
        .replace("/", "_")
    )
    path = os.path.join(output_dir, safe)
    os.makedirs(path, exist_ok=True)
    return path


def _write_text(path: str, content: str):
    with open(path, "w", encoding="utf-8") as f:
        f.write(content)


# ---------------------------------------------------------------------
# DIAGNOSTICS
# ---------------------------------------------------------------------

def _check_vif(df: pd.DataFrame, variables: list[str]) -> pd.DataFrame:
    """
    Calculate VIF for all variables.
    Returns DataFrame with columns: variable, VIF, flag
      flag: 'OK' | 'WARN (VIF>5)' | 'DROP (VIF>10)'
    """
    X = df[variables].copy()
    X = X.replace([np.inf, -np.inf], np.nan).dropna()

    if len(X) == 0:
        raise ValueError("No rows available for VIF calculation.")

    X = X.assign(const=1.0)

    records = []

    for i, col in enumerate(X.columns):
        if col == "const":
            continue

        try:
            vif = variance_inflation_factor(X.values, i)
        except Exception:
            vif = np.nan

        if pd.isna(vif):
            flag = "UNKNOWN"
        elif vif > VIF_DROP_THRESHOLD:
            flag = f"DROP (VIF>{VIF_DROP_THRESHOLD})"
        elif vif > VIF_WARN_THRESHOLD:
            flag = f"WARN (VIF>{VIF_WARN_THRESHOLD})"
        else:
            flag = "OK"

        records.append({"variable": col, "VIF": vif, "flag": flag})

    return pd.DataFrame(records)


def _apply_vif_drop(
    df: pd.DataFrame,
    variables: list[str],
    vif_df: pd.DataFrame,
    label: str,
) -> list[str]:
    """
    FIX #2: Drop variables with VIF > VIF_DROP_THRESHOLD iteratively.
    Returns the cleaned variable list.
    """
    remaining = list(variables)
    iteration = 0

    while True:
        drop_candidates = vif_df[vif_df["flag"].str.startswith("DROP", na=False)]["variable"].tolist()

        if not drop_candidates:
            break

        # Drop the variable with the highest VIF first
        worst = vif_df.loc[
            vif_df["variable"].isin(drop_candidates), ["variable", "VIF"]
        ].sort_values("VIF", ascending=False).iloc[0]["variable"]

        print(
            f"[{label}] VIF DROP iteration {iteration + 1}: "
            f"removing '{worst}' (VIF={vif_df.loc[vif_df['variable']==worst,'VIF'].values[0]:.2f} > {VIF_DROP_THRESHOLD})"
        )
        remaining.remove(worst)
        iteration += 1

        if len(remaining) < 2:
            print(f"[{label}] WARNING: Only {len(remaining)} variable(s) remain after VIF dropping. Stopping.")
            break

        # Recompute VIF on remaining variables
        vif_df = _check_vif(df, remaining)

    return remaining


def _build_knn_weights(gdf: gpd.GeoDataFrame, k: int):
    coords = np.column_stack([
        gdf.geometry.representative_point().x,
        gdf.geometry.representative_point().y,
    ])

    w = ps.weights.KNN.from_array(coords, k=k)
    w.transform = "r"

    neighbor_distances = []

    for i, neighbors in w.neighbors.items():
        xi, yi = coords[i]

        for j in neighbors:
            xj, yj = coords[j]
            neighbor_distances.append(np.sqrt((xi - xj) ** 2 + (yi - yj) ** 2))

    return w, np.array(neighbor_distances)


def _choose_spatial_model(ols_model) -> str:
    """
    Anselin decision rule using LM and robust LM tests.
    """
    lm_lag_p = ols_model.lm_lag[1]
    lm_err_p = ols_model.lm_error[1]
    rlm_lag_p = ols_model.rlm_lag[1]
    rlm_err_p = ols_model.rlm_error[1]

    if lm_lag_p < 0.05 and lm_err_p >= 0.05:
        return "SLM"

    if lm_err_p < 0.05 and lm_lag_p >= 0.05:
        return "SEM"

    if lm_lag_p < 0.05 and lm_err_p < 0.05:
        if rlm_lag_p < 0.05 and rlm_err_p >= 0.05:
            return "SLM"

        if rlm_err_p < 0.05 and rlm_lag_p >= 0.05:
            return "SEM"

        if rlm_lag_p < 0.05 and rlm_err_p < 0.05:
            return "BOTH"

        return "AMBIGUOUS"

    return "OLS"


def _preferred_model_for_reporting(chosen_model: str) -> str:
    """
    Reporting rule:
      - Clear Anselin result → use that model.
      - BOTH: SLM is primary (theoretical: UHI spillovers favour lag formulation);
              SEM retained as sensitivity. See model_choice.txt for rationale.
      - AMBIGUOUS: SEM as conservative primary; SLM as sensitivity.
    """
    if chosen_model == "BOTH":
        return "SLM"

    if chosen_model == "AMBIGUOUS":
        return "SEM"

    return chosen_model


# FIX #3: document the BOTH / AMBIGUOUS decision in model_choice.txt
_MODEL_CHOICE_RATIONALE = {
    "SLM": (
        "Robust LM-Lag significant, Robust LM-Error not significant. "
        "Anselin-Florax rule unambiguously selects SLM."
    ),
    "SEM": (
        "Robust LM-Error significant, Robust LM-Lag not significant. "
        "Anselin-Florax rule unambiguously selects SEM."
    ),
    "BOTH": (
        "Both Robust LM-Lag and Robust LM-Error are significant (the 'BOTH' case). "
        "Anselin-Florax rule does not unambiguously select one model. "
        "SLM is chosen as primary specification on theoretical grounds: "
        "urban heat island effects exhibit documented spatial spillovers, "
        "making the spatial lag formulation more consistent with the underlying "
        "physical process. SEM is retained as sensitivity check."
    ),
    "AMBIGUOUS": (
        "Both Robust LM tests are significant, but neither clearly dominates. "
        "SEM is chosen as the conservative primary specification. "
        "SLM is retained as sensitivity check."
    ),
    "OLS": (
        "Neither LM-Lag nor LM-Error is significant at p<0.05. "
        "OLS is adequate; no spatial model required."
    ),
}


# ---------------------------------------------------------------------
# COEFFICIENT EXTRACTION
# ---------------------------------------------------------------------

def _extract_spreg_coefficients(
    model,
    model_name: str,
    x_names: list[str],
    model_gdf: pd.DataFrame = None,
    y_name: str = None
) -> pd.DataFrame:
    """
    Extract coefficients from spreg model objects.

    Works for:
      - spreg.OLS
      - spreg.GM_Lag
      - spreg.GM_Error

    Saves coefficients in machine-readable CSV for plotting/reporting.
    """
    betas = np.asarray(model.betas).flatten()

    if hasattr(model, 'name_x') and len(model.name_x) == len(betas):
        names = list(model.name_x)
    else:
        names = ["CONSTANT"] + list(x_names)
        if len(betas) > len(names):
            extra_n = len(betas) - len(names)
            for i in range(extra_n):
                if model_name == "SLM" and i == extra_n - 1:
                    names.append("W_lag_dependent")
                elif model_name == "SEM" and i == extra_n - 1:
                    names.append("lambda")
                else:
                    names.append(f"extra_coef_{i + 1}")

    try:
        std_err = np.asarray(model.std_err).flatten()
    except Exception:
        std_err = np.full(len(betas), np.nan)

    try:
        if hasattr(model, 'z_stat'):
            z_stats = model.z_stat
        elif hasattr(model, 't_stat'):
            z_stats = model.t_stat
        else:
            z_stats = [(np.nan, np.nan)] * len(betas)

        stat_values = [z[0] for z in z_stats]
        p_values = [z[1] for z in z_stats]
    except Exception:
        stat_values = [np.nan] * len(betas)
        p_values = [np.nan] * len(betas)

    n = min(len(names), len(betas))

    records = []

    y_sd = np.nan
    if model_gdf is not None and y_name is not None and y_name in model_gdf.columns:
        y_sd = model_gdf[y_name].std()

    for i in range(n):
        se = std_err[i] if i < len(std_err) else np.nan
        stat = stat_values[i] if i < len(stat_values) else np.nan
        p = p_values[i] if i < len(p_values) else np.nan

        if pd.notna(se):
            ci_lower = betas[i] - 1.96 * se
            ci_upper = betas[i] + 1.96 * se
        else:
            ci_lower = np.nan
            ci_upper = np.nan

        var_name = names[i]
        x_sd = np.nan
        if model_gdf is not None and var_name in model_gdf.columns:
            x_sd = model_gdf[var_name].std()

        std_coef = np.nan
        if pd.notna(x_sd) and pd.notna(y_sd) and y_sd != 0:
            std_coef = betas[i] * (x_sd / y_sd)

        records.append(
            {
                "model_name": model_name,
                "variable": var_name,
                "coef": betas[i],
                "std_coef": std_coef,
                "x_sd": x_sd,
                "y_sd": y_sd,
                "std_err": se,
                "statistic": stat,
                "p_value": p,
                "ci_lower": ci_lower,
                "ci_upper": ci_upper,
            }
        )

    return pd.DataFrame(records)


# ---------------------------------------------------------------------
# VARIABLE PREPARATION
# ---------------------------------------------------------------------

def _make_accessible_fraction(gdf: gpd.GeoDataFrame, frac_col: str, candidate_cols: list[str]) -> gpd.GeoDataFrame:
    col = next((c for c in candidate_cols if c in gdf.columns), None)
    if col and "area_m2" in gdf.columns:
        gdf[frac_col] = (
            pd.to_numeric(gdf[col], errors="coerce")
            / pd.to_numeric(gdf["area_m2"], errors="coerce")
        ).clip(0, 1)
    else:
        gdf[frac_col] = np.nan
        print(f"WARNING: Missing accessible area candidates for {frac_col}; set to NaN.")
    return gdf

def _prepare_regression_variables(gdf: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    """
    Create clean regression variables.
    """
    gdf = gdf.copy()

    # Land-cover fractions
    landcover_pct_cols = [
        "tree_pct",
        "grass_pct",
        "crop_pct",
        "water_pct",
        "built_pct",
        "bare_pct",
        "physical_greenblue_pct",
        "urban_green_pct",
        "vegetation_pct",
        "impervious_disturbed_pct",
    ]

    for col in landcover_pct_cols:
        if col in gdf.columns:
            frac_col = col.replace("_pct", "_frac")
            gdf[frac_col] = _to_fraction(gdf[col])

    # JRC water distance
    if "log_dist_jrc_water" not in gdf.columns:
        if "dist_to_jrc_water_m" in gdf.columns:
            gdf["log_dist_jrc_water"] = _safe_log1p(gdf["dist_to_jrc_water_m"])
        elif "log_dist_water" in gdf.columns:
            gdf["log_dist_jrc_water"] = pd.to_numeric(gdf["log_dist_water"], errors="coerce")
            print("WARNING: Using legacy log_dist_water as log_dist_jrc_water.")
        else:
            raise ValueError("Missing log_dist_jrc_water or dist_to_jrc_water_m.")

    # Accessible green-blue fractions from Step 10
    gdf = _make_accessible_fraction(
        gdf,
        "accessible_gb_A_frac",
        [
            "accessible_gb_area_A_final_m2",
            "accessible_gb_area_A_pixel_m2",
            "accessible_gb_area_A_est_m2",
            "accessible_area_A",
        ],
    )

    gdf = _make_accessible_fraction(
        gdf,
        "accessible_gb_AB_frac",
        [
            "accessible_gb_area_AB_final_m2",
            "accessible_gb_area_AB_pixel_m2",
            "accessible_gb_area_AB_est_m2",
            "accessible_area_AB",
        ],
    )

    # Gap variables
    if "relative_gap_A" not in gdf.columns and "relative_gap" in gdf.columns:
        gdf["relative_gap_A"] = pd.to_numeric(gdf["relative_gap"], errors="coerce")

    if "conversion_ratio_A" not in gdf.columns and "conversion_ratio" in gdf.columns:
        gdf["conversion_ratio_A"] = pd.to_numeric(gdf["conversion_ratio"], errors="coerce")

    # Network distance variables
    if "dist_to_accessible_A_m" not in gdf.columns and "dist_to_accessible_A" in gdf.columns:
        gdf["dist_to_accessible_A_m"] = pd.to_numeric(gdf["dist_to_accessible_A"], errors="coerce")

    if "dist_to_accessible_AB_m" not in gdf.columns and "dist_to_accessible_AB" in gdf.columns:
        gdf["dist_to_accessible_AB_m"] = pd.to_numeric(gdf["dist_to_accessible_AB"], errors="coerce")

    if "dist_to_accessible_A_m" in gdf.columns:
        d = pd.to_numeric(gdf["dist_to_accessible_A_m"], errors="coerce")
        d = d.replace([np.inf, -np.inf], np.nan)
        gdf["log_dist_accessible_A"] = np.log1p(d)

    if "dist_to_accessible_AB_m" in gdf.columns:
        d = pd.to_numeric(gdf["dist_to_accessible_AB_m"], errors="coerce")
        d = d.replace([np.inf, -np.inf], np.nan)
        gdf["log_dist_accessible_AB"] = np.log1p(d)

    return gdf


# ---------------------------------------------------------------------
# SINGLE MODEL RUNNER
# ---------------------------------------------------------------------

def run_regression_submodel(
    gdf: gpd.GeoDataFrame,
    dependent_var: str,
    independent_vars: list[str],
    label: str,
    output_dir: str,
    w = None,
) -> dict:
    print("\n=======================================================")
    print(f"MODEL: {label}")
    print(f"Y = {dependent_var}")
    print("=======================================================")

    needed = [dependent_var] + independent_vars + ["geometry"]

    model_gdf = gdf.replace([np.inf, -np.inf], np.nan).dropna(subset=needed).copy()

    if len(model_gdf) < 100:
        raise ValueError(f"Too few observations for {label}: {len(model_gdf)}")

    if w is not None:
        model_gdf = model_gdf.set_index("grid_id").loc[w.id_order].reset_index()

    print(f"Observations: {len(model_gdf)}")
    print(f"Independent variables: {independent_vars}")

    sub_dir = _model_dir(output_dir, label)

    # Save data used
    model_gdf[["grid_id"] + independent_vars + [dependent_var]].to_csv(
        os.path.join(sub_dir, "model_data_used.csv"),
        index=False,
        encoding="utf-8-sig",
    )

    # FIX #2: Two-threshold VIF with iterative drop
    vif_df = _check_vif(model_gdf, independent_vars)

    warn_vif = vif_df[vif_df["flag"].str.startswith("WARN", na=False)]
    drop_vif = vif_df[vif_df["flag"].str.startswith("DROP", na=False)]

    if not warn_vif.empty:
        print(f"[{label}] VIF WARNING (>{VIF_WARN_THRESHOLD}):")
        print(warn_vif[["variable", "VIF"]].to_string(index=False))

    if not drop_vif.empty:
        print(f"[{label}] VIF DROP candidates (>{VIF_DROP_THRESHOLD}):")
        print(drop_vif[["variable", "VIF"]].to_string(index=False))
        independent_vars = _apply_vif_drop(model_gdf, independent_vars, vif_df, label)
        # Recompute final VIF after dropping
        vif_df = _check_vif(model_gdf, independent_vars)
        print(f"[{label}] Variables after VIF drop: {independent_vars}")

    vif_df.to_csv(os.path.join(sub_dir, "vif_results.csv"), index=False, encoding="utf-8-sig")
    print("\nVIF (final):")
    print(vif_df.to_string(index=False))

    if w is None:
        w, neighbor_distances = _build_knn_weights(model_gdf, KNN_K)
    else:
        coords = {
            row.grid_id: (row.geometry.representative_point().x, row.geometry.representative_point().y)
            for row in model_gdf.itertuples()
        }
        neighbor_distances = []
        for i, neighbors in w.neighbors.items():
            xi, yi = coords[i]
            for j in neighbors:
                xj, yj = coords[j]
                neighbor_distances.append(np.sqrt((xi - xj) ** 2 + (yi - yj) ** 2))
        neighbor_distances = np.array(neighbor_distances)

    print(
        f"\nKNN weights: k={len(w.neighbors[w.id_order[0]]) if w.neighbors else 0}, "
        f"mean distance={neighbor_distances.mean():.1f} m, "
        f"max distance={neighbor_distances.max():.1f} m"
    )

    y = model_gdf[dependent_var].to_numpy().reshape(-1, 1)
    X = model_gdf[independent_vars].to_numpy()

    # OLS with spatial diagnostics
    ols_model = OLS(
        y,
        X,
        w=w,
        spat_diag=True,
        moran=True,
        name_y=dependent_var,
        name_x=independent_vars,
    )

    print(f"\nOLS R²: {ols_model.r2:.4f}")

    lm_lag_p = ols_model.lm_lag[1]
    lm_err_p = ols_model.lm_error[1]
    rlm_lag_p = ols_model.rlm_lag[1]
    rlm_err_p = ols_model.rlm_error[1]

    print(
        "LM tests:"
        f" LM-Lag p={lm_lag_p:.5f},"
        f" LM-Error p={lm_err_p:.5f},"
        f" Robust LM-Lag p={rlm_lag_p:.5f},"
        f" Robust LM-Error p={rlm_err_p:.5f}"
    )

    moran_i = np.nan
    moran_p = np.nan

    if hasattr(ols_model, "moran_res"):
        try:
            moran_i = ols_model.moran_res[0]
            moran_z = ols_model.moran_res[1]
            moran_p = ols_model.moran_res[2]
            print(f"OLS residual Moran's I: {moran_i:.5f}, Z={moran_z:.3f}, p={moran_p:.5f}")
        except Exception:
            pass

    chosen_model = _choose_spatial_model(ols_model)
    preferred_model = _preferred_model_for_reporting(chosen_model)

    print(f"Chosen model by Anselin rule: {chosen_model}")
    print(f"Preferred model for reporting: {preferred_model}")

    _write_text(os.path.join(sub_dir, "ols_summary.txt"), ols_model.summary)

    # FIX #3: include human-readable rationale for BOTH/AMBIGUOUS in model_choice.txt
    rationale = _MODEL_CHOICE_RATIONALE.get(chosen_model, "")
    model_choice_text = (
        f"Model: {label}\n"
        f"Dependent variable: {dependent_var}\n"
        f"Observations: {len(model_gdf)}\n"
        f"KNN k: {len(w.neighbors[w.id_order[0]]) if w.neighbors else 0}\n"
        f"KNN mean distance: {neighbor_distances.mean():.3f} m\n"
        f"KNN max distance: {neighbor_distances.max():.3f} m\n"
        f"OLS R2: {ols_model.r2:.6f}\n"
        f"LM Lag p: {lm_lag_p:.8f}\n"
        f"LM Error p: {lm_err_p:.8f}\n"
        f"Robust LM Lag p: {rlm_lag_p:.8f}\n"
        f"Robust LM Error p: {rlm_err_p:.8f}\n"
        f"OLS residual Moran I: {moran_i:.8f}\n"
        f"OLS residual Moran p: {moran_p:.8f}\n"
        f"Chosen model: {chosen_model}\n"
        f"Preferred model for reporting: {preferred_model}\n"
        f"\nRationale:\n{rationale}\n"
    )

    _write_text(os.path.join(sub_dir, "model_choice.txt"), model_choice_text)
    _write_text(
        os.path.join(sub_dir, "preferred_model.txt"),
        f"Chosen model: {chosen_model}\nPreferred model for reporting: {preferred_model}\nDecision Rationale: {rationale}\n",
    )

    # Coefficient tables
    coef_tables = []
    coef_tables.append(_extract_spreg_coefficients(ols_model, "OLS", independent_vars, model_gdf, dependent_var))

    # Spatial models
    spatial_summary = ""

    run_slm = chosen_model in ["SLM", "BOTH", "AMBIGUOUS"]
    run_sem = chosen_model in ["SEM", "BOTH", "AMBIGUOUS"]

    slm_model = None
    sem_model = None

    preferred_moran_i = np.nan
    preferred_moran_p = np.nan

    if run_slm:
        print("Running SLM / Spatial Lag model...")
        slm_model = GM_Lag(
            y,
            X,
            w=w,
            name_y=dependent_var,
            name_x=independent_vars,
        )
        spatial_summary += "=== SPATIAL LAG MODEL (SLM / GM_Lag) ===\n"
        spatial_summary += slm_model.summary + "\n\n"
        coef_tables.append(_extract_spreg_coefficients(slm_model, "SLM", independent_vars, model_gdf, dependent_var))

        # Moran's I of SLM residuals
        if hasattr(slm_model, 'u') and slm_model.u is not None:
            try:
                from esda.moran import Moran
                slm_resid = slm_model.u.flatten()
                mi_slm = Moran(slm_resid, w, permutations=0)
                spatial_summary += (
                    f"SLM residual Moran's I: {mi_slm.I:.5f} "
                    f"(E[I]={mi_slm.EI:.5f})\n\n"
                )
                print(f"SLM residual Moran's I: {mi_slm.I:.5f}")
                if preferred_model == "SLM":
                    preferred_moran_i = mi_slm.I
                    preferred_moran_p = mi_slm.p_norm
            except Exception as e:
                print(f"NOTE: Could not compute SLM residual Moran's I: {e}")

    if run_sem:
        print("Running SEM / Spatial Error model...")
        sem_model = GM_Error(
            y,
            X,
            w=w,
            name_y=dependent_var,
            name_x=independent_vars,
        )
        spatial_summary += "=== SPATIAL ERROR MODEL (SEM / GM_Error) ===\n"
        spatial_summary += sem_model.summary + "\n\n"
        coef_tables.append(_extract_spreg_coefficients(sem_model, "SEM", independent_vars, model_gdf, dependent_var))

        # Moran's I of SEM residuals
        if hasattr(sem_model, 'u') and sem_model.u is not None:
            try:
                from esda.moran import Moran
                sem_resid = sem_model.u.flatten()
                mi_sem = Moran(sem_resid, w, permutations=0)
                spatial_summary += (
                    f"SEM residual Moran's I: {mi_sem.I:.5f} "
                    f"(E[I]={mi_sem.EI:.5f})\n\n"
                )
                print(f"SEM residual Moran's I: {mi_sem.I:.5f}")
                if preferred_model == "SEM":
                    preferred_moran_i = mi_sem.I
                    preferred_moran_p = mi_sem.p_norm
            except Exception as e:
                print(f"NOTE: Could not compute SEM residual Moran's I: {e}")

    if spatial_summary:
        _write_text(os.path.join(sub_dir, "spatial_model_summary.txt"), spatial_summary)

    coef_df = pd.concat(coef_tables, ignore_index=True)
    coef_df["chosen_model"] = chosen_model
    coef_df["preferred_model"] = preferred_model
    coef_df["is_preferred_for_reporting"] = coef_df["model_name"] == preferred_model

    if preferred_model == "OLS":
        coef_df["is_preferred_for_reporting"] = coef_df["model_name"] == "OLS"

    coef_df.to_csv(
        os.path.join(sub_dir, "coefficients.csv"),
        index=False,
        encoding="utf-8-sig",
    )

    print(f"Finished: {label}")
    print(f"Saved to: {sub_dir}")

    # Determine the model object to return
    preferred_model_obj = ols_model
    if preferred_model == "SLM" and slm_model is not None:
        preferred_model_obj = slm_model
    elif preferred_model == "SEM" and sem_model is not None:
        preferred_model_obj = sem_model

    return {
        "label": label,
        "dependent_var": dependent_var,
        "observations": len(model_gdf),
        "independent_vars": ";".join(independent_vars),
        "chosen_model": chosen_model,
        "preferred_model": preferred_model,
        "ols_r2": ols_model.r2,
        "lm_lag_p": lm_lag_p,
        "lm_error_p": lm_err_p,
        "robust_lm_lag_p": rlm_lag_p,
        "robust_lm_error_p": rlm_err_p,
        "moran_i": moran_i,
        "moran_p": moran_p,
        "preferred_moran_i": preferred_moran_i,
        "preferred_moran_p": preferred_moran_p,
        "knn_k": len(w.neighbors[w.id_order[0]]) if w.neighbors else 0,
        "knn_mean_distance_m": neighbor_distances.mean(),
        "knn_max_distance_m": neighbor_distances.max(),
        "output_dir": sub_dir,
        "model_obj": preferred_model_obj,
    }


# ---------------------------------------------------------------------
# MAIN
# ---------------------------------------------------------------------

def run_spatial_regression():
    print("--- STEP 13: SPATIAL REGRESSION ---")

    os.makedirs(OUTPUT_DIR, exist_ok=True)

    print(f"Loading input grid: {INPUT_GRID}")
    gdf = _read_input_grid()
    print(f"Rows loaded: {len(gdf)}")
    print(f"CRS: {gdf.crs}")

    gdf = _prepare_regression_variables(gdf)

    dep_var_obs = _find_first_existing_column(
        gdf,
        ["lst_p75_observed", "lst_p75_mean", "LST_p75_C_mean", "LST_p75"],
        "observed LST variable",
    )

    if "lst_p75_filled" not in gdf.columns:
        raise ValueError("lst_p75_filled missing; run Step 2 IDW filling first.")
    dep_var_fill = "lst_p75_filled"

    print(f"Observed Y: {dep_var_obs}")
    print(f"Filled Y: {dep_var_fill}")

    if dep_var_obs == dep_var_fill:
        print(
            "WARNING: Observed and filled dependent variables resolved to the same column. "
            "Filled sensitivity models will duplicate observed models and should not be interpreted "
            "as an independent IDW-filled robustness check."
        )

    physical_vars = [
        "tree_frac",
        "grass_frac",
        "crop_frac",
        "water_frac",
        "bare_frac",
        "log_dist_jrc_water",
    ]

    if all(col in gdf.columns for col in ["tree_frac", "grass_frac", "crop_frac"]):
        gdf["total_vegetation_frac"] = gdf["tree_frac"] + gdf["grass_frac"] + gdf["crop_frac"]

    robust_veg_vars = [
        "total_vegetation_frac",
        "water_frac",
        "bare_frac",
        "log_dist_jrc_water",
    ]

    model_specs = [
        ("Model-1A-Physical-Observed", dep_var_obs, physical_vars, "observed"),
        ("Model-1B-Physical-Filled", dep_var_fill, physical_vars, "filled"),
        ("Model-1C-Physical-RobustVeg", dep_var_obs, robust_veg_vars, "observed"),
    ]

    base = gdf.copy()

    if "area_m2" in base.columns:
        base = base[base["area_m2"] >= MIN_CELL_AREA_M2].copy()
        print(f"After area filter >= {MIN_CELL_AREA_M2} m²: {len(base)} rows")

    if "edge_cell" in base.columns and EXCLUDE_EDGE_CELLS:
        base = base[base["edge_cell"] == 0].copy()
        print(f"After excluding edge cells: {len(base)} rows")

    observed_base = base.copy()

    if "lst_valid" in observed_base.columns:
        observed_base = observed_base[observed_base["lst_valid"] == 1].copy()
        print(f"Observed-only rows after lst_valid filter: {len(observed_base)}")
    else:
        print("WARNING: lst_valid missing; observed model will use non-null observed LST only.")

    filled_base = base.copy()
    print(f"Filled/sensitivity rows: {len(filled_base)}")

    if "dist_to_accessible_A_m" in observed_base.columns:
        d_a = pd.to_numeric(observed_base["dist_to_accessible_A_m"], errors="coerce")
        reachable_A = np.isfinite(d_a).sum()
        print(f"Reachable observations for network A model: {reachable_A} / {len(observed_base)}")

    if "dist_to_accessible_AB_m" in observed_base.columns:
        d_ab = pd.to_numeric(observed_base["dist_to_accessible_AB_m"], errors="coerce")
        reachable_AB = np.isfinite(d_ab).sum()
        print(f"Reachable observations for network AB model: {reachable_AB} / {len(observed_base)}")

    print("NOTE: Network-access models are estimated only for cells reachable within the Dijkstra cutoff.")
    summary_records = []

    for label, y_var, x_vars, data_type in model_specs:
        missing_x = [v for v in x_vars if v not in gdf.columns]

        if missing_x:
            print(f"\nSkipping {label}; missing variables: {missing_x}")
            continue

        data = observed_base if data_type == "observed" else filled_base

        try:
            result = run_regression_submodel(
                data,
                y_var,
                x_vars,
                label,
                OUTPUT_DIR,
            )
            result["data_type"] = data_type
            summary_records.append(result)

        except Exception as e:
            print(f"ERROR in {label}: {e}")
            summary_records.append(
                {
                    "label": label,
                    "dependent_var": y_var,
                    "independent_vars": ";".join(x_vars),
                    "data_type": data_type,
                    "error": str(e),
                }
            )

    summary_df = pd.DataFrame(summary_records)

    summary_path = os.path.join(OUTPUT_DIR, "model_choice_summary.csv")
    summary_df.to_csv(summary_path, index=False, encoding="utf-8-sig")

    # KNN sensitivity: k=4 (primary k=8 defined at top of file)
    sens_k_values = [4]
    sens_specs = [
        ("Model-1A-Physical-Observed", dep_var_obs, physical_vars, "observed"),
    ]

    knn_sens_records = []

    global KNN_K  # noqa: PLW0603
    original_k = KNN_K

    try:
        for k_sens_val in sens_k_values:
            KNN_K = k_sens_val
            print(f"\n=== KNN SENSITIVITY: k={k_sens_val} ===")

            for label, y_var, x_vars, data_type in sens_specs:
                missing_x = [v for v in x_vars if v not in gdf.columns]
                if missing_x:
                    continue

                data = observed_base if data_type == "observed" else filled_base
                sens_label = f"{label}-KNN{k_sens_val}"

                try:
                    result = run_regression_submodel(
                        data,
                        y_var,
                        x_vars,
                        sens_label,
                        OUTPUT_DIR,
                    )
                    result["data_type"] = data_type
                    result["knn_sensitivity"] = k_sens_val
                    knn_sens_records.append(result)

                except Exception as e:
                    print(f"ERROR in {sens_label} (k={k_sens_val}): {e}")
                    knn_sens_records.append(
                        {
                            "label": sens_label,
                            "knn_sensitivity": k_sens_val,
                            "error": str(e),
                        }
                    )
    finally:
        KNN_K = original_k


    if knn_sens_records:
        knn_sens_df = pd.DataFrame(knn_sens_records)
        knn_sens_path = os.path.join(OUTPUT_DIR, "knn_sensitivity_summary.csv")
        knn_sens_df.to_csv(knn_sens_path, index=False, encoding="utf-8-sig")
        print(f"KNN sensitivity summary saved: {knn_sens_path}")

    print("\nDONE STEP 13.")
    print(f"Summary saved: {summary_path}")
    print("NOTE: Models with data_type='observed' use lst_valid==1 (primary statistical models).")
    print("NOTE: Models with data_type='filled' use IDW-filled LST (sensitivity / continuous surface).")
    print(f"NOTE: Primary KNN k={original_k}. Sensitivity k={sens_k_values}.")


if __name__ == "__main__":
    run_spatial_regression()