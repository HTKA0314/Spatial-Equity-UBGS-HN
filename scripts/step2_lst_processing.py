"""
Step 2: Merge LST data onto base grid and perform IDW interpolation for missing cells.

Inputs:
  config.BASE_GRID_GPKG
  config.LST_MAIN_CSV_PATH

Outputs:
  config.LST_GRID_GPKG
"""

import os
import sys
import warnings

import numpy as np
import pandas as pd
import geopandas as gpd
from scipy.spatial import cKDTree

import config

sys.stdout.reconfigure(encoding="utf-8")
warnings.filterwarnings("ignore", category=UserWarning)


def _idw_interpolation(
    gdf: gpd.GeoDataFrame,
    val_col: str,
    k: int = 8,
    p: float = 2.0,
) -> pd.Series:
    """
    Fill missing values in val_col using Inverse Distance Weighting (IDW).
    """
    known = gdf[gdf[val_col].notna()].copy()
    unknown = gdf[gdf[val_col].isna()].copy()

    if len(unknown) == 0:
        return gdf[val_col]

    if len(known) == 0:
        print(f"WARNING: No known values for {val_col}. Cannot perform IDW.")
        return gdf[val_col]

    known_coords = np.column_stack([
        known.geometry.representative_point().x,
        known.geometry.representative_point().y,
    ])
    unknown_coords = np.column_stack([
        unknown.geometry.representative_point().x,
        unknown.geometry.representative_point().y,
    ])

    tree = cKDTree(known_coords)
    k_actual = min(k, len(known))
    dist, idx = tree.query(unknown_coords, k=k_actual)

    if k_actual == 1:
        dist = dist.reshape(-1, 1)
        idx = idx.reshape(-1, 1)

    # Export diagnostic metric for nearest observation station (for Methodology chapter)
    nn_dist = dist[:, 0]
    print(f"  IDW diagnostic ({val_col}):")
    print(f"    Nearest-neighbour dist — median: {np.median(nn_dist):.0f} m  max: {np.max(nn_dist):.0f} m")
    if np.max(nn_dist) > 2000:
        print(f"    WARNING: Some cells are > 2000 m from any observed LST. IDW quality may be poor.")

    dist = np.maximum(dist, 1e-9)
    weights = 1.0 / (dist ** p)

    known_values = known[val_col].values
    neighbor_values = known_values[idx]
    filled_values = np.sum(weights * neighbor_values, axis=1) / np.sum(weights, axis=1)

    out_series = gdf[val_col].copy()
    out_series.loc[unknown.index] = filled_values
    return out_series


def process_lst_data():
    print("--- STEP 2: LST PROCESSING ---")

    if not os.path.exists(config.BASE_GRID_GPKG):
        raise FileNotFoundError(f"Base grid not found: {config.BASE_GRID_GPKG}")
    if not os.path.exists(config.LST_MAIN_CSV_PATH):
        raise FileNotFoundError(f"LST CSV not found: {config.LST_MAIN_CSV_PATH}")

    print("Loading base grid and LST data...")
    try:
        grid = gpd.read_file(config.BASE_GRID_GPKG, layer="grid_250m")
    except Exception:
        grid = gpd.read_file(config.BASE_GRID_GPKG)

    if grid.crs is None:
        grid = grid.set_crs(config.PROJ_CRS)
    else:
        grid = grid.to_crs(config.PROJ_CRS)

    lst_df = pd.read_csv(config.LST_MAIN_CSV_PATH)
    print(f"  Base grid cells : {len(grid)}")
    print(f"  LST CSV rows    : {len(lst_df)}")

    if "grid_id" not in lst_df.columns:
        raise ValueError("LST CSV must contain 'grid_id' column for merging.")

    # Accurately locate target column: P75 Urban LST
    target_col = "lst_p75_mean"
    if target_col not in lst_df.columns:
        candidates = ["lst_p75_mean", "lst_p75", "LST_p75", "LST_p75_C_mean"]
        target_col = next((c for c in candidates if c in lst_df.columns), None)
        if target_col is None:
            raise ValueError(f"Could not find a P75 LST column. Available: {list(lst_df.columns)}.")
        print(f"  Target column '{target_col}' selected.")

    lst_sub = lst_df[["grid_id", target_col]].copy()
    lst_sub = lst_sub.rename(columns={target_col: "lst_p75_observed"})
    lst_sub["lst_p75_observed"] = pd.to_numeric(lst_sub["lst_p75_observed"], errors="coerce")
    
    print("Merging LST data onto grid...")
    grid = grid.merge(lst_sub, on="grid_id", how="left")

    missing_count = grid["lst_p75_observed"].isna().sum()
    print(f"  Cells with valid LST  : {len(grid) - missing_count}")
    print(f"  Cells missing LST     : {missing_count}")

    print("Performing IDW interpolation for missing cells...")
    grid["lst_p75_filled"] = _idw_interpolation(grid, "lst_p75_observed", k=8, p=2.0)

    grid["lst_valid"] = grid["lst_p75_observed"].notna().astype(int)
    grid["lst_was_filled"] = grid["lst_p75_observed"].isna().astype(int)
    n_interp = int(grid["lst_was_filled"].sum())
    print(f"  Cells interpolated (IDW): {n_interp} ({n_interp / len(grid) * 100:.1f}%)")

    # --- PROCESS SENSITIVITY BLOCK ---
    if hasattr(config, "LST_SENSITIVITY_CSV_PATH") and os.path.exists(config.LST_SENSITIVITY_CSV_PATH):
        print("Processing LST Sensitivity data...")
        sens_df = pd.read_csv(config.LST_SENSITIVITY_CSV_PATH)
        sens_candidates = [target_col, "lst_p75_mean", "lst_p75", "LST_p75", "LST_p75_C_mean"]
        sens_col = next((c for c in sens_candidates if c in sens_df.columns), None)
        
        if "grid_id" in sens_df.columns and sens_col is not None:
            sens_sub = sens_df[["grid_id", sens_col]].rename(columns={sens_col: "lst_p75_sens_observed"})
            grid = grid.merge(sens_sub, on="grid_id", how="left")
            grid["lst_p75_sens_observed"] = pd.to_numeric(grid["lst_p75_sens_observed"], errors="coerce")
            
            grid["lst_p75_sens_filled"] = _idw_interpolation(grid, "lst_p75_sens_observed", k=8, p=2.0)
        else:
            print("  Sensitivity processing SKIPPED due to column mismatch.")

    print(f"Saving to {config.LST_GRID_GPKG}...")
    if "fid" in grid.columns:
        grid = grid.drop(columns=["fid"])
    if grid.index.name == "fid":
        grid.index.name = None
        
    grid.to_file(config.LST_GRID_GPKG, driver="GPKG", layer="grid_250m_lst")
    print("SUCCESS: LST processing complete.")

    return grid


if __name__ == "__main__":
    process_lst_data()