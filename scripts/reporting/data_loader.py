import os
import sys
import numpy as np
import pandas as pd
import geopandas as gpd
from typing import Tuple, Optional

# Add the parent directory to sys.path so we can import config if needed
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))

from scripts.reporting import constants as const

def load_clean_data() -> gpd.GeoDataFrame:
    """
    Robustly loads equity_metrics and walking_access, merges them if needed,
    filters valid cells, and computes standardized metric columns.
    Returns a single, clean GeoDataFrame ready for tables and plots.
    """
    if not os.path.exists(const.EQUITY_GPKG):
        raise FileNotFoundError(f"Missing {const.EQUITY_GPKG}")
        
    try:
        gdf = gpd.read_file(const.EQUITY_GPKG, layer="grid_equity_metrics")
    except Exception:
        gdf = gpd.read_file(const.EQUITY_GPKG)
        
    if "grid_id" not in gdf.columns:
        raise ValueError("Missing 'grid_id' in equity metrics.")
    gdf["grid_id"] = gdf["grid_id"].astype(int)

    # Attempt to merge walking access data to guarantee we have all columns
    if os.path.exists(const.ACCESS_GPKG):
        try:
            try:
                access = gpd.read_file(const.ACCESS_GPKG, layer="grid_walking_access")
            except Exception:
                access = gpd.read_file(const.ACCESS_GPKG)
                
            access["grid_id"] = access["grid_id"].astype(int)
            
            # Select columns from access that are NOT already in gdf
            cols_to_merge = [c for c in access.columns if c not in gdf.columns and c != "geometry"]
            if cols_to_merge:
                cols_to_merge.append("grid_id")
                gdf = gdf.merge(access[cols_to_merge], on="grid_id", how="left")
        except Exception as e:
            print(f"Warning: Failed to merge walking access GPKG: {e}")

    # 1. Base Geometry & Filters
    if "area_m2" not in gdf.columns:
        gdf["area_m2"] = gdf.geometry.area
        
    if "edge_cell" in gdf.columns:
        gdf = gdf[gdf["edge_cell"] == 0].copy()
    
    # Exclude non-full grid cells (assuming 250m grid -> 62500m2, cutoff at 31250)
    gdf = gdf[pd.to_numeric(gdf["area_m2"], errors="coerce") >= 31250].copy()
    
    if "population" not in gdf.columns:
        raise ValueError("Missing 'population' column.")
    gdf["population"] = pd.to_numeric(gdf["population"], errors="coerce").fillna(0)
    
    # 2. Derive Standard Columns
    
    # Density
    density_col = "density" if "density" in gdf.columns else "pop_density_km2"
    if density_col in gdf.columns:
        gdf["density_val"] = pd.to_numeric(gdf[density_col], errors="coerce")
    else:
        gdf["density_val"] = gdf["population"] / (gdf["area_m2"] / 1_000_000)
        
    def _density_group(v: float) -> str:
        if pd.isna(v): return np.nan
        if v < 5000: return "Low-density"
        if v <= 15000: return "Medium-density"
        return "High-density"
        
    gdf["density_group"] = pd.Categorical(
        gdf["density_val"].apply(_density_group),
        categories=const.DENSITY_ORDER,
        ordered=True
    )
    
    # LST (Prioritize P75 filled)
    lst_candidates = ["lst_p75_filled", "lst_p75_observed", "lst_p75_mean", "LST_p75_C_mean"]
    lst_col = next((c for c in lst_candidates if c in gdf.columns), None)
    if lst_col:
        gdf["lst_temp"] = pd.to_numeric(gdf[lst_col], errors="coerce")
    else:
        gdf["lst_temp"] = np.nan
        print("Warning: LST column not found.")

    # Physical Green-Blue Cover
    if "physical_greenblue_pct" in gdf.columns:
        s = pd.to_numeric(gdf["physical_greenblue_pct"], errors="coerce")
        if s.max() <= 1.5: s = s * 100
        gdf["physical_gb_pct"] = s
    elif {"tree_pct", "grass_pct", "water_pct"}.issubset(gdf.columns):
        s = pd.to_numeric(gdf["tree_pct"], errors="coerce") + \
            pd.to_numeric(gdf["grass_pct"], errors="coerce") + \
            pd.to_numeric(gdf["water_pct"], errors="coerce")
        if s.max() <= 1.5: s = s * 100
        gdf["physical_gb_pct"] = s.clip(0, 100)
    else:
        gdf["physical_gb_pct"] = np.nan
        print("Warning: Physical green-blue pct not found.")

    # Accessible Level A Cover
    acc_candidates = ["accessible_gb_area_A_final_m2", "accessible_gb_area_A_pixel_m2", "accessible_area_A", "accessible_gb_A_m2", "accessible_gb_area_A_est_m2"]
    acc_m2_col = next((c for c in acc_candidates if c in gdf.columns), None)
    if acc_m2_col:
        gdf["accessible_gb_A_m2"] = pd.to_numeric(gdf[acc_m2_col], errors="coerce").fillna(0)
        gdf["accessible_A_pct"] = (gdf["accessible_gb_A_m2"] / gdf["area_m2"]) * 100
    else:
        gdf["accessible_A_pct"] = np.nan
        print("Warning: Accessible Level A area not found.")

    # Accessible Level A+B Cover
    acc_ab_candidates = ["accessible_gb_area_AB_final_m2", "accessible_gb_area_AB_pixel_m2", "accessible_area_AB", "accessible_gb_AB_m2", "accessible_gb_area_AB_est_m2"]
    acc_ab_m2_col = next((c for c in acc_ab_candidates if c in gdf.columns), None)
    if acc_ab_m2_col:
        gdf["accessible_gb_AB_m2"] = pd.to_numeric(gdf[acc_ab_m2_col], errors="coerce").fillna(0)
    else:
        gdf["accessible_gb_AB_m2"] = np.nan

    # Physical-Access Gap
    # Equation: Gap = Physical Pct - Accessible Pct
    # Note: Using absolute gap, lower bounded at 0
    gdf["physical_access_gap_pct"] = (gdf["physical_gb_pct"] - gdf["accessible_A_pct"]).clip(lower=0)

    # Conversion Ratio
    # Equation: Accessible Area / Physical Area * 100
    gdf["physical_gb_m2"] = (gdf["physical_gb_pct"] / 100.0) * gdf["area_m2"]
    gdf["conversion_ratio"] = np.where(
        gdf["physical_gb_m2"] > 0,
        (gdf["accessible_gb_A_m2"] / gdf["physical_gb_m2"]) * 100,
        np.nan
    )

    # Network Distance
    dist_candidates = ["dist_to_accessible_A_m", "dist_accessible_A_m", "nearest_accessible_A_m"]
    dist_col = next((c for c in dist_candidates if c in gdf.columns), None)
    if dist_col:
        dist = pd.to_numeric(gdf[dist_col], errors="coerce")
        gdf["dist_to_A_m"] = dist
        
        # Categorize into exact 4 bins
        cls = pd.Series(const.ACCESS_CLASSES[-1], index=gdf.index, dtype="object")
        cls.loc[dist <= 300] = const.ACCESS_CLASSES[0]
        cls.loc[(dist > 300) & (dist <= 500)] = const.ACCESS_CLASSES[1]
        cls.loc[(dist > 500) & (dist <= 1000)] = const.ACCESS_CLASSES[2]
        cls.loc[dist.isna() | (dist > 1000)] = const.ACCESS_CLASSES[3]
        
        gdf["network_access_class"] = pd.Categorical(cls, categories=const.ACCESS_CLASSES, ordered=True)
    else:
        gdf["dist_to_A_m"] = np.nan
        gdf["network_access_class"] = np.nan
        print("Warning: Network distance to Level A not found.")

    # Drop intermediate or non-standard cols to keep memory clean if necessary, 
    # but keeping them is fine for downstream scripts.
    return gdf
