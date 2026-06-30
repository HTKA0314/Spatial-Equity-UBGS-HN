"""
Step 4: JRC Permanent Water Metric Extraction.

Input:
  config.LST_GRID_GPKG
    hanoi_grid_250m_with_lst.gpkg from Step 2/3
  config.WATER_PATH
    JRC permanent water bodies polygon vector layer

Output:
  config.WATER_GRID_GPKG
    hanoi_grid_250m_with_water.gpkg
"""

import os
import sys
import warnings
import numpy as np
import pandas as pd
import geopandas as gpd
import config

# Ensure utf-8 encoding for Windows terminal output
sys.stdout.reconfigure(encoding="utf-8")


def calculate_water_variables():
    print("--- STEP 4: EXTRACT JRC WATER VARIABLES ---")

    if not os.path.exists(config.LST_GRID_GPKG):
        raise FileNotFoundError(f"Input LST grid missing: {config.LST_GRID_GPKG}")
        
    if not os.path.exists(config.WATER_PATH):
        raise FileNotFoundError(f"JRC Water path missing: {config.WATER_PATH}")

    print("Loading datasets...")
    grid = gpd.read_file(config.LST_GRID_GPKG, layer="grid_250m_lst")
    water = gpd.read_file(config.WATER_PATH)

    if grid.crs is None:
        grid = grid.set_crs(config.PROJ_CRS)
    else:
        grid = grid.to_crs(config.PROJ_CRS)

    if water.crs is None:
        raise ValueError("Water layer has no CRS. Please define CRS before running.")
    else:
        water = water.to_crs(config.PROJ_CRS)

    if "grid_id" not in grid.columns:
        raise ValueError("grid_id column not found in grid.")

    if "area_m2" not in grid.columns:
        grid["area_m2"] = grid.geometry.area

    grid["grid_id"] = grid["grid_id"].astype(int)

    # ------------------------------------------------------------------ #
    # 1. Filter JRC water bodies
    # ------------------------------------------------------------------ #
    water = water[water.geometry.notna() & (~water.geometry.is_empty)].copy()
    water["geometry"] = water.geometry.make_valid()

    water = gpd.GeoDataFrame(water).explode(index_parts=False).reset_index(drop=True)
    water = water[water.geometry.geom_type.isin(["Polygon", "MultiPolygon"])].copy()
    water["area_m2"] = water.geometry.area

    min_water_area_m2 = getattr(config, "MIN_WATER_AREA_M2", 1000)
    water_filtered = water[water["area_m2"] >= min_water_area_m2].copy()

    if water_filtered.empty:
        raise ValueError(f"No water polygons remain after filtering area >= {min_water_area_m2} m2.")

    print(f"Valid JRC water bodies >= {min_water_area_m2} m2: {len(water_filtered)}")

    # ------------------------------------------------------------------ #
    # 2. Distance to nearest JRC water polygon (Nearest Neighbor)
    # ------------------------------------------------------------------ #
    print("Calculating nearest distance to JRC water bodies...")
    grid_points = gpd.GeoDataFrame(
        grid[["grid_id"]].copy(),
        geometry=grid.geometry.representative_point(),
        crs=config.PROJ_CRS,
    )

    nearest = gpd.sjoin_nearest(
        grid_points,
        water_filtered[["geometry", "area_m2"]],
        how="left",
        distance_col="dist_to_jrc_water_m",
    )

    # Sort and drop duplicates to preserve matrix stability
    nearest = (
        nearest
        .sort_values(["grid_id", "dist_to_jrc_water_m", "area_m2"], ascending=[True, True, False])
        .drop_duplicates("grid_id")
    )

    nearest = nearest.rename(columns={"area_m2": "nearest_water_area_m2"})

    grid = grid.merge(
        nearest[["grid_id", "dist_to_jrc_water_m", "nearest_water_area_m2"]],
        on="grid_id",
        how="left",
    )

    grid["dist_to_jrc_water_m"] = pd.to_numeric(grid["dist_to_jrc_water_m"], errors="coerce")
    grid["nearest_water_area_m2"] = grid["nearest_water_area_m2"].fillna(0.0)
    grid["log_dist_jrc_water"] = np.log1p(grid["dist_to_jrc_water_m"])

    # ------------------------------------------------------------------ #
    # 3. Calculate local JRC water coverage (Cell-level micro-intersection)
    # ------------------------------------------------------------------ #
    print("Calculating local JRC water coverage share per grid cell...")
    grid["jrc_water_area_m2"] = 0.0

    try:
        joined = gpd.sjoin(
            grid[["grid_id", "geometry"]].reset_index(drop=True),
            water_filtered[["geometry"]].reset_index(drop=True),
            how="inner",
            predicate="intersects",
        )

        if not joined.empty:
            water_geoms_aligned = (
                water_filtered
                .reset_index(drop=True)
                .geometry
                .loc[joined["index_right"]]
            )
            aligned_water = gpd.GeoSeries(water_geoms_aligned.values, index=joined.index, crs=joined.crs)
            joined["inter_geom"] = joined.geometry.intersection(aligned_water)
            joined = joined[~joined["inter_geom"].is_empty & joined["inter_geom"].notna()].copy()

            inter_gdf = gpd.GeoDataFrame(
                joined[["grid_id"]], geometry=joined["inter_geom"], crs=grid.crs
            )
            water_sums = (
                inter_gdf.dissolve(by="grid_id")
                .assign(jrc_water_area_m2=lambda x: x.geometry.area)
                .reset_index()[["grid_id", "jrc_water_area_m2"]]
            )
            grid = grid.drop(columns=["jrc_water_area_m2"], errors="ignore").merge(
                water_sums, on="grid_id", how="left"
            )
            grid["jrc_water_area_m2"] = grid["jrc_water_area_m2"].fillna(0.0)

    except (AttributeError, ValueError, KeyError) as e:
        print(f"Warning: water coverage calculation failed ({type(e).__name__}: {e}). Fallback to 0.0.")
        grid["jrc_water_area_m2"] = 0.0

    grid["jrc_water_frac"] = (grid["jrc_water_area_m2"] / grid["area_m2"]).fillna(0.0).clip(0, 1)
    grid["jrc_water_pct"] = grid["jrc_water_frac"] * 100

    # ------------------------------------------------------------------ #
    # 4. Diagnostics & Validation
    # ------------------------------------------------------------------ #
    print("\nValidation Statistics:")
    print(f"  Avg distance to JRC water   : {grid['dist_to_jrc_water_m'].dropna().mean():.1f} m")
    print(f"  Median distance to JRC water: {grid['dist_to_jrc_water_m'].dropna().median():.1f} m")
    print(f"  Cells with JRC water > 0    : {(grid['jrc_water_pct'] > 0).sum()}")
    print(f"  Mean JRC water pct          : {grid['jrc_water_pct'].mean():.2f}%")
    print(f"  Max JRC water pct           : {grid['jrc_water_pct'].max():.2f}%")

    os.makedirs(os.path.dirname(config.WATER_GRID_GPKG) or ".", exist_ok=True)
    grid.to_file(config.WATER_GRID_GPKG, driver="GPKG", layer="grid_water")
    print(f"\nSaved -> {config.WATER_GRID_GPKG}")

    return grid


if __name__ == "__main__":
    calculate_water_variables()