"""
Step 5: Population Mapping and Commune-level Calibration.

Input:
  config.WATER_GRID_GPKG
    hanoi_grid_250m_with_water.gpkg from Step 4
  config.POP_RASTER_PATH
    WorldPop unconstrained population raster reprojected to config.PROJ_CRS

Output:
  config.POP_GRID_GPKG
    hanoi_grid_250m_with_population.gpkg
"""

import os
import sys
import warnings

import numpy as np
import pandas as pd
import geopandas as gpd
import rasterio
from rasterstats import zonal_stats
import config

# Ensure utf-8 encoding for Windows terminal output
sys.stdout.reconfigure(encoding="utf-8")


def _parse_vn_number(series: pd.Series) -> pd.Series:
    if pd.api.types.is_numeric_dtype(series):
        return series.astype(float)

    s = series.astype(str).str.strip().str.replace(" ", "", regex=False)

    def _parse_one(x: str) -> float:
        if x in ("", "nan", "None", "NULL", "NaN"):
            return float("nan")

        has_comma = "," in x
        has_dot   = "." in x

        if has_comma and has_dot:
            if x.rfind(",") > x.rfind("."):
                x = x.replace(".", "").replace(",", ".")
            else:
                x = x.replace(",", "")
        elif has_comma:
            parts = x.split(",")
            if len(parts[-1]) == 3 and len(parts) > 1:
                x = x.replace(",", "")
            else:
                x = x.replace(",", ".")
        elif has_dot:
            parts = x.split(".")
            if len(parts) > 2:
                x = x.replace(".", "")
            elif len(parts) == 2 and len(parts[-1]) == 3 and len(parts[0]) <= 3:
                x = x.replace(".", "")

        try:
            return float(x)
        except ValueError:
            return float("nan")

    return s.apply(_parse_one)


def _standardize_name(series: pd.Series) -> pd.Series:
    return (
        series.astype(str)
        .str.strip()
        .str.lower()
        .str.replace(r"\s+", " ", regex=True)
    )


def _calibration_factors_by_commune(
    grid: gpd.GeoDataFrame,
    worldpop_total: float,
) -> pd.Series:
    census_total = getattr(config, "CENSUS_TOTAL", None)
    allow_fallback = getattr(config, "ALLOW_POP_FALLBACK", False)

    if census_total is not None:
        print(f"Census total manual override: {census_total:,.0f}")
        global_factor = census_total / worldpop_total
        return pd.Series(global_factor, index=grid.index)

    if not os.path.exists(config.ADMIN_PATH):
        if allow_fallback:
            print("WARNING: Admin file not found — using calibration factor 1.0")
            return pd.Series(1.0, index=grid.index)
        raise FileNotFoundError(f"Admin file not found: {config.ADMIN_PATH}")

    if "commune_name" not in grid.columns:
        if allow_fallback:
            print("WARNING: 'commune_name' missing from grid — using calibration factor 1.0")
            return pd.Series(1.0, index=grid.index)
        raise ValueError("'commune_name' missing from grid — cannot calibrate population by commune.")

    admin = gpd.read_file(config.ADMIN_PATH)

    name_col = getattr(config, "ADMIN_NAME_FIELD", "xaphuong")
    pop_col = getattr(config, "ADMIN_POP_FIELD", "danso")

    if name_col not in admin.columns or pop_col not in admin.columns:
        if allow_fallback:
            print("WARNING: Required admin columns missing — using calibration factor 1.0")
            return pd.Series(1.0, index=grid.index)
        raise ValueError(f"Admin file fields mismatch: name={name_col}, pop={pop_col}")

    admin = admin.copy()
    admin["_name_key"] = _standardize_name(admin[name_col])
    admin["parsed_pop"] = _parse_vn_number(admin[pop_col]).fillna(0)

    census_dict = admin.set_index("_name_key")["parsed_pop"].to_dict()

    grid_temp = grid.copy()
    grid_temp["_name_key"] = _standardize_name(grid_temp["commune_name"])

    worldpop_sums = grid_temp.groupby("_name_key")["pop_raw"].sum().to_dict()
    matched_keys = set(worldpop_sums.keys()) & set(census_dict.keys())

    if len(matched_keys) == 0:
        if allow_fallback:
            print("WARNING: No commune names matched — using calibration factor 1.0")
            return pd.Series(1.0, index=grid.index)
        raise ValueError("No commune names matched between grid and admin census data.")

    total_census_matched = sum(census_dict[k] for k in matched_keys)
    total_worldpop_matched = sum(worldpop_sums[k] for k in matched_keys)

    if total_worldpop_matched <= 0:
        if allow_fallback:
            print("WARNING: Matched WorldPop total is zero — using calibration factor 1.0")
            return pd.Series(1.0, index=grid.index)
        raise ValueError("Matched WorldPop total is zero; cannot calculate calibration factors.")

    global_factor = total_census_matched / total_worldpop_matched

    print(f"Census total matched ({len(matched_keys)} communes): {total_census_matched:,.0f}")
    print(f"WorldPop total matched: {total_worldpop_matched:,.1f}")
    print(f"Global fallback factor: {global_factor:.4f}")

    unmatched_grid = set(worldpop_sums.keys()) - set(census_dict.keys())
    if unmatched_grid:
        print(f"WARNING: {len(unmatched_grid)} commune(s) in grid not found in admin census.")
        print("  First unmatched grid commune(s):", sorted(unmatched_grid)[:10])

    commune_factors = {}
    for key, wp_sum in worldpop_sums.items():
        if key in census_dict and wp_sum > 0:
            commune_factors[key] = census_dict[key] / wp_sum
        else:
            commune_factors[key] = global_factor

    grid_factors = grid_temp["_name_key"].map(commune_factors).fillna(global_factor)
    return grid_factors


def map_population() -> gpd.GeoDataFrame:
    print("--- STEP 5: POPULATION MAPPING ---")

    for path in (config.POP_RASTER_PATH, config.WATER_GRID_GPKG):
        if not os.path.exists(path):
            raise FileNotFoundError(path)

    print("Loading grid with JRC water variables...")
    grid = gpd.read_file(config.WATER_GRID_GPKG, layer="grid_water")

    if grid.crs is None:
        grid = grid.set_crs(config.PROJ_CRS)
    else:
        grid = grid.to_crs(config.PROJ_CRS)

    if "commune_name" not in grid.columns and "commune_na" in grid.columns:
        grid = grid.rename(columns={"commune_na": "commune_name"})

    if "urban_class" not in grid.columns and "urban_clas" in grid.columns:
        grid = grid.rename(columns={"urban_clas": "urban_class"})

    required_cols = ["grid_id", "area_m2", "jrc_water_frac", "commune_name"]
    for col in required_cols:
        if col not in grid.columns:
            raise ValueError(f"Required column missing: '{col}'")

    grid["grid_id"] = grid["grid_id"].astype(int)

    with rasterio.open(config.POP_RASTER_PATH) as src:
        raster_epsg = src.crs.to_epsg() if src.crs else None
        grid_epsg = grid.crs.to_epsg() if grid.crs else None
        if raster_epsg is not None and grid_epsg is not None and raster_epsg != grid_epsg:
            raise ValueError(
                f"CRS mismatch: raster=EPSG:{raster_epsg}, grid=EPSG:{grid_epsg}.\n"
                f"Please reproject the raster to {config.PROJ_CRS} before running Step 5."
            )
        elif raster_epsg is None or grid_epsg is None:
            print("WARNING: Could not resolve EPSG codes — CRS match unverified. Proceeding.")

    print("Calculating WorldPop zonal sums...")
    stats = zonal_stats(
        grid,
        config.POP_RASTER_PATH,
        stats="sum",
        nodata=0,
        all_touched=True,
    )

    grid["pop_raw"] = [
        0.0 if (s["sum"] is None or np.isnan(s["sum"])) else float(s["sum"])
        for s in stats
    ]

    water_threshold = getattr(config, "WATER_CELL_ZERO_POP_THRESHOLD", 0.8)
    water_mask = grid["jrc_water_frac"] > water_threshold
    water_pop_removed = grid.loc[water_mask, "pop_raw"].sum()

    grid.loc[water_mask, "pop_raw"] = 0.0
    grid["pop_zeroed_by_water"] = water_mask.astype(int)

    worldpop_total = grid["pop_raw"].sum()
    if worldpop_total <= 0:
        raise ValueError("WorldPop total is zero after water masking.")

    print(f"Raw WorldPop total after water masking: {worldpop_total:,.1f}")
    print(f"Population removed from water-dominated cells: {water_pop_removed:,.1f}")

    print("Calculating commune-level calibration factors...")
    cal_factors = _calibration_factors_by_commune(grid, worldpop_total)

    grid["cal_factor"] = cal_factors.round(4)
    grid["population"] = (grid["pop_raw"] * cal_factors).round(1)
    grid["pop_density_km2"] = np.where(
        grid["area_m2"] > 0,
        grid["population"] / (grid["area_m2"] / 1e6),
        0.0,
    )

    calibrated_total = grid["population"].sum()
    print("\nPopulation mapping complete:")
    print(f"  Calibrated total         : {calibrated_total:,.0f}")
    print(f"  Max density overall    : {grid['pop_density_km2'].max():,.0f} ppl/km²")
    print(f"  Zero-pop grids         : {(grid['population'] == 0).sum()}")

    grid.to_file(config.POP_GRID_GPKG, driver="GPKG", layer="grid_population")
    print(f"\nSaved -> {config.POP_GRID_GPKG}")

    return grid


if __name__ == "__main__":
    map_population()