"""
Step 8: Merge Sentinel-2 land-cover classification to population grid.

Input:
  config.POP_GRID_GPKG
    hanoi_grid_250m_with_pop.gpkg from Step 5

  config.LANDCOVER_CSV
    GEE-exported CSV with per-cell land-cover class percentages
    (water_pct, tree_pct, grass_pct, crop_pct, built_pct, bare_pct)

Output:
  config.LANDCOVER_GRID_GPKG
    hanoi_grid_250m_with_landcover.gpkg

Derived variables created:
  green_pct               = tree_pct + grass_pct
  urban_green_pct         = tree_pct + grass_pct
  vegetation_pct          = tree_pct + grass_pct + crop_pct
  impervious_disturbed_pct= built_pct + bare_pct
  physical_greenblue_pct  = tree_pct + grass_pct + water_pct

Notes:
  - The CSV must contain .geo column (GeoJSON geometry) for spatial nearest match.
  - LC class columns are expected as 0-100 percentages (not 0-1 fractions).
  - Spatial nearest merge with max_distance=50 m ensures robust alignment.
"""


import os
import sys
import json

import pandas as pd
import geopandas as gpd
from shapely.geometry import shape

import config

sys.stdout.reconfigure(encoding="utf-8")


def _read_population_grid(path: str) -> gpd.GeoDataFrame:
    try:
        grid = gpd.read_file(path, layer="grid_population")
    except Exception:
        grid = gpd.read_file(path)

    if grid.crs is None:
        grid = grid.set_crs(config.PROJ_CRS)
    else:
        grid = grid.to_crs(config.PROJ_CRS)

    if "grid_id" not in grid.columns:
        raise ValueError("grid_id column missing from input population grid.")

    grid["grid_id"] = pd.to_numeric(grid["grid_id"], errors="raise").astype(int)

    if grid["grid_id"].duplicated().any():
        dup_n = int(grid["grid_id"].duplicated().sum())
        raise ValueError(f"Duplicated grid_id in population grid: {dup_n}")

    if "commune_name" not in grid.columns and "commune_na" in grid.columns:
        grid = grid.rename(columns={"commune_na": "commune_name"})

    if "urban_class" not in grid.columns and "urban_clas" in grid.columns:
        grid = grid.rename(columns={"urban_clas": "urban_class"})

    return grid


def _read_landcover_csv(path: str) -> pd.DataFrame:
    lc = pd.read_csv(path)

    if "system:index" in lc.columns:
        lc = lc.drop(columns=["system:index"])
        print("  Dropped GEE column: system:index")

    if ".geo" in lc.columns:
        def parse_geo(x):
            try:
                if pd.isna(x):
                    return None
                return shape(json.loads(x))
            except (ValueError, TypeError, KeyError):
                return None

        print("  Parsing GeoJSON geometries from CSV...")
        lc["geometry"] = lc[".geo"].apply(parse_geo)
        lc = gpd.GeoDataFrame(lc, geometry="geometry", crs="EPSG:4326")
        lc = lc.drop(columns=[".geo"])

    required_base_cols = [
        "grid_id",
        "water_pct",
        "tree_pct",
        "grass_pct",
        "crop_pct",
        "built_pct",
        "bare_pct",
        "lc_valid_pct",
    ]

    missing = [c for c in required_base_cols if c not in lc.columns]
    if missing:
        raise ValueError(f"Missing required land-cover columns: {missing}")

    lc["grid_id"] = pd.to_numeric(lc["grid_id"], errors="raise").astype(int)

    for col in required_base_cols:
        if col != "grid_id":
            lc[col] = pd.to_numeric(lc[col], errors="coerce")


    pct_class_cols = ["water_pct", "tree_pct", "grass_pct", "crop_pct", "built_pct", "bare_pct"]
    
    class_max = lc[[c for c in pct_class_cols if c in lc.columns]].max(skipna=True)
    overall_max = class_max.max(skipna=True)

    if pd.notna(overall_max) and overall_max <= 1.5:
        raise ValueError(
            "Land-cover class columns appear to be 0-1 fractions, but Step 8 expects 0-100 percentages. "
            "Fix the GEE export or multiply LC percentage columns by 100 before merging."
        )

    for col in pct_class_cols:
        if col in lc.columns:
            v = lc[col].dropna()
            if len(v) > 0 and (v.min() < -0.01 or v.max() > 100.01):
                raise ValueError(
                    f"{col} is outside expected 0-100 percentage range: "
                    f"min={v.min():.3f}, max={v.max():.3f}"
                )

    if lc["grid_id"].duplicated().any():
        dup_n = int(lc["grid_id"].duplicated().sum())
        raise ValueError(f"Duplicated grid_id in land-cover CSV: {dup_n}")

    return lc


def _create_landcover_variables(lc: pd.DataFrame) -> pd.DataFrame:
    pct_cols = ["water_pct", "tree_pct", "grass_pct", "crop_pct", "built_pct", "bare_pct"]

    lc = lc.copy()

    lc["lc_total_pct"] = lc[pct_cols].sum(axis=1)
    lc["lc_total_error_pct"] = lc["lc_total_pct"] - 100.0
    lc["lc_quality_pass"] = (lc["lc_valid_pct"] >= 80).astype(int)

    # Core derived variables.
    lc["green_pct"] = lc["tree_pct"].add(lc["grass_pct"], fill_value=0).clip(0, 100)
    lc["urban_green_pct"] = lc["green_pct"]
    lc["vegetation_pct"] = lc["green_pct"].add(lc["crop_pct"], fill_value=0).clip(0, 100)
    lc["impervious_disturbed_pct"] = lc["built_pct"].add(lc["bare_pct"], fill_value=0).clip(0, 100)

    # Crop is intentionally excluded from urban green-blue accessibility.
    lc["physical_greenblue_pct"] = lc["green_pct"].add(lc["water_pct"], fill_value=0).clip(0, 100)

    return lc


def merge_landcover_to_grid():
    print("--- STEP 8: MERGE SENTINEL-2 LAND COVER TO GRID ---")

    lc_csv = getattr(
        config,
        "LANDCOVER_CSV",
        os.path.join(config.DATA_DIR, "hanoi_lc_maroct2024_6class_rededge_grid250m.csv"),
    )

    grid_in = config.POP_GRID_GPKG

    grid_out_path = getattr(
        config,
        "LANDCOVER_GRID_GPKG",
        os.path.join(config.DATA_DIR, "hanoi_grid_250m_with_landcover.gpkg"),
    )

    if not os.path.exists(lc_csv):
        raise FileNotFoundError(f"Land-cover CSV not found: {lc_csv}")

    if not os.path.exists(grid_in):
        raise FileNotFoundError(f"Population grid not found: {grid_in}")

    print("Loading data...")
    lc = _read_landcover_csv(lc_csv)
    grid = _read_population_grid(grid_in)

    print(f"  LC CSV rows : {len(lc)}")
    print(f"  Grid rows   : {len(grid)}")
    print(f"  Grid CRS    : {grid.crs}")

    print("\nCreating land-cover derived variables...")
    lc = _create_landcover_variables(lc)

    bad_total = ((lc["lc_total_pct"] < 90) | (lc["lc_total_pct"] > 110)).sum()
    low_quality = (lc["lc_quality_pass"] == 0).sum()

    print(f"  Cells with lc_valid_pct < 80%: {low_quality}")
    if bad_total > 0:
        print(f"  WARNING: {bad_total} cells have LC total outside 90-110%")

    print("  LC total pct summary:")
    print(lc["lc_total_pct"].describe().to_string())

    print("\nDerived variables:")
    print("  green_pct                 = tree_pct + grass_pct")
    print("  urban_green_pct           = tree_pct + grass_pct")
    print("  vegetation_pct            = tree_pct + grass_pct + crop_pct")
    print("  impervious_disturbed_pct  = built_pct + bare_pct")
    print("  physical_greenblue_pct    = tree_pct + grass_pct + water_pct")
    print("  crop_pct is kept separate and excluded from urban green-blue accessibility")

    lc_output_cols = [
        "water_pct",
        "tree_pct",
        "grass_pct",
        "crop_pct",
        "built_pct",
        "bare_pct",
        "green_pct",
        "physical_greenblue_pct",
        "lc_valid_pct",
        "lc_quality_pass",
        "lc_total_pct",
        "lc_total_error_pct",
        "urban_green_pct",
        "vegetation_pct",
        "impervious_disturbed_pct",
    ]

    drop_existing = [c for c in lc_output_cols if c in grid.columns]
    if drop_existing:
        grid = grid.drop(columns=drop_existing)
        print(f"\nDropped existing LC columns from grid: {drop_existing}")

    print("\nMerging land-cover data spatially using nearest LC centroid...")

    if isinstance(lc, gpd.GeoDataFrame):
        lc_geom = lc.to_crs(grid.crs)

        lc_points = lc_geom.copy()
        lc_points["geometry"] = lc_points.geometry.centroid
        if "grid_id" in lc_points.columns:
            lc_points = lc_points.drop(columns=["grid_id"])

        out = gpd.sjoin_nearest(
            grid,
            lc_points,
            how="left",
            max_distance=50,
            distance_col="lc_join_dist_m"
        )

        if "index_right" in out.columns:
            out = out.drop(columns=["index_right"])

        # Drop extra columns from lc_points that are not needed
        extra_cols = [c for c in out.columns if c not in grid.columns and c not in lc_output_cols and c != "lc_join_dist_m"]
        if extra_cols:
            out = out.drop(columns=extra_cols)

        out = (
            out.sort_values("lc_join_dist_m", na_position="last")
            .drop_duplicates("grid_id", keep="first")
        )

        print("\nSpatial Nearest Match Diagnostics:")
        v = out["lc_join_dist_m"].dropna()
        if not v.empty:
            print(f"  Nearest distance mean: {v.mean():.2f} m")
            print(f"  Nearest distance p95 : {v.quantile(0.95):.2f} m")
            print(f"  Nearest distance max : {v.max():.2f} m")
            print(f"  Cells > 50 m         : {(v > 50).sum()}")
            print(f"  Cells > 100 m        : {(v > 100).sum()}")
        else:
            print("  WARNING: No cells matched within 125 m.")

    else:
        raise ValueError(
            "Land-cover CSV has no geometry column. Spatial nearest merge requires .geo from GEE export."
        )

    n_matched = int(out["tree_pct"].notna().sum())
    n_missing = int(out["tree_pct"].isna().sum())

    print(f"  Matched LC cells : {n_matched} / {len(out)} ({n_matched / len(out) * 100:.1f}%)")
    print(f"  Missing LC cells : {n_missing}")

    if n_missing > 0:
        print("  WARNING: Some grid cells did not match any LC centroid within 125 m.")

    # Derive physical area variables from percentage x cell area
    # These save downstream steps from recomputing pct * area_m2 repeatedly.
    if "area_m2" in out.columns:
        out["physical_gb_area_m2"] = (
            out["physical_greenblue_pct"] / 100.0 * out["area_m2"]
        ).round(2)
        out["urban_green_area_m2"] = (
            out["urban_green_pct"] / 100.0 * out["area_m2"]
        ).round(2)
    else:
        print("  WARNING: 'area_m2' not found in grid — skipping physical area derivation.")

    print("\nScale check for Sentinel-2 LC classes, expected 0-100:")

    for col in ["water_pct", "tree_pct", "grass_pct", "crop_pct", "built_pct", "bare_pct"]:
        v = out[col].dropna()
        if len(v) == 0:
            print(f"  {col:30s}: all missing")
        else:
            print(f"  {col:30s}: min={v.min():.2f}, max={v.max():.2f}, mean={v.mean():.2f}")

    if "jrc_water_pct" in out.columns:
        v = out["jrc_water_pct"].dropna()
        print("\nJRC water variable also present:")
        print(f"  {'jrc_water_pct':30s}: min={v.min():.2f}, max={v.max():.2f}, mean={v.mean():.2f}")
        print("  Note: water_pct = Sentinel-2 RF water; jrc_water_pct = JRC permanent water.")

    print("\nKey variable summary:")

    key_vars = [
        "urban_green_pct",
        "physical_greenblue_pct",
        "vegetation_pct",
        "impervious_disturbed_pct",
        "built_pct",
        "crop_pct",
        "water_pct",
    ]

    for col in key_vars:
        v = out[col].dropna()
        print(f"  {col:30s}: mean={v.mean():.2f}%, median={v.median():.2f}%")

    print("\nLC quality summary:")
    print(out["lc_quality_pass"].value_counts(dropna=False).sort_index().to_string())
    print("  (1 = lc_valid_pct >= 80%, 0 = poor quality / cloud)")

    if "physical_gb_area_m2" in out.columns:
        print(f"  Mean physical GB area : {out['physical_gb_area_m2'].mean():,.1f} m\u00b2/cell")
        print(f"  Mean urban green area : {out['urban_green_area_m2'].mean():,.1f} m\u00b2/cell")

    os.makedirs(os.path.dirname(grid_out_path) or ".", exist_ok=True)

    print(f"\nSaving -> {grid_out_path}")
    out.to_file(grid_out_path, driver="GPKG", layer="grid_with_lc")

    print("\nDone.")
    print(f"  Output rows   : {len(out)}")
    print(f"  Output columns: {len(out.columns)}")
    print("  Next step: run step9_accessible_spaces.py")

    return out


if __name__ == "__main__":
    merge_landcover_to_grid()