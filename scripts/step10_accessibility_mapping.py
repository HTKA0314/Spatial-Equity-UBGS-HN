"""
Step 10: Physical-to-accessible green-blue gap analysis.

Input:
  config.LANDCOVER_GRID_GPKG
    hanoi_grid_250m_with_landcover.gpkg

  config.ACCESSIBLE_SPACES_GPKG
    accessible_spaces.gpkg from Step 9

Optional input:
  config.LANDCOVER_RASTER
    categorical land-cover raster aligned/projected to config.PROJ_CRS.
    Used for pixel-based accessible green-blue area if available.

Output:
  config.GAP_GRID_GPKG
    hanoi_grid_250m_with_gap.gpkg

Main output groups:
  1. Physical green-blue area
     - physical_gb_area_m2

  2. Accessible polygon area
     - accessible_polygon_area_A_m2
     - accessible_polygon_area_AB_m2

  3. Method A: proportional estimate
     - accessible_gb_area_A_est_m2
     - accessible_gb_area_AB_est_m2

  4. Method B: pixel-based estimate, if LANDCOVER_RASTER exists
     - accessible_gb_area_A_pixel_m2
     - accessible_gb_area_AB_pixel_m2

  5. Final accessible green-blue area
     - accessible_gb_area_A_final_m2
     - accessible_gb_area_AB_final_m2
     - accessible_gb_method

  6. Gap and conversion metrics
     - absolute_gap_A_m2
     - absolute_gap_AB_m2
     - relative_gap_A
     - relative_gap_AB
     - conversion_ratio_A
     - conversion_ratio_AB

Definitions:
  Level A  = conservative public accessible green-blue spaces.
  Level AB = Level A + semi-public / de facto accessible spaces.

Important:
  Method A is a linear proportional estimate:
    accessible_gb_area = accessible_polygon_intersection_area x physical_gb_fraction

  It assumes green-blue pixels are uniformly distributed within each 250 m cell.
  This is suitable for relative comparison but not exact square-metre measurement.

  Method B uses raster-vector intersection when a categorical land-cover raster is
  available. It first dissolves accessible polygons to prevent double-counting
  overlapping OSM features, then counts green-blue pixels inside the accessible
  intersection area.

  If Method B is available, final metrics use Method B.
  Otherwise, final metrics fall back to Method A.
"""

import os
import sys
import warnings

import numpy as np
import pandas as pd
import geopandas as gpd

import config

warnings.filterwarnings("ignore", category=UserWarning)
sys.stdout.reconfigure(encoding="utf-8")


# ---------------------------------------------------------------------
# BASIC HELPERS
# ---------------------------------------------------------------------

def _read_grid() -> gpd.GeoDataFrame:
    if not os.path.exists(config.LANDCOVER_GRID_GPKG):
        raise FileNotFoundError(
            f"LANDCOVER_GRID_GPKG not found: {config.LANDCOVER_GRID_GPKG}. "
            "Run Step 8 first."
        )

    try:
        grid = gpd.read_file(config.LANDCOVER_GRID_GPKG, layer="grid_with_lc")
    except Exception:
        grid = gpd.read_file(config.LANDCOVER_GRID_GPKG)

    if grid.empty:
        raise ValueError("Land-cover grid is empty.")

    if grid.crs is None:
        grid = grid.set_crs(config.PROJ_CRS)
    else:
        grid = grid.to_crs(config.PROJ_CRS)

    if "grid_id" not in grid.columns:
        raise ValueError("grid_id missing from land-cover grid.")

    grid["grid_id"] = pd.to_numeric(grid["grid_id"], errors="raise").astype(int)

    if grid["grid_id"].duplicated().any():
        dup_n = int(grid["grid_id"].duplicated().sum())
        raise ValueError(f"Duplicated grid_id in land-cover grid: {dup_n}")

    grid = grid[grid.geometry.notna() & (~grid.geometry.is_empty)].copy()
    grid["geometry"] = grid.geometry.make_valid()

    return grid


def _read_access_layer(layer_name: str, target_crs) -> gpd.GeoDataFrame:
    if not os.path.exists(config.ACCESSIBLE_SPACES_GPKG):
        raise FileNotFoundError(
            f"ACCESSIBLE_SPACES_GPKG not found: {config.ACCESSIBLE_SPACES_GPKG}. "
            "Run Step 9 first."
        )

    try:
        gdf = gpd.read_file(config.ACCESSIBLE_SPACES_GPKG, layer=layer_name)
    except Exception as e:
        raise ValueError(f"Cannot read layer '{layer_name}' from accessible spaces GPKG: {e}")

    if gdf.empty:
        return gpd.GeoDataFrame(columns=["geometry"], geometry="geometry", crs=target_crs)

    if gdf.crs is None:
        gdf = gdf.set_crs(config.PROJ_CRS)
    else:
        gdf = gdf.to_crs(target_crs)

    gdf = gdf[gdf.geometry.notna() & (~gdf.geometry.is_empty)].copy()
    gdf["geometry"] = gdf.geometry.make_valid()

    # Keep only polygonal features.
    gdf = gdf[gdf.geometry.geom_type.isin(["Polygon", "MultiPolygon"])].copy()

    if gdf.empty:
        return gpd.GeoDataFrame(columns=["geometry"], geometry="geometry", crs=target_crs)

    return gdf


def _to_fraction(series: pd.Series) -> pd.Series:
    """
    Convert 0-100 percentage to 0-1 fraction if needed.

    In the current pipeline, Step 8 enforces Sentinel-2 LC variables as 0-100.
    This function is only a defensive fallback for legacy inputs.
    """
    s = pd.to_numeric(series, errors="coerce").astype(float)
    max_val = s.max(skipna=True)

    if pd.notna(max_val) and max_val > 1.5:
        s = s / 100.0

    return s.clip(0, 1)


def _clip_area_to_cell_area(grid: gpd.GeoDataFrame, col: str):
    grid[col] = pd.to_numeric(grid[col], errors="coerce").fillna(0.0)
    grid[col] = grid[col].clip(lower=0)
    grid[col] = np.minimum(grid[col], pd.to_numeric(grid["area_m2"], errors="coerce"))
    return grid


# ---------------------------------------------------------------------
# VECTOR INTERSECTION AREA
# ---------------------------------------------------------------------

def _overlay_area_by_grid(
    grid: gpd.GeoDataFrame,
    spaces: gpd.GeoDataFrame,
    output_col: str,
) -> pd.DataFrame:
    """
    Calculate intersection area between grid cells and accessible-space polygons.
    Optimized: Cell-level spatial join and dissolve to avoid global union bottlenecks.
    Note: cell_intersections.dissolve() internally performs geometric union, which automatically prevents double-counting of overlapping OSM polygons.
    """
    out = grid[["grid_id"]].copy()
    out[output_col] = 0.0

    if spaces.empty:
        return out

    spaces = spaces.to_crs(grid.crs)

    grid_clean = grid[["grid_id", "geometry"]].reset_index(drop=True)
    spaces_clean = spaces[["geometry"]].reset_index(drop=True)

    joined = gpd.sjoin(grid_clean, spaces_clean, how="inner", predicate="intersects")
    if joined.empty:
        return out

    aligned_spaces = gpd.GeoSeries(spaces_clean.geometry.loc[joined["index_right"]].values, index=joined.index, crs=joined.crs)
    joined["inter_geom"] = joined.geometry.intersection(aligned_spaces)
    joined = joined[~joined["inter_geom"].is_empty & joined["inter_geom"].notna()].copy()

    cell_intersections = gpd.GeoDataFrame(joined[["grid_id"]], geometry=joined["inter_geom"], crs=grid.crs)
    dissolved_cells = cell_intersections.dissolve(by="grid_id").reset_index()
    dissolved_cells[output_col] = dissolved_cells.geometry.area

    if output_col in out.columns:
        out = out.drop(columns=[output_col])

    out = out.merge(dissolved_cells[["grid_id", output_col]], on="grid_id", how="left")
    out[output_col] = out[output_col].fillna(0.0)
    return out


# ---------------------------------------------------------------------
# PIXEL-BASED METHOD
# ---------------------------------------------------------------------

def _get_landcover_raster_path() -> str | None:
    path = getattr(config, "LANDCOVER_RASTER", None)
    if path is None:
        return None
    return str(path)


def _get_greenblue_classes() -> list[int]:
    """
    Get categorical raster class IDs representing tree + grass + water.

    Must be defined in config.py to avoid silently using the wrong class mapping.

    Example config:
      LC_GREENBLUE_CLASSES = [0, 1, 2]  # water, tree, grass
      LC_CLASS_MAP = {
          0: "water",
          1: "tree",
          2: "grass",
          3: "crop",
          4: "built",
          5: "bare",
      }
    """
    classes = getattr(config, "LC_GREENBLUE_CLASSES", None)

    if classes is None:
        raise ValueError(
            "LANDCOVER_RASTER exists, but config.LC_GREENBLUE_CLASSES is not defined. "
            "Define the categorical class IDs for tree, grass and water before using Method B."
        )

    if not isinstance(classes, (list, tuple, set)) or len(classes) == 0:
        raise ValueError("config.LC_GREENBLUE_CLASSES must be a non-empty list/tuple/set.")

    return [int(c) for c in classes]


def _pixel_accessible_area_by_grid(
    grid: gpd.GeoDataFrame,
    spaces: gpd.GeoDataFrame,
    lc_raster: str,
    output_col: str,
) -> pd.DataFrame:
    """
    Count green-blue raster pixels inside accessible polygons by grid cell.

    Accessible polygons are dissolved first to avoid double-counting overlapping
    OSM features. Then the dissolved accessible geometry is intersected with the
    grid, and categorical zonal statistics are computed for each intersection.
    """
    from rasterstats import zonal_stats
    import rasterio

    out = grid[["grid_id"]].copy()
    out[output_col] = np.nan

    if spaces.empty:
        out[output_col] = 0.0
        return out

    spaces = spaces.to_crs(grid.crs)
    
    grid_clean = grid[["grid_id", "geometry"]].reset_index(drop=True)
    spaces_clean = spaces[["geometry"]].reset_index(drop=True)

    joined = gpd.sjoin(grid_clean, spaces_clean, how="inner", predicate="intersects")
    if joined.empty:
        out[output_col] = 0.0
        return out

    aligned_spaces = gpd.GeoSeries(spaces_clean.geometry.loc[joined["index_right"]].values, index=joined.index, crs=joined.crs)
    joined["inter_geom"] = joined.geometry.intersection(aligned_spaces)
    joined = joined[~joined["inter_geom"].is_empty & joined["inter_geom"].notna()].copy()
    cell_intersections = gpd.GeoDataFrame(joined[["grid_id"]], geometry=joined["inter_geom"], crs=grid.crs)
    
    ix = cell_intersections.dissolve(by="grid_id").reset_index()

    ix["inter_area_m2"] = ix.geometry.area
    ix = ix[ix["inter_area_m2"] > 0].copy()

    if ix.empty:
        out[output_col] = 0.0
        return out

    with rasterio.open(lc_raster) as src:
        pixel_area_m2 = abs(src.res[0] * src.res[1])
        nodata = src.nodata
        raster_crs = src.crs
        
        # Read a sample to print unique classes
        sample = src.read(1, masked=True)
        vals = np.unique(sample.compressed())
        # Print only on the first call (for output_col ending in A)
        if output_col.endswith("_A_pixel_m2"):
            print(f"  Raster class values: {vals[:20]}")

    if raster_crs is None:
        raise ValueError("LANDCOVER_RASTER has no CRS.")

    if not raster_crs.is_projected:
        raise ValueError(
            "LANDCOVER_RASTER must use a projected CRS with metre units for pixel-area calculation."
        )

    if str(raster_crs) != str(grid.crs):
        if output_col.endswith("_A_pixel_m2"):
            print(f"  WARNING: raster CRS {raster_crs} differs from grid CRS {grid.crs}; reprojecting vectors.")
        ix = ix.to_crs(raster_crs)
        
    if output_col.endswith("_A_pixel_m2"):
        if abs(np.sqrt(pixel_area_m2) - 10) > 0.5:
            print(f"  WARNING: LANDCOVER_RASTER resolution implies pixel area {pixel_area_m2} m2, expected ~100 m2 (10x10).")

    greenblue_classes = _get_greenblue_classes()

    stats = zonal_stats(
        ix,
        lc_raster,
        categorical=True,
        nodata=nodata,
        all_touched=False,
    )

    ix["greenblue_pixel_count"] = [
        sum((s or {}).get(cls, 0) for cls in greenblue_classes)
        for s in stats
    ]

    by_grid = (
        ix.groupby("grid_id", as_index=False)["greenblue_pixel_count"]
        .sum()
    )

    by_grid[output_col] = by_grid["greenblue_pixel_count"] * pixel_area_m2
    by_grid = by_grid[["grid_id", output_col]]

    out = out.drop(columns=[output_col]).merge(by_grid, on="grid_id", how="left")
    out[output_col] = out[output_col].fillna(0.0)

    return out


def _try_pixel_method(
    grid: gpd.GeoDataFrame,
    spaces_A: gpd.GeoDataFrame,
    spaces_AB: gpd.GeoDataFrame,
) -> tuple[gpd.GeoDataFrame, bool]:
    """
    Try Method B pixel-based accessible green-blue estimation.

    Preconditions:
      - physical_gb_area_m2 has already been calculated.
      - accessible_gb_area_A_est_m2 and accessible_gb_area_AB_est_m2 already exist.
      - config.LANDCOVER_RASTER exists.
      - config.LC_GREENBLUE_CLASSES correctly maps tree/grass/water class IDs.
    """
    lc_raster = _get_landcover_raster_path()

    grid["accessible_gb_area_A_pixel_m2"] = np.nan
    grid["accessible_gb_area_AB_pixel_m2"] = np.nan
    grid["accessible_A_method_ratio"] = np.nan
    grid["accessible_AB_method_ratio"] = np.nan

    if lc_raster is None or not os.path.exists(lc_raster):
        print("  Pixel method skipped: LANDCOVER_RASTER not found. Using proportional estimate.")
        return grid, False

    if spaces_A.empty and spaces_AB.empty:
        print("  Pixel method skipped: no accessible spaces.")
        grid["accessible_gb_area_A_pixel_m2"] = 0.0
        grid["accessible_gb_area_AB_pixel_m2"] = 0.0
        return grid, False

    try:
        print("  Attempting pixel-based accessible green-blue area (Method B)...")
        print(f"  LANDCOVER_RASTER: {lc_raster}")
        print(f"  LC_GREENBLUE_CLASSES: {_get_greenblue_classes()}")

        pix_A = _pixel_accessible_area_by_grid(
            grid=grid,
            spaces=spaces_A,
            lc_raster=lc_raster,
            output_col="accessible_gb_area_A_pixel_m2",
        )

        pix_AB = _pixel_accessible_area_by_grid(
            grid=grid,
            spaces=spaces_AB,
            lc_raster=lc_raster,
            output_col="accessible_gb_area_AB_pixel_m2",
        )

        grid = grid.drop(
            columns=[
                "accessible_gb_area_A_pixel_m2",
                "accessible_gb_area_AB_pixel_m2",
            ],
            errors="ignore",
        )

        grid = grid.merge(pix_A, on="grid_id", how="left")
        grid = grid.merge(pix_AB, on="grid_id", how="left")

        grid["accessible_gb_area_A_pixel_m2"] = (
            grid["accessible_gb_area_A_pixel_m2"]
            .fillna(0.0)
            .clip(lower=0)
        )

        grid["accessible_gb_area_AB_pixel_m2"] = (
            grid["accessible_gb_area_AB_pixel_m2"]
            .fillna(0.0)
            .clip(lower=0)
        )

        # Enforce physical upper bound.
        grid["accessible_gb_area_A_pixel_m2"] = np.minimum(
            grid["accessible_gb_area_A_pixel_m2"],
            grid["physical_gb_area_m2"],
        )

        grid["accessible_gb_area_AB_pixel_m2"] = np.minimum(
            grid["accessible_gb_area_AB_pixel_m2"],
            grid["physical_gb_area_m2"],
        )
        
        # Diagnostic check: Pixel GB area shouldn't exceed the accessible polygon area
        bad_A = (grid["accessible_gb_area_A_pixel_m2"] > grid["accessible_polygon_area_A_m2"] + 200).sum()
        bad_AB = (grid["accessible_gb_area_AB_pixel_m2"] > grid["accessible_polygon_area_AB_m2"] + 200).sum()
        print(f"  Pixel GB > accessible polygon area A cells : {bad_A}")
        print(f"  Pixel GB > accessible polygon area AB cells: {bad_AB}")

        grid["accessible_A_method_ratio"] = np.where(
            grid["accessible_gb_area_A_est_m2"] > 0,
            grid["accessible_gb_area_A_pixel_m2"] / grid["accessible_gb_area_A_est_m2"],
            np.nan,
        )

        grid["accessible_AB_method_ratio"] = np.where(
            grid["accessible_gb_area_AB_est_m2"] > 0,
            grid["accessible_gb_area_AB_pixel_m2"] / grid["accessible_gb_area_AB_est_m2"],
            np.nan,
        )

        ratio_A = grid["accessible_A_method_ratio"].replace([np.inf, -np.inf], np.nan)
        ratio_AB = grid["accessible_AB_method_ratio"].replace([np.inf, -np.inf], np.nan)

        print("  Pixel method available.")
        print("  Bias ratio = pixel / proportional")
        print("    Ratio > 1: proportional estimate underestimates accessible green-blue.")
        print("    Ratio < 1: proportional estimate overestimates accessible green-blue.")

        if ratio_A.notna().sum() > 0:
            print(
                f"  Level A ratio: "
                f"mean={ratio_A.mean():.3f}, "
                f"median={ratio_A.median():.3f}, "
                f"P25={ratio_A.quantile(0.25):.3f}, "
                f"P75={ratio_A.quantile(0.75):.3f}"
            )

        if ratio_AB.notna().sum() > 0:
            print(
                f"  Level A+B ratio: "
                f"mean={ratio_AB.mean():.3f}, "
                f"median={ratio_AB.median():.3f}, "
                f"P25={ratio_AB.quantile(0.25):.3f}, "
                f"P75={ratio_AB.quantile(0.75):.3f}"
            )

        return grid, True

    except Exception as e:
        print(f"  Pixel method failed: {e}")
        print("  Falling back to proportional estimate only.")

        grid["accessible_gb_area_A_pixel_m2"] = np.nan
        grid["accessible_gb_area_AB_pixel_m2"] = np.nan
        grid["accessible_A_method_ratio"] = np.nan
        grid["accessible_AB_method_ratio"] = np.nan

        return grid, False


# ---------------------------------------------------------------------
# FINAL METRICS
# ---------------------------------------------------------------------

def _create_final_accessible_area(grid: gpd.GeoDataFrame, pixel_available: bool) -> gpd.GeoDataFrame:
    if pixel_available:
        grid["accessible_gb_area_A_final_m2"] = grid["accessible_gb_area_A_pixel_m2"].fillna(
            grid["accessible_gb_area_A_est_m2"]
        )
        grid["accessible_gb_area_AB_final_m2"] = grid["accessible_gb_area_AB_pixel_m2"].fillna(
            grid["accessible_gb_area_AB_est_m2"]
        )
        grid["accessible_gb_method"] = "pixel"
    else:
        grid["accessible_gb_area_A_final_m2"] = grid["accessible_gb_area_A_est_m2"]
        grid["accessible_gb_area_AB_final_m2"] = grid["accessible_gb_area_AB_est_m2"]
        grid["accessible_gb_method"] = "proportional"

    grid["accessible_gb_area_A_final_m2"] = np.minimum(
        grid["accessible_gb_area_A_final_m2"].clip(lower=0),
        grid["physical_gb_area_m2"],
    )

    grid["accessible_gb_area_AB_final_m2"] = np.minimum(
        grid["accessible_gb_area_AB_final_m2"].clip(lower=0),
        grid["physical_gb_area_m2"],
    )

    return grid


def _calculate_gap_metrics(grid: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    phys = grid["physical_gb_area_m2"]

    acc_A = grid["accessible_gb_area_A_final_m2"]
    acc_AB = grid["accessible_gb_area_AB_final_m2"]

    grid["absolute_gap_A_m2"] = (phys - acc_A).clip(lower=0)
    grid["absolute_gap_AB_m2"] = (phys - acc_AB).clip(lower=0)

    grid["relative_gap_A"] = np.where(grid["area_m2"] > 0, grid["absolute_gap_A_m2"] / grid["area_m2"], 0.0)
    grid["relative_gap_AB"] = np.where(grid["area_m2"] > 0, grid["absolute_gap_AB_m2"] / grid["area_m2"], 0.0)

    grid["gap_share_cell_A"] = grid["relative_gap_A"]
    grid["gap_share_cell_AB"] = grid["relative_gap_AB"]

    grid["gap_share_physical_A"] = np.where(
        phys > 0,
        grid["absolute_gap_A_m2"] / phys,
        np.nan,
    )

    grid["gap_share_physical_AB"] = np.where(
        phys > 0,
        grid["absolute_gap_AB_m2"] / phys,
        np.nan,
    )

    grid["conversion_ratio_A"] = np.where(phys > 0, acc_A / phys, np.nan)
    grid["conversion_ratio_AB"] = np.where(phys > 0, acc_AB / phys, np.nan)

    grid["conversion_ratio_A"] = pd.Series(grid["conversion_ratio_A"]).clip(0, 1).to_numpy()
    grid["conversion_ratio_AB"] = pd.Series(grid["conversion_ratio_AB"]).clip(0, 1).to_numpy()

    # Optional method-specific proportional metrics kept for diagnostics.
    grid["conversion_ratio_A_est"] = np.where(
        phys > 0,
        grid["accessible_gb_area_A_est_m2"] / phys,
        np.nan,
    )

    grid["conversion_ratio_AB_est"] = np.where(
        phys > 0,
        grid["accessible_gb_area_AB_est_m2"] / phys,
        np.nan,
    )

    if "accessible_gb_area_A_pixel_m2" in grid.columns:
        grid["conversion_ratio_A_pixel"] = np.where(
            phys > 0,
            grid["accessible_gb_area_A_pixel_m2"] / phys,
            np.nan,
        )

    if "accessible_gb_area_AB_pixel_m2" in grid.columns:
        grid["conversion_ratio_AB_pixel"] = np.where(
            phys > 0,
            grid["accessible_gb_area_AB_pixel_m2"] / phys,
            np.nan,
        )

    # Backward-compatible aliases use FINAL metrics.
    grid["accessible_area_A"] = grid["accessible_gb_area_A_final_m2"]
    grid["accessible_area_AB"] = grid["accessible_gb_area_AB_final_m2"]

    grid["physical_gb_area"] = grid["physical_gb_area_m2"]

    grid["absolute_gap"] = grid["absolute_gap_A_m2"]
    grid["relative_gap"] = grid["relative_gap_A"]
    grid["conversion_ratio"] = grid["conversion_ratio_A"]

    return grid


def _drop_previous_gap_columns(grid: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    gap_cols = [
        "accessible_polygon_area_A_m2",
        "accessible_polygon_area_AB_m2",

        "accessible_gb_area_A_est_m2",
        "accessible_gb_area_AB_est_m2",

        "accessible_gb_area_A_pixel_m2",
        "accessible_gb_area_AB_pixel_m2",
        "accessible_A_method_ratio",
        "accessible_AB_method_ratio",

        "accessible_gb_area_A_final_m2",
        "accessible_gb_area_AB_final_m2",
        "accessible_gb_method",

        "physical_gb_area_m2",

        "absolute_gap_A_m2",
        "absolute_gap_AB_m2",
        "relative_gap_A",
        "relative_gap_AB",
        "conversion_ratio_A",
        "conversion_ratio_AB",

        "conversion_ratio_A_est",
        "conversion_ratio_AB_est",
        "conversion_ratio_A_pixel",
        "conversion_ratio_AB_pixel",

        "accessible_area_A",
        "accessible_area_AB",
        "physical_gb_area",
        "absolute_gap",
        "relative_gap",
        "conversion_ratio",
    ]

    return grid.drop(columns=[c for c in gap_cols if c in grid.columns], errors="ignore")


def _print_summary(grid: gpd.GeoDataFrame):
    phys = grid["physical_gb_area_m2"]
    acc_A = grid["accessible_gb_area_A_final_m2"]
    acc_AB = grid["accessible_gb_area_AB_final_m2"]

    print("\nGap analysis summary:")
    print(f"  Accessible GB method used       : {grid['accessible_gb_method'].iloc[0]}")
    print(f"  Grids with physical_gb = 0       : {int((phys == 0).sum())}")
    print(f"  Grids with accessible GB A = 0   : {int((acc_A == 0).sum())}")
    print(f"  Grids with accessible GB AB = 0  : {int((acc_AB == 0).sum())}")

    print(f"  Mean physical GB pct             : {grid['physical_greenblue_pct'].mean():.3f}%")

    print(
        f"  Mean accessible polygon A pct    : "
        f"{(grid['accessible_polygon_area_A_m2'] / grid['area_m2'] * 100).mean():.3f}%"
    )

    print(
        f"  Mean accessible polygon AB pct   : "
        f"{(grid['accessible_polygon_area_AB_m2'] / grid['area_m2'] * 100).mean():.3f}%"
    )

    print(
        f"  Mean accessible GB A pct         : "
        f"{(grid['accessible_gb_area_A_final_m2'] / grid['area_m2'] * 100).mean():.3f}%"
    )

    print(
        f"  Mean accessible GB AB pct        : "
        f"{(grid['accessible_gb_area_AB_final_m2'] / grid['area_m2'] * 100).mean():.3f}%"
    )

    print(f"  Mean conversion ratio A          : {grid['conversion_ratio_A'].mean(skipna=True):.3f}")
    print(f"  Mean conversion ratio AB         : {grid['conversion_ratio_AB'].mean(skipna=True):.3f}")
    print(f"  Mean cell gap share A            : {grid['gap_share_cell_A'].mean(skipna=True):.3f}")
    print(f"  Mean cell gap share AB           : {grid['gap_share_cell_AB'].mean(skipna=True):.3f}")

    if "accessible_A_method_ratio" in grid.columns:
        ratio = grid["accessible_A_method_ratio"].replace([np.inf, -np.inf], np.nan)
        if ratio.notna().sum() > 0:
            print(
                f"  Method ratio A pixel/est         : "
                f"mean={ratio.mean():.3f}, median={ratio.median():.3f}"
            )

    mean_cr_a = grid["conversion_ratio_A"].mean(skipna=True)

    if mean_cr_a > 0.9:
        print("WARNING: Mean conversion ratio A is very high. Check accessible-space classification.")

    if mean_cr_a < 0.01:
        print("WARNING: Mean conversion ratio A is very low. Check OSM completeness and access-level rules.")


# ---------------------------------------------------------------------
# MAIN
# ---------------------------------------------------------------------

def calculate_gap_metrics() -> gpd.GeoDataFrame:
    print("--- STEP 10: PHYSICAL-TO-ACCESSIBLE GREEN-BLUE GAP ---")

    print("Loading grid with land-cover variables...")
    grid = _read_grid()

    required_cols = [
        "grid_id",
        "area_m2",
        "tree_pct",
        "grass_pct",
        "water_pct",
        "physical_greenblue_pct",
        "geometry",
    ]

    missing = [c for c in required_cols if c not in grid.columns]
    if missing:
        raise ValueError(f"Missing required grid columns: {missing}")

    print(f"Grid cells: {len(grid)}")
    print(f"Grid CRS: {grid.crs}")

    print("Loading accessible-space layers...")
    spaces_A = _read_access_layer("accessible_spaces_A", grid.crs)
    spaces_AB = _read_access_layer("accessible_spaces_AB", grid.crs)

    print(f"  Level A spaces   : {len(spaces_A)}")
    print(f"  Level A+B spaces : {len(spaces_AB)}")

    grid = _drop_previous_gap_columns(grid)

    print("Calculating accessible polygon area per grid...")

    area_A = _overlay_area_by_grid(
        grid=grid,
        spaces=spaces_A,
        output_col="accessible_polygon_area_A_m2",
    )

    area_AB = _overlay_area_by_grid(
        grid=grid,
        spaces=spaces_AB,
        output_col="accessible_polygon_area_AB_m2",
    )

    grid = grid.merge(area_A, on="grid_id", how="left", validate="one_to_one")
    grid = grid.merge(area_AB, on="grid_id", how="left", validate="one_to_one")

    grid["accessible_polygon_area_A_m2"] = grid["accessible_polygon_area_A_m2"].fillna(0.0)
    grid["accessible_polygon_area_AB_m2"] = grid["accessible_polygon_area_AB_m2"].fillna(0.0)

    grid = _clip_area_to_cell_area(grid, "accessible_polygon_area_A_m2")
    grid = _clip_area_to_cell_area(grid, "accessible_polygon_area_AB_m2")

    print("Calculating physical green-blue area...")

    physical_gb_frac = _to_fraction(grid["physical_greenblue_pct"])

    grid["physical_gb_area_m2"] = (
        physical_gb_frac * pd.to_numeric(grid["area_m2"], errors="coerce")
    ).clip(lower=0)

    print("Calculating Method A: proportional accessible green-blue estimate...")

    grid["accessible_gb_area_A_est_m2"] = (
        grid["accessible_polygon_area_A_m2"] * physical_gb_frac
    ).clip(lower=0)

    grid["accessible_gb_area_AB_est_m2"] = (
        grid["accessible_polygon_area_AB_m2"] * physical_gb_frac
    ).clip(lower=0)

    grid["accessible_gb_area_A_est_m2"] = np.minimum(
        grid["accessible_gb_area_A_est_m2"],
        grid["physical_gb_area_m2"],
    )

    grid["accessible_gb_area_AB_est_m2"] = np.minimum(
        grid["accessible_gb_area_AB_est_m2"],
        grid["physical_gb_area_m2"],
    )

    print("Trying Method B: pixel-based accessible green-blue estimate...")
    grid, pixel_available = _try_pixel_method(grid, spaces_A, spaces_AB)

    grid = _create_final_accessible_area(grid, pixel_available)
    grid = _calculate_gap_metrics(grid)

    _print_summary(grid)

    os.makedirs(os.path.dirname(config.GAP_GRID_GPKG) or ".", exist_ok=True)

    grid.to_file(config.GAP_GRID_GPKG, driver="GPKG", layer="grid_gap")

    print(f"\nSaved -> {config.GAP_GRID_GPKG}")
    
    # Save diagnostic CSV for Method A vs B
    diag_cols = [
        "grid_id",
        "physical_greenblue_pct",
        "physical_gb_area_m2",
        "accessible_polygon_area_A_m2",
        "accessible_gb_area_A_est_m2",
        "accessible_gb_area_A_final_m2",
        "accessible_polygon_area_AB_m2",
        "accessible_gb_area_AB_est_m2",
        "accessible_gb_area_AB_final_m2",
        "gap_share_cell_A",
        "gap_share_cell_AB",
        "gap_share_physical_A",
        "gap_share_physical_AB",
        "conversion_ratio_A",
        "conversion_ratio_AB",
    ]

    if "accessible_gb_area_A_pixel_m2" in grid.columns:
        diag_cols.extend([
            "accessible_gb_area_A_pixel_m2",
            "accessible_A_method_ratio",
            "conversion_ratio_A_pixel",
        ])

    if "accessible_gb_area_AB_pixel_m2" in grid.columns:
        diag_cols.extend([
            "accessible_gb_area_AB_pixel_m2",
            "accessible_AB_method_ratio",
            "conversion_ratio_AB_pixel",
        ])
        
    diag_csv = os.path.join(os.path.dirname(config.GAP_GRID_GPKG), "step10_method_a_b_diagnostics.csv")
    grid[diag_cols].to_csv(diag_csv, index=False)
    print(f"Saved diagnostic CSV -> {diag_csv}")
    print("Next step: run step12_network_accessibility.py")

    return grid


if __name__ == "__main__":
    calculate_gap_metrics()