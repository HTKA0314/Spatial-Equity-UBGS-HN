"""
Step 1: Generate study area communes and 250m base grid.

Input:
  config.ADMIN_PATH
    Hanoi administrative boundaries (126 communes/wards) GeoPackage

Output:
  config.STUDY_AREA_GPKG
    Selected urban commune polygons

  config.BASE_GRID_GPKG
    Clipped 250m x 250m grid with density/commune attributes

Logic:
  - Selects communes that are "Phường" (wards) OR meet the population
    density threshold (config.URBAN_INCLUSION_DENSITY_THRESHOLD).
  - Excludes Sơn Tây and other peri-urban communes via config.EXCLUDED_COMMUNES.
  - Generates a grid of 250m cells clipped to the union of selected communes.
  - Removes geometric slivers (< 1% of expected cell area).
  - Flags edge cells (area < 98% of nominal cell area).
  - Joins commune name and density via centroid spatial join, with
    sjoin_nearest fallback for boundary-clipped cells.
  - Assigns density strata (Low / Medium / High) using config.DENSITY_BINS.
"""

import os
import sys
import warnings
from itertools import product

import numpy as np
import pandas as pd
import geopandas as gpd
from shapely.geometry import box

import config

# Ensure utf-8 encoding for Windows terminal output
sys.stdout.reconfigure(encoding="utf-8")
warnings.filterwarnings("ignore", category=UserWarning)


def generate_study_area():
    if not os.path.exists(config.ADMIN_PATH):
        print(f"Error: {config.ADMIN_PATH} not found.")
        return

    print("--- STEP 1: STUDY AREA & GRID GENERATION ---")
    print("Loading admin boundaries...")
    communes = gpd.read_file(config.ADMIN_PATH)
    communes = communes.to_crs(config.PROJ_CRS)

    name_col = getattr(config, "ADMIN_NAME_FIELD", "xaphuong")
    if name_col not in communes.columns:
        raise ValueError(f"Admin name column '{name_col}' not found in source file.")

    # Clean strings to prevent leading/trailing whitespace issues
    communes[name_col] = communes[name_col].astype(str).str.strip()

    # --- 1. Calculate area in square meters (preserve metric accuracy) ---
    print("Calculating commune population density dynamically from projected geometries...")
    communes['pop_sum'] = pd.to_numeric(communes[getattr(config, "ADMIN_POP_FIELD", "danso")], errors='coerce').fillna(0)
    
    communes['area_km2'] = communes.geometry.area / 1_000_000
    communes['density'] = communes['pop_sum'] / communes['area_km2']

    # --- 2. Apply filter to extract urbanized regions ---
    excluded_list = getattr(config, "EXCLUDED_COMMUNES", [
        "Phường Sơn Tây", "Phường Tùng Thiện", "Phường Trung Hưng", 
        "Phường Viên Sơn", "Phường Xuân Khanh", "Phường Lê Lợi", 
        "Phường Quang Trung", "Phường Ngô Quyền", "Phường Phú Thịnh"
    ])
    
    excluded_list = [str(x).strip() for x in excluded_list]

    is_urban_ward = communes[name_col].str.contains("Phường", na=False, case=False)
    is_high_density = communes['density'] >= config.URBAN_INCLUSION_DENSITY_THRESHOLD
    is_excluded = communes[name_col].isin(excluded_list)
    
    selected = communes[(is_urban_ward | is_high_density) & (~is_excluded)].copy()
    print(f"Selected {len(selected)} / {len(communes)} communes (Urban Wards + Density >= {config.URBAN_INCLUSION_DENSITY_THRESHOLD}, excluding Sơn Tây)")
    
    os.makedirs(os.path.dirname(config.STUDY_AREA_GPKG) or ".", exist_ok=True)
    selected.to_file(config.STUDY_AREA_GPKG, driver='GPKG', layer='study_area')

    # --- 3. Initialize and clip 250m spatial grid ---
    print(f"Generating {config.GRID_SPACING}m base grid...")
    minx, miny, maxx, maxy = selected.total_bounds
    x_coords = np.arange(minx, maxx, config.GRID_SPACING)
    y_coords = np.arange(miny, maxy, config.GRID_SPACING)

    cells = [box(x, y, x + config.GRID_SPACING, y + config.GRID_SPACING)
             for x, y in product(x_coords, y_coords)]

    raw_grid = gpd.GeoDataFrame(geometry=cells, crs=config.PROJ_CRS)
    raw_grid['raw_grid_id'] = raw_grid.index

    try:
        study_union = selected.union_all()
    except AttributeError:
        study_union = selected.unary_union

    print("Clipping grid to exact study area boundaries...")
    clipped_grid = gpd.clip(raw_grid, study_union)
    clipped_grid = clipped_grid.reset_index(drop=True)

    # Calculate clean area and remove geometric slivers
    clipped_grid['area_m2'] = clipped_grid.geometry.area.round(2)
    expected_area = config.GRID_SPACING ** 2
    sliver_threshold = expected_area * 0.01  # Remove slivers < 1% of standard area

    initial_count = len(clipped_grid)
    clipped_grid = clipped_grid[clipped_grid['area_m2'] >= sliver_threshold].copy()
    print(f"  Filtered out {initial_count - len(clipped_grid)} slivers (< {sliver_threshold} m2)")

    clipped_grid['grid_id'] = clipped_grid.index
    # Define Edge Cell (cells with area significantly cropped, < 98%)
    clipped_grid['edge_cell'] = (clipped_grid['area_m2'] < (expected_area * 0.98)).astype(int)

    # --- 4. Synchronize attributes via Spatial Join with boundary patch ---
    print("Joining demographic attributes to grid centroids...")
    grid_centroids = clipped_grid.copy()
    grid_centroids['geometry'] = grid_centroids.geometry.centroid

    joined_centroids = gpd.sjoin(
        grid_centroids[['grid_id', 'geometry']],
        selected[['geometry', 'density', name_col]].rename(columns={name_col: 'commune_name'}),
        how='left',
        predicate='intersects'
    ).drop(columns=['index_right'], errors='ignore')

    joined_centroids = joined_centroids.drop_duplicates(subset=['grid_id'])

    nan_mask = joined_centroids['commune_name'].isna()
    if nan_mask.any():
        nan_grids = joined_centroids[nan_mask].drop(columns=['commune_name', 'density'])
        print(f"  Vá lỗi điểm biên: Phát hiện {len(nan_grids)} ô lưới dính ranh giới. Đang định tuyến bắt cặp láng giềng gần nhất...")
        
        resolved_nan = gpd.sjoin_nearest(
            nan_grids,
            selected[['geometry', 'density', name_col]].rename(columns={name_col: 'commune_name'}),
            how='left',
            max_distance=250
        ).drop(columns=['index_right'], errors='ignore').drop_duplicates(subset=['grid_id'])
        
        # Update processed results back to the main table
        joined_centroids = pd.concat([joined_centroids[~nan_mask], resolved_nan], ignore_index=True)

    grid = clipped_grid.merge(joined_centroids[['grid_id', 'density', 'commune_name']], on='grid_id', how='left')

    # Classify urban density strata
    grid['urban_class'] = pd.cut(
        grid['density'],
        bins=getattr(config, 'DENSITY_GROUP_BINS', config.DENSITY_BINS),
        labels=getattr(config, 'DENSITY_GROUP_LABELS', config.DENSITY_LABELS),
        include_lowest=True,
    )

    # --- 5. Validation & Output ---
    print(f"\nValidation Result:")
    print(f"  Total grids generated  : {len(grid)}")
    print(f"  Edge cells (Boundary)  : {grid['edge_cell'].sum()}")
    print(f"  Unassigned cells (NaN) : {grid['commune_name'].isna().sum()} (Mục tiêu bắt buộc = 0)")
    print(f"  Urban classes distribution:\n{grid['urban_class'].value_counts().to_string()}")

    os.makedirs(os.path.dirname(config.BASE_GRID_GPKG) or ".", exist_ok=True)
    grid.to_file(config.BASE_GRID_GPKG, driver='GPKG', layer='grid_250m')
    print(f"SUCCESS: Saved pristine study area grid to -> {config.BASE_GRID_GPKG}")

    return selected, grid


if __name__ == "__main__":
    generate_study_area()