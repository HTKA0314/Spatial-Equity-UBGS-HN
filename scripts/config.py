import os

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_DIR = os.path.join(BASE_DIR, "data")
OUTPUT_DIR = os.path.join(BASE_DIR, "outputs")

PLOTS_DIR = os.path.join(OUTPUT_DIR, "plots")
MAPS_DIR = os.path.join(OUTPUT_DIR, "maps")

PROJ_CRS = "EPSG:32648"  # UTM Zone 48N for Hanoi
GRID_SPACING = 250       # meters

# Minimum population density (people/km²) for a commune to be included
# in the study area. This is a spatial filter, NOT a grouping boundary.
URBAN_INCLUSION_DENSITY_THRESHOLD = 2000
DENSITY_THRESHOLD = URBAN_INCLUSION_DENSITY_THRESHOLD  # backward-compat alias

# Enclave Communes for Spatial Smoothing
ENCLAVE_COMMUNES = [
    'Phường Lĩnh Nam',   # Vết lõm kẹp ven sông Hồng (Quận Hoàng Mai)
    'Phường Thượng Cát'  # Lỗ thủng ở phía Tây Bắc (Quận Bắc Từ Liêm)
]

# Density strata for grid-cell grouping and plotting (3-class scheme).
# Applied AFTER spatial filtering by URBAN_INCLUSION_DENSITY_THRESHOLD.
DENSITY_BINS   = [0, 5000, 15000, float("inf")]
DENSITY_LABELS = ["Low-density", "Medium-density", "High-density"]

# Legacy aliases kept for backward compatibility with older scripts.
URBAN_BINS   = [URBAN_INCLUSION_DENSITY_THRESHOLD, 5000, 15000, float("inf")]
URBAN_LABELS = ["low_density", "medium_density", "high_density"]

MIN_WATER_AREA_M2 = 1000
NETWORK_BUFFER_M = 2000
MIN_OSM_WALK_NODES = 50000

ADMIN_PATH = os.path.join(DATA_DIR, "hanoi_admin_126.gpkg")
ADMIN_NAME_FIELD = "xaphuong"
ADMIN_POP_FIELD = "danso"

DEMOGRAPHICS_CSV = os.path.join(DATA_DIR, "hanoi_demographics.csv")

STUDY_AREA_GPKG = os.path.join(DATA_DIR, "study_area.gpkg")

# Base grid source file (.gpkg).
BASE_GRID_PATH = os.path.join(DATA_DIR, "hanoi_grid_250m_fn.gpkg")
BASE_GRID_GPKG = BASE_GRID_PATH  # backward-compat alias

LST_MAIN_CSV_PATH = os.path.join(DATA_DIR, "hanoi_grid250m_lst_p75_may_aug_2024.csv")
LST_SENSITIVITY_CSV_PATH = os.path.join(DATA_DIR, "hanoi_grid250m_lst_p75_aug_2024.csv")

LST_GRID_GPKG = os.path.join(DATA_DIR, "hanoi_grid_250m_with_lst.gpkg")

WATER_PATH = os.path.join(DATA_DIR, "hanoi_permanent_water_jrc.shp")
WATER_GRID_GPKG = os.path.join(DATA_DIR, "hanoi_grid_250m_with_water.gpkg")

WATER_CELL_ZERO_POP_THRESHOLD = 0.8
CENSUS_TOTAL = None
ALLOW_POP_FALLBACK = False

POP_RASTER_PATH = os.path.join(DATA_DIR, "hanoi_population_2025.tif")
POP_GRID_GPKG = os.path.join(DATA_DIR, "hanoi_grid_250m_with_pop.gpkg")

# ============================================================
# Cell validity filters
# One shared value used by all steps (Step 14, Step 13, etc.)
# ============================================================
MIN_VALID_CELL_AREA_M2 = 31250
EXCLUDE_EDGE_CELLS     = True

# Aliases for backward compatibility with individual step scripts
RISK_MIN_CELL_AREA_M2        = MIN_VALID_CELL_AREA_M2
RISK_EXCLUDE_EDGE_CELLS      = EXCLUDE_EDGE_CELLS
REGRESSION_MIN_CELL_AREA_M2  = MIN_VALID_CELL_AREA_M2
REGRESSION_EXCLUDE_EDGE_CELLS = EXCLUDE_EDGE_CELLS
HOTSPOT_K_MAIN = 8
HOTSPOT_K_SENSITIVITY = 4
HOTSPOT_PERMUTATIONS = 999
EQUITY_METRICS_GPKG = os.path.join(DATA_DIR, "hanoi_grid_250m_equity_metrics.gpkg")
LANDCOVER_CSV = os.path.join(DATA_DIR, "hanoi_lc_maroct2024_6class_rededge_grid250m.csv")
LANDCOVER_RASTER = os.path.join(DATA_DIR, "hanoi_lc_maroct2024_6class_rededge_raster.tif")
LANDCOVER_GRID_GPKG = os.path.join(DATA_DIR, "hanoi_grid_250m_with_landcover.gpkg")

# ---------------------------------------------------------------------
# Land-cover raster class IDs
# ---------------------------------------------------------------------
# Categorical IDs in LANDCOVER_RASTER corresponding to physical green-blue.
# Must match the exported Sentinel-2 land-cover raster class scheme:
#   0=water, 1=tree, 2=grass, 3=crop, 4=built, 5=bare
LC_GREENBLUE_CLASSES = [0, 1, 2]  # water + tree + grass (physical green-blue)
LC_GREEN_CLASSES     = [1, 2]     # tree + grass only (vegetation without water)

# Legacy alias kept so older scripts that reference LC_GREEN_CLASSES still work.
# Remove once all steps are updated.
_LC_GREEN_CLASSES_LEGACY = LC_GREENBLUE_CLASSES  # was [0,1,2] = water+tree+grass

LC_CLASS_MAP = {
    0: "water",
    1: "tree",
    2: "grass",
    3: "crop",
    4: "built",
    5: "bare",
}

REGRESSION_KNN_K = 8
VIF_WARNING_THRESHOLD = 10.0
# ============================================================
# OSM / Accessibility Data
# ============================================================
OSM_FEATURE_QUERY_BUFFER_M = 500
OSM_NETWORK_GRAPH = os.path.join(DATA_DIR, "hanoi_walking_network.graphml")
ACCESSIBLE_SPACES_GPKG = os.path.join(DATA_DIR, "accessible_spaces.gpkg")
ACCESSIBILITY_GPKG = os.path.join(DATA_DIR, "hanoi_grid_250m_walking_access.gpkg")

DESTINATION_SAMPLE_SPACING_M = 150
MIN_ACCESSIBLE_SPACE_AREA_M2 = 100
DIJKSTRA_CUTOFF_M = 1000
GAP_GRID_GPKG = os.path.join(DATA_DIR, "hanoi_grid_250m_with_gap.gpkg")

# ---------------------------------------------------------------------
# Step 12 network snapping diagnostics
# ---------------------------------------------------------------------
SNAP_WARNING_THRESHOLD_M = 100

# ---------------------------------------------------------------------
# Step 14 output
# ---------------------------------------------------------------------
# Define only once.
FINAL_GRID_GPKG = os.path.join(DATA_DIR, "hanoi_grid_250m_final_vars.gpkg")


for d in [DATA_DIR, OUTPUT_DIR, PLOTS_DIR, MAPS_DIR]:
    os.makedirs(d, exist_ok=True)