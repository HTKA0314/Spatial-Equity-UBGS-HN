# scripts/reporting/constants.py
import os
import sys

# Terminology constants
TERM_PHYSICAL = "Satellite-visible physical green-blue cover"
TERM_ACCESSIBLE = "Publicly accessible UGBS"

# Thresholds
PLANNING_ACCESS_THRESHOLD_M = 1000
ACCESS_CLASSES = [
    "≤300 m",
    "300–500 m",
    "500–1000 m",
    ">1000 m / unreachable"
]

# Density Groups
DENSITY_ORDER = ["Low-density", "Medium-density", "High-density"]
DENSITY_COLORS = {
    "Low-density": "#9ECAE1",     # Light Blue
    "Medium-density": "#4292C6",  # Medium Blue
    "High-density": "#084594",    # Dark Blue
}

DENSITY_CODE_TO_LABEL = {
    "Low_Density": "Low-density",
    "Medium_Density": "Medium-density",
    "High_Density": "High-density",
}

# Directories
PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
DATA_DIR = os.path.join(PROJECT_ROOT, "data")
OUTPUT_DIR = os.path.join(PROJECT_ROOT, "outputs")
PLOTS_DIR = os.path.join(OUTPUT_DIR, "plots")
TABLES_DIR = os.path.join(OUTPUT_DIR, "tables")

# Input files
EQUITY_GPKG = os.path.join(DATA_DIR, "hanoi_grid_250m_equity_metrics.gpkg")
ACCESS_GPKG = os.path.join(DATA_DIR, "hanoi_grid_250m_walking_access.gpkg")

# Ensure output dirs exist
os.makedirs(PLOTS_DIR, exist_ok=True)
os.makedirs(TABLES_DIR, exist_ok=True)
