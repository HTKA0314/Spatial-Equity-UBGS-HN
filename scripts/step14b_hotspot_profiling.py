"""
Step 14B: Profile Network and Area Gap Priority areas.

Input:
  config.EQUITY_METRICS_GPKG
    hanoi_grid_250m_equity_metrics.gpkg from Step 14

  config.ADMIN_PATH
    Hanoi ward/commune administrative boundaries

Output:
  ../outputs/step14_risk_hotspots/hotspots_profile.csv
  ../outputs/step14_risk_hotspots/hotspots_profile.gpkg
  ../outputs/step14_risk_hotspots/hotspots_profile_summary.txt
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


OUTPUT_DIR = os.path.join(config.OUTPUT_DIR, "step14_risk_hotspots")

OUTPUT_CSV  = os.path.join(OUTPUT_DIR, "hotspots_profile.csv")
OUTPUT_GPKG = os.path.join(OUTPUT_DIR, "hotspots_profile.gpkg")
OUTPUT_TXT  = os.path.join(OUTPUT_DIR, "hotspots_profile_summary.txt")


def _read_equity_grid() -> gpd.GeoDataFrame:
    if not os.path.exists(config.EQUITY_METRICS_GPKG):
        raise FileNotFoundError(
            f"Missing: {config.EQUITY_METRICS_GPKG}. Run Step 14 first."
        )

    try:
        grid = gpd.read_file(config.EQUITY_METRICS_GPKG, layer="grid_equity_metrics")
    except Exception:
        grid = gpd.read_file(config.EQUITY_METRICS_GPKG)

    if grid.empty:
        raise ValueError("Equity metrics grid is empty.")

    if grid.crs is None:
        raise ValueError("Equity metrics grid has no CRS.")

    if "grid_id" not in grid.columns:
        raise ValueError("grid_id missing from equity metrics grid.")

    grid["grid_id"] = grid["grid_id"].astype(int)

    return grid


def _read_admin(target_crs) -> gpd.GeoDataFrame:
    if not os.path.exists(config.ADMIN_PATH):
        raise FileNotFoundError(f"Missing admin boundary file: {config.ADMIN_PATH}")

    admin = gpd.read_file(config.ADMIN_PATH)

    if admin.empty:
        raise ValueError("Admin boundary file is empty.")

    if admin.crs is None:
        raise ValueError("Admin boundary file has no CRS.")

    if admin.crs != target_crs:
        admin = admin.to_crs(target_crs)

    admin = admin[admin.geometry.notna() & (~admin.geometry.is_empty)].copy()
    admin["geometry"] = admin.geometry.make_valid()

    return admin


def _standard_filter(grid: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    out = grid.copy()

    if "valid_priority_cell" in out.columns:
        out = out[out["valid_priority_cell"] == 1].copy()
    else:
        if "edge_cell" in out.columns:
            out = out[out["edge_cell"] == 0].copy()

        if "area_m2" in out.columns:
            min_area = getattr(config, "RISK_MIN_CELL_AREA_M2", 31250)
            out = out[
                pd.to_numeric(out["area_m2"], errors="coerce") >= min_area
            ].copy()

    return out


def _ensure_accessible_fractions(grid: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    grid = grid.copy()

    if "area_m2" not in grid.columns:
        return grid

    area = pd.to_numeric(grid["area_m2"], errors="coerce").replace(0, np.nan)

    if (
        "accessible_gb_A_frac" not in grid.columns
        and "accessible_gb_area_A_final_m2" in grid.columns
    ):
        grid["accessible_gb_A_frac"] = (
            pd.to_numeric(grid["accessible_gb_area_A_final_m2"], errors="coerce") / area
        ).clip(0, 1)

    if (
        "accessible_gb_AB_frac" not in grid.columns
        and "accessible_gb_area_AB_final_m2" in grid.columns
    ):
        grid["accessible_gb_AB_frac"] = (
            pd.to_numeric(grid["accessible_gb_area_AB_final_m2"], errors="coerce") / area
        ).clip(0, 1)

    return grid


def _get_priority_areas(grid: gpd.GeoDataFrame) -> tuple[gpd.GeoDataFrame, str]:
    """
    Extract main network-based priority intervention cells defined in Step 14.

    Main priority definition:
      - high-density grid cell
      - more than 1000 m from nearest Level A UGBS, or unreachable
      - located within an FDR-corrected LST hotspot

    Area-gap priority is retained only as a supplementary diagnostic and is
    not included in the main priority profile.
    """
    if "Priority_Area_Network" not in grid.columns:
        raise ValueError("Priority_Area_Network missing. Run Step 14 first.")

    hotspots = grid[grid["Priority_Area_Network"] == 1].copy()

    method = (
        "Main network priority: high-density cells located more than 1000 m "
        "from the nearest Level A UGBS or unreachable, and within an "
        "FDR-corrected LST hotspot"
    )

    return hotspots, method

def _join_admin(hotspots: gpd.GeoDataFrame, admin: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    admin_name_candidates = [
        "commune_name", "xaphuong", "ten_xa", "tenphuong", "ward_name", "name", "NAME_3"
    ]
    admin_code_candidates = ["maxa", "ma_xa", "ward_code", "code", "GID_3"]

    name_col = next((c for c in admin_name_candidates if c in admin.columns), None)
    code_col = next((c for c in admin_code_candidates if c in admin.columns), None)

    keep_cols = [c for c in [name_col, code_col, "geometry"] if c is not None]
    admin_sub = admin[keep_cols].copy()

    rename_map = {}
    if name_col is not None:
        rename_map[name_col] = "ward_name"
    if code_col is not None:
        rename_map[code_col] = "ward_code"

    admin_sub = admin_sub.rename(columns=rename_map)

    hotspot_polygons = hotspots[["grid_id", "geometry"]].copy()

    hotspots_pts = hotspots.copy()
    hotspots_pts["geometry"] = hotspots_pts.geometry.representative_point()

    joined_pts = gpd.sjoin(
        hotspots_pts,
        admin_sub,
        how="left",
        predicate="within",
    )

    joined_pts = (
        joined_pts
        .sort_values("grid_id")
        .drop_duplicates(subset=["grid_id"], keep="first")
    )

    joined_wgs = joined_pts.to_crs("EPSG:4326")
    joined_pts["lon"] = joined_wgs.geometry.x
    joined_pts["lat"] = joined_wgs.geometry.y

    joined_attrs = pd.DataFrame(joined_pts.drop(columns="geometry"))

    joined = hotspot_polygons.merge(
        joined_attrs,
        on="grid_id",
        how="left",
        validate="one_to_one",
    )

    joined = gpd.GeoDataFrame(joined, geometry="geometry", crs=hotspots.crs)

    return joined


def _select_output_columns(joined: gpd.GeoDataFrame) -> pd.DataFrame:
    priority_cols = [
        "Priority_Area_Network",
        "Priority_Area_AreaGap",
    ]

    candidate_cols = [
        "priority_rank",
        "grid_id",
        "ward_name",
        "ward_code",
        "lon",
        "lat",
        "population",
        "lst_p75_observed",
        "lst_p75_filled",
        "lst_p75_mean",
        "LST_p75",
        "tree_pct",
        "grass_pct",
        "built_pct",
        "water_pct",
        "crop_pct",
        "bare_pct",
        "physical_greenblue_pct",
        "accessible_gb_method",
        "accessible_gb_A_frac",
        "accessible_gb_AB_frac",
        "accessible_gb_area_A_final_m2",
        "accessible_gb_area_AB_final_m2",
        "dist_to_accessible_A_m",
        "dist_to_accessible_AB_m",
        "within_300m_A",
        "within_500m_A",
        "within_1000m_A",
        "relative_gap_A",
        "relative_gap_AB",
        "conversion_ratio_A",
        "conversion_ratio_AB",
        "lst_norm",
        "area_gap_norm",
        "network_gap_norm",
        "gb_deficit_norm",
    ] + priority_cols

    out_cols = [c for c in candidate_cols if c in joined.columns]
    return joined[out_cols].copy()


def _write_summary(df_out: pd.DataFrame, method: str):
    lines = []

    lines.append("=== PRIORITY HOTSPOT PROFILE SUMMARY ===")
    lines.append(f"Method: {method}")
    lines.append(f"Total network priority cells: {len(df_out)}")
    lines.append("")

    if "ward_name" in df_out.columns:
        ward_counts = df_out["ward_name"].value_counts(dropna=False).head(15)
        lines.append("Top wards by hotspot cell count:")
        for ward, count in ward_counts.items():
            lines.append(f"  {ward}: {count} cells")
        lines.append("")

    num_cols = {
        "population": "Mean population per cell",
        "built_pct": "Mean built-up (%)",
        "tree_pct": "Mean tree cover (%)",
        "physical_greenblue_pct": "Mean physical green-blue (%)",
        "accessible_gb_A_frac": "Mean accessible GB Level A fraction",
        "accessible_gb_area_A_final_m2": "Mean accessible GB Level A area (m²)",
        "dist_to_accessible_A_m": "Mean distance to Level A (m)",
        "relative_gap_A": "Mean relative gap A",
        "conversion_ratio_A": "Mean conversion ratio A",
    }

    lst_col = next((c for c in ["lst_p75_filled", "lst_p75_observed", "lst_p75_mean", "LST_p75"] if c in df_out.columns), None)

    if lst_col:
        num_cols[lst_col] = "Mean LST P75 (°C)"

    lines.append("Mean characteristics of priority cells:")
    for col, label in num_cols.items():
        if col in df_out.columns:
            val = pd.to_numeric(df_out[col], errors="coerce").replace([np.inf, -np.inf], np.nan).mean()
            lines.append(f"  {label}: {val:.3f}")

    lines.append("")

    if "ward_name" in df_out.columns and "population" in df_out.columns:
        ward_pop = df_out.groupby("ward_name")["population"].sum().sort_values(ascending=False).head(15)
        lines.append("Top wards by priority-cell population:")
        for ward, pop in ward_pop.items():
            lines.append(f"  {ward}: {pop:,.0f} people")
        lines.append("")

    top10_cols = ["priority_rank", "ward_name", "population", "Priority_Area_Network", "Priority_Area_AreaGap", "built_pct", "tree_pct"]
    if lst_col:
        top10_cols.append(lst_col)
    top10_cols = [c for c in top10_cols if c in df_out.columns]

    lines.append("Top 10 priority cells ranked by population:")
    lines.append(df_out[top10_cols].head(10).to_string(index=False))

    summary_text = "\n".join(lines)

    with open(OUTPUT_TXT, "w", encoding="utf-8") as f:
        f.write(summary_text)

    print("\n" + summary_text)
    print(f"\nSaved summary -> {OUTPUT_TXT}")


def profile_hotspots():
    print("--- STEP 14B: PRIORITY HOTSPOT PROFILING ---")

    os.makedirs(OUTPUT_DIR, exist_ok=True)

    print(f"Loading equity metrics: {config.EQUITY_METRICS_GPKG}")
    grid = _read_equity_grid()
    print(f"  Grid rows loaded: {len(grid)}")

    grid_filtered = _standard_filter(grid)
    grid_filtered = _ensure_accessible_fractions(grid_filtered)
    print(f"  Rows after valid sample filtering: {len(grid_filtered)}")

    # Compute priority areas based on Step 14 output
    hotspots, method = _get_priority_areas(grid_filtered)
    print(f"Identified {len(hotspots)} priority cells ({method})")

    if len(hotspots) == 0:
        print("No priority cells found.")
        return pd.DataFrame()

    print(f"Loading admin boundaries: {config.ADMIN_PATH}")
    admin = _read_admin(hotspots.crs)

    joined = _join_admin(hotspots, admin)
    joined = joined.sort_values(by="population", ascending=False).reset_index(drop=True)
    joined["priority_rank"] = joined.index + 1

    df_out = _select_output_columns(joined)
    df_out.to_csv(OUTPUT_CSV, index=False, encoding="utf-8-sig")

    joined.to_file(
        OUTPUT_GPKG,
        driver="GPKG",
        layer="network_priority_profile",
    )

    print(f"Saved profile CSV  -> {OUTPUT_CSV}")
    print(f"Saved profile GPKG -> {OUTPUT_GPKG}")

    _write_summary(df_out, method)

    return df_out


if __name__ == "__main__":
    profile_hotspots()