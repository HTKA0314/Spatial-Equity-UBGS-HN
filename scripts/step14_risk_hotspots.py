import os
import sys
import warnings

import numpy as np
import pandas as pd
import geopandas as gpd

import libpysal as ps
from esda.getisord import G_Local
from esda.moran import Moran_Local_BV
from statsmodels.stats.multitest import multipletests

import config

warnings.filterwarnings("ignore", category=UserWarning)
sys.stdout.reconfigure(encoding="utf-8")


INPUT_GRID = getattr(
    config,
    "ACCESSIBILITY_GPKG",
    os.path.join(config.DATA_DIR, "hanoi_grid_250m_walking_access.gpkg"),
)

OUTPUT_GPKG = getattr(
    config,
    "EQUITY_METRICS_GPKG",
    os.path.join(config.DATA_DIR, "hanoi_grid_250m_equity_metrics.gpkg"),
)

OUTPUT_DIR = os.path.join(config.OUTPUT_DIR, "step14_risk_hotspots")


def scale_min_max(series: pd.Series) -> pd.Series:
    s = pd.to_numeric(series, errors="coerce").astype(float)
    s_min = s.min(skipna=True)
    s_max = s.max(skipna=True)

    if pd.isna(s_min) or pd.isna(s_max):
        return s * np.nan

    if s_max == s_min:
        return s * 0.0

    return ((s - s_min) / (s_max - s_min)).clip(0, 1)


def _read_input_grid() -> gpd.GeoDataFrame:
    if not os.path.exists(INPUT_GRID):
        raise FileNotFoundError(
            f"ACCESSIBILITY_GPKG not found: {INPUT_GRID}. Run Step 12 first."
        )

    try:
        gdf = gpd.read_file(INPUT_GRID, layer="grid_walking_access")
    except Exception:
        gdf = gpd.read_file(INPUT_GRID)

    if gdf.crs is None:
        raise ValueError("Input grid has no CRS.")

    if str(gdf.crs).upper() != str(config.PROJ_CRS).upper():
        gdf = gdf.to_crs(config.PROJ_CRS)

    if "grid_id" not in gdf.columns:
        raise ValueError("grid_id missing from input grid.")

    if "area_m2" not in gdf.columns:
        gdf["area_m2"] = gdf.geometry.area

    gdf["grid_id"] = gdf["grid_id"].astype(int)

    return gdf


def _find_lst_column(gdf: gpd.GeoDataFrame) -> str:
    candidates = [
        "lst_p75_filled",
        "lst_p75_observed",
        "lst_p75_mean",
        "LST_p75",
        "LST_p75_C_mean",
    ]

    col = next((c for c in candidates if c in gdf.columns), None)

    if col is None:
        raise ValueError(f"No LST column found. Available columns: {list(gdf.columns)}")

    return col


def _build_knn_weights(gdf: gpd.GeoDataFrame, k: int):
    coords = np.column_stack([
        gdf.geometry.representative_point().x,
        gdf.geometry.representative_point().y,
    ])

    w = ps.weights.KNN.from_array(coords, k=k)
    w.transform = "r"

    dists = []
    for i, neighs in w.neighbors.items():
        xi, yi = coords[i]
        for j in neighs:
            xj, yj = coords[j]
            dists.append(np.sqrt((xi - xj) ** 2 + (yi - yj) ** 2))

    print(
        f"  Gi* KNN k={k}: mean neighbor distance={np.mean(dists):.1f} m, "
        f"max={np.max(dists):.1f} m"
    )

    return w


def run_hotspot_analysis(
    gdf: gpd.GeoDataFrame,
    val_col: str,
    w,
    k_suffix: str,
    permutations: int = 999,
) -> pd.DataFrame:
    y = pd.to_numeric(gdf[val_col], errors="coerce").astype(float).to_numpy()

    if np.isnan(y).any():
        raise ValueError(f"{val_col} contains NaN before Gi* analysis.")

    g_star = G_Local(
        y,
        w,
        transform="r",
        permutations=permutations,
        star=True,
    )

    p_values = g_star.p_sim
    z_scores = g_star.Zs

    _, p_adjusted, _, _ = multipletests(
        p_values,
        alpha=0.05,
        method="fdr_bh",
    )

    status = np.zeros(len(gdf), dtype=int)
    status[(z_scores > 0) & (p_adjusted < 0.05)] = 1
    status[(z_scores < 0) & (p_adjusted < 0.05)] = -1

    status_raw = np.zeros(len(gdf), dtype=int)
    status_raw[(z_scores > 0) & (p_values < 0.05)] = 1
    status_raw[(z_scores < 0) & (p_values < 0.05)] = -1

    return pd.DataFrame(
        {
            f"{val_col}_gi_z_{k_suffix}": z_scores,
            f"{val_col}_gi_p_{k_suffix}": p_values,
            f"{val_col}_gi_padj_{k_suffix}": p_adjusted,
            f"{val_col}_hotspot_{k_suffix}": status,
            f"{val_col}_hotspot_raw_{k_suffix}": status_raw,
        },
        index=gdf.index,
    )


def run_bivariate_lisa(
    gdf: gpd.GeoDataFrame,
    col_x: str,
    col_y: str,
    w,
    permutations: int = 999,
) -> pd.DataFrame:
    x = pd.to_numeric(gdf[col_x], errors="coerce").astype(float).to_numpy()
    y = pd.to_numeric(gdf[col_y], errors="coerce").astype(float).to_numpy()

    if np.isnan(x).any() or np.isnan(y).any():
        raise ValueError(f"NaNs found in {col_x} or {col_y} before BiLISA.")

    bi_lisa = Moran_Local_BV(x, y, w, permutations=permutations)

    p_values = bi_lisa.p_sim
    q = bi_lisa.q

    _, p_adjusted, _, _ = multipletests(
        p_values,
        alpha=0.05,
        method="fdr_bh",
    )

    status_fdr = np.zeros(len(gdf), dtype=int)
    status_fdr[(q == 1) & (p_adjusted < 0.05)] = 1

    status_raw = np.zeros(len(gdf), dtype=int)
    status_raw[(q == 1) & (p_values < 0.05)] = 1

    return pd.DataFrame(
        {
            f"bilisa_{col_x}_x_{col_y}_q": q,
            f"bilisa_{col_x}_x_{col_y}_p": p_values,
            f"bilisa_{col_x}_x_{col_y}_padj": p_adjusted,
            f"bilisa_{col_x}_x_{col_y}_HH_fdr": status_fdr,
            f"bilisa_{col_x}_x_{col_y}_HH_raw": status_raw,
        },
        index=gdf.index,
    )


def _print_hotspot_stats(grid: gpd.GeoDataFrame, col_name: str, k_suffix: str):
    h_fdr = int((grid[f"{col_name}_hotspot_{k_suffix}"] == 1).sum())
    c_fdr = int((grid[f"{col_name}_hotspot_{k_suffix}"] == -1).sum())
    h_raw = int((grid[f"{col_name}_hotspot_raw_{k_suffix}"] == 1).sum())
    c_raw = int((grid[f"{col_name}_hotspot_raw_{k_suffix}"] == -1).sum())

    print(f"  * {col_name} K={k_suffix}:")
    print(f"    Hotspots FDR     : {h_fdr}")
    print(f"    Hotspots raw     : {h_raw}")
    print(f"    Coldspots FDR    : {c_fdr}")
    print(f"    Coldspots raw    : {c_raw}")


def analyze_risk_and_equity():
    print("--- STEP 14: RISK, EQUITY, AND HOTSPOT ANALYSIS ---")

    os.makedirs(OUTPUT_DIR, exist_ok=True)

    print("Loading Step 12 accessibility grid...")
    grid = _read_input_grid()

    print(f"Rows loaded: {len(grid)}")
    print(f"CRS: {grid.crs}")

    lst_col = _find_lst_column(grid)
    print(f"LST variable used for EBI/TRI: {lst_col}")

    required_cols = [
        lst_col,
        "area_m2",
        "population",
        "relative_gap_A",
        "dist_to_accessible_A_m",
        "geometry",
    ]

    missing = [c for c in required_cols if c not in grid.columns]
    if missing:
        raise ValueError(f"Missing required columns for Step 14: {missing}")

    mask_lst = pd.Series(True, index=grid.index)

    min_area = getattr(config, "RISK_MIN_CELL_AREA_M2", 31250)
    mask_lst = mask_lst & (pd.to_numeric(grid["area_m2"], errors="coerce") >= min_area)

    if "edge_cell" in grid.columns:
        exclude_edge = getattr(config, "RISK_EXCLUDE_EDGE_CELLS", True)
        if exclude_edge:
            mask_lst = mask_lst & (grid["edge_cell"] == 0)

    mask_lst = mask_lst & grid[lst_col].notna()

    population_num = pd.to_numeric(grid["population"], errors="coerce").fillna(0)
    mask_risk = mask_lst & (population_num > 0)

    grid["valid_priority_cell"] = mask_risk.astype(int)

    print(f"Total grid cells: {len(grid)}")
    print(f"Valid LST/risk-surface cells: {int(mask_lst.sum())}")
    print(f"Valid populated risk cells: {int(mask_risk.sum())}")

    print("Calculating burden and deficit indices...")

    grid["lst_norm"] = np.nan
    grid.loc[mask_lst, "lst_norm"] = scale_min_max(grid.loc[mask_lst, lst_col])

    grid["area_gap_norm"] = np.nan
    grid.loc[mask_lst, "area_gap_norm"] = scale_min_max(
        grid.loc[mask_lst, "relative_gap_A"]
    )

    cutoff = getattr(config, "DIJKSTRA_CUTOFF_M", 1000)
    penalty = getattr(config, "UNREACHABLE_PENALTY_M", 3000)

    dist = pd.to_numeric(grid["dist_to_accessible_A_m"], errors="coerce")
    dist_capped = np.where(
        (np.isfinite(dist)) & (dist <= cutoff),
        dist,
        penalty
    )

    grid["dist_accessible_A_capped_m"] = dist_capped

    grid["access_deficit_A"] = dist / 1000.0
    grid["access_deficit_A"] = grid["access_deficit_A"].clip(upper=1.0)
    grid["access_deficit_A"] = grid["access_deficit_A"].fillna(1.0)

    area_km2 = pd.to_numeric(grid["area_m2"], errors="coerce") / 1e6
    grid["pop_density"] = pd.to_numeric(grid["population"], errors="coerce").fillna(0) / area_km2
    
    conditions = [
        (grid["pop_density"] < 5000),
        (grid["pop_density"] >= 5000) & (grid["pop_density"] <= 15000),
        (grid["pop_density"] > 15000)
    ]
    choices = ["Low-density", "Medium-density", "High-density"]
    grid["density_group"] = np.select(conditions, choices, default="Unknown")

    grid["network_gap_norm"] = np.nan
    grid.loc[mask_lst, "network_gap_norm"] = scale_min_max(
        grid.loc[mask_lst, "dist_accessible_A_capped_m"]
    )

    # gb_deficit_norm removed as it is not used in hotspot analysis

    # EBI and TRI removed as requested

    print("Building KNN weights for Gi* hotspot analysis...")

    gdf_clean = grid[mask_lst].copy()
    gdf_clean = gdf_clean.dropna(subset=[lst_col, "dist_accessible_A_capped_m", "relative_gap_A"]).copy()

    k_main = getattr(config, "HOTSPOT_K_MAIN", 8)
    k_sens = getattr(config, "HOTSPOT_K_SENSITIVITY", 4)
    permutations = getattr(config, "HOTSPOT_PERMUTATIONS", 999)

    w_main = _build_knn_weights(gdf_clean, k_main)
    w_sens = _build_knn_weights(gdf_clean, k_sens)

    print(f"Running Getis-Ord Gi* with {permutations} permutations...")

    print("Running Gi* for empirical variables (LST, Network Gap, and Area Gap)...")
    res_lst_main = run_hotspot_analysis(gdf_clean, lst_col, w_main, f"k{k_main}", permutations)
    res_lst_sens = run_hotspot_analysis(gdf_clean, lst_col, w_sens, f"k{k_sens}", permutations)
    res_gap_main = run_hotspot_analysis(gdf_clean, "dist_accessible_A_capped_m", w_main, f"k{k_main}", permutations)
    res_gap_sens = run_hotspot_analysis(gdf_clean, "dist_accessible_A_capped_m", w_sens, f"k{k_sens}", permutations)
    res_area_main = run_hotspot_analysis(gdf_clean, "relative_gap_A", w_main, f"k{k_main}", permutations)
    res_area_sens = run_hotspot_analysis(gdf_clean, "relative_gap_A", w_sens, f"k{k_sens}", permutations)

    print("Running Bivariate LISA (Access Deficit x Pop Density and LST)...")
    res_bilisa_pop = run_bivariate_lisa(gdf_clean, "access_deficit_A", "pop_density", w_main, permutations)
    res_bilisa_lst = run_bivariate_lisa(gdf_clean, "access_deficit_A", lst_col, w_main, permutations)

    gdf_results = gdf_clean.join([
        res_lst_main,
        res_lst_sens,
        res_gap_main,
        res_gap_sens,
        res_area_main,
        res_area_sens,
        res_bilisa_pop,
        res_bilisa_lst,
    ])

    hotspot_cols = (
        list(res_lst_main.columns)
        + list(res_lst_sens.columns)
        + list(res_gap_main.columns)
        + list(res_gap_sens.columns)
        + list(res_area_main.columns)
        + list(res_area_sens.columns)
        + list(res_bilisa_pop.columns)
        + list(res_bilisa_lst.columns)
    )

    index_cols = [
        "grid_id",
        "lst_norm",
        "area_gap_norm",
        "network_gap_norm",
        "valid_priority_cell",
        "dist_accessible_A_capped_m",
        "relative_gap_A",
        "access_deficit_A",
        "pop_density",
        "density_group",
    ]

    result_cols = index_cols + hotspot_cols
    cols_to_drop = [c for c in result_cols if c != "grid_id" and c in grid.columns]
    grid = grid.drop(columns=cols_to_drop, errors="ignore")

    merge_cols = ["grid_id"] + [
        c for c in result_cols
        if c != "grid_id" and c in gdf_results.columns
    ]

    grid = grid.merge(
        gdf_results[merge_cols],
        on="grid_id",
        how="left",
        validate="one_to_one",
    )

    main_lst_hotspot_col = f"{lst_col}_hotspot_k{k_main}"
    main_area_hotspot_col = f"relative_gap_A_hotspot_k{k_main}"

    # Severe access deficit is defined as >1000 m or unreachable.
    # Unreachable cells are encoded as the penalty distance and therefore also satisfy >1000 m.
    grid["Priority_Area_Network"] = np.where(
        ((grid["dist_accessible_A_capped_m"] > 1000) | grid["dist_accessible_A_capped_m"].isna())
        & (grid["density_group"] == "High-density")
        & (grid[main_lst_hotspot_col] == 1),
        1,
        0,
    )
    
    grid["Priority_Area_AreaGap"] = np.where(
        (grid[main_lst_hotspot_col] == 1) & (grid[main_area_hotspot_col] == 1),
        1,
        0,
    )

    valid_mask = mask_lst

    pop_numeric = pd.to_numeric(grid["population"], errors="coerce").fillna(0)
    
    grid["population_in_Priority_Network"] = np.where(
        grid["Priority_Area_Network"] == 1,
        pop_numeric,
        0.0,
    )

    grid["population_in_Priority_AreaGap"] = np.where(
        grid["Priority_Area_AreaGap"] == 1,
        pop_numeric,
        0.0,
    )

    total_population_all = pop_numeric.sum()
    total_population_valid = pop_numeric[valid_mask].sum()
    priority_net_pop = grid["population_in_Priority_Network"].sum()
    priority_area_pop = grid["population_in_Priority_AreaGap"].sum()

    print("\nHotspot summary:")
    _print_hotspot_stats(grid, lst_col, f"k{k_main}")
    _print_hotspot_stats(grid, "dist_accessible_A_capped_m", f"k{k_main}")
    print(f"  * Priority Area (Network): {int((grid['Priority_Area_Network'] == 1).sum())} cells")
    print(f"  * Priority Area (Area Gap): {int((grid['Priority_Area_AreaGap'] == 1).sum())} cells")

    gi_sens_records = []
    for var in [lst_col, "dist_accessible_A_capped_m", "relative_gap_A"]:
        main_col = f"{var}_hotspot_k{k_main}"
        sens_col = f"{var}_hotspot_k{k_sens}"
        if main_col in gdf_results.columns and sens_col in gdf_results.columns:
            agree = (gdf_results[main_col] == gdf_results[sens_col]).mean() * 100
            n_hot_main = (gdf_results[main_col] == 1).sum()
            n_hot_sens = (gdf_results[sens_col] == 1).sum()
            gi_sens_records.append({
                "Variable": var,
                "Agreement_Pct": round(agree, 1),
                "Hotspots_K8": n_hot_main,
                "Hotspots_K4": n_hot_sens,
            })

    if gi_sens_records:
        gi_sens_path = os.path.join(OUTPUT_DIR, f"gi_sensitivity_k{k_main}_vs_k{k_sens}.csv")
        pd.DataFrame(gi_sens_records).to_csv(gi_sens_path, index=False, encoding="utf-8-sig")
        print(f"\nGi* sensitivity table saved: {gi_sens_path}")

    print("\nPopulation in hotspots:")
    if total_population_valid > 0:
        print(
            f"  Population in Priority Area (Network) (% of valid sample): "
            f"{priority_net_pop:,.0f} ({priority_net_pop / total_population_valid * 100:.1f}%)"
        )
        print(
            f"  Population in Priority Area (Area Gap) (% of valid sample): "
            f"{priority_area_pop:,.0f} ({priority_area_pop / total_population_valid * 100:.1f}%)"
        )

    if total_population_all > 0:
        print(
            f"  Population in Priority Area (Network) (% of all grid population): "
            f"{priority_net_pop:,.0f} ({priority_net_pop / total_population_all * 100:.1f}%)"
        )
        print(
            f"  Population in Priority Area (Area Gap) (% of all grid population): "
            f"{priority_area_pop:,.0f} ({priority_area_pop / total_population_all * 100:.1f}%)"
        )

    os.makedirs(os.path.dirname(OUTPUT_GPKG) or ".", exist_ok=True)

    grid.to_file(
        OUTPUT_GPKG,
        driver="GPKG",
        layer="grid_equity_metrics",
    )

    print(f"\nSaved -> {OUTPUT_GPKG}")

    report_path = os.path.join(OUTPUT_DIR, "hotspot_summary.txt")

    with open(report_path, "w", encoding="utf-8") as f:
        f.write("=== HANOI STEP 14: RISK AND HOTSPOT SUMMARY ===\n\n")
        f.write(f"Input: {INPUT_GRID}\n")
        f.write(f"Output: {OUTPUT_GPKG}\n\n")

        f.write("--- BASIC COUNTS ---\n")
        f.write(f"Total grid cells: {len(grid)}\n")
        f.write(f"Valid LST/risk-surface cells: {int(valid_mask.sum())}\n\n")

        f.write("--- VARIABLES ---\n")
        f.write(f"LST variable used: {lst_col}\n")
        f.write(f"Dijkstra cutoff used for network gap: {cutoff} m\n")
        f.write(f"Main Gi* K: {k_main}\n")
        f.write(f"Sensitivity Gi* K: {k_sens}\n")
        f.write(f"Permutations: {permutations}\n\n")

        f.write(f"--- GETIS-ORD GI* HOTSPOTS K={k_main} ---\n")
        f.write(f"Priority Area (Network): {int((grid['Priority_Area_Network'] == 1).sum())}\n")
        f.write(f"Priority Area (Area Gap): {int((grid['Priority_Area_AreaGap'] == 1).sum())}\n")
        f.write(f"LST Hotspots FDR: {(grid[main_lst_hotspot_col] == 1).sum()}\n")
        f.write(f"LST Coldspots FDR: {(grid[main_lst_hotspot_col] == -1).sum()}\n")
        main_gap_hotspot_col = f"dist_accessible_A_capped_m_hotspot_k{k_main}"
        if main_gap_hotspot_col in grid.columns:
            f.write(f"Network Gap Hotspots FDR: {(grid[main_gap_hotspot_col] == 1).sum()}\n")
            f.write(f"Network Gap Coldspots FDR: {(grid[main_gap_hotspot_col] == -1).sum()}\n\n")
            
        f.write(f"--- BIVARIATE LISA ---\n")
        f.write(f"BiLISA (Access Deficit x Pop Density) HH FDR: {(grid['bilisa_access_deficit_A_x_pop_density_HH_fdr'] == 1).sum()}\n")
        f.write(f"BiLISA (Access Deficit x Pop Density) HH raw: {(grid['bilisa_access_deficit_A_x_pop_density_HH_raw'] == 1).sum()}\n")
        f.write(f"BiLISA (Access Deficit x LST) HH FDR: {(grid['bilisa_access_deficit_A_x_' + lst_col + '_HH_fdr'] == 1).sum()}\n")
        f.write(f"BiLISA (Access Deficit x LST) HH raw: {(grid['bilisa_access_deficit_A_x_' + lst_col + '_HH_raw'] == 1).sum()}\n\n")

        f.write("--- POPULATION IN PRIORITY AREAS ---\n")
        if total_population_valid > 0:
            f.write(
                f"Priority Area Network pop (% of valid sample): "
                f"{priority_net_pop:,.0f} ({priority_net_pop / total_population_valid * 100:.1f}%)\n"
            )
            f.write(
                f"Priority Area AreaGap pop (% of valid sample): "
                f"{priority_area_pop:,.0f} ({priority_area_pop / total_population_valid * 100:.1f}%)\n"
            )
            f.write(f"Valid hotspot-analysis population: {total_population_valid:,.0f}\n")

        if total_population_all > 0:
            f.write(
                f"Priority Area Network pop (% of all grid): "
                f"{priority_net_pop:,.0f} ({priority_net_pop / total_population_all * 100:.1f}%)\n"
            )
            f.write(
                f"Priority Area AreaGap pop (% of all grid): "
                f"{priority_area_pop:,.0f} ({priority_area_pop / total_population_all * 100:.1f}%)\n"
            )
            f.write(f"All retained grid population: {total_population_all:,.0f}\n")

    print(f"Report saved -> {report_path}")
    print("\nDONE STEP 14.")

    return grid


if __name__ == "__main__":
    analyze_risk_and_equity()