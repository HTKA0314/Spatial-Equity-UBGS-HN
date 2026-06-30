"""
Generate Q1-style figures for the Hanoi UGBS accessibility paper.

Scientific principle:
    Plot only figures that answer the paper's research questions.
    Do not turn every output into a figure.

Main figures:
    Fig 1. Study area and analytical workflow
    Fig 2. Physical green-blue cover, network accessibility, and physical-accessibility gap
    Fig 3. Population-weighted network accessibility ECDF
    Fig 4. Inequality in accessible environmental cooling resources: Lorenz curve
    Fig 5. Environmental justice bivariate map: access deficit x population density
    Fig 6. Priority intervention areas: heat hotspots with access deficits

Supplementary figures:
    S1. Physical, accessible, and missing green-blue resources by density class
    S2. Conversion ratio of physical UGBS into public accessible UGBS
    S3. Density-dependent land-cover/LST associations from interaction SLM
    S4. Preferred spatial regression coefficients
    S5. LST distribution by density class
    S6. Green illusion binned diagnostic

This script intentionally removes:
    - main "cooling effect" scatterplot
    - main physical-vs-accessible scatterplot
because both can mislead reviewers when accessible cover is strongly zero-inflated.

Author note:
    Regression plots are treated as supplementary association evidence, not causal effects.
"""

from __future__ import annotations

import os
import sys
import glob
from dataclasses import dataclass
from typing import Iterable, Optional, Tuple

import numpy as np
import pandas as pd

import sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))
from scripts.reporting.data_loader import load_clean_data as load_data

import geopandas as gpd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib.ticker import PercentFormatter
from scipy import stats

try:
    import config
except ImportError:
    config = None

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass


# ============================================================
# PATHS
# ============================================================

@dataclass(frozen=True)
class Paths:
    data_dir: str
    output_dir: str
    plots_dir: str
    equity_path: str
    admin_path: str
    step13_dir: str
    step13b_dir: str
    proj_crs: str


def get_paths() -> Paths:
    data_dir = r"d:\KLTN\data"
    output_dir = r"d:\KLTN\outputs"

    plots_dir = r"d:\KLTN\outputs\plots"

    equity_path = (
        getattr(config, "EQUITY_METRICS_GPKG", None)
        if config else None
    ) or os.path.join(data_dir, "hanoi_grid_250m_equity_metrics.gpkg")

    admin_path = os.path.join(data_dir, "hanoi_admin_126.gpkg")
    step13_dir = os.path.join(output_dir, "step13_spatial_regression")
    step13b_dir = os.path.join(output_dir, "step13b_stratified_regression")
    proj_crs = str(getattr(config, "PROJ_CRS", "EPSG:3405")) if config else "EPSG:3405"

    return Paths(
        data_dir=data_dir,
        output_dir=output_dir,
        plots_dir=plots_dir,
        equity_path=equity_path,
        admin_path=admin_path,
        step13_dir=step13_dir,
        step13b_dir=step13b_dir,
        proj_crs=proj_crs,
    )


PATHS = get_paths()
os.makedirs(PATHS.plots_dir, exist_ok=True)


# ============================================================
# CONSTANTS AND STYLE
# ============================================================

DEFAULT_DIJKSTRA_CUTOFF_M = 1500

DENSITY_ORDER = ["Low-density", "Medium-density", "High-density"]

DENSITY_CODE_TO_LABEL = {
    "Low_Density": "Low-density",
    "Medium_Density": "Medium-density",
    "High_Density": "High-density",
}

DENSITY_COLORS = {
    "Low-density": "#4C78A8",
    "Medium-density": "#F58518",
    "High-density": "#E45756",
}

BOUNDARY_COLOR = "#333333"


def set_style() -> None:
    import seaborn as sns
    sns.set_theme(
        style="whitegrid",
        palette="muted",
        font="Arial",
        rc={
            "font.size": 11,
            "figure.dpi": 300,
            "savefig.dpi": 400,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.labelsize": 12,
            "axes.titlesize": 14,
            "axes.titleweight": "bold",
            "xtick.labelsize": 10,
            "ytick.labelsize": 10,
        }
    )


def save_fig(fig: plt.Figure, filename: str) -> None:
    out = os.path.join(PATHS.plots_dir, filename)
    fig.savefig(out, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"Saved: {out}")


# ============================================================
# DATA HELPERS
# ============================================================

def find_col(
    gdf: gpd.GeoDataFrame,
    candidates: Iterable[str],
    label: str,
    required: bool = True,
) -> Optional[str]:
    for col in candidates:
        if col in gdf.columns:
            return col

    if required:
        raise ValueError(f"Missing {label}. Tried columns: {list(candidates)}")
    return None


def to_pct(s: pd.Series) -> pd.Series:
    x = pd.to_numeric(s, errors="coerce").astype(float)
    if x.max(skipna=True) <= 1.5:
        x = x * 100.0
    return x


def density_group(v: float) -> str | float:
    if pd.isna(v):
        return np.nan
    if v < 5000:
        return "Low-density"
    if v <= 15000:
        return "Medium-density"
    return "High-density"


def pop_weighted_mean(df: pd.DataFrame, val_col: str, weight_col: str = "population") -> float:
    valid = df.dropna(subset=[val_col, weight_col])
    if valid.empty or valid[weight_col].sum() == 0:
        return np.nan
    return np.average(valid[val_col], weights=valid[weight_col])



def quantile_vmax(
    series: pd.Series,
    q: float = 0.99,
    minimum: float = 1.0,
    positive_only: bool = False,
) -> float:
    x = pd.to_numeric(series, errors="coerce").replace([np.inf, -np.inf], np.nan).dropna()
    if positive_only:
        x = x[x > 0]
    if x.empty:
        return minimum
    return max(float(np.nanquantile(x, q)), minimum)


PLANNING_ACCESS_THRESHOLD_M = 1000
DIJKSTRA_SEARCH_LIMIT_M = 1500


def get_dijkstra_search_limit() -> int:
    return int(getattr(config, "DIJKSTRA_CUTOFF_M", DIJKSTRA_SEARCH_LIMIT_M)) if config else DIJKSTRA_SEARCH_LIMIT_M


def get_planning_access_threshold() -> int:
    return PLANNING_ACCESS_THRESHOLD_M


def get_physical_gb_area(gdf: gpd.GeoDataFrame) -> pd.Series:
    if "physical_gb_area_m2" in gdf.columns:
        return pd.to_numeric(gdf["physical_gb_area_m2"], errors="coerce")

    if "physical_greenblue_pct" in gdf.columns:
        return to_pct(gdf["physical_greenblue_pct"]) / 100.0 * pd.to_numeric(gdf["area_m2"], errors="coerce")

    if {"tree_pct", "grass_pct", "water_pct"}.issubset(gdf.columns):
        gb = to_pct(gdf["tree_pct"]) + to_pct(gdf["grass_pct"]) + to_pct(gdf["water_pct"])
        return gb.clip(0, 100) / 100.0 * pd.to_numeric(gdf["area_m2"], errors="coerce")

    raise ValueError("Cannot derive physical green-blue area.")


def get_accessible_area(gdf: gpd.GeoDataFrame, level: str = "A") -> pd.Series:
    candidates = [
        f"accessible_gb_area_{level}_final_m2",
        f"accessible_gb_area_{level}_pixel_m2",
        f"accessible_gb_area_{level}_est_m2",
        f"accessible_gb_{level}_m2",
        f"accessible_area_{level}",
    ]

    col = next((c for c in candidates if c in gdf.columns), None)
    if col is None:
        return pd.Series(0.0, index=gdf.index)

    return pd.to_numeric(gdf[col], errors="coerce").fillna(0)


def load_data_old() -> gpd.GeoDataFrame:
    if not os.path.exists(PATHS.equity_path):
        raise FileNotFoundError(PATHS.equity_path)

    try:
        gdf = gpd.read_file(PATHS.equity_path, layer="grid_equity_metrics")
    except Exception:
        gdf = gpd.read_file(PATHS.equity_path)

    if gdf.crs is None:
        raise ValueError("Input grid has no CRS.")

    if str(gdf.crs).upper() != PATHS.proj_crs.upper():
        gdf = gdf.to_crs(PATHS.proj_crs)

    if "edge_cell" in gdf.columns:
        gdf = gdf[gdf["edge_cell"] == 0].copy()

    if "area_m2" in gdf.columns:
        gdf = gdf[pd.to_numeric(gdf["area_m2"], errors="coerce") >= 31250].copy()

    for required_col in ["grid_id", "population", "area_m2"]:
        if required_col not in gdf.columns:
            raise ValueError(f"Missing required column: {required_col}")

    density_col = find_col(gdf, ["density", "pop_density_km2"], "population density")
    gdf["density_val"] = pd.to_numeric(gdf[density_col], errors="coerce")

    gdf["density_group"] = pd.Categorical(
        gdf["density_val"].apply(density_group),
        categories=DENSITY_ORDER,
        ordered=True,
    )

    lst_col = find_col(
        gdf,
        ["lst_p75_filled", "lst_p75_observed", "lst_p75_mean", "LST_p75", "LST_p75_C_mean"],
        "LST variable",
    )
    gdf["lst_temp"] = pd.to_numeric(gdf[lst_col], errors="coerce")
    gdf.attrs["lst_col"] = lst_col

    if "physical_greenblue_pct" in gdf.columns:
        gdf["physical_greenblue_pct"] = to_pct(gdf["physical_greenblue_pct"])
    elif {"tree_pct", "grass_pct", "water_pct"}.issubset(gdf.columns):
        gdf["physical_greenblue_pct"] = (
            to_pct(gdf["tree_pct"])
            + to_pct(gdf["grass_pct"])
            + to_pct(gdf["water_pct"])
        ).clip(0, 100)
    else:
        raise ValueError("Cannot derive physical_greenblue_pct.")

    gdf["physical_gb_m2"] = get_physical_gb_area(gdf)

    gdf["accessible_gb_A_m2"] = get_accessible_area(gdf, "A")
    gdf["accessible_gb_AB_m2"] = get_accessible_area(gdf, "AB")

    gdf["accessible_A_pct"] = (
        gdf["accessible_gb_A_m2"] / pd.to_numeric(gdf["area_m2"], errors="coerce") * 100
    ).clip(0, 100)

    gdf["accessible_AB_pct"] = (
        gdf["accessible_gb_AB_m2"] / pd.to_numeric(gdf["area_m2"], errors="coerce") * 100
    ).clip(0, 100)

    gdf["physical_access_gap_pct"] = (
        gdf["physical_greenblue_pct"] - gdf["accessible_A_pct"]
    ).clip(lower=0)

    gdf["population"] = pd.to_numeric(gdf["population"], errors="coerce").fillna(0)

    print(f"Loaded grid cells: {len(gdf):,}")
    print(f"LST column used: {lst_col}")
    print(gdf["density_group"].value_counts(dropna=False).sort_index().to_string())

    return gdf




# ============================================================
# MAP HELPERS
# ============================================================









def add_panel_label(ax, label: str) -> None:
    ax.text(
        0.03,
        0.97,
        label,
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=14,
        fontweight="bold",
    )


def create_distance_classes(gdf: gpd.GeoDataFrame) -> Tuple[gpd.GeoDataFrame, Optional[str]]:
    dist_col = find_col(
        gdf,
        ["dist_to_accessible_A_m", "dist_accessible_A_m", "nearest_accessible_A_m"],
        "distance to accessible Level A UGBS",
        required=False,
    )

    if dist_col is None:
        return gdf.copy(), None

    search_limit = get_dijkstra_search_limit()
    df = gdf.copy()
    dist = pd.to_numeric(df[dist_col], errors="coerce")

    # Enforce standard 1000m maximum limit for plotting classes
    labels = [
        "≤300 m",
        "300–500 m",
        "500–1000 m",
        ">1000 m /\nunreachable",
    ]
    cls = pd.Series(labels[-1], index=df.index, dtype="object")
    cls.loc[dist <= 300] = labels[0]
    cls.loc[(dist > 300) & (dist <= 500)] = labels[1]
    cls.loc[(dist > 500) & (dist <= 1000)] = labels[2]
    cls.loc[dist.isna() | (dist > 1000)] = labels[3]
    
    df["_network_access_class"] = pd.Categorical(cls, categories=labels, ordered=True)
    return df, dist_col


def safe_tercile(series: pd.Series) -> pd.Series:
    s = pd.to_numeric(series, errors="coerce")
    r = s.rank(method="average", na_option="keep")
    out = pd.qcut(r, 3, labels=False, duplicates="drop")
    return out.astype("float")


def normalize_hotspot(v) -> int:
    if pd.isna(v):
        return 0

    if isinstance(v, (int, float, np.integer, np.floating)):
        iv = int(round(float(v)))
        return iv if iv in (-1, 0, 1) else 0

    s = str(v).strip().upper()
    if s in ("HH", "HOT", "HOTSPOT", "1"):
        return 1
    if s in ("LL", "COLD", "COLDSPOT", "-1"):
        return -1
    return 0


def write_summary_tables(gdf: gpd.GeoDataFrame) -> None:
    total_pop = float(gdf["population"].sum())

    rows = {
        "n_grid_cells": len(gdf),
        "total_population": total_pop,
        "mean_physical_gb_pct": pop_weighted_mean(gdf, "physical_greenblue_pct"),
        "mean_accessible_A_pct": pop_weighted_mean(gdf, "accessible_A_pct"),
        "mean_physical_access_gap_pct": pop_weighted_mean(gdf, "physical_access_gap_pct"),
        "zero_accessible_A_pct_cells": (gdf["accessible_A_pct"] <= 0).mean() * 100,
        "zero_accessible_A_pct_population": (
            gdf.loc[gdf["accessible_A_pct"] <= 0, "population"].sum() / total_pop * 100
            if total_pop > 0 else np.nan
        ),
    }

    pd.DataFrame([rows]).to_csv(
        os.path.join(PATHS.plots_dir, "figure_summary_statistics.csv"),
        index=False,
    )


# ============================================================
# FIGURE 1. STUDY AREA AND WORKFLOW
# ============================================================



# ============================================================
# FIGURE 2. PHYSICAL COVER + NETWORK ACCESSIBILITY + GAP
# ============================================================


def fig3_conversion_ratio_by_density(gdf: gpd.GeoDataFrame) -> None:
    df = gdf.dropna(subset=["density_group"]).copy()

    summary = df.groupby("density_group", observed=True).apply(
        lambda x: pd.Series({
            "mean_physical": pop_weighted_mean(x, "physical_greenblue_pct", "population"),
            "mean_accessible": pop_weighted_mean(x, "accessible_A_pct", "population"),
            "mean_gap": pop_weighted_mean(x, "physical_access_gap_pct", "population"),
        })
    ).reset_index()

    summary["conversion_ratio"] = (summary["mean_accessible"] / summary["mean_physical"]) * 100

    fig, ax = plt.subplots(figsize=(5, 5))
    x = np.arange(len(DENSITY_ORDER))

    bars = ax.bar(
        x,
        summary["conversion_ratio"],
        color=[DENSITY_COLORS[g] for g in DENSITY_ORDER],
        edgecolor="black",
        linewidth=1,
    )

    ax.set_xticks(x)
    ax.set_xticklabels(DENSITY_ORDER)
    ax.set_ylabel("Conversion Ratio (%)\n(Accessible Cover / Physical Cover)")
    ax.set_title("Physical to Accessible Conversion Ratio", fontsize=11, fontweight="bold")
    ax.grid(axis="y", linestyle="--", alpha=0.3)

    for bar, ratio in zip(bars, summary["conversion_ratio"]):
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            ratio + 0.2,
            f"{ratio:.1f}%",
            ha="center",
            fontsize=9,
            fontweight="bold"
        )

    fig.tight_layout()
    save_fig(fig, "fig3_conversion_ratio_by_density.png")

    summary.to_csv(os.path.join(PATHS.plots_dir, "fig3_conversion_ratio_by_density_summary.csv"), index=False)


# ============================================================
# FIGURE 4. POPULATION-WEIGHTED ECDF
# ============================================================

def fig4_population_weighted_accessibility(gdf: gpd.GeoDataFrame) -> None:
    dist_col = find_col(
        gdf,
        ["dist_to_accessible_A", "dist_to_accessible_A_m", "dist_accessible_A_m", "nearest_accessible_A_m"],
        "distance to accessible public UGBS",
        required=False,
    )

    if dist_col is None:
        print("Skipping Fig 4: missing distance column.")
        return

    dist_col_ab = find_col(
        gdf,
        ["dist_to_accessible_AB", "dist_to_accessible_AB_m", "dist_accessible_AB_m", "nearest_accessible_AB_m"],
        "distance to public + semi-public UGBS",
        required=False,
    )

    cutoff = get_dijkstra_search_limit()
    data = gdf[gdf["population"] > 0].copy()
    total_pop = data["population"].sum()

    data[dist_col] = pd.to_numeric(data[dist_col], errors="coerce")
    reachable = data[dist_col].notna() & (data[dist_col] <= cutoff)

    d = data.loc[reachable, [dist_col, "population"]].sort_values(dist_col)

    if d.empty or total_pop <= 0:
        print("Skipping Fig 4: no reachable population.")
        return

    x = np.insert(d[dist_col].to_numpy(), 0, 0)
    y = np.insert(np.cumsum(d["population"].to_numpy()) / total_pop, 0, 0)

    fig, ax = plt.subplots(figsize=(8, 6.2))
    ax.plot(x, y, color="#1F77B4", linewidth=2.8, label="Level A UGBS")

    if dist_col_ab is not None:
        data[dist_col_ab] = pd.to_numeric(data[dist_col_ab], errors="coerce")
        reachable_ab = data[dist_col_ab].notna() & (data[dist_col_ab] <= cutoff)
        d_ab = data.loc[reachable_ab, [dist_col_ab, "population"]].sort_values(dist_col_ab)
        if not d_ab.empty:
            x_ab = np.insert(d_ab[dist_col_ab].to_numpy(), 0, 0)
            y_ab = np.insert(np.cumsum(d_ab["population"].to_numpy()) / total_pop, 0, 0)
            ax.plot(x_ab, y_ab, color="#B2182B", linewidth=2.3, linestyle="--", label="Level A+B UGBS")
        ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.16), ncol=2, frameon=False, fontsize=10)

    thresholds = [300, 500, 1000]
    summary = []

    for thr in thresholds:
        if thr <= cutoff:
            pct = data.loc[reachable & (data[dist_col] <= thr), "population"].sum() / total_pop * 100
            summary.append((thr, pct))
            ax.axvline(thr, color="#777777", linestyle=":" if thr != 500 else "--", linewidth=1.0)
            
            ax.text(
                thr - 15,
                0.045,
                f"{thr:,} m",
                fontsize=8.5,
                color="#333333",
                ha="right",
            )

    planning_threshold = get_planning_access_threshold()
    within_threshold_pop = data.loc[reachable & (data[dist_col] <= planning_threshold), "population"].sum()
    beyond_threshold_pct = (1.0 - (within_threshold_pop / total_pop)) * 100

    idx_threshold = np.searchsorted(d[dist_col].to_numpy(), planning_threshold, side='right')
    y_threshold = y[idx_threshold] if idx_threshold < len(y) else y[-1]

    # ax.annotate(
    #     f"Level A: {beyond_threshold_pct:.1f}% beyond {planning_threshold} m\nor unreachable",
    #     xy=(min(planning_threshold, cutoff), y_threshold),
    #     xytext=(min(planning_threshold, cutoff) * 0.50, max(0.08, y_threshold - 0.15)),
    #     fontsize=9,
    #     color="#1F77B4",
    #     arrowprops=dict(arrowstyle="->", color="#1F77B4", linewidth=1.0),
    # )

    ax.set_xlim(0, 1050)
    ax.set_ylim(0, 1.03)
    ax.yaxis.set_major_formatter(PercentFormatter(1))
    ax.set_xlabel("Network walking distance to nearest UGBS (m)")
    ax.set_ylabel("Cumulative share of residential population")
    # ax.set_title("Population-weighted accessibility to Level A and Level A+B UGBS", fontsize=12, fontweight="bold")
    ax.grid(True, linestyle="--", alpha=0.25)
    fig.tight_layout(rect=[0, 0.15, 1, 0.95])
    save_fig(fig, "fig4_population_weighted_accessibility.png")

    pd.DataFrame(summary, columns=["threshold_m", "population_within_pct"]).to_csv(
        os.path.join(PATHS.plots_dir, "fig4_accessibility_thresholds.csv"),
        index=False,
    )


# ============================================================
# FIGURE 5. LORENZ CURVE AND GINI
# ============================================================

def fig5_lorenz_inequality(gdf: gpd.GeoDataFrame) -> None:
    x_a, y_a, g_a = lorenz_xy(gdf, "accessible_gb_A_m2")
    x_ab, y_ab, g_ab = lorenz_xy(gdf, "accessible_gb_AB_m2")

    fig, ax = plt.subplots(figsize=(6.8, 6.8))

    ax.plot(
        [0, 1],
        [0, 1],
        color="#888888",
        linestyle="--",
        linewidth=1.2,
        label="Perfect equality",
    )

    if len(x_a):
        ax.plot(
            x_a,
            y_a,
            color="#1F77B4",
            linewidth=2.7,
            label=f"Level A UGBS (Gini = {g_a:.3f})",
        )
        ax.fill_between(x_a, x_a, y_a, color="#1F77B4", alpha=0.10)

    if len(x_ab):
        ax.plot(
            x_ab,
            y_ab,
            color="#B2182B",
            linewidth=2.4,
            linestyle="-.",
            label=f"Level A+B UGBS (Gini = {g_ab:.3f})",
        )
        # ax.fill_between(x_ab, x_ab, y_ab, color="#B2182B", alpha=0.05)

    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_aspect("equal")
    ax.xaxis.set_major_formatter(PercentFormatter(1))
    ax.yaxis.set_major_formatter(PercentFormatter(1))
    ax.set_xlabel("Cumulative share of population")
    ax.set_ylabel("Cumulative share of accessible UGBS area")
    # ax.set_title("Lorenz curves of population-weighted accessible UGBS", fontsize=12.7, fontweight="bold")
    ax.grid(True, linestyle="--", alpha=0.25)
    ax.legend(loc="upper left", fontsize=9)

    fig.tight_layout(rect=[0, 0.04, 1, 1])
    save_fig(fig, "fig5_lorenz_inequality.png")




# ============================================================
# FIGURE 6. LST P75 BY DENSITY (BOXPLOT)
# ============================================================

def fig6_lst_distribution_by_density(gdf: gpd.GeoDataFrame) -> None:
    data = gdf.dropna(subset=["density_group", "lst_temp", "population"]).copy()

    plot_data = []
    labels = []

    for group in DENSITY_ORDER:
        sub = data[data["density_group"] == group]
        if sub.empty:
            continue
        plot_data.append(sub["lst_temp"].values)
        labels.append(group)

    if not plot_data:
        print("Skipping Fig 6: no LST data.")
        return

    fig, ax = plt.subplots(figsize=(6.0, 4.0))

    bplot = ax.boxplot(
        plot_data,
        patch_artist=True,
        tick_labels=labels,
        showfliers=True,
        widths=0.4,
    )

    for patch in bplot["boxes"]:
        patch.set_facecolor("#1F77B4")
        patch.set_edgecolor("black")
        patch.set_alpha(1.0)

    for median in bplot["medians"]:
        median.set(color="#FF7F0E", linewidth=1.5)

    ax.set_ylabel("Land Surface Temperature P75 (°C)")
    # ax.set_title("Fig 6. Grid-cell LST P75 by population-density class", fontsize=11, fontweight="bold")
    ax.grid(axis="y", linestyle="-", alpha=0.2)
    ax.grid(axis="x", linestyle="-", alpha=0.2)

    fig.tight_layout()
    save_fig(fig, "fig6_lst_p75_by_density.png")


# ============================================================
# FIGURE 7. PRIORITY INTERVENTION AREAS
# ============================================================
# FIGURE 5. BIVARIATE LISA CLUSTERS
# ============================================================



# ============================================================
# FIGURE 6. PRIORITY INTERVENTION AREAS
# ============================================================



# ============================================================
# LORENZ CURVE HELPER
# ============================================================

def lorenz_xy(gdf: gpd.GeoDataFrame, asset_col: str) -> Tuple[np.ndarray, np.ndarray, float]:
    d = gdf[["population", asset_col]].copy()
    d["population"] = pd.to_numeric(d["population"], errors="coerce").fillna(0)
    d[asset_col] = pd.to_numeric(d[asset_col], errors="coerce").fillna(0)

    d = d[(d["population"] > 0) & (d[asset_col] >= 0)].copy()

    if d.empty or d[asset_col].sum() <= 0:
        return np.array([]), np.array([]), np.nan

    d["asset_pc"] = d[asset_col] / d["population"]
    d = d.sort_values("asset_pc")

    x = np.insert(np.cumsum(d["population"]) / d["population"].sum(), 0, 0)
    y = np.insert(np.cumsum(d[asset_col]) / d[asset_col].sum(), 0, 0)

    gini = 1.0 - 2.0 * (np.trapezoid(y, x) if hasattr(np, "trapezoid") else np.trapz(y, x))
    return x, y, gini


def to_bool_priority(s: pd.Series) -> pd.Series:
    """
    Convert common boolean / numeric / text priority columns to boolean.
    """
    def one(v):
        if pd.isna(v):
            return False
        if isinstance(v, (bool, np.bool_)):
            return bool(v)
        if isinstance(v, (int, float, np.integer, np.floating)):
            return float(v) != 0
        txt = str(v).strip().lower()
        return txt in {"true", "t", "yes", "y", "1", "priority", "hotspot"}

    return s.apply(one)


def find_bool_col(
    gdf: gpd.GeoDataFrame,
    candidates: Iterable[str],
    label: str,
) -> Optional[str]:
    for col in candidates:
        if col in gdf.columns:
            print(f"Using {label}: {col}")
            return col
    return None


def get_network_access_deficit(gdf: gpd.GeoDataFrame) -> Tuple[pd.Series, Optional[str], str]:
    """
    Network-accessibility deficit used for environmental-justice mapping.
    Based on the planning access threshold (e.g., 1000 m).

    0 = very close to public accessible UGBS
    1 = beyond threshold or unreachable
    """
    threshold = get_planning_access_threshold()

    dist_col = find_col(
        gdf,
        ["dist_to_accessible_A_m", "dist_accessible_A_m", "nearest_accessible_A_m"],
        "distance to public accessible UGBS",
        required=False,
    )

    if dist_col is not None:
        dist = pd.to_numeric(gdf[dist_col], errors="coerce")
        deficit = (dist / threshold).clip(lower=0, upper=1)
        deficit.loc[dist.isna() | (dist > threshold)] = 1.0
        note = f"Network-accessibility deficit is based on walking distance capped at planning threshold ({threshold:,} m); areas beyond this or unreachable are assigned maximum deficit (1.0)."
        return deficit, dist_col, note

    deficit = 1.0 - (pd.to_numeric(gdf["accessible_A_pct"], errors="coerce") / 100.0).clip(0, 1)
    note = "Network distance was unavailable; deficit is approximated as 1 − public accessible cover."
    return deficit, None, note


def classify_access_deficit_for_bivariate(gdf: gpd.GeoDataFrame) -> Tuple[pd.Series, str]:
    """
    Classify access deficit using absolute thresholds for proper representation.
    """
    dist_col = find_col(
        gdf,
        ["dist_to_accessible_A_m", "dist_accessible_A_m", "nearest_accessible_A_m"],
        "distance to public accessible UGBS",
        required=False,
    )

    if dist_col is not None:
        dist = pd.to_numeric(gdf[dist_col], errors="coerce")
        cls = pd.Series(2, index=gdf.index, dtype=float)
        cls.loc[dist <= 300] = 0
        cls.loc[(dist > 300) & (dist <= 1000)] = 1
        label = "Access deficit (Absolute: <=300m, 300-1000m, >1000m/Unreachable)"
        return cls, label

    acc = pd.to_numeric(gdf["accessible_A_pct"], errors="coerce")
    cls = pd.Series(2, index=gdf.index, dtype=float)
    cls.loc[acc >= 5] = 0
    cls.loc[(acc > 0) & (acc < 5)] = 1
    label = "Access deficit (Inverse accessible cover)"
    return cls, label
# ============================================================
# FIGURE 4. LORENZ CURVE
# ============================================================

def lorenz_xy(gdf: gpd.GeoDataFrame, asset_col: str) -> Tuple[np.ndarray, np.ndarray, float]:
    d = gdf[["population", asset_col]].copy()
    d["population"] = pd.to_numeric(d["population"], errors="coerce").fillna(0)
    d[asset_col] = pd.to_numeric(d[asset_col], errors="coerce").fillna(0)

    d = d[(d["population"] > 0) & (d[asset_col] >= 0)].copy()

    if d.empty or d[asset_col].sum() <= 0:
        return np.array([]), np.array([]), np.nan

    d["asset_pc"] = d[asset_col] / d["population"]
    d = d.sort_values("asset_pc")

    x = np.insert(np.cumsum(d["population"]) / d["population"].sum(), 0, 0)
    y = np.insert(np.cumsum(d[asset_col]) / d[asset_col].sum(), 0, 0)

    gini = 1.0 - 2.0 * (np.trapezoid(y, x) if hasattr(np, "trapezoid") else np.trapz(y, x))
    return x, y, gini


# ============================================================
# SUPPLEMENTARY 4. BIVARIATE NETWORK ACCESS DEFICIT & POPULATION
# ============================================================





# ============================================================
# SUPPLEMENTARY 5. PREFERRED SPATIAL REGRESSION COEFFICIENTS
# ============================================================

def supp5_preferred_spatial_regression_coefficients() -> None:
    pass

# ============================================================
# SUPPLEMENTARY FIGURE S1
# ============================================================



# ============================================================
# SUPPLEMENTARY FIGURE S2
# ============================================================
def supp2_conversion_ratio_by_density(gdf: gpd.GeoDataFrame) -> None:
    """
    Conversion ratio should be calculated as:
        sum(accessible area) / sum(physical area)

    Do not average cell-level ratios, because cells with tiny physical
    green-blue area can distort the mean ratio.
    """
    rows = []

    for group in DENSITY_ORDER:
        sub = gdf[gdf["density_group"] == group].copy()

        physical_sum = float(sub["physical_gb_m2"].sum())
        accessible_sum = float(sub["accessible_gb_A_m2"].sum())

        ratio = accessible_sum / physical_sum * 100 if physical_sum > 0 else np.nan

        rows.append({
            "density_group": group,
            "physical_gb_m2": physical_sum,
            "accessible_gb_A_m2": accessible_sum,
            "conversion_ratio_pct": ratio,
        })

    df = pd.DataFrame(rows)

    fig, ax = plt.subplots(figsize=(7, 4.8))
    x = np.arange(len(df))

    bars = ax.bar(
        x,
        df["conversion_ratio_pct"],
        color=[DENSITY_COLORS[g] for g in df["density_group"]],
        edgecolor="white",
    )

    for bar in bars:
        h = bar.get_height()
        if pd.notna(h):
            ax.text(
                bar.get_x() + bar.get_width() / 2,
                h + 0.4,
                f"{h:.1f}%",
                ha="center",
                fontsize=9,
            )

    ax.set_xticks(x)
    ax.set_xticklabels(df["density_group"])
    ax.set_ylabel("Accessible share of physical green-blue resources (%)")
    ax.set_title(
        "The Green Illusion: Very Little Green Space is Publicly Accessible",
        fontsize=13,
        fontweight="bold",
    )
    ax.yaxis.set_major_formatter(PercentFormatter(100))
    ax.grid(axis="y", linestyle="--", alpha=0.25)

    fig.tight_layout()
    save_fig(fig, "supp2_conversion_ratio_by_density.png")

    df.to_csv(
        os.path.join(PATHS.plots_dir, "supp2_conversion_ratio_by_density_summary.csv"),
        index=False,
    )

# ============================================================
# SUPPLEMENTARY FIGURE S3
# ============================================================

def supp3_density_dependent_lst_associations() -> None:
    path = os.path.join(PATHS.step13b_dir, "marginal_effects_by_density.csv")
    sig_path = os.path.join(PATHS.step13b_dir, "interaction_significance.csv")

    if not os.path.exists(path):
        print("Skipping S3: missing marginal_effects_by_density.csv")
        return

    df = pd.read_csv(path)
    df = df[
        (df["model_label"] == "Model-1A-Physical-Observed")
        & (df["model_name"] == "SLM")
    ].copy()

    if df.empty:
        print("Skipping S3: no Model-1A SLM marginal effects.")
        return

    sig = pd.read_csv(sig_path) if os.path.exists(sig_path) else pd.DataFrame()

    var_order = ["tree_frac", "grass_frac", "crop_frac", "water_frac"]
    var_labels = {
        "tree_frac": "Tree cover",
        "grass_frac": "Grass cover",
        "crop_frac": "Cropland",
        "water_frac": "Water cover",
    }

    strata = ["Low_Density", "Medium_Density", "High_Density"]

    fig, axes = plt.subplots(
        len(var_order),
        1,
        figsize=(8.4, 2.35 * len(var_order)),
        sharex=True,
    )

    if len(var_order) == 1:
        axes = [axes]

    for ax, var in zip(axes, var_order):
        d = df[df["variable"] == var].copy()
        d["stratum"] = pd.Categorical(d["stratum"], categories=strata, ordered=True)
        d = d.sort_values("stratum")

        labels = []

        for i, row in enumerate(d.itertuples()):
            label = DENSITY_CODE_TO_LABEL.get(row.stratum, row.stratum)

            if row.stratum != "Low_Density" and not sig.empty:
                match = sig[
                    (sig["model_label"] == "Model-1A-Physical-Observed")
                    & (sig["variable"] == var)
                ]

                if not match.empty:
                    p_col = "Medium_p_value" if row.stratum == "Medium_Density" else "High_p_value"
                    p = match.iloc[0].get(p_col, np.nan)

                    if pd.notna(p) and p < 0.05:
                        label += " *"

            labels.append(label)

            color = DENSITY_COLORS.get(
                DENSITY_CODE_TO_LABEL.get(row.stratum, ""),
                "#555555",
            )

            ax.errorbar(
                row.coef,
                i,
                xerr=[[row.coef - row.ci_lower], [row.ci_upper - row.coef]],
                fmt="o",
                color=color,
                ecolor=color,
                markersize=7,
                capsize=4,
                linewidth=1.5,
            )

        ax.axvline(0, color="#666666", linestyle="--", linewidth=1)
        ax.set_yticks(range(len(d)))
        ax.set_yticklabels(labels)
        ax.set_title(var_labels[var], fontsize=10.8, fontweight="bold")
        ax.grid(axis="x", linestyle="--", alpha=0.25)

    fig.supxlabel(
        "SLM coefficient (°C per unit land-cover fraction; not decomposed direct impact)",
        fontsize=10
    )

    fig.suptitle(
        "Density-Dependent Land-Cover/LST Associations",
        fontsize=13,
        fontweight="bold",
        y=0.99,
    )
    fig.text(
        0.5,
        0.01,
        "Error bars = 95% CI. * indicates the interaction term differs from low-density baseline at p < 0.05.",
        ha="center",
        fontsize=8.2,
        color="#555555",
        style="italic",
    )

    fig.tight_layout(rect=[0, 0.04, 1, 0.95])
    save_fig(fig, "supp3_density_dependent_lst_associations.png")


# ============================================================
# SUPPLEMENTARY FIGURE S4
# ============================================================

def _sig_color(p) -> str:
    if pd.isna(p):
        return "#AAAAAA"
    if p < 0.01:
        return "#C0392B"
    if p < 0.05:
        return "#E67E22"
    return "#AAAAAA"


def supp4_preferred_spatial_regression_coefficients() -> None:
    files = glob.glob(
        os.path.join(PATHS.step13_dir, "**", "coefficients.csv"),
        recursive=True,
    )

    if not files:
        print("Skipping S4: no coefficients.csv found.")
        return

    records = []

    for f in files:
        try:
            d = pd.read_csv(f)
            if {"variable", "coef", "ci_lower", "ci_upper"}.issubset(d.columns):
                d["_src"] = os.path.basename(os.path.dirname(f))
                records.append(d)
        except Exception:
            pass

    if not records:
        return

    df = pd.concat(records, ignore_index=True)

    if "is_preferred_for_reporting" in df.columns:
        pref_mask = (
            df["is_preferred_for_reporting"]
            .astype(str)
            .str.strip()
            .str.lower()
            .isin(["true", "1", "yes", "t"])
        )
        preferred = df[pref_mask].copy()
        if not preferred.empty:
            df = preferred

    wanted = [
        "tree_frac",
        "grass_frac",
        "crop_frac",
        "water_frac",
        "bare_frac",
        "log_dist_jrc_water",
    ]

    labels = {
        "tree_frac": "Tree cover",
        "grass_frac": "Grass cover",
        "crop_frac": "Cropland",
        "water_frac": "Water cover",
        "bare_frac": "Bare soil",
        "log_dist_jrc_water": "Log distance to water",
    }

    df = df[df["variable"].isin(wanted)].copy()
    df = df.drop_duplicates(subset=["variable"], keep="first")

    if df.empty:
        return

    if "p_value" not in df.columns:
        df["p_value"] = np.nan

    df["_color"] = df["p_value"].apply(_sig_color)

    y_order = [v for v in wanted if v in set(df["variable"])]
    y_pos = {v: i for i, v in enumerate(reversed(y_order))}

    fig, ax = plt.subplots(figsize=(8, 5.4))
    ax.axvline(0, color="#666666", linestyle="--", linewidth=1)

    for _, row in df.iterrows():
        v = row["variable"]

        if v not in y_pos:
            continue

        ax.errorbar(
            row["coef"],
            y_pos[v],
            xerr=[[row["coef"] - row["ci_lower"]], [row["ci_upper"] - row["coef"]]],
            fmt="o",
            color=row["_color"],
            ecolor=row["_color"],
            markersize=7,
            capsize=4,
            linewidth=1.5,
        )

    ax.set_yticks([y_pos[v] for v in reversed(y_order)])
    ax.set_yticklabels([labels[v] for v in reversed(y_order)])
    ax.set_xlabel("Estimated coefficient (°C per unit change in predictor)")
    ax.set_title("Preferred Spatial Regression Coefficients", fontsize=12, fontweight="bold")
    ax.grid(axis="x", linestyle="--", alpha=0.25)

    handles = [
        mpatches.Patch(color="#C0392B", label="p < 0.01"),
        mpatches.Patch(color="#E67E22", label="p < 0.05"),
        mpatches.Patch(color="#AAAAAA", label="Not significant"),
    ]

    ax.legend(
        handles=handles,
        loc="lower right",
        title="Significance",
        fontsize=8.5,
    )

    fig.text(
        0.5,
        0.01,
        "Built-up is the reference class. Coefficients are spatial model parameters and not decomposed direct, indirect, or total impacts.",
        ha="center",
        fontsize=8,
        color="#555555",
        style="italic",
    )

    fig.tight_layout(rect=[0, 0.04, 1, 1])
    save_fig(fig, "supp5_preferred_spatial_regression_coefficients.png")


# ============================================================
# SUPPLEMENTARY FIGURE S6
# ============================================================

    fig.tight_layout(rect=[0, 0.04, 1, 1])
    save_fig(fig, "supp5_preferred_spatial_regression_coefficients.png")


# ============================================================
# SUPPLEMENTARY FIGURE S6
# ============================================================

def supp6_green_illusion_binned(gdf: gpd.GeoDataFrame) -> None:
    """
    Diagnostic figure replacing the misleading raw scatter.

    It answers:
        Even among cells with high physical green-blue cover,
        how much public accessible cover exists?
    """
    data = gdf[
        (gdf["population"] > 0)
        & gdf["physical_greenblue_pct"].notna()
        & gdf["accessible_A_pct"].notna()
    ].copy()

    if data.empty:
        print("Skipping S6: no valid data.")
        return

    bins = [0, 10, 25, 50, 75, 100]
    labels = ["0–10", "10–25", "25–50", "50–75", "75–100"]

    data["_phys_bin"] = pd.cut(
        data["physical_greenblue_pct"],
        bins=bins,
        labels=labels,
        include_lowest=True,
        right=True,
    )

    rows = []

    for lab in labels:
        sub = data[data["_phys_bin"] == lab]

        if sub.empty:
            rows.append({
                "physical_bin": lab,
                "n": 0,
                "median_accessible_pct": np.nan,
                "zero_accessible_cells_pct": np.nan,
                "population": 0,
                "zero_accessible_population_pct": np.nan,
            })
            continue

        pop_sum = sub["population"].sum()

        rows.append({
            "physical_bin": lab,
            "n": len(sub),
            "mean_accessible_pct": sub["accessible_A_pct"].mean(),
            "zero_accessible_cells_pct": (sub["accessible_A_pct"] <= 0).mean() * 100,
            "population": pop_sum,
            "zero_accessible_population_pct": (
                sub.loc[sub["accessible_A_pct"] <= 0, "population"].sum() / pop_sum * 100
                if pop_sum > 0 else np.nan
            ),
        })

    df = pd.DataFrame(rows)

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8.2, 7.2), sharex=True, gridspec_kw={'height_ratios': [1, 1]})
    x = np.arange(len(df))

    # Panel A: Zero accessible cells
    bars = ax1.bar(
        x,
        df["zero_accessible_cells_pct"],
        color="#C7C7C7",
        edgecolor="white",
        label="Cells with zero public accessible cover",
    )

    ax1.set_ylim(0, 100)
    ax1.yaxis.set_major_formatter(PercentFormatter(100))
    ax1.set_ylabel("Cells with zero\naccessible cover (%)")
    ax1.grid(axis="y", linestyle="--", alpha=0.25)

    for bar, val in zip(bars, df["zero_accessible_cells_pct"]):
        if pd.notna(val):
            ax1.text(
                bar.get_x() + bar.get_width() / 2,
                val + 1.5,
                f"{val:.0f}%",
                ha="center",
                fontsize=8,
            )

    add_panel_label(ax1, "A")
    ax1.set_title(
        "Panel A: Cells with zero publicly accessible UGBS",
        fontsize=10,
        fontweight="bold"
    )

    fig.suptitle(
        "High Physical Green-Blue Cover Still Often Lacks Public Accessibility",
        fontsize=12,
        fontweight="bold",
        y=0.98
    )

    # Panel B: Mean accessible cover
    ax2.plot(
        x,
        df["mean_accessible_pct"],
        color="#1F77B4",
        marker="o",
        linewidth=2,
        label="Mean public accessible park cover",
    )
    ax2.set_ylabel("Mean publicly\naccessible UGBS (%)", color="#1F77B4")
    ax2.tick_params(axis="y", labelcolor="#1F77B4")
    ax2.grid(axis="y", linestyle="--", alpha=0.5)
    
    ax2.set_xticks(x)
    ax2.set_xticklabels(df["physical_bin"])
    ax2.set_xlabel("Satellite-visible physical green-blue cover (%)")
    add_panel_label(ax2, "B")

    fig.tight_layout(rect=[0, 0.04, 1, 0.95])
    save_fig(fig, "supp6_green_blue_illusion_binned.png")

    df.to_csv(os.path.join(PATHS.plots_dir, "supp6_green_illusion_binned_summary.csv"), index=False)


# ============================================================
# RUN
# ============================================================

def export_qgis_ready_data(gdf: gpd.GeoDataFrame) -> None:
    print("\n--- Exporting QGIS Ready GeoPackage ---")
    df = gdf.copy()
    
    # Priority Class
    threshold = get_planning_access_threshold()
    lst_col = find_col(df, ["lst_p75_filled_hotspot_k8", "lst_p75_observed_hotspot_k8", "lst_p75_mean_hotspot_k8", "LST_p75_hotspot_k8"], "LST hotspot", required=False)
    area_gap_col = find_col(df, ["physical_access_gap_pct_hotspot_k8", "relative_gap_A_hotspot_k8", "area_gap_A_hotspot_k8"], "area gap hotspot", required=False)
    dist_col = find_col(df, ["dist_to_accessible_A_m", "dist_accessible_A_m", "nearest_accessible_A_m"], "distance to public accessible UGBS", required=False)
    
    if "Priority_Area_Network" in df.columns or "Priority_Area_AreaGap" in df.columns:
        is_network_priority = df.get('Priority_Area_Network', 0) == 1
        is_areagap_priority = df.get('Priority_Area_AreaGap', 0) == 1
        
        df["priority_class"] = "Not priority"
        df.loc[is_network_priority, "priority_class"] = "Network priority"
        df.loc[is_areagap_priority, "priority_class"] = "Area-gap priority"
        df.loc[is_network_priority & is_areagap_priority, "priority_class"] = "Both priority types"
    
    # Bivariate Key
    access_cls, _ = classify_access_deficit_for_bivariate(df)
    pop = df["density_val"]
    pop_cls = pd.Series(0, index=df.index, dtype=float)
    pop_cls.loc[(pop >= 5000) & (pop <= 15000)] = 1
    pop_cls.loc[pop > 15000] = 2
    
    df["_bivar_key_str"] = access_cls.astype(str) + "_" + pop_cls.astype(str)
    
    # Network access class
    df_dist, _ = create_distance_classes(df)
    if "_network_access_class" in df_dist.columns:
        df["_network_access_class_label"] = df_dist["_network_access_class"].astype(str)
    
    # Drop complex object types before saving
    cols_to_drop = [c for c in df.columns if df[c].dtype == 'object' and c not in ['grid_id', 'priority_class', '_bivar_key_str', '_network_access_class_label']]
    df = df.drop(columns=cols_to_drop, errors="ignore")
    
    out_path = os.path.join(PATHS.output_dir, "hanoi_grid_250m_qgis_ready.gpkg")
    df.to_file(out_path, driver="GPKG")
    print(f"Exported QGIS data to: {out_path}")

def run() -> None:
    set_style()
    gdf = load_data()
    write_summary_tables(gdf)
    export_qgis_ready_data(gdf)

    print("\n=== MAIN FIGURES ===")
    fig3_conversion_ratio_by_density(gdf)
    fig4_population_weighted_accessibility(gdf)
    fig5_lorenz_inequality(gdf)
    fig6_lst_distribution_by_density(gdf)

    print("\n=== SUPPLEMENTARY FIGURES ===")
    supp2_conversion_ratio_by_density(gdf)
    supp5_preferred_spatial_regression_coefficients()
    supp6_green_illusion_binned(gdf)

    print(f"\nAll figures saved to: {PATHS.plots_dir}")


if __name__ == "__main__":
    run()
