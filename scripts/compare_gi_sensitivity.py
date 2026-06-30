"""
compare_gi_sensitivity.py
--------------------------
Compare Getis-Ord Gi* hotspot results between KNN k=8 (main) and k=4 (sensitivity).

Input:
  ../data/hanoi_grid_250m_equity_metrics.gpkg  (Step 14 output)

Outputs:
  ../outputs/step14_risk_hotspots/gi_sensitivity_table.csv
  ../outputs/step14_risk_hotspots/gi_sensitivity_table.md
  ../outputs/step14_risk_hotspots/gi_sensitivity_maps_[var].png

Purpose:
  Provide a transparent comparison table and side-by-side maps for the
  thesis defense, showing that thermal and gap hotspot patterns are stable
  across k=4 and k=8 neighbourhood specifications.
"""

import os
import sys

import geopandas as gpd
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

sys.stdout.reconfigure(encoding="utf-8")

# ------------------------------------------------------------------
# Paths
# ------------------------------------------------------------------

DATA_DIR   = os.path.join(os.path.dirname(__file__), "..", "data")
OUTPUT_DIR = os.path.join(os.path.dirname(__file__), "..", "outputs",
                          "step14_risk_hotspots")

EQUITY_PATH = os.path.join(DATA_DIR, "hanoi_grid_250m_equity_metrics.gpkg")

os.makedirs(OUTPUT_DIR, exist_ok=True)


# ------------------------------------------------------------------
# Load
# ------------------------------------------------------------------

def load_grid() -> gpd.GeoDataFrame:
    if not os.path.exists(EQUITY_PATH):
        raise FileNotFoundError(
            f"Equity metrics GPKG not found: {EQUITY_PATH}\n"
            "Run step14_risk_hotspots.py first."
        )
    try:
        gdf = gpd.read_file(EQUITY_PATH, layer="grid_equity_metrics")
    except Exception:
        gdf = gpd.read_file(EQUITY_PATH)

    print(f"Loaded {len(gdf)} cells from equity metrics.")
    return gdf


# ------------------------------------------------------------------
# Hotspot count table
# ------------------------------------------------------------------

HOTSPOT_LABEL = {1: "Hotspot", -1: "Coldspot", 0: "Not significant"}
HOTSPOT_COLOUR = {1: "#d73027", -1: "#4575b4", 0: "#f7f7f7"}


def count_table(gdf: gpd.GeoDataFrame) -> pd.DataFrame:
    """Build comparison table of hotspot counts for k=8 (main) and k=4 (sensitivity)."""
    records = []

    for var in ("lst_p75_filled", "dist_accessible_A_capped_m", "relative_gap_A"):
        for k in (8, 4):
            col = f"{var}_hotspot_k{k}"
            if col not in gdf.columns:
                print(f"WARNING: column {col} not found — skipping.")
                continue

            s = pd.to_numeric(gdf[col], errors="coerce").fillna(0)
            n_total   = int(s.notna().sum())
            n_hot     = int((s == 1).sum())
            n_cold    = int((s == -1).sum())
            n_neutral = int((s == 0).sum())

            records.append({
                "Variable" : var,
                "k"        : k,
                "Role"     : "Primary" if k == 8 else "Sensitivity",
                "Total cells": n_total,
                "Hotspots (FDR)" : n_hot,
                "Coldspots (FDR)": n_cold,
                "Not significant": n_neutral,
                "Hotspot %"      : round(n_hot  / n_total * 100, 1) if n_total else np.nan,
                "Coldspot %"     : round(n_cold / n_total * 100, 1) if n_total else np.nan,
            })

    df = pd.DataFrame(records)
    return df


def spearman_z_scores(gdf: gpd.GeoDataFrame) -> pd.DataFrame:
    """Compute Spearman rank correlation of Gi* z-scores between k=8 and k=4."""
    records = []
    for var in ("lst_p75_filled", "dist_accessible_A_capped_m", "relative_gap_A"):
        z8_col = f"{var}_gi_z_k8"
        z4_col = f"{var}_gi_z_k4"
        if z8_col not in gdf.columns or z4_col not in gdf.columns:
            print(f"WARNING: z-score columns for {var} not found.")
            continue

        z8 = pd.to_numeric(gdf[z8_col], errors="coerce")
        z4 = pd.to_numeric(gdf[z4_col], errors="coerce")

        mask = z8.notna() & z4.notna()
        rho = z8[mask].corr(z4[mask], method="spearman")
        records.append({
            "Variable": var,
            "Spearman rho (z-score k8 vs k4)": round(rho, 4),
            "Stable (rho >= 0.90)": "YES" if rho >= 0.90 else "NO",
        })

    return pd.DataFrame(records)


# ------------------------------------------------------------------
# Maps
# ------------------------------------------------------------------

def _hotspot_series(gdf: gpd.GeoDataFrame, col: str) -> pd.Series:
    if col not in gdf.columns:
        return pd.Series(0, index=gdf.index)
    return pd.to_numeric(gdf[col], errors="coerce").fillna(0)


def plot_side_by_side(
    gdf: gpd.GeoDataFrame,
    variable: str,
    out_path: str,
) -> None:
    """Plot k=8 (main) and k=4 (sensitivity) hotspot maps side by side."""

    col_k8 = f"{variable}_hotspot_k8"
    col_k4 = f"{variable}_hotspot_k4"

    s8 = _hotspot_series(gdf, col_k8)
    s4 = _hotspot_series(gdf, col_k4)

    colour_map = {1: "#d73027", -1: "#4575b4", 0: "#e8e8e8"}

    gdf = gdf.copy()
    gdf["_c8"] = s8.map(colour_map).fillna("#e8e8e8")
    gdf["_c4"] = s4.map(colour_map).fillna("#e8e8e8")

    fig, axes = plt.subplots(1, 2, figsize=(14, 7), constrained_layout=True)
    fig.suptitle(
        f"Getis-Ord Gi* Hotspot Comparison — {variable}\n"
        "Left: KNN k=8 (primary)   Right: KNN k=4 (sensitivity check)",
        fontsize=12, fontweight="bold",
    )

    for ax, col_name, colour_col, k_label in [
        (axes[0], col_k8, "_c8", "k=8 (Primary)"),
        (axes[1], col_k4, "_c4", "k=4 (Sensitivity)"),
    ]:
        # Add ultra-thin white border to separate 250m urban grid structure
        gdf.plot(
            ax=ax,
            color=gdf[colour_col],
            linewidth=0.05,
            edgecolor="#ffffff",
        )
        ax.set_title(k_label, fontsize=11)
        ax.set_axis_off()

        # Counts
        s = _hotspot_series(gdf, col_name)
        n_h = int((s == 1).sum())
        n_c = int((s == -1).sum())
        ax.text(
            0.02, 0.02,
            f"Hotspots: {n_h:,}\nColdspots: {n_c:,}",
            transform=ax.transAxes,
            fontsize=9,
            verticalalignment="bottom",
            bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8),
        )

    legend_patches = [
        mpatches.Patch(color="#d73027", label="Hotspot (FDR p<0.05)"),
        mpatches.Patch(color="#4575b4", label="Coldspot (FDR p<0.05)"),
        mpatches.Patch(color="#e8e8e8", label="Not significant"),
    ]
    fig.legend(
        handles=legend_patches,
        loc="lower center",
        ncol=3,
        fontsize=9,
        framealpha=0.9,
    )

    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Saved map: {out_path}")


# ------------------------------------------------------------------
# Markdown table helper
# ------------------------------------------------------------------

def df_to_md(df: pd.DataFrame) -> str:
    headers = list(df.columns)
    lines = [
        "| " + " | ".join(str(h) for h in headers) + " |",
        "| " + " | ".join(["---"] * len(headers)) + " |",
    ]
    for _, row in df.iterrows():
        lines.append("| " + " | ".join(str(row[c]) for c in headers) + " |")
    return "\n".join(lines)


# ------------------------------------------------------------------
# Main
# ------------------------------------------------------------------

def main() -> None:
    print("=== Gi* Sensitivity Comparison ===\n")

    gdf = load_grid()

    # -- Count table
    counts = count_table(gdf)
    print("\nHotspot count comparison (k=8 vs k=4):")
    print(counts.to_string(index=False))

    # -- Spearman z-score stability
    rho_df = spearman_z_scores(gdf)
    print("\nSpearman rho of Gi* z-scores (k=8 vs k=4):")
    print(rho_df.to_string(index=False))

    # -- Save table
    counts_path = os.path.join(OUTPUT_DIR, "gi_sensitivity_table.csv")
    counts.to_csv(counts_path, index=False, encoding="utf-8-sig")
    print(f"\nSaved: {counts_path}")

    rho_path = os.path.join(OUTPUT_DIR, "gi_sensitivity_zscore_rho.csv")
    rho_df.to_csv(rho_path, index=False, encoding="utf-8-sig")
    print(f"Saved: {rho_path}")

    # -- Markdown
    md_path = os.path.join(OUTPUT_DIR, "gi_sensitivity_table.md")
    with open(md_path, "w", encoding="utf-8") as f:
        f.write("# Getis-Ord Gi* Sensitivity: k=8 (Primary) vs k=4 (Sensitivity)\n\n")
        f.write(
            "> Overlap statistic is computed among cells in the Step 14 valid sample\n"
            "> (non-edge, area ≥ 31,250 m², valid LST/Gap variables). Not the full study grid.\n\n"
        )
        f.write("## Hotspot / Coldspot Counts\n\n")
        f.write(df_to_md(counts))
        f.write("\n\n")
        f.write("## Spearman Rank Correlation of Gi* Z-scores\n\n")
        f.write(df_to_md(rho_df))
        f.write("\n\n")
        f.write(
            "**Interpretation:** If Spearman rho ≥ 0.90 between k=8 and k=4 z-scores,\n"
            "the spatial pattern of high/low clustering is stable across neighbourhood\n"
            "size choices, and conclusions do not depend on the specific k value.\n"
        )
    print(f"Saved: {md_path}")

    # -- Maps (only if geometry present)
    if hasattr(gdf, "geometry") and gdf.geometry.notna().any():
        for var in ("lst_p75_filled", "dist_accessible_A_capped_m", "relative_gap_A"):
            map_path = os.path.join(OUTPUT_DIR, f"gi_sensitivity_maps_{var.lower()}.png")
            try:
                plot_side_by_side(gdf, var, map_path)
            except Exception as e:
                print(f"WARNING: Could not generate map for {var}: {e}")
    else:
        print("No geometry found; skipping maps.")

    print("\nDONE.")


if __name__ == "__main__":
    main()
