import os
import sys
import numpy as np
import pandas as pd

# Add the parent directory to sys.path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))

from scripts.reporting import constants as const
from scripts.reporting.data_loader import load_clean_data

def compile_table_2():
    print("Loading clean data...")
    gdf = load_clean_data()
    
    # Identify Top 20% LST hotspot threshold
    lst_valid = gdf.loc[gdf["population"] > 0, "lst_temp"].dropna()
    if not lst_valid.empty:
        lst_threshold = lst_valid.quantile(0.80)
        gdf["top_20_lst"] = ((gdf["population"] > 0) & (gdf["lst_temp"] >= lst_threshold)).astype(int)
    else:
        gdf["top_20_lst"] = 0
        
    summary_rows = []
    
    for group in const.DENSITY_ORDER:
        sub = gdf[gdf["density_group"] == group].copy()
        if len(sub) == 0:
            continue
            
        total_pop = float(sub["population"].sum())
        pop_hotspot = sub.loc[sub["top_20_lst"] == 1, "population"].sum()
        pct_hotspot = (pop_hotspot / total_pop * 100) if total_pop > 0 else np.nan
        
        # Calculate Conversion Ratio for the stratum as an aggregate
        total_phys_m2 = sub["physical_gb_m2"].sum()
        total_acc_m2 = sub["accessible_gb_A_m2"].sum()
        agg_conversion = (total_acc_m2 / total_phys_m2 * 100) if total_phys_m2 > 0 else np.nan
        
        # Network distance unreachable
        valid_dists = sub["dist_to_A_m"].replace([np.inf, -np.inf], np.nan)
        unreachable_pct = valid_dists.isna().mean() * 100
        
        row = {
            "Stratum": group,
            "Cells": len(sub),
            "Population": f"{total_pop:,.0f}",
            "Mean LST P75 (°C)": f"{sub['lst_temp'].mean():.2f}",
            "Mean Physical GB (%)": f"{sub['physical_gb_pct'].mean():.2f}%",
            "Mean Accessible UGBS (%)": f"{sub['accessible_A_pct'].mean():.2f}%",
            "Mean Gap (%)": f"{sub['physical_access_gap_pct'].mean():.2f}%",
            "Conversion Ratio (%)": f"{agg_conversion:.2f}%",
            "Mean Reachable Distance (m)": f"{valid_dists.mean():.1f}",
            "Unreachable Cells (%)": f"{unreachable_pct:.1f}%",
            "Pop in LST Hotspot (%)": f"{pop_hotspot:,.0f} ({pct_hotspot:.1f}%)"
        }
        summary_rows.append(row)
        
    df_summary = pd.DataFrame(summary_rows)
    
    out_csv = os.path.join(const.TABLES_DIR, "Table2_summary_by_density.csv")
    df_summary.to_csv(out_csv, index=False, encoding="utf-8-sig")
    
    print("\n=== Table 2: Summary by Density Class ===")
    print(df_summary.to_markdown(index=False))
    print(f"\nSaved to: {out_csv}")

if __name__ == "__main__":
    compile_table_2()
