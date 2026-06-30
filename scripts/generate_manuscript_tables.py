"""
Generate CSV tables for the manuscript to accompany or replace figures.
"""

import os
import pandas as pd
import geopandas as gpd
import config

OUTPUT_DIR = os.path.join(config.OUTPUT_DIR, "tables")
os.makedirs(OUTPUT_DIR, exist_ok=True)

def generate_table2(gdf: gpd.GeoDataFrame):
    """Table 2: Summary statistics by density class"""
    print("Generating Table 2: Summary Statistics...")
    
    # Calculate means grouped by density group
    grouped = gdf.groupby('density_group')
    
    table = pd.DataFrame()
    table['Physical Green-Blue Cover (%)'] = grouped['physical_greenblue_pct'].mean()
    table['Accessible Level A Cover (%)'] = grouped['accessible_A_pct'].mean() if 'accessible_A_pct' in gdf.columns else grouped['accessible_area_A'].sum() / grouped['area_m2'].sum() * 100
    table['Physical-Access Gap (%)'] = grouped['physical_access_gap_pct'].mean()
    
    # Conversion ratio: Total accessible A area / Total physical area
    table['Conversion Ratio (%)'] = (grouped['accessible_area_A'].sum() / grouped['physical_gb_area'].sum()) * 100
    
    # Population weighted mean LST
    def pop_weighted_mean(df, val_col, pop_col):
        v = df[val_col]
        p = df[pop_col].fillna(0)
        return (v * p).sum() / p.sum() if p.sum() > 0 else v.mean()

    table['Mean LST (°C)'] = grouped.apply(lambda x: pop_weighted_mean(x, 'lst_p75_filled', 'population'))
    
    # Sort index appropriately
    order = ['Low-density', 'Medium-density', 'High-density']
    table = table.reindex(order)
    
    out_path = os.path.join(OUTPUT_DIR, "table2_summary_statistics_by_density.csv")
    table.to_csv(out_path)
    print(f"Saved: {out_path}")

def generate_fig3_table(gdf: gpd.GeoDataFrame):
    """Inline table for Fig 3: % population within distance thresholds"""
    print("Generating Fig 3 Inline Table...")
    
    total_pop = gdf['population'].sum()
    if total_pop == 0:
        return
        
    results = []
    thresholds = [300, 500, 1000]
    
    for level, dist_col in [("Level A", "dist_to_accessible_A"), ("Level A+B", "dist_to_accessible_AB")]:
        row = {'Level': level}
        for t in thresholds:
            pop_within = gdf.loc[gdf[dist_col] <= t, 'population'].sum()
            row[f'<= {t}m (%)'] = (pop_within / total_pop) * 100
        results.append(row)
        
    df = pd.DataFrame(results)
    out_path = os.path.join(OUTPUT_DIR, "table_fig3_threshold_percentages.csv")
    df.to_csv(out_path, index=False)
    print(f"Saved: {out_path}")

def generate_fig6_table(gdf: gpd.GeoDataFrame):
    """Inline table for Fig 6: Top wards by priority cell count"""
    print("Generating Fig 6 Inline Table: Top Wards by Priority...")
    
    if 'priority_class' not in gdf.columns:
        cls = pd.Series("Not priority", index=gdf.index, dtype="object")
        is_network_priority = gdf.get('Priority_Area_Network', 0) == 1
        is_areagap_priority = gdf.get('Priority_Area_AreaGap', 0) == 1
        
        cls.loc[is_network_priority] = "Network priority"
        cls.loc[is_areagap_priority] = "Area-gap priority"
        cls.loc[is_network_priority & is_areagap_priority] = "Both priority types"
        gdf['priority_class'] = cls
        
    priority_only = gdf[gdf['Priority_Area_Network'] == 1].copy()
    
    # Count cells and sum population per commune
    grouped = priority_only.groupby('commune_name').agg({
        'grid_id': 'count',
        'population': 'sum',
    }).rename(columns={'grid_id': 'Priority Cells', 'population': 'Affected Population'})
    
    # Sort by number of priority cells
    grouped = grouped.sort_values('Priority Cells', ascending=False).head(10)
    
    out_path = os.path.join(OUTPUT_DIR, "table_fig6_top_wards.csv")
    grouped.to_csv(out_path)
    print(f"Saved: {out_path}")

if __name__ == "__main__":
    print("Loading clean data...")
    import sys
    sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
    from scripts.reporting.data_loader import load_clean_data
    gdf = load_clean_data()
    
    # Map cleaned columns to the names expected by the legacy functions
    if 'physical_greenblue_pct' not in gdf.columns and 'physical_gb_pct' in gdf.columns:
        gdf['physical_greenblue_pct'] = gdf['physical_gb_pct']
    if 'physical_gb_area' not in gdf.columns and 'physical_gb_m2' in gdf.columns:
        gdf['physical_gb_area'] = gdf['physical_gb_m2']
    if 'accessible_area_A' not in gdf.columns and 'accessible_gb_A_m2' in gdf.columns:
        gdf['accessible_area_A'] = gdf['accessible_gb_A_m2']
    if 'dist_to_accessible_A' not in gdf.columns and 'dist_to_A_m' in gdf.columns:
        gdf['dist_to_accessible_A'] = gdf['dist_to_A_m']
        
    # dist_to_accessible_AB mapping
    if 'dist_to_accessible_AB' not in gdf.columns:
        dist_candidates_ab = ["dist_to_accessible_AB_m", "dist_accessible_AB_m", "nearest_accessible_AB_m"]
        col_ab = next((c for c in dist_candidates_ab if c in gdf.columns), None)
        if col_ab:
            gdf['dist_to_accessible_AB'] = pd.to_numeric(gdf[col_ab], errors='coerce')
        
    generate_table2(gdf)
    generate_fig3_table(gdf)
    generate_fig6_table(gdf)
    print("Done!")
