import os
import sys
import geopandas as gpd
import pandas as pd
import config

sys.stdout.reconfigure(encoding='utf-8')

def update_admin_data():
    print("--- STEP 0: UPDATING ADMIN DATA WITH PRECISE CENSUS (dshn.csv) ---")
    
    csv_path = os.path.join(config.DATA_DIR, "dshn.csv")
    old_gpkg = os.path.join(config.DATA_DIR, "hanoi_admin_126.gpkg")
    new_gpkg = os.path.join(config.DATA_DIR, "hanoi_admin_126_updated.gpkg")
    
    if not os.path.exists(csv_path):
        raise FileNotFoundError(f"Missing {csv_path}")
    if not os.path.exists(old_gpkg):
        raise FileNotFoundError(f"Missing {old_gpkg}")
        
    print("Loading datasets...")
    df_csv = pd.read_csv(csv_path, dtype={'maxa': str})
    gdf_admin = gpd.read_file(old_gpkg)
    
    # Ensure string match by padding commune codes (e.g., '00331')
    if 'maxa' in gdf_admin.columns:
        gdf_admin['maxa'] = gdf_admin['maxa'].astype(str).str.zfill(5)
    else:
        print("WARNING: 'maxa' column not found in GPKG. Trying to join by name as fallback.")
    
    df_csv['maxa'] = df_csv['maxa'].astype(str).str.zfill(5)
    
    # Merge data
    print("Merging data...")
    if 'maxa' in gdf_admin.columns:
        gdf_merged = gdf_admin.merge(df_csv[['maxa', 'danso', 'st_area_sh']], on='maxa', how='left', suffixes=('_old', ''))
    else:
        # Fallback if no maxa in GPKG
        gdf_admin['xaphuong_clean'] = gdf_admin[config.ADMIN_NAME_FIELD].str.lower().str.strip()
        df_csv['xaphuong_clean'] = df_csv['xaphuong'].str.lower().str.strip()
        gdf_merged = gdf_admin.merge(df_csv[['xaphuong_clean', 'danso', 'st_area_sh']], on='xaphuong_clean', how='left', suffixes=('_old', ''))
        gdf_merged.drop(columns=['xaphuong_clean'], inplace=True)
        
    # Check for missing joins
    missing = gdf_merged['danso'].isna().sum()
    if missing > 0:
        print(f"WARNING: {missing} communes could not be matched with dshn.csv.")
    
    # Remove old columns if they exist
    if 'danso_old' in gdf_merged.columns:
        gdf_merged.drop(columns=['danso_old'], inplace=True)
    if 'st_area_sh_old' in gdf_merged.columns:
        gdf_merged.drop(columns=['st_area_sh_old'], inplace=True)
        
    print(f"New Total Population in GPKG: {gdf_merged['danso'].sum():,.0f}")
    
    gdf_merged.to_file(new_gpkg, driver="GPKG")
    print(f"SUCCESS: Saved updated admin boundaries to {new_gpkg}")

if __name__ == "__main__":
    update_admin_data()
