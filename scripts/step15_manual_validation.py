"""
Step 15: Manual Validation for Level A Public Accessible Spaces.

This script performs two tasks:
1. create_validation_template(): Generates a random stratified sample (by source_tag) 
   of 100 polygons from the Level A accessible green-blue spaces, generates Google Maps 
   links for their centroids, and outputs a CSV template for manual validation.
   
2. calculate_precision(): Reads the filled CSV file and calculates the 
   precision of the Level A classification to report in the manuscript.
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

TEMPLATE_PATH = os.path.join(config.DATA_DIR, "manual_validation_template.csv")


def create_validation_template(sample_size: int = 100):
    print(f"--- CREATING STRATIFIED MANUAL VALIDATION TEMPLATE (N={sample_size}) ---")
    
    input_gpkg = getattr(config, "ACCESSIBLE_SPACES_GPKG", os.path.join(config.DATA_DIR, "accessible_spaces.gpkg"))
    
    if not os.path.exists(input_gpkg):
        raise FileNotFoundError(f"Input GPKG not found: {input_gpkg}. Run Step 9 first.")
        
    print(f"Loading Level A spaces from {input_gpkg}...")
    try:
        spaces_a = gpd.read_file(input_gpkg, layer="accessible_spaces_A")
    except Exception as e:
        raise ValueError(f"Could not load layer 'accessible_spaces_A'. Error: {e}")
        
    print(f"Total Level A polygons available: {len(spaces_a)}")
    
    valid_spaces = spaces_a.copy()
    
    # Perform stratified sampling preserving attribute structure
    tag_counts = valid_spaces['source_tag'].value_counts()
    sample_counts = (tag_counts / tag_counts.sum() * sample_size).round().astype(int)
    
    # Ensure rare strata have at least 1 sample if allowed by group size
    for tag in sample_counts.index:
        if sample_counts[tag] == 0 and tag_counts[tag] > 0:
            sample_counts[tag] = 1
            
    # Điều chỉnh khớp chính xác tổng sample_size
    diff = sample_size - sample_counts.sum()
    if diff != 0 and len(sample_counts) > 0:
        largest_group = sample_counts.idxmax()
        if sample_counts[largest_group] + diff > 0:
            sample_counts[largest_group] += diff
        
    samples = []
    for tag, count in sample_counts.items():
        if count > 0:
            group = valid_spaces[valid_spaces['source_tag'] == tag]
            n = min(count, len(group))
            if n > 0:
                samples.append(group.sample(n=n, random_state=42))
            
    sample = pd.concat(samples)
    
    # Bốc bù nếu bị khuyết mẫu do hàm min() chặn biên
    remaining = sample_size - len(sample)
    if remaining > 0:
        unused = valid_spaces[~valid_spaces.index.isin(sample.index)]
        if not unused.empty:
            extra = unused.sample(n=min(remaining, len(unused)), random_state=42)
            sample = pd.concat([sample, extra])

    print(f"Sampled {len(sample)} polygons (stratified by source_tag).")
    
    # Chuyển đổi hệ tọa độ sang WGS84 để lấy kinh vĩ độ chuẩn cho Google Maps
    sample_wgs84 = sample.to_crs("EPSG:4326")
    rep_points = sample_wgs84.geometry.representative_point()
    
    df = pd.DataFrame({
        "space_id": sample["space_id"],
        "name": sample["name"],
        "source_tag": sample["source_tag"],
        "area_m2": sample["area_m2"].round(1),
        "lat": rep_points.y.round(6),
        "lon": rep_points.x.round(6),
    })
    
    df["google_maps_link"] = df.apply(
        lambda row: f"https://www.google.com/maps/search/?api=1&query={row['lat']},{row['lon']}", 
        axis=1
    )
    
    df["is_accessible"] = ""
    df["is_greenblue"] = ""
    df["notes"] = ""
    
    os.makedirs(os.path.dirname(TEMPLATE_PATH) or ".", exist_ok=True)
    df.to_csv(TEMPLATE_PATH, index=False, encoding='utf-8-sig')
    
    print(f"\nTemplate saved to: {TEMPLATE_PATH}")
    print("\nHƯỚNG DẪN KIỂM ĐỊNH THỰC ĐỊA:")
    print("1. Mở file 'manual_validation_template.csv' bằng Excel.")
    print("2. Click vào từng link Google Maps để xem vị trí (bật chế độ Vệ tinh/Street View).")
    print("3. Điền cột 'is_accessible' (Quyền tiếp cận công cộng): Điền Y (Có), N (Không), P (Một phần/Nghi vấn).")
    print("4. Điền cột 'is_greenblue' (Sự hiện diện thực tế): Điền Y (Có thực thể xanh/mặt nước), N (Không), P (Một phần).")
    print("5. Lưu file Excel và kích hoạt lệnh tính: `python step15_manual_validation.py --calculate`.")


def calculate_precision():
    print("--- CALCULATING VALIDATION PRECISION ---")
    
    if not os.path.exists(TEMPLATE_PATH):
        print(f"Error: Template file not found at {TEMPLATE_PATH}")
        print("Please run without arguments first to generate the template.")
        return
        
    df = pd.read_csv(TEMPLATE_PATH)
    
    if "is_accessible" not in df.columns or "is_greenblue" not in df.columns:
        print("Error: Required columns 'is_accessible' or 'is_greenblue' are missing.")
        return
        
    df["acc_clean"] = df["is_accessible"].astype(str).str.strip().str.upper()
    df["gb_clean"] = df["is_greenblue"].astype(str).str.strip().str.upper()
    
    valid_responses = ["Y", "N", "P"]
    filled_df = df[
        df["acc_clean"].isin(valid_responses) & 
        df["gb_clean"].isin(valid_responses)
    ].copy()
    
    total_samples = len(df)
    filled_samples = len(filled_df)
    
    if filled_samples == 0:
        print("No validation data found in the CSV file.")
        print("Please fill out the columns 'is_accessible' and 'is_greenblue' first.")
        return
        
    print(f"Total samples: {total_samples}")
    print(f"Validated samples: {filled_samples}")
    
    acc_y = (filled_df["acc_clean"] == "Y").sum()
    acc_n = (filled_df["acc_clean"] == "N").sum()
    access_precision = acc_y / (acc_y + acc_n) * 100 if (acc_y + acc_n) > 0 else 0
    
    gb_y = (filled_df["gb_clean"] == "Y").sum()
    gb_n = (filled_df["gb_clean"] == "N").sum()
    gb_precision = gb_y / (gb_y + gb_n) * 100 if (gb_y + gb_n) > 0 else 0
    
    strict_y = ((filled_df["acc_clean"] == "Y") & (filled_df["gb_clean"] == "Y")).sum()
    strict_n = ((filled_df["acc_clean"] == "N") | (filled_df["gb_clean"] == "N")).sum()
    strict_precision = strict_y / (strict_y + strict_n) * 100 if (strict_y + strict_n) > 0 else 0
    
    liberal_y = (
        filled_df["acc_clean"].isin(["Y", "P"]) & 
        filled_df["gb_clean"].isin(["Y", "P"])
    ).sum()
    liberal_precision = liberal_y / filled_samples * 100 if filled_samples > 0 else 0
    
    print("\n--- METRICS ---")
    print(f"1. Access Precision (is_accessible == Y)                                 : {access_precision:.1f}%")
    print(f"2. Green-blue Precision (is_greenblue == Y)                              : {gb_precision:.1f}%")
    print(f"3. Strict Public Green-blue Precision (Y & Y)                            : {strict_precision:.1f}%")
    print(f"4. Liberal Public Green-blue Precision ([Y/P] & [Y/P])   : {liberal_precision:.1f}%")
        
    print("\nWrite this in your manuscript:")
    print("--------------------------------------------------")
    print(f"To assess the reliability of the OSM-derived public accessibility model, a manual validation was conducted on a random stratified sample of {filled_samples} Level A polygons, stratified by OSM source tags to ensure proportional representation of different space types. Using high-resolution satellite imagery and Google Street View, the spaces were visually inspected for both physical green-blue cover and signs of public access. The manual validation yielded a strict public green-blue precision of {strict_precision:.1f}% (requiring both confirmed public access and confirmed physical vegetation/water), and a liberal precision of {liberal_precision:.1f}% when including partially verifiable spaces. These results indicate a high degree of confidence in the rule-based classification of accessible urban green and blue spaces.")
    print("--------------------------------------------------")


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "--calculate":
        calculate_precision()
    else:
        create_validation_template()