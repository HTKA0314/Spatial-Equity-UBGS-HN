# PHƯƠNG PHÁP NGHIÊN CỨU — GHI CHÉP CHÍNH XÁC
# Đối chiếu 100% với source code (không bịa, không nói quá)
# Cập nhật: 2026-06-17

---

## A. TỔNG QUAN PIPELINE

Nghiên cứu thực hiện tuần tự qua 15 bước (step1 → step15), chia thành 5 nhóm phân tích:

1. Xây dựng lưới không gian + LST (Step 1, 2)
2. Dữ liệu phụ trợ: nước, dân số, mạng lưới đường (Step 4, 5, 6)
3. Phân loại lớp phủ + không gian xanh tiếp cận được (Step 8, 9, 10)
4. Khả năng tiếp cận mạng lưới (Step 12)
5. Hồi quy không gian + Hotspot + Ưu tiên (Step 13, 13b, 14)

---

## B. CÁC PHƯƠNG PHÁP — CHI TIẾT CHÍNH XÁC

---

### B1. PHÂN LOẠI LỚP PHỦ ĐẤT (LAND COVER CLASSIFICATION)

**Phương pháp:** Random Forest Classification trên Google Earth Engine

**Source:** `gee_export_lst_builtup.js`

**Dữ liệu đầu vào:**
- Sentinel-2 Level-2A Harmonized (COPERNICUS/S2_SR_HARMONIZED)
- Thời gian: 2024-05-01 → 2024-08-31
- Lọc mây: CLOUDY_PIXEL_PERCENTAGE < 45%
- Che mây pixel: SCL mask (loại SCL = 0,1,2,3,8,9,10,11)
- Composite: Per-pixel **median** (không phải mean)

**Features đầu vào (inputBands):**
- 10 kênh phổ: B2, B3, B4, B5, B6, B7, B8, B8A, B11, B12
- 7 chỉ số: NDVI, NDWI, MNDWI, NDBI, NDRE, SAVI, BSI
- Tổng cộng: **17 bands**

**6 lớp phủ:**
| Code | Nhãn |
|:---:|:---|
| 0 | Water |
| 1 | Tree/Forest |
| 2 | Grass |
| 3 | Cropland |
| 4 | Built-up/Impervious |
| 5 | Bare soil |

**Tham số RF:**
- numberOfTrees: **300**
- bagFraction: **0.7**
- seed: **42**

**Training/Validation:**
- Split tại cấp **polygon** (70% train / 30% test), random seed = 42
- Cân bằng mẫu: tối đa **2,000 pixel/lớp train**, **700 pixel/lớp test**
- Tổng ≈ 12,000 train pixels, ≈ 4,200 test pixels

**Đánh giá độ chính xác:** Overall Accuracy + Kappa từ error matrix trên held-out test set
(giá trị thực tế lấy từ GEE console — **chưa ghi vào file nào trong project**)

**Đầu ra Python:**
- `reduceRegions` với `ee.Reducer.mean()` tại scale = 10m → phần trăm từng lớp/ô 250m
- Export CSV → `config.LANDCOVER_CSV`
- Merge vào grid: Step 8 (`step8_merge_landcover.py`)

---

### B2. LAND SURFACE TEMPERATURE (LST)

**Phương pháp:** Landsat Level-2 Surface Temperature + Temporal Aggregation + IDW Fill

**Source:** `gee_export_lst_builtup.js` (export từ GEE) + `step2_lst_processing.py`

**Dữ liệu:** Landsat 8/9 Collection 2, Level-2 Surface Temperature
**Thời gian:** 2024-05-01 → 2024-08-31

**Metric:** **75th percentile** (`lst_p75_observed`) của tất cả pixel LST hợp lệ trong ô lưới 250m trong khoảng thời gian trên

**Hai biến LST được tạo ra:**

| Biến | Mô tả | Dùng cho |
|:---|:---|:---|
| `lst_p75_observed` | Giá trị thực đo, ô nào không đủ ảnh = NaN | **Hồi quy không gian** |
| `lst_p75_filled` | IDW fill NaN từ `lst_p75_observed` | **Gi* hotspot**, bản đồ trực quan |

**IDW Interpolation:**
- Thư viện: `scipy.spatial.cKDTree`
- k = **8** nearest neighbours
- power p = **2.0**
- Chỉ fill ô có `lst_p75_observed = NaN`

**Lý do phân tách:** Hồi quy chỉ dùng `lst_p75_observed` để tránh suy luận vòng; Gi* dùng `lst_p75_filled` để đảm bảo tính liên tục không gian.

---

### B3. DỮ LIỆU DÂN SỐ

**Phương pháp:** WorldPop + Commune-level Calibration

**Source:** `step5_population_mapping.py`

**Dữ liệu gốc:** WorldPop unconstrained, Vietnam 2025, ~100m resolution

**Quy trình:**
1. Zonal summation: tổng pop từ raster vào từng ô lưới 250m
2. **Commune-level calibration:** Với mỗi xã/phường, tính factor = (dân số census chính thức) / (tổng WorldPop raw trong xã đó) → nhân đều cho mọi ô trong xã
3. Ô lưới có > **80%** diện tích là JRC water → gán pop = 0

**Dữ liệu census:** File admin `config.ADMIN_PATH` với hai cột: `commune_name` và `pop`

**Calibration fallback:** Nếu không match được tên xã → dùng global scale factor hoặc giữ nguyên = 1.0 (xem `_calibration_factors_by_commune()`)

---

### B4. PHÂN LOẠI TIẾP CẬN KHÔNG GIAN XANH-XAnh (UGBS ACCESSIBILITY CLASSIFICATION)

**Phương pháp:** OSM Tag-based Classification (3-tier)

**Source:** `step9_accessible_spaces.py`

**Dữ liệu:** OpenStreetMap features (download qua OSMnx)
**Query area:** Buffered bounding box của toàn bộ grid + **500m buffer** → clip về ranh giới chính xác sau khi download

**Tags được query (OSM_TAGS):**
- leisure: park, garden, recreation_ground, playground, nature_reserve, sports_centre, pitch
- landuse: grass, recreation_ground, village_green, forest, cemetery, reservoir, basin
- natural: wood, grassland, scrub, heath, wetland, water
- water: lake, pond, reservoir, river, basin, canal
- waterway: river, canal, stream, drain, ditch
- amenity: grave_yard
- boundary: protected_area
- tourism: picnic_site

**Ba mức phân loại:**

**Level A — Công cộng rõ ràng (de jure):**
- leisure ∈ {park, playground, recreation_ground}
- landuse ∈ {village_green, recreation_ground}
- tourism = picnic_site
- access ∈ {yes, public, designated} hoặc foot ∈ {yes, public, designated}

**Level B — Bán công cộng (de facto):**
- leisure ∈ {garden, nature_reserve, sports_centre, pitch}
- landuse ∈ {cemetery, grass, forest}
- natural ∈ {wood, grassland, scrub, heath, wetland}
- amenity = grave_yard
- boundary = protected_area
- access = permissive hoặc foot = permissive

**Level C — Hạn chế/riêng tư (không tính):**
- access ∈ {private, no, customers, permit, delivery, destination}
- foot ∈ {no, private}
- natural = water (mặc định) — không có access tag rõ ràng
- landuse ∈ {reservoir, basin} không có access tag
- **Toàn bộ JRC water polygon** → Level C mặc định

**Ưu tiên rule:** Explicit restriction override tất cả (dù leisure=park, nếu access=private → Level C)

**Bộ lọc diện tích:** Polygon < **100 m²** bị loại

**Waterline buffer:** Linestring waterway → buffer **5.0m** → chuyển thành polygon trước khi filter

**Dissolve trước khi tính diện tích giao với ô lưới:** Tránh double-counting (một ô lưới có thể bị tính nhiều lần nếu có 2 polygon overlapping)

**Đầu ra (GPKG layers):**
- `accessible_spaces_all`
- `accessible_spaces_A`
- `accessible_spaces_AB`
- `jrc_water_reference`
- `accessible_spaces_A_dissolved`
- `accessible_spaces_AB_dissolved`

---

### B5. ĐO LƯỜNG DIỆN TÍCH TIẾP CẬN TRONG Ô LƯỚI

**Phương pháp:** Geometric Intersection + Area Capping

**Source:** `step10_accessibility_mapping.py`

**Quy trình:**
1. Dissolve toàn bộ Level A polygon → single geometry
2. Tính intersection area giữa dissolved geometry và từng ô lưới 250m
3. Cap: `accessible_area ≤ actual_cell_area`
4. Tính fraction: `accessible_gb_A_frac = accessible_area_m2 / cell_area_m2`

**Physical GB fraction:**
`Physical_GB_frac = tree_frac + grass_frac + water_frac`
*(Cropland không tính vào physical GB)*

**Gap metrics:**
- `relative_gap_A = Physical_GB_frac - accessible_gb_A_frac` (đơn vị: percentage points, pp)
- `conversion_ratio_A = (accessible_gb_A_frac / Physical_GB_frac) × 100%`

---

### B6. MẠNG LƯỚI GIAO THÔNG PEDESTRIAN (OSM NETWORK)

**Phương pháp:** OSMnx Graph Download + Undirected Walk Graph

**Source:** `step6_osm_network.py`, `step12_network_accessibility.py`

**Tham số:**
- `network_type = "walk"`
- Buffer: **2,000m** ngoài ranh giới study area
- CRS: EPSG:32648 (UTM Zone 48N)
- Graph type: Undirected (`ox.convert.to_undirected()`)

**Lý do undirected:** Hà Nội không encode one-way restriction cho đường pedestrian trong OSM → undirected defensible hơn

---

### B7. KHOẢNG CÁCH ĐI BỘ MẠNG LƯỚI (NETWORK WALKING ACCESSIBILITY)

**Phương pháp:** Multi-Source Dijkstra Algorithm

**Source:** `step12_network_accessibility.py`

**Thuật toán:** `networkx.multi_source_dijkstra_path_length()`
- Nguồn: Tất cả destination nodes (park entry points)
- Đích: Tất cả origin nodes (ô lưới)
- Cutoff: **1,000m**
- Weight: edge `length` (metres)

**Origin:**
- Mỗi ô lưới: `geometry.representative_point()` (điểm đảm bảo nằm trong ô)
- Snap sang nearest network node: `osmnx.nearest_nodes()`
- Snap warning nếu > **100m** (ghi log, không loại)

**Destination sampling:**
- Mỗi accessible polygon (sau dissolve):
  - 1 representative point
  - Boundary points cứ mỗi **150m** (`DESTINATION_SAMPLE_SPACING_M = 150`)
- Snap sang network node, loại nếu snap > **150m** (`DEST_SNAP_MAX_M = 150`)

**Kết quả:**
- `dist_to_accessible_A_m` = khoảng cách ngắn nhất từ ô lưới đến Level A park
- Ô không reach được trong 1,000m → `np.inf`
- Ô có `accessible_gb_area_A_final_m2 > 0` → gán **0m** (overlapping)

**Threshold flags:**
- `within_300m_A` = (dist ≤ 300) → 1, else 0
- `within_500m_A` = (dist ≤ 500) → 1, else 0
- `within_1000m_A` = (dist ≤ 1000) → 1, else 0
- Ô `np.inf`: tất cả flags = **0**

**Population-weighted accessibility:**
`Acc(d) = Σ(pop_i | dist_i ≤ d) / Σ(pop_i)` — chỉ tính ô có pop > 0

**Access deficit:**
`access_deficit_A = min(dist_to_accessible_A_m / 1000, 1.0)`
- Ô `np.inf` → access_deficit_A = **1.0**
- (Có issue #5 trong audit: có 1,094 ô bị NaN thay vì 1.0 — cần rerun Step 14)

---

### B8. BÌNH ĐẲNG PHÂN PHỐI (EQUITY METRICS)

**Phương pháp:** Lorenz Curve + Gini Coefficient

**Source:** tính từ `accessible_gb_A_frac` và `population`

**Lorenz curve:** Sort ô lưới theo per-capita accessible area tăng dần → tính cumulative pop share vs cumulative accessible area share

**Gini:**
`G = 1 - Σ(x_i - x_{i-1})(y_i + y_{i-1})`
(trapezoidal approximation)

**Reported value:** Gini = **0.979** (từ ALL_RESULTS_COMPILED.md)

---

### B9. HỒI QUY KHÔNG GIAN (SUPPLEMENTARY SPATIAL REGRESSION)

**Phương pháp:** OLS → Anselin-Florax Decision → SLM hoặc SEM

**Source:** `step13_spatial_regression.py`, `step13b_stratified_regression.py`

**Biến phụ thuộc (Y):**
- **Primary:** `lst_p75_observed` (non-interpolated)
- **Sensitivity:** `lst_p75_filled` (Model-1B)

**Ba đặc tả mô hình:**

| Tên | Y | Biến X |
|:---|:---|:---|
| Model-1A-Physical-Observed | lst_p75_observed | tree_frac, grass_frac, crop_frac, water_frac, bare_frac, log_dist_jrc_water |
| Model-1B-Physical-Filled | lst_p75_filled | (same) |
| Model-1C-Physical-RobustVeg | lst_p75_observed | total_vegetation_frac, water_frac, bare_frac, log_dist_jrc_water |

**Biến reference:** `built_frac` bị omit (reference category)

**`total_vegetation_frac` = tree_frac + grass_frac + crop_frac**

**VIF screening:**
- Warn: VIF > **5.0**
- Drop (iterative, highest-first): VIF > **10.0**

**Spatial weights:**
- KNN, k = **8** (primary), k = **4** (sensitivity)
- `libpysal.weights.KNN.from_array()`
- Row-standardised (`transform = "r"`)

**Anselin-Florax Decision Rules** (implemented in `_choose_spatial_model()`):

| Điều kiện | Kết quả | Mô hình báo cáo |
|:---|:---|:---|
| LM-Lag sig, LM-Error không | SLM | SLM |
| LM-Error sig, LM-Lag không | SEM | SEM |
| Cả hai sig: RLM-Lag sig, RLM-Error không | SLM | SLM |
| Cả hai sig: RLM-Error sig, RLM-Lag không | SEM | SEM |
| Cả hai RLM sig ("BOTH") | BOTH → **SLM** (lý do: UHI spillover) | SLM (primary) + SEM (sensitivity) |
| Cả hai không sig | OLS | OLS |

**Ước lượng:**
- SLM: `spreg.GM_Lag` (Kelejian-Prucha GMM 2SLS)
- SEM: `spreg.GM_Error`
- OLS: `spreg.OLS`

**Kiểm tra phần dư:** Moran's I của residuals (`esda.moran.Moran`, permutations=0)

**Loại bỏ ô:**
- `area_m2 < 31,250 m²` (edge cells)
- `EXCLUDE_EDGE_CELLS = True` (xem config)

**Phân tầng (Step 13b):**
- Low-density: pop_density < 5,000 người/km²
- Medium-density: 5,000–15,000 người/km²
- High-density: > 15,000 người/km²
- Cùng KNN k=8; cảnh báo nếu mean KNN dist > 1,000m (vì strata fragmented)
- **Mục đích:** Heterogeneity check, KHÔNG phải thay thế global model

---

### B10. PHÂN TÍCH ĐIỂM NÓNG (GETIS-ORD Gi*)

**Phương pháp:** Local Spatial Statistics

**Source:** `step14_risk_hotspots.py`

**Thư viện:** `esda.getisord.G_Local`

**Tham số:**
- `star = True` (Gi* — bao gồm chính điểm i trong tính toán)
- `transform = "r"` (row-standardised)
- Permutations: **999**

**Spatial weights:** KNN
- Primary: k = **8** (`config.HOTSPOT_K_MAIN`)
- Sensitivity: k = **4** (`config.HOTSPOT_K_SENSITIVITY`)

**Biến input:**
- Gi* cho LST: `lst_p75_filled` (bắt buộc phải continuous, không NaN)
- Gi* cho access deficit: `dist_accessible_A_capped_m` (= min(dist, 1000), unreachable → 1000)
- Gi* cho area gap: `relative_gap_A`

**Hiệu chỉnh multiple testing:**
- `statsmodels.stats.multitest.multipletests(method="fdr_bh")` — Benjamini-Hochberg FDR
- Significance: FDR-adjusted p < **0.05**

**Phân loại:**
- Hotspot (1): z > 0 AND p_adj < 0.05
- Coldspot (-1): z < 0 AND p_adj < 0.05
- Không đáng kể (0): p_adj ≥ 0.05

**Sensitivity check:** So sánh hotspot classification k=8 vs k=4 → tính % đồng thuận (agreement)

---

### B11. BIVARIATE LISA

**Phương pháp:** Bivariate Local Moran's I

**Source:** `step14_risk_hotspots.py`

**Thư viện:** `esda.moran.Moran_Local_BV`

**Tham số:** KNN k=8, permutations=999, FDR correction

**Hai phân tích:**

| BiLISA | X (Access) | Y (Exposure) |
|:---|:---|:---|
| Equity co-location | `access_deficit_A` | `pop_density` (người/km²) |
| Heat co-location | `access_deficit_A` | `lst_p75_filled` |

**Output được dùng:**
- `bilisa_access_deficit_A_x_pop_density_HH_fdr`: Ô có HH cluster (đồng thời access deficit cao VÀ pop density cao), FDR-corrected
- `bilisa_access_deficit_A_x_{lst_col}_HH_fdr`: HH cho access deficit x LST

**Vai trò trong bài:** **Bằng chứng không gian bổ sung** — KHÔNG dùng để định nghĩa Priority Area chính

---

### B12. VÙNG ƯU TIÊN (PRIORITY AREA DEFINITION)

**Source:** `step14_risk_hotspots.py`, lines 396–410

**Định nghĩa chính (`Priority_Area_Network`):**

```python
Priority_Area_Network = (
    ((dist_accessible_A_capped_m > 1000) | dist_accessible_A_capped_m.isna())
    & (density_group == "High-density")      # > 15,000 người/km²
    & (lst_col_hotspot_k8 == 1)              # LST hotspot FDR-corrected
)
```

Ba điều kiện phải **đồng thời** thỏa mãn:
1. `dist_accessible_A_capped_m > 1000` hoặc là `NaN` (unreachable/severe deficit)
2. `density_group == "High-density"` (> 15,000 người/km²)
3. `{lst_col}_hotspot_k8 == 1` (FDR LST hotspot k=8)

**Định nghĩa phụ (`Priority_Area_AreaGap`):**
```python
Priority_Area_AreaGap = (lst_hotspot_k8 == 1) & (relative_gap_A_hotspot_k8 == 1)
```
→ Chỉ dùng như **supplementary diagnostic**, không dùng cho policy recommendation chính

---

## C. NHỮNG GÌ KHÔNG LÀM / KHÔNG CÓ TRONG CODE

**Danh sách các phương pháp KHÔNG có trong bài** (để tránh đề cập trong paper):

- ❌ Scenario simulation / counterfactual analysis — Không thấy step15 làm điều này (step15 là manual validation)
- ❌ Kernel density estimation (KDE) cho park demand
- ❌ Buffer-based accessibility (Euclidean buffer) — Bài chỉ dùng network distance
- ❌ Gravity model accessibility
- ❌ GWR (Geographically Weighted Regression)
- ❌ Spatial interpolation cho land cover
- ❌ Time-series analysis của LST
- ❌ Accuracy assessment của step chính thức từ GEE (OA, Kappa chưa lưu vào file nào trong project)
- ❌ Thống kê Moran's I Global (chỉ có local: Gi* và LISA)

---

## D. ĐỐI CHIẾU NHANH: CLAIM vs CODE

| Claim trong bài | Đúng trong code? | Ghi chú |
|:---|:---:|:---|
| "Sentinel-2 period: May–Aug 2024" | ✅ | `startDate = '2024-05-01'` |
| "Cloud threshold 45%" | ✅ | `CLOUDY_PIXEL_PERCENTAGE < 45` |
| "Composite = median" | ✅ | `s2.median()` |
| "RF: 300 trees" | ✅ | `numberOfTrees: 300` |
| "RF: 70/30 split at polygon level" | ✅ | `randomColumn('poly_random', 42)` |
| "17 input features (bands + indices)" | ✅ | inputBands có 17 items |
| "IDW k=8, p=2.0" | ✅ | `step2: line 129` |
| "lst_p75_filled dùng cho hotspot" | ✅ | `step14: _find_lst_column()` |
| "lst_p75_observed dùng cho regression" | ✅ | `step13: dep_var_obs` |
| "Dijkstra cutoff = 1,000m" | ✅ | `config: DIJKSTRA_CUTOFF_M = 1000` |
| "Unreachable = np.inf (không phải 1000m)" | ✅ | `step12: dist_A.get(node, np.inf)` |
| "Destination sample spacing = 150m" | ✅ | `DESTINATION_SAMPLE_SPACING_M = 150` |
| "Dest snap max = 150m" | ✅ | `DEST_SNAP_MAX_M = 150` |
| "Origin snap warning = 100m" | ✅ | `SNAP_WARNING_THRESHOLD_M = 100` |
| "VIF warn = 5, drop = 10" | ✅ | `step13: lines 76-77` |
| "SLM = GM_Lag (spreg)" | ✅ | `step13: line 46` |
| "Gi* star=True, permutations=999" | ✅ | `step14: G_Local(star=True, permutations=999)` |
| "FDR = Benjamini-Hochberg" | ✅ | `step14: multipletests(method='fdr_bh')` |
| "BiLISA X = access_deficit_A" | ✅ | `step14: Moran_Local_BV` block |
| "Priority: 3 điều kiện AND" | ✅ | `step14: lines 398-404` |
| "Priority: high-density > 15,000" | ✅ | `step14: density_group == 'High-density'` |
| "BiLISA chỉ là supplementary" | ✅ | audit_report.md: Issue #8 |
| "JRC min water area = 1,000 m²" | ✅ | `step9: MIN_WATER_AREA_M2 = 1000` |
| "Min accessible space area = 100 m²" | ✅ | `step9: MIN_ACCESSIBLE_SPACE_AREA_M2 = 100` |
| "Network buffer = 2,000m" | ✅ | `config: NETWORK_BUFFER_M = 2000` |
| "OSM query buffer = 500m" | ✅ | `config: OSM_FEATURE_QUERY_BUFFER_M = 500` |
| "CRS = EPSG:32648" | ✅ | `config: PROJ_CRS` |
| "Undirected walk graph" | ✅ | `step12: _to_undirected_walk_graph()` |
| "Pop calibration commune-level" | ✅ | `step5: _calibration_factors_by_commune()` |
| "Pop zero if >80% water" | ✅ | `config: WATER_CELL_ZERO_POP_THRESHOLD = 0.8` |
| "13,956 valid cells" | ✅ | ALL_RESULTS_COMPILED.md |
| "5.57 triệu dân" | ✅ | audit_report.md, Issue #7 |
