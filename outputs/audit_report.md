# Báo cáo Tổng rà soát Dự án (End-to-End Audit Report) - Revised V3

Báo cáo này liệt kê 10 hạng mục rà soát cốt lõi nhất, nhằm bảo vệ tính toàn vẹn (integrity) của toàn bộ số liệu trước khi xuất bản bản thảo chính thức.

| Issue | Severity | Location | Evidence & Analysis | Recommended fix |
| :--- | :--- | :--- | :--- | :--- |
| **1. CRS Mismatch** | **Major** | `config.py` vs Manuscript | Output `.gpkg` và `config.py` thực tế đang dùng `EPSG:32648` (UTM 48N). Bản thảo ghi `EPSG:3405` là sai lệch. | Sửa toàn bộ Manuscript, Captions thành: *"All spatial datasets were projected to UTM Zone 48N (EPSG:32648)"*. |
| **2. Edge-cell Area Anomalies** | Minor | `.gpkg` (`area_m2`) | Các ô ở rìa lưới bị clip theo ranh giới nghiên cứu. Việc tính mật độ và tỷ lệ phần trăm đã dùng diện tích thực tế. | Cập nhật gọn vào Methods: *"Grid cells along the study boundary were clipped to the administrative boundary, and all area-normalized indicators were calculated using actual cell area"*. |
| **3. LST Outliers / IDW Artifacts** | Minor | `.gpkg` (`lst_p75_filled`) | Max LST đạt `65.8°C`. Hơi cao nhưng phản ánh đúng hiện tượng đảo nhiệt khốc liệt tại các siêu bề mặt nhân tạo. | Bổ sung câu vào bài: *"IDW-filled LST was used for continuous mapping and hotspot detection, while observed-only LST was used for regression models"*. |
| **4. Spatial-weights Sensitivity** | **Moderate / Major** | Spatial Regression Pipeline | Kết luận về liên hệ land-cover–LST cần được chứng minh là vững khi thay đổi cấu trúc không gian (KNN). | Primary k=8; sensitivity k=4 and k=24. Báo cáo xem dấu và thứ hạng tương đối của các biến cốt lõi có duy trì ổn định không. Không thêm distance band. |
| **5. Unreachable/access_deficit Mismatch** | **Major** *(Data Consistency)* | `.gpkg` (`access_deficit_A`) | Có 1094 ô lưới >1000m bị dính `NaN` ở biến `access_deficit_A` thay vì tuyệt đối bằng $1.0$. | All cells with `dist_to_accessible_A_m` > 1000 or `NaN` must have `access_deficit_A` = 1.0. Yêu cầu rerun lại Step 14. |
| **6. Gini and Threshold Consistency** | **Major** | Figure Annotations vs Tables | ECDF và Gini annotations trên Figure phải khớp tuyệt đối 100% với Table và Text sau khi lọc Edge cells. | Thống nhất rounding rule trên toàn bộ bài. Nếu Gini final là `0.979` ở bảng, Figure ghi `0.979`. |
| **7. Population/Admin Consistency** | Moderate | `dshn.csv` & Boundary | Join rate ổn định với vùng lõi đạt 5.57 triệu dân (valid non-edge cells). | Ghi rõ trong Study Area: Phân tích chỉ tập trung vào *"selected urbanized/peri-urban administrative units"* (71 phường xã). |
| **8. Priority Definition Consistency** | **Major** | Step 14, 14B, Methods, Fig 6, Tab 4 | The final manuscript should use Network-based Priority Areas as the main intervention definition. These cells are defined as high-density grid cells located more than 1000 m from the nearest Level A UGBS, or unreachable, and falling within an FDR-corrected LST hotspot. AreaGap priority cells are retained only as a supplementary diagnostic. Bivariate LISA results provide supporting spatial evidence of co-location between access deficits, population density, and LST, but they are not used to define the main priority map. | Remove all wording that describes the main priority areas as “dual-hotspots” or as BiLISA High–High cells. Ensure Fig. 6 and Table 4 are based only on Priority_Area_Network. |
| **9. Gap Unit Consistency** | Minor | Reporting Text / Figures | Physical–access gap mang bản chất là chênh lệch giữa hai số phần trăm. | Đổi đơn vị của Gap từ `%` thành **`percentage points`** (hoặc `pp`) trên toàn bộ văn bản và trục tọa độ. |
| **10. Accessible Area Capping** | **Passed** | `.gpkg` (Area logic) | Đã kiểm tra: `accessible_gb_A_m2 <= physical_greenblue_m2` đạt tỷ lệ tuân thủ 100% ở mọi cell. | Status: Passed. Severity: Not an issue. |

---

## Kế hoạch Hành động Tức thời (Immediate Action Plan)
1. **Rerun Step 14**: Step 14 must be re-run, followed by Step 14B and reporting scripts.
2. **Verify Output**: After rerun, the audit should verify that all >1000 m/unreachable cells have `access_deficit_A = 1.0` and that BiLISA/priority outputs are regenerated from the corrected layer.
3. **Synchronize Reporting**: Sau khi xác minh, chạy lại toàn bộ script reporting (Make Tables, Make Plots) để áp dụng mẫu số "sạch" (5.57 triệu dân) cho mọi file đầu ra.
