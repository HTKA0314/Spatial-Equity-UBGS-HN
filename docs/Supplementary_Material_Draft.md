# Supplementary Material Draft

## Supplementary Table S2: OpenStreetMap (OSM) Semantic Tag Classification Scheme for Urban Green-Blue Spaces (UGBS)

| Accessibility Level | Definition | OSM Tag Conditions (Key = Value) |
| :--- | :--- | :--- |
| **Level A** | **Strictly Public UGBS:** Spaces with designated or officially recognized public access. | `leisure` = park, playground, recreation_ground<br>`landuse` = village_green, recreation_ground<br>`tourism` = picnic_site<br>*(Any space with explicit `access` = yes, public, designated)* |
| **Level B** | **Semi-Public / Permissive UGBS:** Spaces that are physically open or conditionally accessible, but lack guaranteed permanent public recreational status. | `leisure` = garden, nature_reserve, sports_centre, pitch<br>`landuse` = grass, forest, cemetery<br>`natural` = wood, grassland, scrub, heath, wetland<br>`amenity` = grave_yard<br>`boundary` = protected_area<br>*(Any space with explicit `access` = permissive)* |
| **Level C** | **Restricted / Private / Physical Blue Reference:** Inaccessible spaces or water bodies lacking explicit pedestrian waterfront access. | Explicit restrictions: `access` / `foot` = no, private, customers, permit, delivery, destination.<br>Physical blue reference: `natural` = water, `landuse` = reservoir, basin.<br>All JRC permanent water reference polygons. |

*(Note: Polygons smaller than 100 m² were excluded. Explicit access restrictions (`access=private/no`) strictly overrode all other tags).*

---

## Phụ lục Hồi quy Không gian (Methods S4)
Bài báo yêu cầu chèn kết quả các model đầy đủ vào S4. Các bảng kết quả đã được trích xuất sẵn dưới dạng file CSV trong thư mục `outputs/`:

*   **Bảng hệ số Hồi quy (Model Coefficients & p-values)**: Nằm tại file `outputs/step13b_stratified_regression/stratified_model_summary.csv`. Bảng này đã gom sẵn hệ số của cả Global model và 3 model theo Low, Medium, High density.
*   **Bảng kiểm định (LM Diagnostics & Moran's I)**: Nằm tại file `outputs/step13b_stratified_regression/global_model_diagnostics.csv`. Trong đó chứa đầy đủ các chỉ số chứng minh tại sao Spatial Lag Model (SLM) lại tốt hơn OLS thông thường.

*(Bạn có thể mở các file CSV này bằng Excel, format lại viền bảng và dán thẳng vào bản Word hoàn chỉnh).*
