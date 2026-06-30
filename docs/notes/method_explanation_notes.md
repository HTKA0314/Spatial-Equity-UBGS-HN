# Chi Tiết Phương Pháp Đo Lường (Methodology Notes)

Tài liệu này giải thích các khái niệm cốt lõi trong thuật toán đo lường Khả năng tiếp cận không gian xanh (Network-based accessibility) để dùng làm tư liệu viết bài.

## 1. Cơ chế tính biểu đồ ECDF gia quyền dân số (Population-weighted accessibility)
*(Nguồn tham chiếu: `script_2_make_plots_full.py` - hàm `fig4_population_weighted_accessibility`)*

Biểu đồ ECDF trả lời câu hỏi: **"Bao nhiêu phần trăm dân số thực sự tiếp cận được công viên ở các khoảng cách nhất định?"**

**Cách tính toán:**
1. **Lọc dữ liệu:** Loại bỏ các ô grid (250m) không có người ở (`population > 0`). Tính tổng dân số toàn thành phố.
2. **Sắp xếp theo khoảng cách (Sorting):** Lấy tất cả các ô grid có người ở, gán với khoảng cách mạng lưới (tới công viên gần nhất) và dân số của ô đó. Sau đó, sắp xếp tất cả các ô theo thứ tự khoảng cách tăng dần (từ ô gần công viên nhất đến ô xa nhất).
3. **Cộng dồn dân số (Cumulative Sum):** 
   - Trục hoành (X-axis) là khoảng cách đi bộ.
   - Trục tung (Y-axis) là tỷ lệ phần trăm cộng dồn của dân số. Tại một mốc khoảng cách $x$ (ví dụ 500m), thuật toán quét tất cả các ô grid có khoảng cách $\le 500m$, cộng tổng dân số lại rồi chia cho tổng dân số thành phố.
4. **Tại sao lại dùng Population-Weighted?** Nếu chỉ đếm "số ô grid", kết quả sẽ bị sai lệch vì một ô ở ngoại thành không có người ở sẽ có trọng số ngang bằng một ô ở quận trung tâm đông đúc. Dùng dân số làm trọng số (Y-axis là tích lũy dân số) giúp biểu đồ phản ánh chính xác **trải nghiệm thực tế của con người** và sự bất bình đẳng.

## 2. Làm sao biết đâu là "Công viên" công cộng?
*(Nguồn tham chiếu: `step9_accessible_spaces.py`)*

Nghiên cứu sử dụng dữ liệu mã nguồn mở từ OpenStreetMap (OSM), phân loại các nhãn (tags) do cộng đồng đóng góp để xác định tính chất công cộng của không gian:

*   **Level A (Hoàn toàn công cộng):**
    *   Các không gian có nhãn `leisure` là `park`, `playground`, `recreation_ground`.
    *   Hoặc có gắn nhãn truy cập rõ ràng: `access="yes"` hoặc `access="public"`.
*   **Level B (Bán công cộng / Có điều kiện):**
    *   Các không gian như `cemetery` (nghĩa trang), `nature_reserve` (khu bảo tồn), `sports_centre` (trung tâm thể thao), `garden` (vườn).
    *   Các mảng cỏ (`grass`) hoặc rừng (`forest`) không ghi rõ quyền truy cập, vì không loại trừ khả năng có rào chắn.
*   **Level C (Tư nhân / Không thể tiếp cận):**
    *   Các mảng xanh/hồ nước bị gắn nhãn `access="private"` hoặc `access="no"`.
    *   Mặt nước (hồ/sông) mặc định thuộc Level C nếu không có bờ kè/đường đi dạo ven hồ (để tránh lỗi "ảo giác xanh").

## 3. Điểm xuất phát và "Điểm đại diện" (Representative Point) là gì?
*(Nguồn tham chiếu: `step12_network_accessibility.py`)*

Nghiên cứu chia bản đồ thành một mạng lưới các ô vuông (Grid) 250m x 250m (tương đương với quy mô một khối phố).
*   **Điểm đại diện (`representative_point()`):** Là tọa độ tâm chính xác của ô vuông 250m đó. 
*   Điểm này đóng vai trò là "vạch xuất phát" chung cho toàn bộ người dân sống trong khối phố đó.

## 4. "Ngõ/Con phố gần nhất" và cơ chế "Bắt dính" (Snapping)
*(Nguồn tham chiếu: `step6_osm_network.py` và `step12_network_accessibility.py`)*

*   **Mạng lưới đi bộ (Pedestrian Network):** Tải từ OSM với thuộc tính `network_type="walk"`. Hệ thống này bao gồm vỉa hè, đường nội bộ, ngõ, ngách, hẻm và các con phố cho phép đi bộ (loại trừ đường cao tốc). Mạng lưới này tồn tại dưới dạng một **Đồ thị (Graph)** với các Nút (Ngã ba, ngã tư, ngõ cụt) và Cạnh (Đoạn đường có chiều dài thực tế).
*   **Cơ chế Bắt dính (Snapping):** Điểm tâm của ô lưới 250m có thể nằm trên nóc nhà hoặc giữa hồ nước. Thuật toán `osmnx.nearest_nodes()` sẽ tìm kiếm "nút giao thông" (đầu ngõ, ngã ba) gần nhất trên đồ thị, và kéo điểm xuất phát xuống nút đó.
*   **Thuật toán Dijkstra:** Từ nút xuất phát thực tế trên đường, thuật toán Dijkstra sẽ men theo mạng lưới các con ngõ, ngách để tìm **quãng đường đi bộ thực tế ngắn nhất** tới điểm đích (cũng được bắt dính từ đường viền của công viên). Nếu quãng đường ngoằn ngoèo này vượt quá ngưỡng cho phép (1000m), khu vực đó bị coi là "không thể tiếp cận" (unreachable). Cơ chế này loại bỏ hoàn toàn sai số của khoảng cách đường chim bay (Euclidean distance).

## 5. Mô hình Hồi quy Không gian (Spatial Regression: OLS vs. SEM/SLM)
*(Nguồn tham chiếu: `step13_spatial_regression.py`)*

*Mục tiêu:* Chứng minh bằng thống kê rằng: "Ngay cả khi kiểm soát lượng cây xanh vật lý hiện có, khu vực nào có khả năng tiếp cận (khoảng cách đi bộ tới công viên) tốt hơn thì sẽ có hiệu ứng làm mát (LST) hiệu quả hơn."

**Cách mô hình hoạt động:**
1. **Biến số (Variables):** 
   - Biến phụ thuộc (Y): Nhiệt độ bề mặt (LST).
   - Biến độc lập (X): Tỷ lệ bao phủ cây xanh (`tree_frac`), mặt nước, và quan trọng nhất là **khoảng cách đi bộ tới công viên (`log_dist_accessible_A`)**.
2. **Loại bỏ Đa cộng tuyến (VIF):** Thuật toán tự động kiểm tra hệ số VIF (Variance Inflation Factor). Bất kỳ biến nào có VIF > 10 sẽ bị loại bỏ để mô hình không bị nhiễu do các biến có ý nghĩa trùng lặp nhau.
3. **Tại sao không dùng hồi quy OLS thông thường?** Vì nhiệt độ có tính "lan truyền" (Spatial Autocorrelation) - một khu vực nóng sẽ tỏa nhiệt làm khu vực bên cạnh nóng theo. Mô hình OLS truyền thống không xử lý được tương tác không gian này.
4. **Cơ chế ra quyết định Anselin-Florax:** Hệ thống dựa trên bài test Lagrange Multiplier (LM) để tự động chọn giữa hai mô hình:
   - **SLM (Spatial Lag Model):** Xử lý hiệu ứng lan truyền trực tiếp của biến phụ thuộc (nhiệt độ) giữa các ô lưới lân cận. Đây là mô hình được ưu tiên (Primary) vì nó phản ánh đúng bản chất vật lý lan tỏa của Đảo nhiệt đô thị (UHI spillover effects).
   - **SEM (Spatial Error Model):** Xử lý các yếu tố chưa quan sát được nhưng có tính phân bố không gian (vd: hướng gió, bóng râm tòa nhà). Được dùng làm mô hình độ nhạy (Sensitivity).

## 6. Nhận diện Điểm nóng Rủi ro và Bất công (Hotspot & Bivariate LISA)
*(Nguồn tham chiếu: `step14_risk_hotspots.py`)*

*Mục tiêu:* Tìm ra các "Vùng Ưu Tiên" (Priority Areas) - nơi cư dân đông đúc phải chịu đựng sức nóng tồi tệ nhất nhưng lại không có công viên để giải nhiệt, làm cơ sở cho quy hoạch đô thị.

**Cách mô hình hoạt động:**
1. **Trọng số Không gian KNN (k=8):** Định nghĩa "khu vực lân cận" là 8 ô vuông bao quanh một ô trung tâm. Thuật toán sẽ phân tích một ô vuông trong bối cảnh của các ô xung quanh nó.
2. **Phân tích Getis-Ord Gi* (Điểm nóng LST):** Tìm kiếm các "cụm" (clusters) có nhiệt độ cao bất thường so với phần còn lại của thành phố. Hệ thống sử dụng hiệu chỉnh False Discovery Rate (FDR) cực kỳ khắt khe để đảm bảo kết quả có ý nghĩa thống kê (p < 0.05), loại bỏ các khu vực chỉ nóng ngẫu nhiên.
3. **Bivariate LISA (Moran's I hai biến):** Đo lường sự đồng xuất hiện (chồng chéo không gian) giữa hai biến độc lập: "Sự thiếu hụt khả năng tiếp cận" (Access Deficit) và "Mật độ dân số" (Pop Density). Mô hình tìm ra các cụm **High-High (HH)**: Nơi dân cư RẤT ĐÔNG nhưng công viên lại RẤT XA.
4. **Vùng Ưu tiên Mạng lưới (Priority Area - Network):** Là khu vực đáp ứng ĐỒNG THỜI 3 điều kiện khắt khe:
   - Cư dân đông đúc: Thuộc nhóm mật độ dân số cực cao (> 15,000 người/km²).
   - Thiếu hụt trầm trọng: Cách công viên > 1000m (hoặc hoàn toàn không có mạng lưới đường đi bộ tới).
   - Nguy cơ nhiệt độ cao: Là điểm nóng nhiệt độ (Gi* Hotspot) đã được chứng minh thống kê.
   *=> Phân tích này chuyển đổi từ câu hỏi "Chỗ nào trồng nhiều cây" sang câu hỏi sắc bén hơn: "Chính xác khu dân cư nào cần được nhà nước ưu tiên xây dựng công viên bỏ túi (pocket parks) ngay lập tức để giải quyết bất bình đẳng môi trường (Urban cooling justice)."*
