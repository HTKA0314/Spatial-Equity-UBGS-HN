# Methodological Framework

Dưới đây là sơ đồ Mermaid mô phỏng khung phương pháp luận (Methodological Framework) cho bài báo của bạn. Bạn có thể copy mã nguồn này vào [Mermaid Live Editor](https://mermaid.live/) hoặc Draw.io để vẽ lại, hoặc xem trực tiếp trên các trình đọc Markdown.

```mermaid
graph TD
    %% Define Styles
    classDef input fill:#e1f5fe,stroke:#01579b,stroke-width:2px,color:#000
    classDef process fill:#f3e5f5,stroke:#4a148c,stroke-width:2px,color:#000
    classDef output fill:#e8f5e9,stroke:#1b5e20,stroke-width:2px,color:#000
    classDef core fill:#fff3e0,stroke:#e65100,stroke-width:2px,color:#000

    %% 1. Data Inputs
    subgraph Data ["1. Data Acquisition"]
        A1[Sentinel-2 <br/> 10m Resolution]:::input
        A2[Landsat-8 <br/> 30m Resolution]:::input
        A3[OpenStreetMap <br/> Network & Tags]:::input
        A4[WorldPop <br/> 100m Population]:::input
    end

    %% 2. Physical Layer
    subgraph Physical ["2. Physical Environment Layer"]
        B1[Land Cover Classification <br/> Trees, Grass, Water]:::process
        B2[Land Surface Temp. LST <br/> Cooling metrics]:::process
        A1 --> B1
        A2 --> B2
    end

    %% 3. Accessibility Layer
    subgraph Access ["3. Functional & Network Accessibility Layer"]
        C1[Functional Classification <br/> Level A: Strict Public <br/> Level B: Semi-public]:::process
        C2[Network Routing <br/> Dijkstra Algorithm <br/> Snapping to 250m Grid]:::process
        A3 --> C1
        A3 --> C2
        B1 -.-> C1
    end

    %% 4. Statistical Layer
    subgraph Stats ["4. Cooling Effect Validation"]
        D1[Spatial Regression <br/> SLM / SEM vs OLS]:::core
        D2[Variance Inflation Factor VIF <br/> Multicollinearity Check]:::process
        B1 --> D2
        C2 --> D2
        B2 --> D1
        D2 --> D1
    end

    %% 5. Equity Layer
    subgraph Equity ["5. Spatial Equity & Priority Mapping"]
        E1[Getis-Ord Gi* <br/> LST Hotspots Identification]:::core
        E2[Bivariate LISA <br/> High Pop & Low Access HH]:::core
        E3[Priority Area Identification <br/> High Density + Severe Deficit + Hotspot]:::output
        
        B2 --> E1
        A4 --> E2
        C2 --> E2
        E1 --> E3
        E2 --> E3
    end

    %% Connections across subgraphs
    Access --> Stats
    Stats --> Equity
```

## Gợi ý cách thiết kế trên file thực tế (Draw.io / Visio / PowerPoint):
Để biểu đồ nhìn chuyên nghiệp nhất trong bài báo quốc tế, bạn nên phân bổ theo **chiều dọc** (từ trên xuống dưới) với 5 "tầng" (Layers) rõ rệt:
1. **Tầng trên cùng (Inputs):** Hiển thị các hộp dữ liệu đầu vào. Hãy chèn thêm các icon nhỏ minh họa (vệ tinh cho Sentinel/Landsat, logo bản đồ cho OSM, biểu tượng hình người cho WorldPop) để nhìn trực quan và đỡ nhàm chán.
2. **Tầng thứ hai (Physical):** Thể hiện các thuật toán trích xuất dữ liệu vật lý cơ bản (LST, Land Cover).
3. **Tầng thứ ba (Accessibility - Trái tim của bài):** Ở tầng này, bạn vẽ hai luồng chạy song song: 
    - Một luồng phân loại tính chất công cộng (Functional Classification: Level A, B).
    - Một luồng tính toán khoảng cách đi bộ mạng lưới (Network Routing: Dijkstra).
4. **Tầng thứ tư (Analysis/Validation):** Khối Hồi quy Không gian (Spatial Regression) để chứng minh giả thuyết.
5. **Tầng dưới cùng (Outputs/Equity):** Các phân tích công bằng không gian (Hotspots, Bivariate LISA) và Kết quả cuối cùng - Bản đồ Vùng Ưu tiên (Priority Areas).

**Mẹo về màu sắc:** 
Nên dùng màu sắc thống nhất và có ý nghĩa đồ họa. Ví dụ: Dùng viền xanh lá cho khối dữ liệu mảng xanh, xanh dương cho nước, màu cam/đỏ cho các khối liên quan đến nhiệt độ (LST/Hotspot), và màu đen/xám cho dân số. Điều này giúp người đọc lướt qua là bắt được logic luồng dữ liệu ngay lập tức.
