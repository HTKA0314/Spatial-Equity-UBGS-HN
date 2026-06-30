# Methodological Framework

The following diagram illustrates the 6-layer spatial data science pipeline used in this study. You can copy this Mermaid code into tools like [Mermaid Live Editor](https://mermaid.live/) or draw.io to generate a high-resolution figure for your manuscript.

```mermaid
graph TD
    %% Define Styles
    classDef input fill:#e1f5fe,stroke:#0288d1,stroke-width:2px,color:#000
    classDef process fill:#fff3e0,stroke:#f57c00,stroke-width:2px,color:#000
    classDef output fill:#e8f5e9,stroke:#388e3c,stroke-width:2px,color:#000
    classDef layer fill:#f3e5f5,stroke:#8e24aa,stroke-width:2px,stroke-dasharray: 5 5,color:#000

    %% Inputs
    subgraph Inputs ["Data Sources"]
        S2["Sentinel-2 Level-2A<br/>(Earth Engine)"]:::input
        ESA["ESA WorldCover v200"]:::input
        JRC["JRC Surface Water"]:::input
        OSM["OpenStreetMap<br/>(POIs & Streets)"]:::input
        LST["Landsat 8<br/>(Land Surface Temp)"]:::input
        POP["WorldPop<br/>(Population Data)"]:::input
    end

    %% Layer 1 & 2: Land Cover & Classification
    subgraph L1_2 ["Layer 1 & 2: Physical vs. Functional Spaces"]
        RF["Random Forest<br/>Classification"]:::process
        Mask["Physical Green-Blue<br/>Footprint"]:::output
        Filter["POI Filtering &<br/>Tag Matching"]:::process
        LevelA["Strictly Public<br/>(Level A)"]:::output
        LevelAB["Semi-Public<br/>(Level A+B)"]:::output

        S2 --> RF
        ESA --> RF
        JRC --> RF
        RF --> Mask
        Mask --> Filter
        OSM --> Filter
        Filter --> LevelA
        Filter --> LevelAB
    end

    %% Layer 3: Network Accessibility
    subgraph L3 ["Layer 3: Network-Based Accessibility"]
        Network["Pedestrian Network<br/>Construction (OSMnx)"]:::process
        Routing["Shortest-path Routing<br/>(300m, 500m, 1000m)"]:::process
        AccMetrics["Population-weighted<br/>Accessibility"]:::output

        OSM --> Network
        LevelA --> Routing
        LevelAB --> Routing
        Network --> Routing
        POP --> Routing
        Routing --> AccMetrics
    end

    %% Layer 4: Regression
    subgraph L4 ["Layer 4: Spatial Regression Modeling"]
        VIF["Multicollinearity Check<br/>(VIF < 10)"]:::process
        SEM["Spatial Error Model (SEM)<br/>& Spatial Lag Model (SLM)"]:::process
        Cooling["Cooling Effect Evaluation<br/>(Physical vs Accessible)"]:::output

        LST --> VIF
        Mask --> VIF
        AccMetrics --> VIF
        VIF --> SEM
        SEM --> Cooling
    end

    %% Layer 5: Spatial Equity
    subgraph L5 ["Layer 5: Spatial Equity & Risk"]
        LISA["Bivariate LISA<br/>(Moran's I)"]:::process
        Hotspots["High Heat - Low Access<br/>Risk Hotspots"]:::output

        LST --> LISA
        AccMetrics --> LISA
        LISA --> Hotspots
    end

    %% Layer 6: Validation
    subgraph L6 ["Layer 6: Manual Validation"]
        Sampling["Stratified Random<br/>Sampling"]:::process
        GSV["Google Street View<br/>Verification"]:::process
        Precision["Classification<br/>Precision Metrics"]:::output

        LevelA --> Sampling
        Sampling --> GSV
        GSV --> Precision
    end

    %% Flow connections between main blocks
    AccMetrics -.-> L4
    AccMetrics -.-> L5
    LevelA -.-> L6
```
