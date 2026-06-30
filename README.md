# Understanding Green-Blue Space Inequalities in a Rapidly Urbanizing City: Evidence from Hanoi, Vietnam

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

This repository contains the complete analytical pipeline, codebase, and reproducible scripts for the paper: **"Understanding Green-Blue Space Inequalities in a Rapidly Urbanizing City: Evidence from Hanoi, Vietnam."**

The study investigates the relationship between Urban Heat Island (UHI) effects, actual walking accessibility to functional green/blue spaces via street networks, and spatial environmental equity in Hanoi, Vietnam.

## Overview of the Pipeline

The project implements a robust 6-layer spatial data science pipeline, processing high-resolution spatial data at a 250m x 250m grid. The methodology replaces flawed Euclidean (straight-line) distance estimates with actual walking network accessibility using Dijkstra's algorithm.

### Analytical Flow

1. **Layer 1: Physical Green-Blue Extraction:** Extracts physical footprint of tree canopy and water bodies using Sentinel-2, ESA WorldCover, and JRC Water data.
2. **Layer 2: Functional Accessibility Classification:** Filters physical footprints into strictly public (Level A) and semi-public (Level AB) spaces based on OpenStreetMap POIs.
3. **Layer 3: Network-Based Walking Accessibility:** Calculates actual pedestrian routing distance across the urban street network.
4. **Layer 4: Spatial Regression Modeling:** Employs Ordinary Least Squares (OLS) and Spatial Error Models (SEM) to isolate the cooling effect of accessible vs. physical green spaces.
5. **Layer 5: Spatial Equity & Risk Hotspots:** Uses Local Indicators of Spatial Association (LISA) and Moran's I to pinpoint thermal risk hotspots suffering from access deficits.
6. **Layer 6: Manual Validation:** Evaluates the classification accuracy of accessible spaces via a stratified random sampling approach and Google Street View verification.

---

## Contact

- hothikimanh17@gmail.com
- For questions regarding the methodology or codebase, please open an issue in this repository.

