# Getis-Ord Gi* Sensitivity: k=8 (Primary) vs k=4 (Sensitivity)

> Overlap statistic is computed among cells in the Step 14 valid sample
> (non-edge, area ≥ 31,250 m², valid LST/Gap variables). Not the full study grid.

## Hotspot / Coldspot Counts

| Variable | k | Role | Total cells | Hotspots (FDR) | Coldspots (FDR) | Not significant | Hotspot % | Coldspot % |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| lst_p75_filled | 8 | Primary | 15077 | 3315 | 2939 | 8823 | 22.0 | 19.5 |
| lst_p75_filled | 4 | Sensitivity | 15077 | 1323 | 774 | 12980 | 8.8 | 5.1 |
| dist_accessible_A_capped_m | 8 | Primary | 15077 | 0 | 3631 | 11446 | 0.0 | 24.1 |
| dist_accessible_A_capped_m | 4 | Sensitivity | 15077 | 0 | 2647 | 12430 | 0.0 | 17.6 |
| relative_gap_A | 8 | Primary | 15077 | 1364 | 1743 | 11970 | 9.0 | 11.6 |
| relative_gap_A | 4 | Sensitivity | 15077 | 603 | 671 | 13803 | 4.0 | 4.5 |

## Spearman Rank Correlation of Gi* Z-scores

| Variable | Spearman rho (z-score k8 vs k4) | Stable (rho >= 0.90) |
| --- | --- | --- |
| lst_p75_filled | 0.9922 | YES |
| dist_accessible_A_capped_m | 0.976 | YES |
| relative_gap_A | 0.9634 | YES |

**Interpretation:** If Spearman rho ≥ 0.90 between k=8 and k=4 z-scores,
the spatial pattern of high/low clustering is stable across neighbourhood
size choices, and conclusions do not depend on the specific k value.
