# ALL CORE RESULTS COMPILED (NEWEST LAND COVER DATA)

This file contains all the statistical tables, regression results, and hotspot analyses generated from the pipeline.

## Table 2: Summary Statistics by Density

| Stratum        |   Cells |   Population |   Mean LST P75 (°C) | Mean Physical GB (%)   | Mean Accessible UGBS (%)   | Mean Gap (%)   | Conversion Ratio (%)   |   Mean Reachable Distance (m) | Unreachable Cells (%)   | Pop in LST Hotspot (%)   |
|:---------------|--------:|-------------:|--------------------:|:-----------------------|:---------------------------|:---------------|:-----------------------|------------------------------:|:------------------------|:-------------------------|
| Low-density    |   10281 |    1,964,953 |               42.81 | 37.81%                 | 0.12%                      | 37.70%         | 0.31%                  |                         534.8 | 83.5%                   | 425,628 (21.7%)          |
| Medium-density |    2647 |    1,604,008 |               45.41 | 38.92%                 | 2.23%                      | 36.69%         | 5.72%                  |                         461.8 | 28.3%                   | 668,681 (41.7%)          |
| High-density   |    1028 |    2,004,339 |               49.69 | 23.85%                 | 2.87%                      | 20.98%         | 12.04%                 |                         381.9 | 5.8%                    | 1,570,389 (78.3%)        |

## Figure 3: Threshold Percentages

| Level     |   <= 300m (%) |   <= 500m (%) |   <= 1000m (%) |
|:----------|--------------:|--------------:|---------------:|
| Level A   |       22.0618 |       37.8674 |        64.1672 |
| Level A+B |       48.9035 |       68.2288 |        86.2806 |

## Table: Top Priority Wards

| commune_name       |   Priority Cells |   Affected Population |
|:-------------------|-----------------:|----------------------:|
| Phường Hà Đông     |               15 |               22466.1 |
| Phường Khương Đình |               12 |               21086.8 |
| Phường Xuân Phương |               12 |               14194.1 |
| Phường Tây Hồ      |                6 |                6370.4 |
| Phường Phú Diễn    |                4 |                4138.8 |
| Phường Từ Liêm     |                3 |                3462.4 |
| Phường Hồng Hà     |                3 |                4367.2 |
| Xã Đan Phượng      |                3 |                3004.3 |
| Phường Bồ Đề       |                2 |                2708.8 |
| Phường Giảng Võ    |                2 |                4689.1 |

## Risk Hotspots Overview

```text
=== HANOI STEP 14: RISK AND HOTSPOT SUMMARY ===

Input: d:\KLTN\data\hanoi_grid_250m_walking_access.gpkg
Output: d:\KLTN\data\hanoi_grid_250m_equity_metrics.gpkg

--- BASIC COUNTS ---
Total grid cells: 15077
Valid LST/risk-surface cells: 13956

--- VARIABLES ---
LST variable used: lst_p75_filled
Dijkstra cutoff used for network gap: 1000 m
Main Gi* K: 8
Sensitivity Gi* K: 4
Permutations: 999

--- GETIS-ORD GI* HOTSPOTS K=8 ---
Priority Area (Network): 81
Priority Area (Area Gap): 29
LST Hotspots FDR: 3249
LST Coldspots FDR: 2989
Network Gap Hotspots FDR: 0
Network Gap Coldspots FDR: 3428

--- BIVARIATE LISA ---
BiLISA (Access Deficit x Pop Density) HH FDR: 180
BiLISA (Access Deficit x Pop Density) HH raw: 247
BiLISA (Access Deficit x LST) HH FDR: 1581
BiLISA (Access Deficit x LST) HH raw: 1864

--- POPULATION IN PRIORITY AREAS ---
Priority Area Network pop (% of valid sample): 124,120 (2.2%)
Priority Area AreaGap pop (% of valid sample): 7,692 (0.1%)
Valid hotspot-analysis population: 5,573,301
Priority Area Network pop (% of all grid): 124,120 (2.2%)
Priority Area AreaGap pop (% of all grid): 7,692 (0.1%)
All retained grid population: 5,638,640

```

## SEM Regression: High Density

```text
=== SPATIAL LAG MODEL (SLM / GM_Lag) ===
REGRESSION RESULTS
------------------

SUMMARY OF OUTPUT: SPATIAL TWO STAGE LEAST SQUARES
------------------------------------------------------------------------------------
Data set            :     unknown
Weights matrix      :     unknown
Dependent Variable  :lst_p75_observed                Number of Observations:        1188
Mean dependent var  :     49.7391                Number of Variables   :           8
S.D. dependent var  :      2.5102                Degrees of Freedom    :        1180
Pseudo R-squared    :      0.6445
Spatial Pseudo R-squared:  0.4247

------------------------------------------------------------------------------------
            Variable     Coefficient       Std.Error     z-Statistic     Probability
------------------------------------------------------------------------------------
            CONSTANT        35.28803         2.21793        15.91035         0.00000
           tree_frac        -3.21405         0.44815        -7.17176         0.00000
          grass_frac        -6.52798         0.99457        -6.56360         0.00000
           crop_frac       -17.29088         1.84596        -9.36690         0.00000
          water_frac       -13.18604         0.82335       -16.01516         0.00000
           bare_frac        -5.87939         0.83394        -7.05015         0.00000
  log_dist_jrc_water         0.00015         0.05889         0.00260         0.99793
  W_lst_p75_observed         0.32131         0.04312         7.45153         0.00000
------------------------------------------------------------------------------------
Instrumented: W_lst_p75_observed
Instruments: W_bare_frac, W_crop_frac, W_grass_frac, W_log_dist_jrc_water,
             W_tree_frac, W_water_frac

DIAGNOSTICS FOR SPATIAL DEPENDENCE
TEST                              DF         VALUE           PROB
Anselin-Kelejian Test             1         50.314           0.0000

SPATIAL LAG MODEL IMPACTS
Impacts computed using the 'simple' method.
            Variable         Direct        Indirect          Total
           tree_frac        -3.2140         -1.5216         -4.7357
          grass_frac        -6.5280         -3.0906         -9.6186
           crop_frac       -17.2909         -8.1861        -25.4770
          water_frac       -13.1860         -6.2428        -19.4288
           bare_frac        -5.8794         -2.7835         -8.6629
  log_dist_jrc_water         0.0002          0.0001          0.0002
================================ END OF REPORT =====================================

SLM residual Moran's I: 0.46894 (E[I]=-0.00084)

=== SPATIAL ERROR MODEL (SEM / GM_Error) ===
REGRESSION RESULTS
------------------

SUMMARY OF OUTPUT: GM SPATIALLY WEIGHTED LEAST SQUARES
------------------------------------------------------------------------------------
Data set            :     unknown
Weights matrix      :     unknown
Dependent Variable  :lst_p75_observed                Number of Observations:        1188
Mean dependent var  :     49.7391                Number of Variables   :           7
S.D. dependent var  :      2.5102                Degrees of Freedom    :        1181
Pseudo R-squared    :      0.3982

------------------------------------------------------------------------------------
            Variable     Coefficient       Std.Error     z-Statistic     Probability
------------------------------------------------------------------------------------
            CONSTANT        51.28787         0.64506        79.50845         0.00000
           tree_frac        -5.74585         0.40130       -14.31813         0.00000
          grass_frac        -6.98026         0.73154        -9.54189         0.00000
           crop_frac       -13.13668         1.39905        -9.38969         0.00000
          water_frac       -12.79263         0.71509       -17.88958         0.00000
           bare_frac        -5.62964         0.67028        -8.39895         0.00000
  log_dist_jrc_water         0.04187         0.09527         0.43949         0.66031
              lambda         0.77579    
------------------------------------------------------------------------------------
================================ END OF REPORT =====================================

SEM residual Moran's I: 0.69839 (E[I]=-0.00084)


```

## SEM Regression: Medium Density

```text
=== SPATIAL LAG MODEL (SLM / GM_Lag) ===
REGRESSION RESULTS
------------------

SUMMARY OF OUTPUT: SPATIAL TWO STAGE LEAST SQUARES
------------------------------------------------------------------------------------
Data set            :     unknown
Weights matrix      :     unknown
Dependent Variable  :lst_p75_observed                Number of Observations:        3013
Mean dependent var  :     47.0962                Number of Variables   :           8
S.D. dependent var  :      4.3036                Degrees of Freedom    :        3005
Pseudo R-squared    :      0.7461
Spatial Pseudo R-squared:  0.5332

------------------------------------------------------------------------------------
            Variable     Coefficient       Std.Error     z-Statistic     Probability
------------------------------------------------------------------------------------
            CONSTANT        34.31904         1.02407        33.51231         0.00000
           tree_frac        -5.25394         0.37710       -13.93252         0.00000
          grass_frac        -7.44273         0.35041       -21.24029         0.00000
           crop_frac       -12.61095         0.41353       -30.49551         0.00000
          water_frac       -13.00913         0.52289       -24.87914         0.00000
           bare_frac        -5.13890         0.37972       -13.53335         0.00000
  log_dist_jrc_water         0.14148         0.04373         3.23505         0.00122
  W_lst_p75_observed         0.33256         0.02144        15.51091         0.00000
------------------------------------------------------------------------------------
Instrumented: W_lst_p75_observed
Instruments: W_bare_frac, W_crop_frac, W_grass_frac, W_log_dist_jrc_water,
             W_tree_frac, W_water_frac

DIAGNOSTICS FOR SPATIAL DEPENDENCE
TEST                              DF         VALUE           PROB
Anselin-Kelejian Test             1        188.729           0.0000

SPATIAL LAG MODEL IMPACTS
Impacts computed using the 'simple' method.
            Variable         Direct        Indirect          Total
           tree_frac        -5.2539         -2.6178         -7.8718
          grass_frac        -7.4427         -3.7084        -11.1511
           crop_frac       -12.6109         -6.2835        -18.8945
          water_frac       -13.0091         -6.4819        -19.4911
           bare_frac        -5.1389         -2.5605         -7.6994
  log_dist_jrc_water         0.1415          0.0705          0.2120
================================ END OF REPORT =====================================

SLM residual Moran's I: 0.60307 (E[I]=-0.00033)

=== SPATIAL ERROR MODEL (SEM / GM_Error) ===
REGRESSION RESULTS
------------------

SUMMARY OF OUTPUT: GM SPATIALLY WEIGHTED LEAST SQUARES
------------------------------------------------------------------------------------
Data set            :     unknown
Weights matrix      :     unknown
Dependent Variable  :lst_p75_observed                Number of Observations:        3013
Mean dependent var  :     47.0962                Number of Variables   :           7
S.D. dependent var  :      4.3036                Degrees of Freedom    :        3006
Pseudo R-squared    :      0.5153

------------------------------------------------------------------------------------
            Variable     Coefficient       Std.Error     z-Statistic     Probability
------------------------------------------------------------------------------------
            CONSTANT        50.03062         0.45721       109.42643         0.00000
           tree_frac        -6.38248         0.29360       -21.73870         0.00000
          grass_frac        -7.55891         0.25164       -30.03866         0.00000
           crop_frac       -11.80089         0.29667       -39.77824         0.00000
          water_frac       -12.37712         0.41244       -30.00921         0.00000
           bare_frac        -6.01676         0.28569       -21.06070         0.00000
  log_dist_jrc_water         0.16365         0.06145         2.66315         0.00774
              lambda         0.78743    
------------------------------------------------------------------------------------
================================ END OF REPORT =====================================

SEM residual Moran's I: 0.82038 (E[I]=-0.00033)


```

## SEM Regression: Low Density

```text
=== SPATIAL LAG MODEL (SLM / GM_Lag) ===
REGRESSION RESULTS
------------------

SUMMARY OF OUTPUT: SPATIAL TWO STAGE LEAST SQUARES
------------------------------------------------------------------------------------
Data set            :     unknown
Weights matrix      :     unknown
Dependent Variable  :lst_p75_observed                Number of Observations:        9069
Mean dependent var  :     41.7541                Number of Variables   :           8
S.D. dependent var  :      4.8470                Degrees of Freedom    :        9061
Pseudo R-squared    :      0.7991
Spatial Pseudo R-squared:  0.6246

------------------------------------------------------------------------------------
            Variable     Coefficient       Std.Error     z-Statistic     Probability
------------------------------------------------------------------------------------
            CONSTANT        34.69942         0.61432        56.48443         0.00000
           tree_frac        -5.88545         0.22840       -25.76771         0.00000
          grass_frac        -8.16764         0.18227       -44.81109         0.00000
           crop_frac       -10.73752         0.15151       -70.86871         0.00000
          water_frac       -13.83217         0.22366       -61.84487         0.00000
           bare_frac        -7.16338         0.18429       -38.87048         0.00000
  log_dist_jrc_water         0.25403         0.02040        12.45242         0.00000
  W_lst_p75_observed         0.31151         0.01412        22.06735         0.00000
------------------------------------------------------------------------------------
Instrumented: W_lst_p75_observed
Instruments: W_bare_frac, W_crop_frac, W_grass_frac, W_log_dist_jrc_water,
             W_tree_frac, W_water_frac

DIAGNOSTICS FOR SPATIAL DEPENDENCE
TEST                              DF         VALUE           PROB
Anselin-Kelejian Test             1        587.580           0.0000

SPATIAL LAG MODEL IMPACTS
Impacts computed using the 'simple' method.
            Variable         Direct        Indirect          Total
           tree_frac        -5.8855         -2.6629         -8.5484
          grass_frac        -8.1676         -3.6955        -11.8632
           crop_frac       -10.7375         -4.8583        -15.5958
          water_frac       -13.8322         -6.2585        -20.0906
           bare_frac        -7.1634         -3.2411        -10.4045
  log_dist_jrc_water         0.2540          0.1149          0.3690
================================ END OF REPORT =====================================

SLM residual Moran's I: 0.73340 (E[I]=-0.00011)

=== SPATIAL ERROR MODEL (SEM / GM_Error) ===
REGRESSION RESULTS
------------------

SUMMARY OF OUTPUT: GM SPATIALLY WEIGHTED LEAST SQUARES
------------------------------------------------------------------------------------
Data set            :     unknown
Weights matrix      :     unknown
Dependent Variable  :lst_p75_observed                Number of Observations:        9069
Mean dependent var  :     41.7541                Number of Variables   :           7
S.D. dependent var  :      4.8470                Degrees of Freedom    :        9062
Pseudo R-squared    :      0.6108

------------------------------------------------------------------------------------
            Variable     Coefficient       Std.Error     z-Statistic     Probability
------------------------------------------------------------------------------------
            CONSTANT        47.06474         0.21113       222.91753         0.00000
           tree_frac        -7.82543         0.16414       -47.67546         0.00000
          grass_frac        -8.13439         0.11446       -71.06912         0.00000
           crop_frac       -10.87636         0.08816      -123.37062         0.00000
          water_frac       -13.89018         0.14860       -93.47376         0.00000
           bare_frac        -7.27072         0.12610       -57.65877         0.00000
  log_dist_jrc_water         0.41656         0.02359        17.65565         0.00000
              lambda         0.88959    
------------------------------------------------------------------------------------
================================ END OF REPORT =====================================

SEM residual Moran's I: 0.86302 (E[I]=-0.00011)


```

## SEM Regression: Global Model (All densities)

```text
=== SPATIAL LAG MODEL (SLM / GM_Lag) ===
REGRESSION RESULTS
------------------

SUMMARY OF OUTPUT: SPATIAL TWO STAGE LEAST SQUARES
------------------------------------------------------------------------------------
Data set            :     unknown
Weights matrix      :     unknown
Dependent Variable  :lst_p75_observed                Number of Observations:       13373
Mean dependent var  :     43.7016                Number of Variables   :           8
S.D. dependent var  :      5.4111                Degrees of Freedom    :       13365
Pseudo R-squared    :      0.8892
Spatial Pseudo R-squared:  0.7454

------------------------------------------------------------------------------------
            Variable     Coefficient       Std.Error     z-Statistic     Probability
------------------------------------------------------------------------------------
            CONSTANT        29.22511         0.35443        82.45710         0.00000
           tree_frac        -5.71616         0.12971       -44.06924         0.00000
          grass_frac        -8.78335         0.11627       -75.54136         0.00000
           crop_frac        -9.38257         0.11129       -84.31002         0.00000
          water_frac       -12.52278         0.15430       -81.15712         0.00000
           bare_frac        -4.14819         0.15742       -26.35183         0.00000
  log_dist_jrc_water         0.17782         0.01398        12.72248         0.00000
  W_lst_p75_observed         0.43018         0.00730        58.94409         0.00000
------------------------------------------------------------------------------------
Instrumented: W_lst_p75_observed
Instruments: W_bare_frac, W_crop_frac, W_grass_frac, W_log_dist_jrc_water,
             W_tree_frac, W_water_frac

DIAGNOSTICS FOR SPATIAL DEPENDENCE
TEST                              DF         VALUE           PROB
Anselin-Kelejian Test             1       1441.073           0.0000

SPATIAL LAG MODEL IMPACTS
Impacts computed using the 'simple' method.
            Variable         Direct        Indirect          Total
           tree_frac        -5.7162         -4.3153        -10.0315
          grass_frac        -8.7834         -6.6309        -15.4142
           crop_frac        -9.3826         -7.0832        -16.4658
          water_frac       -12.5228         -9.4539        -21.9767
           bare_frac        -4.1482         -3.1316         -7.2798
  log_dist_jrc_water         0.1778          0.1342          0.3121
================================ END OF REPORT =====================================

SLM residual Moran's I: 0.60547 (E[I]=-0.00007)

=== SPATIAL ERROR MODEL (SEM / GM_Error) ===
REGRESSION RESULTS
------------------

SUMMARY OF OUTPUT: GM SPATIALLY WEIGHTED LEAST SQUARES
------------------------------------------------------------------------------------
Data set            :     unknown
Weights matrix      :     unknown
Dependent Variable  :lst_p75_observed                Number of Observations:       13373
Mean dependent var  :     43.7016                Number of Variables   :           7
S.D. dependent var  :      5.4111                Degrees of Freedom    :       13366
Pseudo R-squared    :      0.7062

------------------------------------------------------------------------------------
            Variable     Coefficient       Std.Error     z-Statistic     Probability
------------------------------------------------------------------------------------
            CONSTANT        47.41655         0.18643       254.34036         0.00000
           tree_frac        -8.02399         0.11554       -69.44713         0.00000
          grass_frac        -8.30027         0.08827       -94.03159         0.00000
           crop_frac       -10.37232         0.08690      -119.36353         0.00000
          water_frac       -13.31739         0.13844       -96.19865         0.00000
           bare_frac        -6.02887         0.13503       -44.64856         0.00000
  log_dist_jrc_water         0.36402         0.02272        16.02331         0.00000
              lambda         0.87517    
------------------------------------------------------------------------------------
================================ END OF REPORT =====================================

SEM residual Moran's I: 0.85216 (E[I]=-0.00007)


```

