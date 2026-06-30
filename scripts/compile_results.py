import os
import pandas as pd

output_file = r'd:\KLTN\outputs\ALL_RESULTS_COMPILED.md'
files_to_read = [
    (r'd:\KLTN\outputs\tables\Table2_summary_by_density.csv', 'Table 2: Summary Statistics by Density'),
    (r'd:\KLTN\outputs\tables\table_fig3_threshold_percentages.csv', 'Figure 3: Threshold Percentages'),
    (r'd:\KLTN\outputs\tables\table_fig6_top_wards.csv', 'Table: Top Priority Wards'),
    (r'd:\KLTN\outputs\step14_risk_hotspots\hotspot_summary.txt', 'Risk Hotspots Overview'),
    (r'd:\KLTN\outputs\step13b_stratified_regression\model_1a_physical_observed_high_density\spatial_model_summary.txt', 'SEM Regression: High Density'),
    (r'd:\KLTN\outputs\step13b_stratified_regression\model_1a_physical_observed_medium_density\spatial_model_summary.txt', 'SEM Regression: Medium Density'),
    (r'd:\KLTN\outputs\step13b_stratified_regression\model_1a_physical_observed_low_density\spatial_model_summary.txt', 'SEM Regression: Low Density'),
    (r'd:\KLTN\outputs\step13_spatial_regression\model_1a_physical_observed\spatial_model_summary.txt', 'SEM Regression: Global Model (All densities)'),
]

with open(output_file, 'w', encoding='utf-8') as f_out:
    f_out.write('# ALL CORE RESULTS COMPILED (NEWEST LAND COVER DATA)\n\n')
    f_out.write('This file contains all the statistical tables, regression results, and hotspot analyses generated from the pipeline.\n\n')
    
    for file_path, title in files_to_read:
        f_out.write(f'## {title}\n\n')
        if not os.path.exists(file_path):
            f_out.write(f'*File not found: {file_path}*\n\n')
            continue
            
        if file_path.endswith('.csv'):
            try:
                df = pd.read_csv(file_path)
                f_out.write(df.to_markdown(index=False) + '\n\n')
            except Exception as e:
                f_out.write(f'*Error reading CSV: {e}*\n\n')
        else:
            try:
                with open(file_path, 'r', encoding='utf-8') as f_in:
                    f_out.write('```text\n')
                    f_out.write(f_in.read() + '\n')
                    f_out.write('```\n\n')
            except Exception as e:
                f_out.write(f'*Error reading TXT: {e}*\n\n')

print(f"Successfully compiled all results to {output_file}")
