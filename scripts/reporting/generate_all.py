import os
import sys

# Add the parent directory to sys.path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))

from scripts.reporting.script_1_make_tables import compile_table_2
from scripts.reporting.script_2_make_plots_full import run as make_all_plots
from scripts.reporting.script_3_make_regression import make_regression_figures

def run_all():
    print("="*50)
    print("Starting Q1 Thesis Reporting Pipeline (Full Cartography)")
    print("="*50)
    
    compile_table_2()
    
    print("\n" + "="*50)
    make_all_plots()
    
    print("\n" + "="*50)
    make_regression_figures()
    
    print("\n" + "="*50)
    print("ALL REPORTS GENERATED SUCCESSFULLY.")

if __name__ == "__main__":
    run_all()
