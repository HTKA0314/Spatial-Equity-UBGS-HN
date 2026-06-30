import os
import sys
import time
import subprocess

sys.stdout.reconfigure(encoding="utf-8")

def run_script(script_name):
    print(f"\n{'='*80}")
    print(f"🚀 RUNNING: {script_name}")
    print(f"{'='*80}")
    
    start_time = time.time()
    
    # Run the script using the current python executable
    result = subprocess.run([sys.executable, script_name], capture_output=False)
    
    elapsed = time.time() - start_time
    
    if result.returncode != 0:
        print(f"\n❌ ERROR: {script_name} failed with exit code {result.returncode} after {elapsed:.1f} seconds.")
        sys.exit(result.returncode)
    else:
        print(f"\n✅ SUCCESS: {script_name} completed in {elapsed:.1f} seconds.")

def main():
    scripts_to_run = [
        "step2_lst_processing.py",
        "step4_water_processing.py",
        "step5_population_mapping.py",
        "step6_osm_network.py",
        "step8_merge_landcover.py",
        "step9_accessible_spaces.py",
        "step10_accessibility_mapping.py",
        "step12_network_accessibility.py",
        "step13_spatial_regression.py",
        "step13b_stratified_regression.py",
        "step14_risk_hotspots.py",
        "step14b_hotspot_profiling.py",
        "generate_thesis_plots.py",
        "compile_results_table.py"
    ]
    
    # Ensure we are in the scripts directory
    current_dir = os.path.basename(os.getcwd())
    if current_dir != "scripts":
        scripts_dir = os.path.join(os.getcwd(), "scripts")
        if os.path.exists(scripts_dir):
            os.chdir(scripts_dir)
        else:
            print("❌ ERROR: Must run this script from the project root or the 'scripts' directory.")
            sys.exit(1)
            
    print("\n" + "#"*80)
    print("🌟 HANOI URBAN HEAT AND GREEN-BLUE ACCESSIBILITY PIPELINE 🌟")
    print("#"*80)
    
    total_start_time = time.time()
    
    for script in scripts_to_run:
        if not os.path.exists(script):
            print(f"\n❌ ERROR: Cannot find script '{script}'. Aborting pipeline.")
            sys.exit(1)
        run_script(script)
        
    total_elapsed = time.time() - total_start_time
    hours, rem = divmod(total_elapsed, 3600)
    minutes, seconds = divmod(rem, 60)
    
    print("\n" + "#"*80)
    print(f"🎉 PIPELINE COMPLETED SUCCESSFULLY! 🎉")
    print(f"⏱️ Total Execution Time: {int(hours)}h {int(minutes)}m {seconds:.1f}s")
    print("#"*80 + "\n")

if __name__ == "__main__":
    main()
