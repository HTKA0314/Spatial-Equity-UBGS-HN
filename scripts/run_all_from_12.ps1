$ErrorActionPreference = "Stop"
$python = "d:\KLTN\.venv\Scripts\python.exe"

$scripts = @(
    "step12_network_accessibility.py",
    "step13_spatial_regression.py",
    "step13b_stratified_regression.py",
    "step14_risk_hotspots.py",
    "step14b_hotspot_profiling.py",
    "compare_gi_sensitivity.py",
    "generate_q1_figures.py"
)

foreach ($script in $scripts) {
    Write-Host "`n========================================================"
    Write-Host "Executing $script at $(Get-Date)"
    Write-Host "========================================================`n"
    
    $scriptPath = Join-Path "d:\KLTN\scripts" $script
    & $python $scriptPath
    
    if ($LASTEXITCODE -ne 0) {
        throw "Script $script failed with exit code $LASTEXITCODE"
    }
}

Write-Host "`n========================================================"
Write-Host "ALL PIPELINE SCRIPTS COMPLETED SUCCESSFULLY!"
Write-Host "========================================================"
