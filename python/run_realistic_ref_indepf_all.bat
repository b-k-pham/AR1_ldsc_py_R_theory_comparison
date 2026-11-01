@echo off
setlocal enabledelayedexpansion
REM Runs do_realistic_sim_Sam_ref_indepf.py (reference panel present, but its f is NOT
REM shared with the main replicates) for every Fst value in the grid 0 -> 0.5. Reuses the
REM existing ref_panel_npy_realistic/ files as-is - does NOT re-run
REM make_ref_panel_realistic.py or touch python_realistic_refpanel_fixedf/.
REM
REM Usage: run_realistic_ref_indepf_all.bat

set FST_LIST=0 0.01 0.02 0.03 0.04 0.05 0.1 0.2 0.3 0.4 0.5
set FAILED=

for %%F in (%FST_LIST%) do (
    echo.
    echo ##################################################
    echo # Fst=%%F
    echo ##################################################
    python do_realistic_sim_Sam_ref_indepf.py %%F
    if errorlevel 1 (
        echo *** Fst=%%F FAILED ***
        set FAILED=!FAILED! %%F
    )
)

echo.
echo ================ Summary ================
if defined FAILED (
    echo Failed Fst values:!FAILED!
    exit /b 1
) else (
    echo All Fst values completed successfully: %FST_LIST%
)
endlocal
