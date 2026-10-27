@echo off
setlocal enabledelayedexpansion
REM Runs do_AR1_sim_Sam_ref_indepf.py (independent reference-panel f, not shared with the
REM main replicates) for every Fst value in the grid 0 -> 0.5, writing output into
REM with_ref_panel_indepf/ - does NOT touch with_ref_panel/ (the current fixed-f results)
REM or re-run make_ref_panel.py (the existing reference-panel .npy files are reused as-is).
REM
REM Usage: run_ref_indepf_all.bat

set FST_LIST=0 0.01 0.02 0.03 0.04 0.05 0.1 0.2 0.3 0.4 0.5
set FAILED=

for %%F in (%FST_LIST%) do (
    echo.
    echo ##################################################
    echo # Fst=%%F
    echo ##################################################
    python do_AR1_sim_Sam_ref_indepf.py %%F
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
