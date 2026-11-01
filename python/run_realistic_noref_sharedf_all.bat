@echo off
setlocal enabledelayedexpansion
REM Runs do_realistic_sim_Sam_noref_sharedf.py (no reference panel, f shared across all 100
REM replicates, real chr1 LD structure) for every Fst value in the grid 0 -> 0.5.
REM
REM Usage: run_realistic_noref_sharedf_all.bat

set FST_LIST=0 0.01 0.02 0.03 0.04 0.05 0.1 0.2 0.3 0.4 0.5
set FAILED=

for %%F in (%FST_LIST%) do (
    echo.
    echo ##################################################
    echo # Fst=%%F
    echo ##################################################
    python do_realistic_sim_Sam_noref_sharedf.py %%F
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
