@echo off
setlocal enabledelayedexpansion
REM Runs run_ref_pipeline.bat (make_ref_panel.py + do_AR1_sim_Sam_ref.py) for every
REM Fst value in the grid 0 -> 0.5 used throughout this project. Keeps going even if
REM one Fst fails, and prints a summary of any failures at the end.
REM
REM Usage: run_ref_pipeline_all.bat

set FST_LIST=0 0.01 0.02 0.03 0.04 0.05 0.1 0.2 0.3 0.4 0.5
set FAILED=

for %%F in (%FST_LIST%) do (
    echo.
    echo ##################################################
    echo # Fst=%%F
    echo ##################################################
    call run_ref_pipeline.bat %%F
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
