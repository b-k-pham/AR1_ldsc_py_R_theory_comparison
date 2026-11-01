@echo off
REM Runs make_ref_panel_realistic.py then do_realistic_sim_Sam_ref_fixedf.py with the same
REM Fst value (real chr1 LD structure, reference panel's f shared with main replicates).
REM Usage: run_realistic_ref_pipeline.bat <Fst>
REM   e.g. run_realistic_ref_pipeline.bat 0.3

if "%~1"=="" (
    echo Usage: %~nx0 ^<Fst^>
    echo   e.g. %~nx0 0.3
    exit /b 1
)

set FST=%~1

echo === make_ref_panel_realistic.py %FST% ===
python make_ref_panel_realistic.py %FST%
if errorlevel 1 (
    echo make_ref_panel_realistic.py failed for Fst=%FST% - aborting.
    exit /b 1
)

echo.
echo === do_realistic_sim_Sam_ref_fixedf.py %FST% ===
python do_realistic_sim_Sam_ref_fixedf.py %FST%
if errorlevel 1 (
    echo do_realistic_sim_Sam_ref_fixedf.py failed for Fst=%FST%.
    exit /b 1
)

echo.
echo Done: Fst=%FST%
