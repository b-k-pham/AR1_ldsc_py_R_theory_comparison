@echo off
REM Runs make_ref_panel.py then do_AR1_sim_Sam_ref.py with the same Fst value.
REM Usage: run_ref_pipeline.bat <Fst>
REM   e.g. run_ref_pipeline.bat 0.3

if "%~1"=="" (
    echo Usage: %~nx0 ^<Fst^>
    echo   e.g. %~nx0 0.3
    exit /b 1
)

set FST=%~1

echo === make_ref_panel.py %FST% ===
python make_ref_panel.py %FST%
if errorlevel 1 (
    echo make_ref_panel.py failed for Fst=%FST% - aborting.
    exit /b 1
)

echo.
echo === do_AR1_sim_Sam_ref.py %FST% ===
python do_AR1_sim_Sam_ref.py %FST%
if errorlevel 1 (
    echo do_AR1_sim_Sam_ref.py failed for Fst=%FST%.
    exit /b 1
)

echo.
echo Done: Fst=%FST%
