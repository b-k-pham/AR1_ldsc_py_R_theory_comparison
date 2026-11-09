<#
.SYNOPSIS
  Resumable, crash-safe runner for the n=10000/m=20000 chr1 "realistic LD" sweep across
  the three scenarios that had zero output before the last crash: ref_indepf,
  noref_sharedf, noref_indepf. (ref_fixedf was already complete for all 20 Fst values and
  is intentionally NOT re-run here.)

.DESCRIPTION
  Runs each (script, Fst) combination ONE AT A TIME (never in parallel) as its own Python
  subprocess, so memory from one run is fully released before the next starts. Before
  running a combination, it checks whether the output CSV already exists and skips it if
  so - safe to re-launch this script after a crash without redoing finished work. Each
  run's stdout/stderr is redirected to its own log file under python/logs/n10000_other3/,
  and the subprocess priority is lowered to BelowNormal so the desktop stays responsive
  and you can tell "still working" apart from "actually frozen". A failure in one
  combination is logged and the sweep continues to the next rather than aborting.

.PARAMETER DryRun
  Print the plan (what would run, what would be skipped) without launching anything.

.PARAMETER FstList
  Override the default 20-value Fst grid (strings, so formatting is exact).

.PARAMETER Scripts
  Override the default list of driver scripts to run.

.PARAMETER ThreadLimit
  Caps the BLAS thread count (OMP/MKL/OPENBLAS/NUMEXPR_NUM_THREADS) used by each simulation
  subprocess. 0 (default) leaves it unset - the job will use every core it can get, which is
  fastest but leaves nothing for anything else (e.g. gaming) running at the same time. Set
  this to something like 8-12 (out of 24 cores) to keep the machine usable for other
  foreground work while the sweep runs in the background; the sweep itself just runs slower.

.EXAMPLE
  .\run_n10000_other3.ps1 -DryRun
  .\run_n10000_other3.ps1
  .\run_n10000_other3.ps1 -ThreadLimit 10   # leave ~14 cores free, e.g. for gaming
#>

param(
    [switch]$DryRun,
    [int]$ThreadLimit = 0,
    [string[]]$FstList = @("0.0","0.001","0.002","0.003","0.004","0.005","0.006","0.007","0.008","0.009",
                            "0.01","0.02","0.03","0.04","0.05","0.06","0.07","0.08","0.09","0.1"),
    [string[]]$Scripts = @("do_realistic_sim_Sam_ref_indepf.py",
                           "do_realistic_sim_Sam_noref_sharedf.py",
                           "do_realistic_sim_Sam_noref_indepf.py")
)

if ($ThreadLimit -gt 0) {
    $env:OMP_NUM_THREADS = "$ThreadLimit"
    $env:MKL_NUM_THREADS = "$ThreadLimit"
    $env:OPENBLAS_NUM_THREADS = "$ThreadLimit"
    $env:NUMEXPR_NUM_THREADS = "$ThreadLimit"
    Write-Host "Thread cap: $ThreadLimit cores per job (OMP/MKL/OPENBLAS/NUMEXPR_NUM_THREADS)"
}

$N = 10000
$M = 20000

# script -> (output dir under data/, filename template)
$OutputSpec = @{
    "do_realistic_sim_Sam_ref_indepf.py"    = @{ Dir = "python_realistic_refpanel_indepf"; Needs = "refpanel" }
    "do_realistic_sim_Sam_noref_sharedf.py" = @{ Dir = "python_realistic_noref_sharedf";   Needs = $null }
    "do_realistic_sim_Sam_noref_indepf.py"  = @{ Dir = "python_realistic_noref_indepf";    Needs = $null }
}

$RepoRoot = Split-Path -Parent $PSScriptRoot   # python/.. == repo root
$PythonDir = $PSScriptRoot                     # this script lives in python/
$DataDir = Join-Path $RepoRoot "data"
$RefPanelDir = Join-Path $DataDir "ref_panel_npy_realistic"
$LogDir = Join-Path $PythonDir "logs\n10000_other3"

New-Item -ItemType Directory -Force -Path $LogDir | Out-Null
$SummaryLog = Join-Path $LogDir "_summary.log"

function Write-Summary($msg) {
    $line = "[{0}] {1}" -f (Get-Date -Format "yyyy-MM-dd HH:mm:ss"), $msg
    Write-Host $line
    Add-Content -Path $SummaryLog -Value $line
}

$planned = 0
$skipped = 0
$succeeded = 0
$failed = 0
$failedCombos = @()

Write-Summary "=== run_n10000_other3.ps1 starting (DryRun=$DryRun) ==="

foreach ($script in $Scripts) {
    $spec = $OutputSpec[$script]
    if (-not $spec) {
        Write-Summary "SKIP (unknown script, no output spec): $script"
        continue
    }
    $outDir = Join-Path $DataDir $spec.Dir

    foreach ($fst in $FstList) {
        $label = $fst.Replace(".", "")
        $csvName = "realisticLD_n{0}_m{1}_Fst{2}_sigmas03_h2pop03.csv" -f $N, $M, $label
        $csvPath = Join-Path $outDir $csvName

        $planned++

        if (Test-Path $csvPath) {
            Write-Summary "SKIP (already done): $script Fst=$fst -> $csvName"
            $skipped++
            continue
        }

        if ($spec.Needs -eq "refpanel") {
            $refName = "n{0}_m{1}_ref_panel_Fst{2}_sigma_s_03_realisticLD_chr1_seed123.npy" -f $N, $M, $label
            $refPath = Join-Path $RefPanelDir $refName
            if (-not (Test-Path $refPath)) {
                Write-Summary "SKIP (missing prerequisite ref panel): $script Fst=$fst -> expected $refName"
                $failed++
                $failedCombos += "$script Fst=$fst (missing ref panel)"
                continue
            }
        }

        if ($DryRun) {
            Write-Summary "WOULD RUN: $script Fst=$fst n=$N -> $csvName"
            continue
        }

        $timestamp = Get-Date -Format "yyyyMMdd_HHmmss"
        $jobTag = "{0}_Fst{1}_n{2}" -f ($script -replace '\.py$',''), $label, $N
        $outLog = Join-Path $LogDir "$jobTag.out.log"
        $errLog = Join-Path $LogDir "$jobTag.err.log"

        Write-Summary "RUN: $script Fst=$fst n=$N (logging to $jobTag.out.log / .err.log)"

        try {
            $proc = Start-Process -FilePath "python" `
                -ArgumentList @($script, $fst, "$N") `
                -WorkingDirectory $PythonDir `
                -NoNewWindow -PassThru `
                -RedirectStandardOutput $outLog `
                -RedirectStandardError $errLog

            try { $proc.PriorityClass = [System.Diagnostics.ProcessPriorityClass]::BelowNormal } catch {}

            $proc.WaitForExit()
            $exitCode = $proc.ExitCode
        }
        catch {
            Write-Summary "ERROR launching ${script} Fst=${fst}: $($_.Exception.Message)"
            $exitCode = -1
        }

        # Note: $proc.ExitCode is unreliable under -NoNewWindow -PassThru (often comes back
        # blank even on a clean exit), so success is judged by the actual output file existing,
        # not the exit code.
        if (Test-Path $csvPath) {
            Write-Summary "OK: $script Fst=$fst completed -> $csvName"
            $succeeded++
        }
        else {
            Write-Summary "FAILED: $script Fst=$fst (exit code $exitCode) - see $jobTag.err.log"
            $failed++
            $failedCombos += "$script Fst=$fst (exit $exitCode)"
        }
    }
}

Write-Summary "=== Summary: planned=$planned skipped=$skipped succeeded=$succeeded failed=$failed ==="
if ($failedCombos.Count -gt 0) {
    Write-Summary "Failed combos:"
    foreach ($c in $failedCombos) { Write-Summary "  - $c" }
}
