param(
    [int]$Trials = 600,
    [int]$ThreadsPerProcess = 4,
    [ValidateSet("gtpsa", "central-difference")]
    [string]$ResponseMethod = "gtpsa",
    [string]$PythonExecutable = "",
    [switch]$RecomputeResponse,
    [switch]$BuildOnly
)

$ErrorActionPreference = "Stop"
$env:JULIA_PKG_PRECOMPILE_AUTO = "0"
$env:PYTHONDONTWRITEBYTECODE = "1"

if ($Trials -lt 2) {
    throw "Trials must be at least 2"
}
if ($ThreadsPerProcess -lt 1) {
    throw "ThreadsPerProcess must be positive"
}

function Resolve-PythonExecutable([string]$Requested) {
    if (-not [string]::IsNullOrWhiteSpace($Requested)) {
        return $Requested
    }
    $profiles = @($env:USERPROFILE, [Environment]::GetFolderPath("UserProfile")) |
        Where-Object { -not [string]::IsNullOrWhiteSpace($_) } |
        Select-Object -Unique
    foreach ($profile in $profiles) {
        $bundled = Join-Path $profile ".cache/codex-runtimes/codex-primary-runtime/dependencies/python/python.exe"
        if (Test-Path -LiteralPath $bundled) {
            return $bundled
        }
    }
    foreach ($candidate in @("python", "python3", "py")) {
        $command = Get-Command $candidate -ErrorAction SilentlyContinue
        if ($null -ne $command) {
            return $command.Source
        }
    }
    throw "No Python interpreter found; pass -PythonExecutable explicitly"
}

function Invoke-External([string]$Executable, [string[]]$Arguments) {
    Write-Host "`n> $Executable $($Arguments -join ' ')" -ForegroundColor Cyan
    & $Executable @Arguments
    if ($LASTEXITCODE -ne 0) {
        throw "$Executable exited with status $LASTEXITCODE"
    }
}

$python = Resolve-PythonExecutable $PythonExecutable
$chapter = $PSScriptRoot
$projectRoot = (Resolve-Path (Join-Path $chapter "../..")).Path
$errorAnalysis = Join-Path $projectRoot "orbit/error_analysis"
$results = Join-Path $chapter "results/latest_cesr"
$chunks = Join-Path $results "chunks"
$merged = Join-Path $results "merged"
$figures = Join-Path $chapter "figures"
$tables = Join-Path $chapter "tables"

New-Item -ItemType Directory -Path $chunks, $merged, $figures, $tables -Force | Out-Null

$previousLocation = Get-Location
try {
    Set-Location $projectRoot
    if (-not $BuildOnly) {
        $parallel = Join-Path $errorAnalysis "run_response_rho_sweep_parallel.ps1"
        $arguments = @{
            OutputRoot = $chunks
            ThreadsPerProcess = $ThreadsPerProcess
            Trials = $Trials
            Ring = "latest"
            ResponseMethod = $ResponseMethod
        }
        if ($RecomputeResponse) {
            $arguments["RecomputeResponse"] = $true
        }
        & $parallel @arguments
        if (-not $?) {
            throw "The rho-sweep production launcher failed"
        }

        Invoke-External $python @(
            (Join-Path $errorAnalysis "merge_response_rho_sweep_chunks.py"),
            "--root", $chunks,
            "--output-dir", $merged
        )
    }

    $summary = Join-Path $merged "rho_sweep_summary.csv"
    $metadata = Join-Path $merged "rho_sweep_metadata.json"
    if (-not (Test-Path -LiteralPath $summary) -or -not (Test-Path -LiteralPath $metadata)) {
        throw "Merged production output is missing; run without -BuildOnly first"
    }

    $renderer = Join-Path $errorAnalysis "render_response_rho_sweep_svg.py"
    Invoke-External $python @(
        $renderer,
        "--summary", $summary,
        "--layout", "stacked",
        "--output", (Join-Path $figures "scibmad_orbit_response_error_stacked.svg")
    )
    Invoke-External $python @(
        $renderer,
        "--summary", $summary,
        "--layout", "stacked",
        "--normalize-rho-squared",
        "--output", (Join-Path $figures "scibmad_orbit_response_error_rho2_normalized_stacked.svg")
    )
    Invoke-External $python @(
        (Join-Path $chapter "build_paper_outputs.py"),
        "--summary", $summary,
        "--metadata", $metadata,
        "--figures-dir", $figures,
        "--tables-dir", $tables,
        "--report", (Join-Path $chapter "RESULTS.md")
    )
} finally {
    Set-Location $previousLocation
}

Write-Host "`nChapter 1 production and paper outputs are complete." -ForegroundColor Green
