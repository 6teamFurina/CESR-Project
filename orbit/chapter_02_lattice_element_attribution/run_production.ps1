param(
    [int]$Trials = 100,
    [int]$Threads = 4,
    [int]$Seed = 20260804,
    [double]$BaseKickRad = 5e-6,
    [string]$PythonExecutable = "",
    [switch]$BuildOnly
)

$ErrorActionPreference = "Stop"
$env:JULIA_PKG_PRECOMPILE_AUTO = "0"
$env:PYTHONDONTWRITEBYTECODE = "1"

if ($Trials -lt 1) {
    throw "Trials must be positive"
}
if ($Threads -lt 1) {
    throw "Threads must be positive"
}
if (-not [double]::IsFinite($BaseKickRad) -or $BaseKickRad -le 0) {
    throw "BaseKickRad must be finite and positive"
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
$study = Join-Path $projectRoot "orbit/error_analysis/thick_element_sextupole_sourcing"
$runner = Join-Path $study "run_thick_element_sourcing.jl"
$analyzer = Join-Path $study "analyze_thick_element_sourcing.py"
$results = Join-Path $chapter "results/latest_cesr"
$horizontal = Join-Path $results "horizontal"
$vertical = Join-Path $results "vertical"
$figures = Join-Path $chapter "figures"
$tables = Join-Path $chapter "tables"

New-Item -ItemType Directory -Path $horizontal, $vertical, $figures, $tables -Force | Out-Null

$previousLocation = Get-Location
try {
    Set-Location $projectRoot
    if (-not $BuildOnly) {
        foreach ($plane in @("x", "y")) {
            $output = if ($plane -eq "x") { $horizontal } else { $vertical }
            Invoke-External "julia" @(
                "--threads=$Threads",
                "--project=.",
                $runner,
                "--ring=latest",
                "--output-plane=$plane",
                "--trials=$Trials",
                "--seed=$Seed",
                "--base-kick-rad=$BaseKickRad",
                "--output-dir=$output"
            )
            Invoke-External $python @($analyzer, $output)
        }
    }

    foreach ($required in @(
        (Join-Path $horizontal "metadata.toml"),
        (Join-Path $horizontal "element_contribution_summary.csv"),
        (Join-Path $vertical "metadata.toml"),
        (Join-Path $vertical "element_contribution_summary.csv")
    )) {
        if (-not (Test-Path -LiteralPath $required)) {
            throw "Required production artifact is missing: $required"
        }
    }

    Invoke-External $python @(
        $analyzer,
        $horizontal,
        "--paired-with", $vertical,
        "--paired-output", (Join-Path $figures "all_element_signed_contributions_paired.svg")
    )
    Invoke-External $python @(
        (Join-Path $chapter "build_paper_outputs.py"),
        "--horizontal", $horizontal,
        "--vertical", $vertical,
        "--figures-dir", $figures,
        "--tables-dir", $tables,
        "--report", (Join-Path $chapter "RESULTS.md")
    )
} finally {
    Set-Location $previousLocation
}

Write-Host "`nChapter 2 production and paper outputs are complete." -ForegroundColor Green
