param(
    [ValidateSet("all", "rho", "element-family", "sextupole")]
    [string]$Experiment = "all",
    [ValidateSet("smoke", "production")]
    [string]$Mode = "smoke",
    [string]$Ring = "latest",
    [ValidateSet("gtpsa", "central-difference")]
    [string]$ResponseMethod = "gtpsa",
    [int]$SmokeDirections = 2,
    [int]$ProductionDirections = 100,
    [int]$RhoProductionTrials = 600,
    [int]$ThreadsPerRhoProcess = 4,
    [string]$JuliaExecutable = "julia",
    [string]$PythonExecutable = "",
    [switch]$RecomputeResponse
)

$ErrorActionPreference = "Stop"
$env:JULIA_PKG_PRECOMPILE_AUTO = "0"

if ($SmokeDirections -lt 2) {
    throw "SmokeDirections must be at least 2 (the rho sweep requires two directions)"
}
if ($ProductionDirections -lt 1) {
    throw "ProductionDirections must be positive"
}
if ($RhoProductionTrials -lt 1) {
    throw "RhoProductionTrials must be positive"
}

$ringAliases = @(
    "latest", "latest_cesr", "repaired_latest",
    "legacy", "legacy_cesr", "historical"
)
if ($Ring -notin $ringAliases) {
    throw "Ring must be one of: $($ringAliases -join ', ')"
}

function Resolve-PythonExecutable([string]$Requested) {
    if (-not [string]::IsNullOrWhiteSpace($Requested)) {
        return $Requested
    }
    foreach ($candidate in @("python", "python3", "py")) {
        $command = Get-Command $candidate -ErrorAction SilentlyContinue
        if ($null -ne $command) {
            return $command.Source
        }
    }
    $profileCandidates = @($env:USERPROFILE, [Environment]::GetFolderPath("UserProfile")) |
        Where-Object { -not [string]::IsNullOrWhiteSpace($_) } |
        Select-Object -Unique
    foreach ($profilePath in $profileCandidates) {
        $bundled = Join-Path $profilePath ".cache/codex-runtimes/codex-primary-runtime/dependencies/python/python.exe"
        if (Test-Path -LiteralPath $bundled) {
            return $bundled
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
$projectRoot = (Resolve-Path (Join-Path $PSScriptRoot "../..")).Path
$errorAnalysis = Join-Path $projectRoot "orbit/error_analysis"
$artifact = if ($Ring -in @("legacy", "legacy_cesr", "historical")) {
    "legacy"
} else {
    "latest_cesr"
}
$directions = if ($Mode -eq "smoke") { $SmokeDirections } else { $ProductionDirections }
$runRho = $Experiment -in @("all", "rho")
$runElementFamily = $Experiment -in @("all", "element-family")
$runSextupole = $Experiment -in @("all", "sextupole")

$previousLocation = Get-Location
try {
    Set-Location $projectRoot

    if ($runRho) {
        $rhoRunner = Join-Path $errorAnalysis "run_response_rho_sweep.jl"
        $rhoRenderer = Join-Path $errorAnalysis "render_response_rho_sweep_svg.py"
        $rhoStudy = Join-Path $errorAnalysis "response_rho_sweep_600/$artifact"
        if ($Mode -eq "smoke") {
            $rhoOutput = Join-Path $rhoStudy "smoke_runner"
            $rhoArguments = @(
                "--project=.", $rhoRunner,
                "--ring=$Ring", "--rhos=0,0.1", "--trials=$directions",
                "--response-method=$ResponseMethod", "--output-dir=$rhoOutput"
            )
            if ($RecomputeResponse) {
                $rhoArguments += "--recompute-response=true"
            }
            Invoke-External $JuliaExecutable $rhoArguments
            Invoke-External $python @(
                $rhoRenderer,
                "--summary", (Join-Path $rhoOutput "rho_sweep_summary.csv"),
                "--output", (Join-Path $rhoOutput "figures/scibmad_orbit_response_error.svg")
            )
        } else {
            $rhoChunks = Join-Path $rhoStudy "chunks"
            $rhoMerged = Join-Path $rhoStudy "merged"
            $launcher = Join-Path $errorAnalysis "run_response_rho_sweep_parallel.ps1"
            $launcherArguments = @{
                Ring = $Ring
                Trials = $RhoProductionTrials
                ThreadsPerProcess = $ThreadsPerRhoProcess
                ResponseMethod = $ResponseMethod
                OutputRoot = $rhoChunks
            }
            if ($RecomputeResponse) {
                $launcherArguments["RecomputeResponse"] = $true
            }
            & $launcher @launcherArguments
            if ($LASTEXITCODE -ne 0) {
                throw "rho production launcher exited with status $LASTEXITCODE"
            }
            Invoke-External $python @(
                (Join-Path $errorAnalysis "merge_response_rho_sweep_chunks.py"),
                "--root", $rhoChunks, "--output-dir", $rhoMerged
            )
            Invoke-External $python @(
                $rhoRenderer,
                "--summary", (Join-Path $rhoMerged "rho_sweep_summary.csv"),
                "--output", (Join-Path $rhoMerged "figures/scibmad_orbit_response_error.svg")
            )
        }
    }

    if ($runElementFamily) {
        $study = Join-Path $errorAnalysis "thick_element_sextupole_sourcing"
        $runner = Join-Path $study "run_thick_element_sourcing.jl"
        $analyzer = Join-Path $study "analyze_thick_element_sourcing.py"
        $suffix = if ($Mode -eq "smoke") { "${artifact}_smoke_runner" } else { $artifact }
        $horizontal = Join-Path $study "horizontal_results/$suffix"
        $vertical = Join-Path $study "vertical_results/$suffix"
        foreach ($plane in @("x", "y")) {
            $output = if ($plane -eq "x") { $horizontal } else { $vertical }
            Invoke-External $JuliaExecutable @(
                "--project=.", $runner,
                "--ring=$Ring", "--output-plane=$plane", "--trials=$directions",
                "--output-dir=$output"
            )
            Invoke-External $python @($analyzer, $output)
        }
        $paired = Join-Path $study "paired_results/$suffix/element_signed_contributions_paired.svg"
        Invoke-External $python @(
            $analyzer, $horizontal,
            "--paired-with", $vertical,
            "--paired-output", $paired
        )
    }

    if ($runSextupole) {
        $contributionStudy = Join-Path $errorAnalysis "sextupole_detector_contributions"
        $opticsStudy = Join-Path $errorAnalysis "sextupole_beta_phase_correlation"
        $scope = if ($Mode -eq "smoke") { "smoke/$artifact" } else { "results/$artifact" }
        $contributionOutput = Join-Path $contributionStudy $scope
        $opticsOutput = Join-Path $opticsStudy $scope
        Invoke-External $JuliaExecutable @(
            "--project=.",
            (Join-Path $contributionStudy "run_sextupole_detector_contributions.jl"),
            "--ring=$Ring", "--trials=$directions", "--output-dir=$contributionOutput"
        )
        Invoke-External $JuliaExecutable @(
            "--project=.",
            (Join-Path $opticsStudy "export_nominal_optics.jl"),
            "--ring=$Ring", "--output-dir=$opticsOutput"
        )
        Invoke-External $JuliaExecutable @(
            "--project=.",
            (Join-Path $opticsStudy "export_direction_optics.jl"),
            "--ring=$Ring", "--trials=$directions", "--output-dir=$opticsOutput"
        )
        Invoke-External $python @(
            (Join-Path $opticsStudy "analyze_beta_phase_correlation.py"),
            "--optics", (Join-Path $opticsOutput "nominal_optics_points.csv"),
            "--optics-metadata", (Join-Path $opticsOutput "nominal_optics_metadata.toml"),
            "--direction-optics", (Join-Path $opticsOutput "direction_optics_points.csv"),
            "--direction-tunes", (Join-Path $opticsOutput "direction_optics_tunes.csv"),
            "--contributions", (Join-Path $contributionOutput "sextupole_direction_contributions.csv"),
            "--closure", (Join-Path $contributionOutput "direction_closure.csv"),
            "--output-dir", $opticsOutput
        )
    }
} finally {
    Set-Location $previousLocation
}

Write-Host "`nCompleted $Experiment ($Mode) for ring '$Ring'." -ForegroundColor Green
