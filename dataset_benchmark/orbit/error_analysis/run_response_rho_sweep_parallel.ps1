param(
    [string]$OutputRoot = "",
    [int]$ThreadsPerProcess = 4,
    [int]$Trials = 600,
    [string]$Ring = "latest",
    [string]$Inputs = "",
    [string]$DetectorResponse = "",
    [string]$ClosedOrbitResponse = "",
    [ValidateSet("gtpsa", "central-difference")]
    [string]$ResponseMethod = "gtpsa",
    [switch]$RecomputeResponse
)

$ErrorActionPreference = "Stop"
$env:JULIA_PKG_PRECOMPILE_AUTO = "0"
$projectRoot = Resolve-Path (Join-Path $PSScriptRoot "../../..")
$runner = Join-Path $PSScriptRoot "run_response_rho_sweep.jl"
$ringAliases = @("latest", "latest_cesr", "repaired_latest", "legacy", "legacy_cesr", "historical")
if ($Ring -notin $ringAliases) {
    throw "Ring must be one of: $($ringAliases -join ', ')"
}
$outputRootPath = if ([System.IO.Path]::IsPathRooted($OutputRoot)) {
    [System.IO.Path]::GetFullPath($OutputRoot)
} else {
    $relativeOutput = if ([string]::IsNullOrWhiteSpace($OutputRoot)) {
        $artifact = if ($Ring -in @("legacy", "legacy_cesr", "historical")) {
            "legacy"
        } else {
            "latest_cesr"
        }
        Join-Path "dataset_benchmark/orbit/error_analysis" (
            "response_rho_sweep_600/{0}/chunks" -f $artifact
        )
    } else {
        $OutputRoot
    }
    [System.IO.Path]::GetFullPath((Join-Path $projectRoot $relativeOutput))
}
$groups = @(
    @{ Name = "chunk_1"; Rhos = "0,0.1,0.14,0.2,0.28" },
    @{ Name = "chunk_2"; Rhos = "0,0.4,0.57,0.8,1.13" },
    @{ Name = "chunk_3"; Rhos = "0,1.6,2.26,3.2,4.53" },
    @{ Name = "chunk_4"; Rhos = "0,6.4,9.05,12.8,18.1" },
    @{ Name = "chunk_5"; Rhos = "0,25.6,36.2,51.2,64" }
)

New-Item -ItemType Directory -Path $outputRootPath -Force | Out-Null

function Start-RhoChunk($group) {
    $directory = Join-Path $outputRootPath $group.Name
    New-Item -ItemType Directory -Path $directory -Force | Out-Null
    $stdout = Join-Path $directory "run.stdout.log"
    $recomputeForChunk = [bool]$RecomputeResponse -and ($group.Name -eq $groups[0].Name)
    $job = Start-Job `
        -Name $group.Name `
        -ArgumentList @(
            $projectRoot,
            $runner,
            $group.Rhos,
            $directory,
            $ThreadsPerProcess,
            $stdout,
            $Trials,
            $Ring,
            $Inputs,
            $DetectorResponse,
            $ClosedOrbitResponse,
            $ResponseMethod,
            $recomputeForChunk
        ) `
        -ScriptBlock {
            param(
                $ProjectRoot,
                $Runner,
                $Rhos,
                $Directory,
                $Threads,
                $Log,
                $Trials,
                $Ring,
                $Inputs,
                $DetectorResponse,
                $ClosedOrbitResponse,
                $ResponseMethod,
                $RecomputeResponse
            )
            $env:JULIA_PKG_PRECOMPILE_AUTO = "0"
            Set-Location $ProjectRoot
            $runnerArguments = @(
                "--rhos=$Rhos",
                "--trials=$Trials",
                "--ring=$Ring",
                "--output-dir=$Directory",
                "--response-method=$ResponseMethod"
            )
            if ($RecomputeResponse) {
                $runnerArguments += "--recompute-response=true"
            }
            if (-not [string]::IsNullOrWhiteSpace($Inputs)) {
                $runnerArguments += "--inputs=$Inputs"
            }
            if (-not [string]::IsNullOrWhiteSpace($DetectorResponse)) {
                $runnerArguments += "--detector-response=$DetectorResponse"
            }
            if (-not [string]::IsNullOrWhiteSpace($ClosedOrbitResponse)) {
                $runnerArguments += "--closed-orbit-response=$ClosedOrbitResponse"
            }
            & julia `
                "--threads=$Threads" `
                "--project=." `
                $Runner `
                @runnerArguments 2>&1 | Tee-Object -FilePath $Log
            if ($LASTEXITCODE -ne 0) {
                throw "Julia exited with status $LASTEXITCODE"
            }
        }
    return [pscustomobject]@{
        Name = $group.Name
        Job = $job
        Directory = $directory
        Stdout = $stdout
    }
}

# Run the first chunk alone so it can create (or explicitly rebuild) the
# shared ring-scoped detector and closed-orbit response caches.  The remaining
# chunks then read complete labeled caches instead of racing to truncate/write
# the same CSV files.  A method is passed explicitly so a stale FD/GTPSA cache
# cannot silently change the production experiment.
$first = Start-RhoChunk $groups[0]
$first.Job | Wait-Job | Receive-Job
if ($first.Job.State -ne "Completed") {
    throw "$($first.Name) failed; inspect $($first.Stdout)"
}
Remove-Job -Job $first.Job -Force

$active = @()
foreach ($group in ($groups | Select-Object -Skip 1)) {
    $active += Start-RhoChunk $group
}

foreach ($run in $active) {
    $run.Job | Wait-Job | Receive-Job
    if ($run.Job.State -ne "Completed") {
        throw "$($run.Name) failed; inspect $($run.Stdout)"
    }
    Remove-Job -Job $run.Job -Force
}

Write-Output "All five rho chunks completed: $outputRootPath"
