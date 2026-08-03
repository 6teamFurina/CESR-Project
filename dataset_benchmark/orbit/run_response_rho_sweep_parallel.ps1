param(
    [string]$OutputRoot = "dataset_benchmark/orbit/error_analysis/response_rho_sweep_600/chunks",
    [int]$ThreadsPerProcess = 4
)

$ErrorActionPreference = "Stop"
$env:JULIA_PKG_PRECOMPILE_AUTO = "0"
$projectRoot = Resolve-Path (Join-Path $PSScriptRoot "../..")
$runner = Join-Path $PSScriptRoot "run_response_rho_sweep.jl"
$outputRootPath = if ([System.IO.Path]::IsPathRooted($OutputRoot)) {
    [System.IO.Path]::GetFullPath($OutputRoot)
} else {
    [System.IO.Path]::GetFullPath((Join-Path $projectRoot $OutputRoot))
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
    $job = Start-Job `
        -Name $group.Name `
        -ArgumentList @(
            $projectRoot,
            $runner,
            $group.Rhos,
            $directory,
            $ThreadsPerProcess,
            $stdout
        ) `
        -ScriptBlock {
            param($ProjectRoot, $Runner, $Rhos, $Directory, $Threads, $Log)
            $env:JULIA_PKG_PRECOMPILE_AUTO = "0"
            Set-Location $ProjectRoot
            & julia `
                "--threads=$Threads" `
                "--project=." `
                $Runner `
                "--rhos=$Rhos" `
                "--trials=600" `
                "--output-dir=$Directory" 2>&1 | Tee-Object -FilePath $Log
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

$active = @()
foreach ($group in $groups[0..3]) {
    $active += Start-RhoChunk $group
}

# The low-rho first chunk normally finishes first. Starting the fifth chunk as
# soon as it exits keeps at most four four-thread Julia processes active.
$active[0].Job | Wait-Job | Receive-Job
if ($active[0].Job.State -ne "Completed") {
    throw "$($active[0].Name) failed; inspect $($active[0].Stdout)"
}
$fifth = Start-RhoChunk $groups[4]

foreach ($run in $active[1..3] + @($fifth)) {
    $run.Job | Wait-Job | Receive-Job
    if ($run.Job.State -ne "Completed") {
        throw "$($run.Name) failed; inspect $($run.Stdout)"
    }
}

Write-Output "All five rho chunks completed: $outputRootPath"
