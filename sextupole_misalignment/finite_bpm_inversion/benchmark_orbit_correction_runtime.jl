#!/usr/bin/env julia
"""Instrument correction calls without editing production source or outputs."""

using LinearAlgebra, Statistics, TOML
BLAS.set_num_threads(1)
const CORRECTION_SAMPLES = Float64[]
const CORRECTION_ITERATIONS = Int[]
const BENCH_SOURCE = abspath(joinpath(@__DIR__, "..", "quadrupole_orbit_correction",
                                    "generate_corrected_joint_machine_scans.jl"))
const BENCH_OUTPUT = joinpath(@__DIR__, "results", "runtime_breakdown_20260905")
source = read(BENCH_SOURCE, String)
needle = "result = solve_noisy_corrected_machine!("
count(needle, source) == 1 || error("Production correction call changed")
include_string(Main, replace(source, needle => "result = timed_noisy_correction!("), BENCH_SOURCE)

function timed_noisy_correction!(args...; kwargs...)
    timed = @timed solve_noisy_corrected_machine!(args...; kwargs...)
    push!(CORRECTION_SAMPLES, timed.time)
    push!(CORRECTION_ITERATIONS, timed.value.iterations)
    println("Correction timing machine $(length(CORRECTION_SAMPLES)): $(timed.time) s")
    return timed.value
end

main_corrected_scans([
    "--machines=16", "--target-limit=1", "--target-parallelism=serial",
    "--thread-equivalence-check=false", "--output-root=$(joinpath(BENCH_OUTPUT, "correction_probe"))",
])
report = Dict(
    "date" => "2026-09-05", "julia" => string(VERSION),
    "cpu" => Sys.CPU_NAME, "blas_threads" => BLAS.get_num_threads(),
    "semantics" => "one correction per latent machine; includes SVD, three updates, exact SciBmad trial orbits, readout simulation and validation; excludes nominal ORM setup",
    "seconds_by_machine" => CORRECTION_SAMPLES,
    "accepted_iterations_by_machine" => CORRECTION_ITERATIONS,
    "first_call_seconds" => first(CORRECTION_SAMPLES),
    "warm_median_seconds" => median(CORRECTION_SAMPLES[2:end]),
    "warm_min_seconds" => minimum(CORRECTION_SAMPLES[2:end]),
    "warm_max_seconds" => maximum(CORRECTION_SAMPLES[2:end]),
)
open(joinpath(BENCH_OUTPUT, "orbit_correction_timing.toml"), "w") do io
    TOML.print(io, report)
end
println(report)
