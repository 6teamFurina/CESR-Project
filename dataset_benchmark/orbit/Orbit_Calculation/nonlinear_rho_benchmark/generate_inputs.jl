#!/usr/bin/env julia

using Printf
using Random
using TOML

const HERE = @__DIR__
const ORBIT_ROOT = normpath(joinpath(HERE, "..", ".."))
const RESPONSE_PATH = joinpath(ORBIT_ROOT, "reference", "closed_orbit_response_6x119.csv")
const OUTPUT_DIR = joinpath(HERE, "shared_input")
const INPUT_PATH = joinpath(OUTPUT_DIR, "nonlinear_rho_correctors.csv")
const MANIFEST_PATH = joinpath(OUTPUT_DIR, "sample_manifest.csv")
const METADATA_PATH = joinpath(OUTPUT_DIR, "input_metadata.toml")
const SCENARIOS = ("all", "horizontal", "vertical")
const RHOS = [1.13, 3.2, 4.53, 6.4, 9.05]
const TRIALS = 600
const SEED = 20260803
const BASE_KICK = 5.0e-6

function control_names(path)
    header = split(first(readlines(path)), ',')
    first(header) == "coordinate" || error("Unexpected response header: $path")
    names = String.(header[2:end])
    length(names) == 119 || error("Expected 119 controls, found $(length(names))")
    return names
end

function active_indices(names, scenario)
    scenario == "all" && return collect(eachindex(names))
    prefix = scenario == "horizontal" ? "H" : "V"
    indices = findall(name -> startswith(name, prefix), names)
    expected = scenario == "horizontal" ? 58 : 61
    length(indices) == expected || error("Expected $expected $scenario controls")
    return indices
end

function unit_rms_directions(rng, names, scenario)
    indices = active_indices(names, scenario)
    directions = zeros(TRIALS, length(names))
    for trial in 1:TRIALS
        active = randn(rng, length(indices))
        active ./= sqrt(sum(abs2, active) / length(active))
        directions[trial, indices] .= active
    end
    return directions, indices
end

function main()
    names = control_names(RESPONSE_PATH)
    mkpath(OUTPUT_DIR)
    open(INPUT_PATH, "w") do input_io
        open(MANIFEST_PATH, "w") do manifest_io
            println(input_io, join(vcat("sample_id", names), ','))
            println(manifest_io, "sample_id,scenario,rho,trial_id,active_rms_rad,global_rms_rad,maximum_abs_kick_rad")
            println(input_io, join(vcat("0", fill("0", length(names))), ','))
            println(manifest_io, "0,baseline,0,0,0,0,0")
            sample_id = 0
            for (scenario_index, scenario) in enumerate(SCENARIOS)
                directions, indices = unit_rms_directions(
                    MersenneTwister(SEED + scenario_index - 1), names, scenario,
                )
                for rho in RHOS, trial in 1:TRIALS
                    sample_id += 1
                    values = rho * BASE_KICK .* view(directions, trial, :)
                    active_rms = sqrt(sum(abs2, view(values, indices)) / length(indices))
                    global_rms = sqrt(sum(abs2, values) / length(values))
                    max_kick = maximum(abs, values)
                    fields = [string(sample_id); [@sprintf("%.17g", value) for value in values]]
                    println(input_io, join(fields, ','))
                    println(
                        manifest_io,
                        join(
                            (
                                sample_id,
                                scenario,
                                @sprintf("%.17g", rho),
                                trial,
                                @sprintf("%.17g", active_rms),
                                @sprintf("%.17g", global_rms),
                                @sprintf("%.17g", max_kick),
                            ),
                            ',',
                        ),
                    )
                end
            end
            sample_id == length(SCENARIOS) * length(RHOS) * TRIALS ||
                error("Sample count mismatch")
        end
    end

    metadata = Dict(
        "format" => "cesr-nonlinear-rho-shared-input-v1",
        "scenarios" => collect(SCENARIOS),
        "rhos" => RHOS,
        "trials_per_scenario_rho" => TRIALS,
        "positive_sample_count" => length(SCENARIOS) * length(RHOS) * TRIALS,
        "total_sample_count" => 1 + length(SCENARIOS) * length(RHOS) * TRIALS,
        "seed" => SEED,
        "scenario_seeds" => Dict(
            scenario => SEED + index - 1 for (index, scenario) in enumerate(SCENARIOS)
        ),
        "base_kick_rad" => BASE_KICK,
        "direction_distribution" => "Gaussian direction normalized to exact unit RMS over active controls",
        "direction_reuse" => "same 600 directions reused at every rho within each scenario",
        "control_order_source" => RESPONSE_PATH,
        "input_csv" => INPUT_PATH,
        "manifest_csv" => MANIFEST_PATH,
    )
    open(METADATA_PATH, "w") do io
        TOML.print(io, metadata; sorted=true)
    end
    total_sample_count = metadata["total_sample_count"]
    println("Wrote $total_sample_count shared states to $INPUT_PATH")
end

main()
