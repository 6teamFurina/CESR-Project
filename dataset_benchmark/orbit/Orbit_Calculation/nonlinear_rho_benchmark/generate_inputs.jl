#!/usr/bin/env julia

using Printf
using Random
using TOML

const HERE = @__DIR__
const ORBIT_ROOT = normpath(joinpath(HERE, "..", ".."))
const PROJECT_ROOT = normpath(joinpath(HERE, "..", "..", "..", ".."))

function selected_ring()
    ring = :latest
    for argument in ARGS
        startswith(argument, "--ring=") || continue
        ring = Symbol(lowercase(split(argument, "="; limit=2)[2]))
    end
    ring in (:latest, :latest_cesr, :repaired_latest, :legacy, :legacy_cesr, :historical) ||
        error("--ring must be latest or explicitly legacy")
    return ring in (:latest, :latest_cesr, :repaired_latest) ? :latest_cesr : :legacy
end

const ARTIFACT_RING = selected_ring()
const CONTROL_SOURCE = ARTIFACT_RING == :latest_cesr ?
    joinpath(PROJECT_ROOT, "Latest_Lattice", "bmad_reference", "control_tracking", "controls.csv") :
    joinpath(ORBIT_ROOT, "reference", "closed_orbit_response_6x119.csv")
const OUTPUT_DIR = joinpath(HERE, "shared_input", String(ARTIFACT_RING))
const INPUT_PATH = joinpath(OUTPUT_DIR, "nonlinear_rho_correctors.csv")
const MANIFEST_PATH = joinpath(OUTPUT_DIR, "sample_manifest.csv")
const METADATA_PATH = joinpath(OUTPUT_DIR, "input_metadata.toml")
const SCENARIOS = ("all", "horizontal", "vertical")
const RHOS = [1.13, 3.2, 4.53, 6.4, 9.05]
const TRIALS = 600
const SEED = 20260803
const BASE_KICK = 5.0e-6

function control_inventory(path)
    lines = readlines(path)
    isempty(lines) && error("Empty control reference: $path")
    header = split(first(lines), ',')
    names, planes = if first(header) == "lord_id"
        lord_name = findfirst(==("lord_name"), header)
        lord_key = findfirst(==("lord_key"), header)
        variable = findfirst(==("variable"), header)
        (isnothing(lord_name) || isnothing(lord_key) || isnothing(variable)) &&
            error("Latest control metadata has an incomplete header: $path")
        selected = [
            fields for fields in (split(line, ','; keepempty=true) for line in lines[2:end])
            if length(fields) == length(header) &&
               uppercase(fields[lord_key]) == "OVERLAY" &&
               uppercase(fields[variable]) in ("HKICK", "VKICK")
        ]
        (
            [fields[lord_name] for fields in selected],
            [uppercase(fields[variable]) == "HKICK" ? :horizontal : :vertical for fields in selected],
        )
    elseif first(header) in ("coordinate", "observable")
        names = String.(header[2:end])
        (names, [startswith(name, "H") ? :horizontal : :vertical for name in names])
    else
        error("Unexpected control reference header: $path")
    end
    isempty(names) && error("No steering controls found in $path")
    length(unique(names)) == length(names) || error("Control names are not unique")
    return names, planes
end

function active_indices(names, planes, scenario)
    scenario == "all" && return collect(eachindex(names))
    selected_plane = scenario == "horizontal" ? :horizontal : :vertical
    indices = findall(==(selected_plane), planes)
    isempty(indices) && error("No $scenario controls found in selected ring")
    return indices
end

function unit_rms_directions(rng, names, planes, scenario)
    indices = active_indices(names, planes, scenario)
    directions = zeros(TRIALS, length(names))
    for trial in 1:TRIALS
        active = randn(rng, length(indices))
        active ./= sqrt(sum(abs2, active) / length(active))
        directions[trial, indices] .= active
    end
    return directions, indices
end

function main()
    names, planes = control_inventory(CONTROL_SOURCE)
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
                    MersenneTwister(SEED + scenario_index - 1), names, planes, scenario,
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
        "format" => "ring-nonlinear-rho-shared-input-v2",
        "ring" => String(ARTIFACT_RING),
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
        "control_order_source" => CONTROL_SOURCE,
        "control_planes" => String.(planes),
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
