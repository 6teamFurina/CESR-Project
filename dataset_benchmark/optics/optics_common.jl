using Beamlines
using GTPSA
using Printf
using SciBmad
using Statistics
using TOML

const OPTICS_DIR = @__DIR__
const DATASET_DIR = normpath(joinpath(OPTICS_DIR, ".."))
const ORBIT_DIR = joinpath(DATASET_DIR, "orbit")
const ORBIT_CALCULATION_DIR = joinpath(ORBIT_DIR, "Orbit_Calculation")
const PROJECT_DIR = normpath(joinpath(DATASET_DIR, ".."))

include(joinpath(ORBIT_CALCULATION_DIR, "benchmark_scibmad.jl"))

function parse_optics_args(args)
    options = Dict{String,String}(
        "inputs" => joinpath(ORBIT_CALCULATION_DIR, "inputs", "cesr_corrector_samples_1000.csv"),
        "sample-count" => "10",
        "output-dir" => joinpath(OPTICS_DIR, "results", "chromatic_test_10"),
        "reltol" => "1e-8",
        "abstol" => "1e-10",
        "maxiter" => "100",
        "warmup" => "true",
        "twiss-mode" => "fresh",
        "response-matrix-cache" => joinpath(ORBIT_DIR, "reference", "closed_orbit_response_6x119.csv"),
    )
    for argument in args
        startswith(argument, "--") ||
            error("Arguments must have --name=value form: $argument")
        fields = split(argument[3:end], "="; limit=2)
        length(fields) == 2 || error("Missing value in argument: $argument")
        haskey(options, fields[1]) || error("Unknown option: --$(fields[1])")
        options[fields[1]] = fields[2]
    end
    lowercase(options["warmup"]) in ("true", "false") ||
        error("--warmup must be true or false")
    options["twiss-mode"] in ("fresh", "reuse") ||
        error("--twiss-mode must be fresh or reuse")
    return options
end

detector_elements(ring) = [
    element for element in ring.line
    if startswith(uppercase(String(element.name)), "DET_")
]

function apply_sample!(model, names, values)
    length(names) == length(values) || error("Control row has the wrong width")
    for (column, name) in enumerate(names)
        model.controls[name] = values[column]
    end
    return model
end

function solve_input_closed_orbits(
    names,
    values;
    response_matrix_cache,
    reltol,
    abstol,
    maxiter,
)
    guess = prepare_initial_guess(
        names,
        values,
        "response-linear";
        response_matrix_cache,
        recompute_response=false,
        reltol,
        abstol,
        maxiter,
    )
    model = prepare_batch_model(names, values)
    result = frozen_solve_and_track(
        model,
        size(values, 1),
        guess.nominal_jacobian;
        initial_v0=guess.v0,
        reltol,
        abstol,
        maxiter,
    )
    result = apply_full_newton_fallback(
        result,
        names,
        values;
        reltol,
        abstol,
        maxiter,
    )
    all(result.converged) || error(
        "Closed orbit failed for sample rows $(findall(.!result.converged))",
    )
    return (; guess, result)
end

function write_start_orbits(path, sample_ids, closed_orbits, closure_norms)
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, "sample_id,x,px,y,py,z,pz,closure_norm")
        for row in eachindex(sample_ids)
            fields = Any[sample_ids[row]]
            append!(fields, closed_orbits[row, :])
            push!(fields, closure_norms[row])
            println(io, join(fields, ','))
        end
    end
    return path
end
