#!/usr/bin/env julia

"""Export RF-on Twiss data for the 100 maintained corrector-direction pairs."""

using Beamlines
using Dates
using GTPSA
using Printf
using Random
using SciBmad
using Statistics
using TOML

const HERE = @__DIR__
const PROJECT_DIR = normpath(joinpath(HERE, "..", "..", "..", ".."))
const ORBIT_CALCULATION_DIR = normpath(joinpath(HERE, "..", "..", "Orbit_Calculation"))
include(joinpath(PROJECT_DIR, "cesr_model.jl"))

constant_term(value) = Float64(GTPSA.scalar(value))

function parse_options(args)
    options = Dict{String,String}(
        "trials" => "100",
        "seed" => "20260804",
        "base-kick-rad" => "5e-6",
        "inputs" => joinpath(ORBIT_CALCULATION_DIR, "inputs", "cesr_corrector_samples_1000.csv"),
        "nominal-metadata" => joinpath(HERE, "results", "nominal_optics_metadata.toml"),
        "output-dir" => joinpath(HERE, "results"),
    )
    for argument in args
        startswith(argument, "--") || error("Arguments must have --name=value form: $argument")
        fields = split(argument[3:end], '='; limit=2)
        length(fields) == 2 || error("Missing value in argument: $argument")
        haskey(options, fields[1]) || error("Unknown option: --$(fields[1])")
        options[fields[1]] = fields[2]
    end
    return options
end

function control_names(path)
    header = split(readline(path), ',')
    first(header) == "sample_id" || error("Input CSV must begin with sample_id")
    names = String.(header[2:end])
    length(names) == 119 || error("Expected 119 controls, found $(length(names))")
    return names
end

function gaussian_unit_rms_directions(generator, trials, n_controls, active_indices)
    directions = zeros(trials, n_controls)
    for trial in 1:trials
        active = randn(generator, length(active_indices))
        active ./= sqrt(mean(abs2, active))
        directions[trial, active_indices] .= active
    end
    return directions
end

function direction_controls(names, trials, seed, base_kick)
    horizontal = findall(name -> startswith(name, "H"), names)
    vertical = findall(name -> startswith(name, "V"), names)
    length(horizontal) == 58 || error("Expected 58 horizontal correctors")
    length(vertical) == 61 || error("Expected 61 vertical correctors")
    h = gaussian_unit_rms_directions(MersenneTwister(seed), trials, length(names), horizontal)
    v = gaussian_unit_rms_directions(MersenneTwister(seed + 1), trials, length(names), vertical)
    return base_kick .* (h + v)
end

function integrated_normal_sextupole_strength(element)
    kn2 = constant_term(Beamlines.deval(element.Kn2))
    kn2l = constant_term(Beamlines.deval(element.Kn2L))
    length_m = constant_term(Beamlines.deval(element.L))
    return iszero(kn2) ? kn2l : kn2 * length_m
end

function point_inventory(ring)
    s_m = 0.0
    specifications = NamedTuple[]
    target_indices = Int[]
    for (index, element) in enumerate(ring.line)
        name = String(element.name)
        length_m = constant_term(Beamlines.deval(element.L))
        exit_s_m = s_m + length_m
        k2l = integrated_normal_sextupole_strength(element)
        if !iszero(k2l)
            reference_index = index + 1
            push!(specifications, (;
                point_type="sextupole_exit", element_index=index, element_name=name,
                s_m=exit_s_m, k2l_m2=k2l, reference_index,
            ))
            push!(target_indices, reference_index)
        end
        if startswith(uppercase(name), "DET_")
            push!(specifications, (;
                point_type="detector", element_index=index, element_name=name,
                s_m, k2l_m2=0.0, reference_index=index,
            ))
            push!(target_indices, index)
        end
        s_m = exit_s_m
    end
    count(row -> row.point_type == "sextupole_exit", specifications) == 76 ||
        error("Expected 76 active normal sextupoles")
    count(row -> row.point_type == "detector", specifications) == 99 ||
        error("Expected 99 detectors")
    return specifications, sort!(unique(target_indices))
end

function write_rows(path, rows)
    isempty(rows) && error("Refusing to write an empty CSV")
    mkpath(dirname(path))
    columns = propertynames(first(rows))
    open(path, "w") do io
        println(io, join(columns, ','))
        for row in rows
            println(io, join((getproperty(row, column) for column in columns), ','))
        end
    end
    return path
end

function solved_closed_orbit(ring, initial)
    solution = find_closed_orbit(
        ring;
        v0=reshape(copy(initial), 1, 6),
        coasting_beam=false,
        batch=Val{false}(),
        reltol=1e-12,
        abstol=1e-13,
        maxiter=100,
        warn=false,
    )
    all(solution.sol.retcode .== SciBmad.BatchSolve.RETCODE_SUCCESS) ||
        error("Direction closed-orbit solve failed: $(solution.sol.retcode)")
    return vec(Float64.(solution.v0[1, :]))
end

function main(args=ARGS)
    options = parse_options(args)
    trials = parse(Int, options["trials"])
    seed = parse(Int, options["seed"])
    base_kick = parse(Float64, options["base-kick-rad"])
    output_dir = abspath(options["output-dir"])
    names = control_names(options["inputs"])
    values = direction_controls(names, trials, seed, base_kick)
    nominal_metadata = TOML.parsefile(options["nominal-metadata"])
    nominal_full_tunes = (
        Float64(nominal_metadata["full_tune_1_turn"]),
        Float64(nominal_metadata["full_tune_2_turn"]),
    )

    model = load_cesr_model(zero_value=0.0, rf_on=true)
    specifications, target_indices = point_inventory(model.ring)
    target_elements = [model.ring.line[index] for index in target_indices]
    descriptor = Descriptor(6, 1)
    nominal = solved_closed_orbit(model.ring, zeros(6))
    point_rows = NamedTuple[]
    tune_rows = NamedTuple[]
    elapsed = @elapsed begin
        for trial in 1:trials
            @printf("Direction optics %d/%d\n", trial, trials)
            for (column, name) in enumerate(names)
                model.controls[name] = values[trial, column]
            end
            closed = solved_closed_orbit(model.ring, nominal)
            optics = twiss(
                model.ring;
                GTPSA_descriptor=descriptor,
                at=target_elements,
                v0=reshape(closed, 1, 6),
                v0_and_coast=(reshape(closed, 1, 6), false),
                spin=false,
                RDTs=false,
                normalizing_map=false,
            )
            table = optics.table
            row_by_index = Dict(Int(table.beamline_index[row]) => row for row in eachindex(table.beamline_index))
            full_tune_1 = floor(nominal_full_tunes[1]) + mod(constant_term(optics.tunes[1]), 1)
            full_tune_2 = floor(nominal_full_tunes[2]) + mod(constant_term(optics.tunes[2]), 1)
            push!(tune_rows, (; trial, full_tune_1_turn=full_tune_1, full_tune_2_turn=full_tune_2))
            for specification in specifications
                table_row = row_by_index[specification.reference_index]
                push!(point_rows, (;
                    trial,
                    point_type=specification.point_type,
                    element_index=specification.element_index,
                    element_name=specification.element_name,
                    s_m=specification.s_m,
                    k2l_m2=specification.k2l_m2,
                    beta_1_m=constant_term(table.beta_1[table_row]),
                    beta_2_m=constant_term(table.beta_2[table_row]),
                    phi_1_turn=constant_term(table.phi_1[table_row]),
                    phi_2_turn=constant_term(table.phi_2[table_row]),
                    twiss_reference_index=specification.reference_index,
                    twiss_reference_name=String(table.name[table_row]),
                    twiss_reference_s_m=constant_term(table.s[table_row]),
                ))
            end
        end
    end
    points_path = write_rows(joinpath(output_dir, "direction_optics_points.csv"), point_rows)
    tunes_path = write_rows(joinpath(output_dir, "direction_optics_tunes.csv"), tune_rows)
    metadata = Dict(
        "format" => "cesr-sextupole-beta-phase-direction-optics-v1",
        "date" => string(Dates.today()),
        "trials" => trials,
        "seed" => seed,
        "base_kick_rad" => base_kick,
        "direction_state" => "simultaneous h+v direction at unit rho",
        "lattice_mode" => "RF-on CESR",
        "phase_units" => "turn",
        "solve_seconds" => elapsed,
        "points_csv" => points_path,
        "tunes_csv" => tunes_path,
    )
    open(joinpath(output_dir, "direction_optics_metadata.toml"), "w") do io
        TOML.print(io, metadata; sorted=true)
    end
    @printf("Direction optics time: %.3f s\n", elapsed)
    println("Direction optics: $points_path")
    return 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main())
end
