#!/usr/bin/env julia

"""
Generate sextupole--quadrupole affinity response dictionaries with the repaired
CHESS-U 6 GeV SciBmad lattice.

The calculation is split into three resumable processes because the nominal
nuisance dictionary and candidate-condition dictionaries use different GTPSA
descriptors:

    --stage=screen
    --stage=nominal
    --stage=candidates

`screen` ranks active quadrupoles from exact scalar SciBmad optics at
`Kn1 = Kn1_0 +/- fraction*abs(Kn1_0)`. `nominal` uses batched GTPSA parameters
to calculate every sextupole's own mixed Kn2--offset response. For a given
target, the other 75 own responses form the 150 nuisance columns. The response dictionary is
built from directly measurable one-turn launch experiments: the horizontal and
vertical BPM trajectories for a fixed launch, plus the trajectory changes from
four small x, px, y, and py launch probes. This avoids a non-analytic coupled-
Twiss gauge at zero coupling. `candidates` calculates the target mixed responses
at each retained quadrupole's positive and negative condition. Candidate-
condition nuisances deliberately reuse the nominal dictionary, matching the
maintained affinity definition.
"""

using Beamlines
using Dates
using GTPSA
using Printf
using SciBmad
using Statistics
using TOML

const HERE = @__DIR__
const PROJECT_DIR = normpath(joinpath(HERE, "..", "..", "..", ".."))
const LATEST_LATTICE = joinpath(
    PROJECT_DIR,
    "Latest_Lattice",
    "latest_cesr_scibmad_repaired.jl",
)
include(LATEST_LATTICE)

const DETECTOR_COLUMNS = (
    :trajectory_x,
    :trajectory_y,
    :x_from_x_probe,
    :x_from_px_probe,
    :x_from_y_probe,
    :x_from_py_probe,
    :y_from_x_probe,
    :y_from_px_probe,
    :y_from_y_probe,
    :y_from_py_probe,
)
const TRANSVERSE_PROBE = (1.0e-3, 1.0e-4, 1.0e-3, 1.0e-4)

constant_term(value) = Float64(GTPSA.scalar(Beamlines.deval(value)))
base_name(element) = first(split(uppercase(String(element.name)), '!'))

function parse_options(args)
    options = Dict{String,String}(
        "stage" => "screen",
        "output-dir" => joinpath(HERE, "results", "scibmad_latest", "responses"),
        "top-k" => "15",
        "quadrupole-fraction" => "0.001",
        "max-tune-shift" => "0.01",
        "max-beta-beating" => "0.20",
        "targets" => "all",
        "target-batch" => "4",
        "candidate-target-batch" => "4",
        "phase-reference" => "auto",
        "overwrite" => "false",
    )
    for argument in args
        startswith(argument, "--") || error("Arguments must have --name=value form: $argument")
        fields = split(argument[3:end], '='; limit=2)
        length(fields) == 2 || error("Missing value in argument: $argument")
        haskey(options, fields[1]) || error("Unknown option: --$(fields[1])")
        options[fields[1]] = fields[2]
    end
    options["stage"] in ("screen", "nominal", "candidates") ||
        error("--stage must be screen, nominal, or candidates")
    lowercase(options["overwrite"]) in ("true", "false") ||
        error("--overwrite must be true or false")
    top_k = parse(Int, options["top-k"])
    10 <= top_k <= 20 || error("--top-k must be between 10 and 20")
    for key in ("target-batch", "candidate-target-batch")
        parse(Int, options[key]) > 0 || error("--$key must be positive")
    end
    for key in ("quadrupole-fraction", "max-tune-shift", "max-beta-beating")
        value = parse(Float64, options[key])
        isfinite(value) && value > 0 || error("--$key must be positive and finite")
    end
    return options
end

function active_sextupole_inventory(ring)
    entries = NamedTuple[]
    s_m = 0.0
    for (index, element) in enumerate(ring.line)
        length_m = constant_term(element.L)
        if string(element.kind) == "Sextupole"
            kn2_m3 = constant_term(element.Kn2)
            if !iszero(kn2_m3)
                push!(entries, (;
                    index,
                    name=base_name(element),
                    runtime_name=String(element.name),
                    s_m=s_m + length_m / 2,
                    kn2_m3,
                    kn1_m2=constant_term(element.Kn1),
                    ks1_m2=constant_term(element.Ks1),
                    length_m,
                    x_offset_m=constant_term(element.x_offset),
                    y_offset_m=constant_term(element.y_offset),
                    original_kn2=element.Kn2,
                    original_x_offset=element.x_offset,
                    original_y_offset=element.y_offset,
                    original_tracking_method=element.tracking_method,
                ))
            end
        end
        s_m += length_m
    end
    names = getproperty.(entries, :name)
    length(unique(names)) == length(names) || error("Active sextupole names are not unique")
    length(entries) == 76 || error("Expected 76 active sextupoles, found $(length(entries))")
    return entries
end

function active_quadrupole_inventory(ring)
    slices_by_name = Dict{String,Vector{NamedTuple}}()
    order = String[]
    s_m = 0.0
    for (index, element) in enumerate(ring.line)
        length_m = constant_term(element.L)
        if string(element.kind) == "Quadrupole"
            kn1_m2 = constant_term(element.Kn1)
            if !iszero(kn1_m2)
                name = base_name(element)
                haskey(slices_by_name, name) || push!(order, name)
                push!(get!(slices_by_name, name, NamedTuple[]), (;
                    index,
                    runtime_name=String(element.name),
                    s_m=s_m + length_m / 2,
                    kn1_m2,
                    length_m,
                    original_kn1=element.Kn1,
                ))
            end
        end
        s_m += length_m
    end
    entries = NamedTuple[]
    for name in order
        slices = slices_by_name[name]
        strengths = getproperty.(slices, :kn1_m2)
        maximum(abs.(strengths .- first(strengths))) <= 1.0e-12 * max(abs(first(strengths)), 1.0) ||
            error("Tracking slices for $name do not share one Kn1 value")
        lengths = getproperty.(slices, :length_m)
        centers = getproperty.(slices, :s_m)
        push!(entries, (;
            name,
            index=first(slices).index,
            indices=getproperty.(slices, :index),
            runtime_names=getproperty.(slices, :runtime_name),
            s_m=sum(centers .* lengths) / sum(lengths),
            kn1_m2=first(strengths),
            length_m=sum(lengths),
            original_kn1=getproperty.(slices, :original_kn1),
        ))
    end
    length(entries) == 113 || error("Expected 113 independent active quadrupoles, found $(length(entries))")
    return entries
end

function measurable_bpms(ring)
    detectors = [
        element for element in ring.line
        if startswith(uppercase(String(element.name)), "DET_") &&
           occursin("BPM DETECTOR", uppercase(String(element.label)))
    ]
    length(detectors) == 111 || error("Expected 111 labeled BPMs, found $(length(detectors))")
    return detectors
end

function phase_reference_name(requested, detectors)
    names = base_name.(detectors)
    requested_name = uppercase(strip(requested))
    reference = requested_name == "AUTO" ? first(names) : requested_name
    reference in names || error("Phase-reference BPM $reference is not in the measurable BPM inventory")
    return reference
end

function requested_targets(requested, sextupoles)
    all_names = getproperty.(sextupoles, :name)
    lowercase(strip(requested)) == "all" && return all_names
    selected = [uppercase(strip(name)) for name in split(requested, ',') if !isempty(strip(name))]
    unknown = setdiff(selected, all_names)
    isempty(unknown) || error("Unknown target sextupoles: $(join(unknown, ", "))")
    return [name for name in all_names if name in selected]
end

function solve_closed_orbit(ring)
    solution = find_closed_orbit(
        ring;
        coasting_beam=false,
        batch=Val{false}(),
        warn=true,
    )
    all(solution.sol.retcode .== SciBmad.BatchSolve.RETCODE_SUCCESS) ||
        error("RF-on closed-orbit solve failed: $(solution.sol.retcode)")
    return solution
end

function calculate_twiss(ring, at, closed; descriptor=Descriptor(6, 3))
    return twiss(
        ring;
        GTPSA_descriptor=descriptor,
        at=at,
        v0=closed.v0,
        v0_and_coast=(closed.v0, closed.coasting_beam),
        spin=false,
        RDTs=false,
        normalizing_map=false,
    )
end

function calculate_direct_transport(ring, at, closed; descriptor)
    _, names, _, step_save = SciBmad._twiss_1(ring, at)
    _, eye, _, _, _, _ = SciBmad._twiss_2(
        step_save,
        (closed.v0, closed.coasting_beam),
        descriptor,
        Val{false}(),
        Val{false}(),
    )
    maps = [zero(eye) for _ in step_save]
    first_saved = 1
    if !isempty(step_save) && first(step_save) == 0
        maps[1] = one(eye)
        SciBmad.NNF.setscalar!(maps[1], eye.v0)
        first_saved = 2
    end
    let
        saved_steps = step_save
        saved_maps = maps
        current_step = Ref{Int}(0)
        saved_index = Ref{Int}(first_saved)
        callback = (
            i,
            coordinates,
            current_s,
            current_t_ref,
            last_ds_step,
            last_g,
            transforms_out!,
            transforms_in!,
        ) -> begin
            current_step[] += 1
            if saved_index[] <= length(saved_steps) &&
               current_step[] == saved_steps[saved_index[]]
                transforms_out!(i, coordinates, current_s, current_t_ref)
                SciBmad.NNF.setray!(
                    saved_maps[saved_index[]].v;
                    v=reshape(coordinates.v, :),
                )
                transforms_in!(i, coordinates, current_s, current_t_ref)
                saved_index[] += 1
            end
            nothing
        end
        SciBmad._twiss_4(eye, callback, ring)
    end
    return (; names, maps)
end

function table_index_by_name(optics)
    result = Dict{String,Int}()
    for (index, name) in enumerate(optics.table.name)
        normalized = first(split(uppercase(String(name)), '!'))
        haskey(result, normalized) && error("Duplicate Twiss table name: $normalized")
        result[normalized] = index
    end
    return result
end

wrap_turn(value) = mod(value + 0.5, 1.0) - 0.5

function csv_value(value)
    text = string(value)
    occursin(r"[,\"\n]", text) || return text
    return "\"" * replace(text, "\"" => "\"\"") * "\""
end

function write_rows(path, rows)
    isempty(rows) && error("Refusing to write empty CSV: $path")
    mkpath(dirname(path))
    columns = propertynames(first(rows))
    open(path, "w") do io
        println(io, join(string.(columns), ','))
        for row in rows
            println(io, join((csv_value(getproperty(row, column)) for column in columns), ','))
        end
    end
    return path
end

function read_simple_csv(path)
    lines = readlines(path)
    isempty(lines) && error("Empty CSV: $path")
    header = split(first(lines), ',')
    return [
        Dict(header[index] => fields[index] for index in eachindex(header))
        for line in lines[2:end]
        for fields in (split(line, ','),)
    ]
end

function write_lines(path, values)
    mkpath(dirname(path))
    open(path, "w") do io
        for value in values
            println(io, value)
        end
    end
    return path
end

function write_npy(path, values)
    array = Array{Float64}(values)
    shape_text = length(size(array)) == 1 ? "$(length(array))," : join(size(array), ", ")
    header = "{'descr': '<f8', 'fortran_order': True, 'shape': ($shape_text), }"
    preamble_bytes = 10
    padding = mod(-(preamble_bytes + ncodeunits(header) + 1), 64)
    padded_header = header * repeat(" ", padding) * "\n"
    header_length = ncodeunits(padded_header)
    header_length <= typemax(UInt16) || error("NPY header is too long")
    mkpath(dirname(path))
    open(path, "w") do io
        write(io, UInt8(0x93))
        write(io, codeunits("NUMPY"))
        write(io, UInt8(1), UInt8(0))
        write(io, UInt8(header_length & 0xff), UInt8((header_length >> 8) & 0xff))
        write(io, codeunits(padded_header))
        write(io, reinterpret(UInt8, vec(array)))
    end
    return path
end

function write_metadata(path, metadata)
    mkpath(dirname(path))
    open(path, "w") do io
        TOML.print(io, metadata; sorted=true)
    end
    return path
end

function observation_labels(detectors)
    labels = String[]
    for detector in detectors
        detector_name = lowercase(base_name(detector))
        for observable in DETECTOR_COLUMNS
            push!(labels, "$detector_name:$(String(observable))")
        end
    end
    return labels
end

function direct_observable_values(transport, detectors)
    index_by_name = Dict(
        first(split(uppercase(String(name)), '!')) => index
        for (index, name) in enumerate(transport.names)
    )
    values = Any[]
    for detector in detectors
        map = transport.maps[index_by_name[base_name(detector)]]
        push!(values, map.v[1], map.v[3])
        for coordinate in 1:4
            push!(values, TRANSVERSE_PROBE[coordinate] * GTPSA.deriv(map.v[1], coordinate))
        end
        for coordinate in 1:4
            push!(values, TRANSVERSE_PROBE[coordinate] * GTPSA.deriv(map.v[3], coordinate))
        end
    end
    return values
end

function restore_sextupole!(ring, entry)
    element = ring.line[entry.index]
    element.Kn2 = entry.original_kn2
    element.x_offset = entry.original_x_offset
    element.y_offset = entry.original_y_offset
    element.tracking_method = entry.original_tracking_method
    return nothing
end

function set_target_kn2_parameter!(element, entry, parameter)
    # SciBmadStandard classifies a sextupole carrying zero-valued Kn1/Ks1
    # control fields as a general quadrupole. Once Kn2 is a TPS, the standard
    # MatrixKick setup evaluates sqrt(zero(TPS)^2 + zero(TPS)^2), which is
    # non-analytic. For the scalar lattice the same branch immediately falls
    # back to DriftKick, so selecting it explicitly is algebraically identical
    # at the nominal point and keeps the GTPSA parameterization analytic.
    if iszero(entry.kn1_m2) && iszero(entry.ks1_m2)
        element.tracking_method = BeamTracking.DriftKick()
    end
    element.Kn2 = entry.kn2_m3 + parameter
    return nothing
end

function restore_quadrupole!(ring, entry)
    for (index, original) in zip(entry.indices, entry.original_kn1)
        ring.line[index].Kn1 = original
    end
    return nothing
end

function set_quadrupole!(ring, entry, value)
    for index in entry.indices
        ring.line[index].Kn1 = value
    end
    return nothing
end

function screen_stage(options, ring, sextupoles, quadrupoles, detectors, reference_name, output_dir)
    top_k = parse(Int, options["top-k"])
    fraction = parse(Float64, options["quadrupole-fraction"])
    max_tune_shift = parse(Float64, options["max-tune-shift"])
    max_beta_beating = parse(Float64, options["max-beta-beating"])
    at = vcat(detectors, [ring.line[entry.index] for entry in sextupoles])

    baseline_closed_timed = @timed solve_closed_orbit(ring)
    baseline_twiss_timed = @timed calculate_twiss(ring, at, baseline_closed_timed.value)
    baseline = baseline_twiss_timed.value
    baseline_index = table_index_by_name(baseline)
    reference_index = baseline_index[reference_name]
    baseline_beta_1 = [
        constant_term(baseline.table.beta_1[baseline_index[base_name(detector)]])
        for detector in detectors
    ]
    baseline_beta_2 = [
        constant_term(baseline.table.beta_2[baseline_index[base_name(detector)]])
        for detector in detectors
    ]
    baseline_tunes = constant_term.(baseline.tunes[1:2])

    rows = NamedTuple[]
    screening_start = time()
    for (quadrupole_inventory_index, quadrupole) in enumerate(quadrupoles)
        step = abs(quadrupole.kn1_m2) * fraction
        conditioned = NamedTuple[]
        for sign in (1.0, -1.0)
            set_quadrupole!(ring, quadrupole, quadrupole.kn1_m2 + sign * step)
            closed_timed = @timed solve_closed_orbit(ring)
            optics_timed = @timed calculate_twiss(ring, at, closed_timed.value)
            push!(conditioned, (;
                optics=optics_timed.value,
                closed_seconds=closed_timed.time,
                twiss_seconds=optics_timed.time,
            ))
        end
        restore_quadrupole!(ring, quadrupole)
        plus, minus = conditioned
        plus_index = table_index_by_name(plus.optics)
        minus_index = table_index_by_name(minus.optics)
        plus_reference = plus_index[reference_name]
        minus_reference = minus_index[reference_name]

        tune_shift = maximum(vcat(
            abs.(constant_term.(plus.optics.tunes[1:2]) .- baseline_tunes),
            abs.(constant_term.(minus.optics.tunes[1:2]) .- baseline_tunes),
        ))
        beta_beating = 0.0
        for (detector_index, detector) in enumerate(detectors)
            name = base_name(detector)
            beta_beating = max(
                beta_beating,
                abs(constant_term(plus.optics.table.beta_1[plus_index[name]]) / baseline_beta_1[detector_index] - 1),
                abs(constant_term(plus.optics.table.beta_2[plus_index[name]]) / baseline_beta_2[detector_index] - 1),
                abs(constant_term(minus.optics.table.beta_1[minus_index[name]]) / baseline_beta_1[detector_index] - 1),
                abs(constant_term(minus.optics.table.beta_2[minus_index[name]]) / baseline_beta_2[detector_index] - 1),
            )
        end
        allowed = tune_shift <= max_tune_shift && beta_beating <= max_beta_beating

        for sextupole in sextupoles
            name = sextupole.name
            plus_row = plus_index[name]
            minus_row = minus_index[name]
            delta_log_beta_1 = 0.5 * log(
                constant_term(plus.optics.table.beta_1[plus_row]) /
                constant_term(minus.optics.table.beta_1[minus_row])
            )
            delta_log_beta_2 = 0.5 * log(
                constant_term(plus.optics.table.beta_2[plus_row]) /
                constant_term(minus.optics.table.beta_2[minus_row])
            )
            plus_phi_1 = constant_term(plus.optics.table.phi_1[plus_row]) -
                         constant_term(plus.optics.table.phi_1[plus_reference])
            minus_phi_1 = constant_term(minus.optics.table.phi_1[minus_row]) -
                          constant_term(minus.optics.table.phi_1[minus_reference])
            plus_phi_2 = constant_term(plus.optics.table.phi_2[plus_row]) -
                         constant_term(plus.optics.table.phi_2[plus_reference])
            minus_phi_2 = constant_term(minus.optics.table.phi_2[minus_row]) -
                          constant_term(minus.optics.table.phi_2[minus_reference])
            delta_phi_1_rad = pi * wrap_turn(plus_phi_1 - minus_phi_1)
            delta_phi_2_rad = pi * wrap_turn(plus_phi_2 - minus_phi_2)
            leverage = sqrt(
                delta_log_beta_1^2 + delta_log_beta_2^2 +
                delta_phi_1_rad^2 + delta_phi_2_rad^2
            )
            push!(rows, (;
                sextupole=name,
                sextupole_s_m=sextupole.s_m,
                quadrupole=quadrupole.name,
                quadrupole_inventory_index,
                quadrupole_ring_index=quadrupole.index,
                quadrupole_s_m=quadrupole.s_m,
                quadrupole_k1_m2=quadrupole.kn1_m2,
                delta_k1_m2=step,
                optics_leverage=leverage,
                delta_log_beta_x=delta_log_beta_1,
                delta_log_beta_y=delta_log_beta_2,
                delta_phi_x_rad=delta_phi_1_rad,
                delta_phi_y_rad=delta_phi_2_rad,
                max_abs_tune_shift=tune_shift,
                max_detector_beta_beating=beta_beating,
                allowed=Int(allowed),
                closed_orbit_seconds=plus.closed_seconds + minus.closed_seconds,
                twiss_seconds=plus.twiss_seconds + minus.twiss_seconds,
            ))
        end
        @printf(
            "Optics screen %d/%d %-10s tune %.3e beta %.3e\n",
            quadrupole_inventory_index,
            length(quadrupoles),
            quadrupole.name,
            tune_shift,
            beta_beating,
        )
        flush(stdout)
    end

    selected_rank = Dict{Tuple{String,String},Int}()
    selected_rows = NamedTuple[]
    for sextupole in sextupoles
        target_rows = [row for row in rows if row.sextupole == sextupole.name && row.allowed == 1]
        length(target_rows) >= top_k ||
            error("Only $(length(target_rows)) allowed quadrupoles for $(sextupole.name)")
        sort!(target_rows; by=row -> row.optics_leverage, rev=true)
        for (rank, row) in enumerate(target_rows[1:top_k])
            selected_rank[(row.sextupole, row.quadrupole)] = rank
            push!(selected_rows, (;
                sextupole=row.sextupole,
                sextupole_s_m=row.sextupole_s_m,
                selected_rank=rank,
                quadrupole=row.quadrupole,
                quadrupole_s_m=row.quadrupole_s_m,
                quadrupole_k1_m2=row.quadrupole_k1_m2,
                delta_k1_m2=row.delta_k1_m2,
                optics_leverage=row.optics_leverage,
            ))
        end
    end
    final_rows = [
        (;
            row...,
            selected=Int(haskey(selected_rank, (row.sextupole, row.quadrupole))),
            selected_rank=get(selected_rank, (row.sextupole, row.quadrupole), ""),
        )
        for row in rows
    ]
    write_rows(joinpath(output_dir, "quadrupole_optics_screen.csv"), final_rows)
    write_rows(joinpath(output_dir, "selected_candidates.csv"), selected_rows)
    write_metadata(
        joinpath(output_dir, "screen_metadata.toml"),
        Dict(
            "format" => "cesr-scibmad-quadrupole-affinity-screen-v1",
            "engine" => "SciBmad exact scalar RF-on closed orbit and Twiss",
            "date" => string(Dates.today()),
            "lattice" => LATEST_LATTICE,
            "tracking_elements" => length(ring.line),
            "sextupole_count" => length(sextupoles),
            "quadrupole_count" => length(quadrupoles),
            "detector_count" => length(detectors),
            "phase_reference" => reference_name,
            "candidate_count_per_target" => top_k,
            "quadrupole_fraction" => fraction,
            "max_tune_shift" => max_tune_shift,
            "max_beta_beating" => max_beta_beating,
            "baseline_closed_orbit_seconds" => baseline_closed_timed.time,
            "baseline_twiss_seconds" => baseline_twiss_timed.time,
            "screen_seconds" => time() - screening_start,
            "julia_version" => string(VERSION),
        ),
    )
    return nothing
end

function selected_candidates(output_dir, sextupoles, quadrupoles)
    path = joinpath(output_dir, "selected_candidates.csv")
    isfile(path) || error("Missing screen output: $path")
    rows = read_simple_csv(path)
    quadrupole_by_name = Dict(entry.name => entry for entry in quadrupoles)
    result = Dict{String,Vector{NamedTuple}}(entry.name => NamedTuple[] for entry in sextupoles)
    for row in rows
        quadrupole = quadrupole_by_name[row["quadrupole"]]
        push!(result[row["sextupole"]], (;
            rank=parse(Int, row["selected_rank"]),
            quadrupole=quadrupole.name,
            delta_k1_m2=parse(Float64, row["delta_k1_m2"]),
        ))
    end
    for values in values(result)
        sort!(values; by=value -> value.rank)
    end
    return result
end

function target_directory(output_dir, target_name)
    return joinpath(output_dir, "targets", "$(lowercase(target_name))_responses")
end

function nuisance_labels(target, sextupoles)
    labels = String[]
    for sextupole in sextupoles
        sextupole.name == target.name && continue
        push!(labels, "$(sextupole.name):x_offset", "$(sextupole.name):y_offset")
    end
    return labels
end

function nominal_stage(options, ring, sextupoles, detectors, reference_name, output_dir)
    target_batch_size = parse(Int, options["target-batch"])
    overwrite = lowercase(options["overwrite"]) == "true"
    selected_names = requested_targets(options["targets"], sextupoles)
    selected = [entry for entry in sextupoles if entry.name in selected_names]
    labels = observation_labels(detectors)
    write_lines(joinpath(output_dir, "observation_labels.txt"), labels)

    complete = all(
        isfile(joinpath(target_directory(output_dir, target.name), "nuisance_response_nominal.npy"))
        for target in selected
    )
    if complete && !overwrite
        @printf("All %d requested nominal response bundles reused\n", length(selected))
        flush(stdout)
        return nothing
    end

    # Every nuisance column is another sextupole's *own* Kn2--offset response,
    # not the much smaller cross derivative of the target Kn2 with that other
    # offset. Calculate all 76 own responses once, then omit the target pair
    # when assembling its 150-column nuisance dictionary.
    parameter_count = 3 * target_batch_size
    descriptor = Descriptor(6, 3, parameter_count, 2)
    parameters = params(descriptor)
    closed_timed = @timed solve_closed_orbit(ring)
    closed = closed_timed.value
    stage_start = time()
    timing_rows = NamedTuple[]
    own_response = Dict(
        sextupole.name => fill(NaN, length(labels), 2)
        for sextupole in sextupoles
    )

    for batch_start in 1:target_batch_size:length(sextupoles)
        batch_stop = min(batch_start + target_batch_size - 1, length(sextupoles))
        batch = sextupoles[batch_start:batch_stop]
        for (slot, target) in enumerate(batch)
            element = ring.line[target.index]
            set_target_kn2_parameter!(element, target, parameters[3 * slot - 2])
            element.x_offset = target.x_offset_m + parameters[3 * slot - 1]
            element.y_offset = target.y_offset_m + parameters[3 * slot]
        end

        transport_timed = @timed calculate_direct_transport(
            ring,
            detectors,
            closed;
            descriptor=descriptor,
        )
        observable_tps = direct_observable_values(transport_timed.value, detectors)
        length(observable_tps) == length(labels) || error("Observation count changed")
        for observation_index in eachindex(observable_tps)
            hessian = GTPSA.hessian(observable_tps[observation_index]; include_params=true)
            for (slot, target) in enumerate(batch)
                k2_parameter = 6 + 3 * slot - 2
                own_response[target.name][observation_index, 1] =
                    Float64(hessian[k2_parameter, 6 + 3 * slot - 1])
                own_response[target.name][observation_index, 2] =
                    Float64(hessian[k2_parameter, 6 + 3 * slot])
            end
        end
        for target in batch
            restore_sextupole!(ring, target)
        end
        push!(timing_rows, (;
            target_batch_start=batch_start,
            target_batch_stop=batch_stop,
            targets=join(getproperty.(batch, :name), ';'),
            transport_seconds=transport_timed.time,
        ))
        @printf(
            "Nominal own-response batch %d:%d transport %.3f s\n",
            batch_start,
            batch_stop,
            transport_timed.time,
        )
        flush(stdout)
        GC.gc()
    end

    for target in selected
        target_response = own_response[target.name]
        nuisance_response = hcat((
            own_response[nuisance.name]
            for nuisance in sextupoles
            if nuisance.name != target.name
        )...)
        all(isfinite, target_response) || error("Incomplete target response for $(target.name)")
        all(isfinite, nuisance_response) || error("Incomplete nuisance response for $(target.name)")
        directory = target_directory(output_dir, target.name)
        mkpath(directory)
        write_npy(joinpath(directory, "target_response_nominal.npy"), target_response)
        write_npy(joinpath(directory, "nuisance_response_nominal.npy"), nuisance_response)
        write_lines(joinpath(directory, "nuisance_labels.txt"), nuisance_labels(target, sextupoles))
        write_metadata(
            joinpath(directory, "nominal_metadata.toml"),
            Dict(
                "target" => target.name,
                "sextupole_s_m" => target.s_m,
                "observation_count" => length(labels),
                "nuisance_count" => 2 * (length(sextupoles) - 1),
                "nuisance_definition" => "Own Kn2--offset responses of all other active sextupoles",
            ),
        )
    end

    write_rows(joinpath(output_dir, "nominal_target_batch_timings.csv"), timing_rows)
    write_metadata(
        joinpath(output_dir, "nominal_metadata.toml"),
        Dict(
            "format" => "cesr-scibmad-quadrupole-affinity-nominal-v2",
            "engine" => "SciBmad/GTPSA direct launch-trajectory mixed derivatives",
            "date" => string(Dates.today()),
            "lattice" => LATEST_LATTICE,
            "descriptor" => "Descriptor(6, 3, $parameter_count, 2)",
            "target_batch_size" => target_batch_size,
            "target_count" => length(selected),
            "nuisance_columns_per_target" => 2 * (length(sextupoles) - 1),
            "nuisance_definition" => "Own Kn2--offset responses of all other active sextupoles",
            "detector_count" => length(detectors),
            "observation_count" => length(labels),
            "phase_reference" => reference_name,
            "probe_x_m" => TRANSVERSE_PROBE[1],
            "probe_px_rad" => TRANSVERSE_PROBE[2],
            "probe_y_m" => TRANSVERSE_PROBE[3],
            "probe_py_rad" => TRANSVERSE_PROBE[4],
            "closed_orbit_seconds" => closed_timed.time,
            "stage_seconds" => time() - stage_start,
            "julia_version" => string(VERSION),
        ),
    )
    return nothing
end

safe_name(name) = lowercase(replace(name, r"[^A-Za-z0-9]+" => "_"))

function candidate_path(output_dir, target_name, quadrupole_name, sign_name)
    return joinpath(
        target_directory(output_dir, target_name),
        "candidate_$(safe_name(quadrupole_name))_$(sign_name).npy",
    )
end

function candidate_stage(options, ring, sextupoles, quadrupoles, detectors, reference_name, output_dir)
    target_batch_size = parse(Int, options["candidate-target-batch"])
    overwrite = lowercase(options["overwrite"]) == "true"
    selected_names = requested_targets(options["targets"], sextupoles)
    selected = selected_candidates(output_dir, sextupoles, quadrupoles)
    selected = Dict(name => selected[name] for name in selected_names)
    quadrupole_by_name = Dict(entry.name => entry for entry in quadrupoles)
    target_by_name = Dict(entry.name => entry for entry in sextupoles)
    labels = observation_labels(detectors)
    isfile(joinpath(output_dir, "observation_labels.txt")) ||
        write_lines(joinpath(output_dir, "observation_labels.txt"), labels)

    for target_name in selected_names
        directory = target_directory(output_dir, target_name)
        isfile(joinpath(directory, "target_response_nominal.npy")) ||
            error("Nominal response is missing for $target_name")
        candidates = selected[target_name]
        write_lines(joinpath(directory, "candidate_names.txt"), getproperty.(candidates, :quadrupole))
        write_npy(
            joinpath(directory, "candidate_delta_k1_m2.npy"),
            getproperty.(candidates, :delta_k1_m2),
        )
    end

    targets_by_quadrupole = Dict{String,Vector{String}}()
    for (target_name, candidates) in selected
        for candidate in candidates
            push!(get!(targets_by_quadrupole, candidate.quadrupole, String[]), target_name)
        end
    end
    for names in values(targets_by_quadrupole)
        sort!(names; by=name -> target_by_name[name].s_m)
    end

    parameter_count = 3 * target_batch_size
    descriptor = Descriptor(6, 3, parameter_count, 2)
    parameters = params(descriptor)
    stage_start = time()
    timing_rows = NamedTuple[]
    ordered_quadrupoles = [
        entry for entry in quadrupoles if haskey(targets_by_quadrupole, entry.name)
    ]
    for (quadrupole_counter, quadrupole) in enumerate(ordered_quadrupoles)
        target_names = targets_by_quadrupole[quadrupole.name]
        step = abs(quadrupole.kn1_m2) * parse(Float64, options["quadrupole-fraction"])
        for (sign, sign_name) in ((1.0, "plus"), (-1.0, "minus"))
            pending_names = [
                name for name in target_names
                if overwrite || !isfile(candidate_path(output_dir, name, quadrupole.name, sign_name))
            ]
            isempty(pending_names) && continue
            set_quadrupole!(ring, quadrupole, quadrupole.kn1_m2 + sign * step)
            closed_timed = @timed solve_closed_orbit(ring)
            condition_twiss_seconds = 0.0
            for batch_start in 1:target_batch_size:length(pending_names)
                batch_stop = min(batch_start + target_batch_size - 1, length(pending_names))
                batch_names = pending_names[batch_start:batch_stop]
                batch = [target_by_name[name] for name in batch_names]
                for (slot, target) in enumerate(batch)
                    element = ring.line[target.index]
                    set_target_kn2_parameter!(element, target, parameters[3 * slot - 2])
                    element.x_offset = target.x_offset_m + parameters[3 * slot - 1]
                    element.y_offset = target.y_offset_m + parameters[3 * slot]
                end
                transport_timed = @timed calculate_direct_transport(
                    ring,
                    detectors,
                    closed_timed.value;
                    descriptor=descriptor,
                )
                condition_twiss_seconds += transport_timed.time
                observable_tps = direct_observable_values(transport_timed.value, detectors)
                responses = Dict(name => fill(NaN, length(labels), 2) for name in batch_names)
                for observation_index in eachindex(observable_tps)
                    hessian = GTPSA.hessian(observable_tps[observation_index]; include_params=true)
                    for (slot, target) in enumerate(batch)
                        k2_parameter = 6 + 3 * slot - 2
                        responses[target.name][observation_index, 1] = Float64(hessian[k2_parameter, 6 + 3 * slot - 1])
                        responses[target.name][observation_index, 2] = Float64(hessian[k2_parameter, 6 + 3 * slot])
                    end
                end
                for target in batch
                    all(isfinite, responses[target.name]) || error("Incomplete candidate response for $(target.name)")
                    write_npy(
                        candidate_path(output_dir, target.name, quadrupole.name, sign_name),
                        responses[target.name],
                    )
                    restore_sextupole!(ring, target)
                end
                @printf(
                    "Candidate %d/%d %-10s %s target %d:%d transport %.3f s\n",
                    quadrupole_counter,
                    length(ordered_quadrupoles),
                    quadrupole.name,
                    sign_name,
                    batch_start,
                    batch_stop,
                    transport_timed.time,
                )
                flush(stdout)
            end
            restore_quadrupole!(ring, quadrupole)
            push!(timing_rows, (;
                quadrupole=quadrupole.name,
                sign=sign_name,
                target_count=length(pending_names),
                closed_orbit_seconds=closed_timed.time,
                transport_seconds=condition_twiss_seconds,
            ))
            GC.gc()
        end
    end

    isempty(timing_rows) || write_rows(joinpath(output_dir, "candidate_condition_timings.csv"), timing_rows)
    write_metadata(
        joinpath(output_dir, "candidate_metadata.toml"),
        Dict(
            "format" => "cesr-scibmad-quadrupole-affinity-candidates-v2",
            "engine" => "SciBmad/GTPSA direct launch-trajectory mixed derivatives at fixed quadrupole conditions",
            "date" => string(Dates.today()),
            "lattice" => LATEST_LATTICE,
            "descriptor" => "Descriptor(6, 3, $parameter_count, 2)",
            "candidate_target_batch_size" => target_batch_size,
            "target_count" => length(selected_names),
            "quadrupole_condition_count" => length(ordered_quadrupoles),
            "candidate_count_per_target" => parse(Int, options["top-k"]),
            "detector_count" => length(detectors),
            "observation_count" => length(labels),
            "phase_reference" => reference_name,
            "quadrupole_fraction" => parse(Float64, options["quadrupole-fraction"]),
            "nuisance_policy" => "Nominal 150-column nuisance dictionary is reused in candidate +/- blocks",
            "stage_seconds" => time() - stage_start,
            "julia_version" => string(VERSION),
        ),
    )
    return nothing
end

function main(args=ARGS)
    options = parse_options(args)
    output_dir = abspath(options["output-dir"])
    mkpath(output_dir)
    ring = cesr
    sextupoles = active_sextupole_inventory(ring)
    quadrupoles = active_quadrupole_inventory(ring)
    detectors = measurable_bpms(ring)
    reference_name = phase_reference_name(options["phase-reference"], detectors)
    stage = options["stage"]
    @printf(
        "SciBmad affinity stage=%s sextupoles=%d quadrupoles=%d BPMs=%d reference=%s\n",
        stage,
        length(sextupoles),
        length(quadrupoles),
        length(detectors),
        reference_name,
    )
    if stage == "screen"
        screen_stage(options, ring, sextupoles, quadrupoles, detectors, reference_name, output_dir)
    elseif stage == "nominal"
        nominal_stage(options, ring, sextupoles, detectors, reference_name, output_dir)
    else
        candidate_stage(options, ring, sextupoles, quadrupoles, detectors, reference_name, output_dir)
    end
    return 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main())
end
