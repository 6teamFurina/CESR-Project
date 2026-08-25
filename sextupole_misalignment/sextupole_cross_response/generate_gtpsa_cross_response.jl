#!/usr/bin/env julia

"""Build compact all-sextupole orbit-propagation response matrices.

The calculation deliberately avoids a descriptor containing every sextupole's
K2, bump, and alignment parameters.  It uses two small GTPSA calculations:

1. `Descriptor(6, 1)` gives nominal cumulative maps and the periodic response
   from a unit local px/py kick at each active sextupole to the orbit at every
   active sextupole.
2. `Descriptor(6, 2, Ncorrector, 1)` gives only the first derivatives of the
   sextupole orbits with respect to the independent correctors.  Multiplying
   by the saved two-plane bump knobs gives the bump cross-response.

The exact local normal-sextupole polynomial is then composed analytically with
the periodic kick response.  This produces the selected K2--bump--center
derivative needed by the alignment inverse without constructing unused global
Hessian or third-derivative entries.
"""

include(joinpath(
    @__DIR__, "..", "quadrupole_affinity", "exact_11_triplet_validation", "common.jl",
))

const DEFAULT_BUMP_KNOBS = joinpath(
    @__DIR__, "..", "quadrupole_affinity", "exact_11_triplet_validation",
    "results", "bump_knobs", "local_bump_knobs.csv",
)

function map_jacobian(map)
    return [
        Float64(constant_term(GTPSA.deriv(map.v[i], j)))
        for i in 1:6, j in 1:6
    ]
end

function bump_knob_array(path, sextupoles, controls)
    rows = read_simple_csv(path)
    array = fill(NaN, length(sextupoles), 2, length(controls))
    target_index = Dict(target.name => index for (index, target) in enumerate(sextupoles))
    control_index = Dict(
        (control.name, String(control.axis)) => index
        for (index, control) in enumerate(controls)
    )
    seen = Set{Tuple{Int,Int}}()
    for row in rows
        target = uppercase(row["target_sextupole"])
        key = (row["corrector"], row["field"])
        haskey(target_index, target) || error("Unknown bump-knob target: $target")
        haskey(control_index, key) || error("Unknown bump-knob control: $key")
        ti, ci = target_index[target], control_index[key]
        (ti, ci) in seen && error("Duplicate bump-knob row for $target and $key")
        push!(seen, (ti, ci))
        array[ti, 1, ci] = parse(Float64, row["field_per_x_bump_m"])
        array[ti, 2, ci] = parse(Float64, row["field_per_y_bump_m"])
    end
    expected = length(sextupoles) * length(controls)
    length(seen) == expected || error("Expected $expected bump-knob rows, found $(length(seen))")
    all(isfinite, array) || error("Bump-knob array contains non-finite values")
    return array
end

function periodic_kick_response(target_maps, one_turn_map, line_indices)
    target_count = size(target_maps, 1)
    response = zeros(target_count, target_count, 2, 2)
    closure = inv(I - one_turn_map)
    kicks = (Float64[0, 1, 0, 0, 0, 0], Float64[0, 0, 0, 1, 0, 0])
    for source in 1:target_count
        inverse_source = inv(target_maps[source, :, :])
        for kick_plane in 1:2
            kick = kicks[kick_plane]
            one_turn_source = one_turn_map * inverse_source * kick
            start_closed = closure * one_turn_source
            for observation in 1:target_count
                propagated = target_maps[observation, :, :] * start_closed
                if line_indices[observation] > line_indices[source]
                    propagated += target_maps[observation, :, :] * inverse_source * kick
                end
                response[source, observation, 1, kick_plane] = propagated[1]
                response[source, observation, 2, kick_plane] = propagated[3]
            end
        end
    end
    return response
end

function target_control_response(ring, sextupoles, controls, closed)
    descriptor = Descriptor(6, 2, length(controls), 1)
    parameters = params(descriptor)
    for (control_index, control) in enumerate(controls)
        baseline = Float64(constant_term(first(control.originals)))
        for element_index in control.indices
            if control.axis == :Kn0
                ring.line[element_index].Kn0 = baseline + parameters[control_index]
            else
                ring.line[element_index].Ks0 = baseline + parameters[control_index]
            end
        end
    end

    at = [ring.line[target.index] for target in sextupoles]
    timed = @timed calculate_twiss(ring, at, closed; descriptor)
    optics = timed.value
    index_by_name = table_index_by_name(optics)
    result = zeros(length(sextupoles), 2, length(controls))
    for (target_index, target) in enumerate(sextupoles)
        table_index = index_by_name[target.name]
        values = [optics.table.orbit_x[table_index], optics.table.orbit_y[table_index]]
        jacobian = Matrix(GTPSA.jacobian(values; include_params=true))
        result[target_index, :, :] .= jacobian[:, 7:end]
    end
    return result, timed.time, descriptor
end

function main(args=ARGS)
    defaults = Dict(
        "bump-knobs-csv" => DEFAULT_BUMP_KNOBS,
        "output-dir" => joinpath(@__DIR__, "results", "raw"),
        "overwrite" => "false",
    )
    options = parse_exact11_options(defaults, args)
    lowercase(options["overwrite"]) in ("true", "false") ||
        error("--overwrite must be true or false")
    output_dir = abspath(options["output-dir"])
    metadata_path = joinpath(output_dir, "response_metadata.toml")
    if isfile(metadata_path) && lowercase(options["overwrite"]) != "true"
        println("Existing response retained: $metadata_path")
        return 0
    end

    ring = cesr
    sextupoles = active_sextupole_inventory(ring)
    controls = independent_corrector_inventory(ring)
    target_names = String.(getproperty.(sextupoles, :name))
    line_indices = Int.(getproperty.(sextupoles, :index))
    bump_knobs = bump_knob_array(abspath(options["bump-knobs-csv"]), sextupoles, controls)

    closed_timed = @timed solve_closed_orbit(ring)
    closed = closed_timed.value
    nominal = track_orbits_at_names(ring, closed, target_names)

    transport_at = vcat(
        [ring.line[target.index] for target in sextupoles],
        [ring.line[end]],
    )
    transport_descriptor = Descriptor(6, 1)
    transport_timed = @timed calculate_direct_transport(
        ring, transport_at, closed; descriptor=transport_descriptor,
    )
    transport = transport_timed.value
    transport_index = Dict(
        first(split(uppercase(String(name)), '!')) => index
        for (index, name) in enumerate(transport.names)
    )
    target_maps = zeros(length(sextupoles), 6, 6)
    for (index, name) in enumerate(target_names)
        target_maps[index, :, :] .= map_jacobian(transport.maps[transport_index[name]])
    end
    one_turn_map = map_jacobian(transport.maps[end])
    kick_response = periodic_kick_response(target_maps, one_turn_map, line_indices)

    corrector_response, corrector_seconds, corrector_descriptor =
        target_control_response(ring, sextupoles, controls, closed)
    bump_response = zeros(length(sextupoles), length(sextupoles), 2, 2)
    for target in eachindex(sextupoles), observation in eachindex(sextupoles), axis in 1:2
        bump_response[target, observation, :, axis] .=
            corrector_response[observation, :, :] * bump_knobs[target, axis, :]
    end

    # Source 1 multiplies q_n = 0.5*(x^2-y^2) and gives -L*q_n in px.
    # Source 2 multiplies q_s = x*y and gives +L*q_s in py.
    source_response = zeros(size(kick_response))
    for target in eachindex(sextupoles)
        length_m = sextupoles[target].length_m
        source_response[target, :, :, 1] .= -length_m .* kick_response[target, :, :, 1]
        source_response[target, :, :, 2] .= +length_m .* kick_response[target, :, :, 2]
    end

    local_bump_jacobian = zeros(length(sextupoles), 2, 2)
    alignment_design = zeros(length(sextupoles), 2, length(sextupoles), 2, 2)
    for target in eachindex(sextupoles)
        local_bump_jacobian[target, :, :] .= bump_response[target, target, :, :]
        normal = source_response[target, :, :, 1]
        skew = source_response[target, :, :, 2]
        for bump_axis in 1:2
            dx = local_bump_jacobian[target, 1, bump_axis]
            dy = local_bump_jacobian[target, 2, bump_axis]
            alignment_design[target, bump_axis, :, :, 1] .= -dx .* normal .- dy .* skew
            alignment_design[target, bump_axis, :, :, 2] .= +dy .* normal .- dx .* skew
        end
    end

    all(isfinite, kick_response) || error("Non-finite periodic kick response")
    all(isfinite, bump_response) || error("Non-finite bump response")
    all(isfinite, source_response) || error("Non-finite sextupole source response")
    all(isfinite, alignment_design) || error("Non-finite alignment design")
    identity_error = maximum(abs, local_bump_jacobian .- reshape(Matrix{Float64}(I, 2, 2), 1, 2, 2))
    identity_error < 1.0e-8 || error("Target-local bump Jacobian is not identity: $identity_error")

    target_rows = NamedTuple[]
    for (index, target) in enumerate(sextupoles)
        push!(target_rows, (;
            target_index=index,
            target=target.name,
            line_index=target.index,
            s_m=target.s_m,
            length_m=target.length_m,
            nominal_kn2_m3=target.kn2_m3,
            nominal_orbit_x_m=nominal.horizontal[target.name][1],
            nominal_orbit_y_m=nominal.vertical[target.name][1],
        ))
    end
    control_rows = [
        (;
            control_index=index,
            corrector=control.name,
            field=String(control.axis),
            total_length_m=control.total_length_m,
        )
        for (index, control) in enumerate(controls)
    ]
    response_rows = NamedTuple[]
    plane_names = ("x", "y")
    for target in eachindex(sextupoles), observation in eachindex(sextupoles), plane in 1:2
        push!(response_rows, (;
            source_target_index=target,
            source_target=sextupoles[target].name,
            observation_index=observation,
            observation_sextupole=sextupoles[observation].name,
            observation_plane=plane_names[plane],
            kick_px_response=kick_response[target, observation, plane, 1],
            kick_py_response=kick_response[target, observation, plane, 2],
            bump_x_response_m_per_m=bump_response[target, observation, plane, 1],
            bump_y_response_m_per_m=bump_response[target, observation, plane, 2],
            normal_source_response=source_response[target, observation, plane, 1],
            skew_source_response=source_response[target, observation, plane, 2],
            k2_bx_center_x_response=alignment_design[target, 1, observation, plane, 1],
            k2_bx_center_y_response=alignment_design[target, 1, observation, plane, 2],
            k2_by_center_x_response=alignment_design[target, 2, observation, plane, 1],
            k2_by_center_y_response=alignment_design[target, 2, observation, plane, 2],
        ))
    end

    mkpath(output_dir)
    write_npy(joinpath(output_dir, "target_cumulative_maps.npy"), target_maps)
    write_npy(joinpath(output_dir, "one_turn_map.npy"), one_turn_map)
    write_npy(joinpath(output_dir, "target_control_response.npy"), corrector_response)
    write_npy(joinpath(output_dir, "bump_knobs.npy"), bump_knobs)
    write_npy(joinpath(output_dir, "periodic_kick_response.npy"), kick_response)
    write_npy(joinpath(output_dir, "bump_response.npy"), bump_response)
    write_npy(joinpath(output_dir, "sextupole_source_response.npy"), source_response)
    write_npy(joinpath(output_dir, "local_bump_jacobian.npy"), local_bump_jacobian)
    write_npy(joinpath(output_dir, "alignment_design.npy"), alignment_design)
    write_rows(joinpath(output_dir, "target_inventory.csv"), target_rows)
    write_rows(joinpath(output_dir, "control_inventory.csv"), control_rows)
    write_rows(joinpath(output_dir, "cross_response_long.csv"), response_rows)
    write_metadata(metadata_path, Dict(
        "format" => "cesr-sextupole-cross-response-v1",
        "date" => string(Dates.today()),
        "engine" => "SciBmad/GTPSA nominal RF-on periodic response with analytic normal-sextupole source composition",
        "lattice" => LATEST_LATTICE,
        "target_count" => length(sextupoles),
        "observation_count" => length(sextupoles),
        "corrector_count" => length(controls),
        "coordinate_order" => ["x", "y"],
        "bump_axis_order" => ["x", "y"],
        "kick_axis_order" => ["px", "py"],
        "source_order" => ["normal_qn", "skew_qs"],
        "center_axis_order" => ["x", "y"],
        "periodic_kick_response_axes" => ["source_target", "observation_sextupole", "output_plane", "kick_plane"],
        "bump_response_axes" => ["source_target", "observation_sextupole", "output_plane", "bump_axis"],
        "source_response_axes" => ["source_target", "observation_sextupole", "output_plane", "source"],
        "alignment_design_axes" => ["source_target", "bump_axis", "observation_sextupole", "output_plane", "center_axis"],
        "transport_descriptor" => string(transport_descriptor),
        "corrector_response_descriptor" => string(corrector_descriptor),
        "derivative_policy" => "No global Hessian is formed. The selected K2-bump-center derivative is the exact local sextupole polynomial composed with order-1 periodic GTPSA transport and the first-order target-local bump Jacobian.",
        "location_convention" => "first element-entry occurrence by base name; a source kick is applied immediately after the source entry",
        "thin_source_qualification" => "The source composition applies the integrated normal-sextupole kick at the element entry; finite-length and nonlinear propagation require exact-scan validation.",
        "bump_knobs_csv" => abspath(options["bump-knobs-csv"]),
        "target_local_bump_identity_max_abs_error" => identity_error,
        "nominal_closed_orbit_seconds" => closed_timed.time,
        "transport_seconds" => transport_timed.time,
        "corrector_response_seconds" => corrector_seconds,
        "julia_version" => string(VERSION),
    ))
    @printf(
        "Cross response complete: %d targets, %d observations, %d correctors\n",
        length(sextupoles), length(sextupoles), length(controls),
    )
    @printf(
        "Timings: closed %.3f s, transport %.3f s, corrector response %.3f s\n",
        closed_timed.time, transport_timed.time, corrector_seconds,
    )
    println("Output: $output_dir")
    return 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main())
end
