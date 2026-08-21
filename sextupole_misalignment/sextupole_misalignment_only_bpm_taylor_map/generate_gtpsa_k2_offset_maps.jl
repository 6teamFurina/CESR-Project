#!/usr/bin/env julia

"""Generate latest-lattice high-order GTPSA K2--offset BPM-orbit maps.

For each requested target, the three parameters are delta Kn2, delta target
x_offset, and delta target y_offset.  The saved tensor contains true mixed
derivatives

    d^(1+a+b) O / dKn2 dx_offset^a dy_offset^b

for every a+b <= 3 and every measurable BPM x/y channel.  The corresponding
Taylor coefficient is the derivative divided by a! b!.
"""

include(joinpath(
    @__DIR__, "..", "quadrupole_affinity", "exact_11_triplet_validation", "common.jl",
))

function selected_target_entries(requested, sextupoles)
    text = lowercase(strip(requested))
    text == "all" && return collect(enumerate(sextupoles))
    wanted = Set(uppercase(strip(name)) for name in split(requested, ',') if !isempty(strip(name)))
    known = Set(getproperty.(sextupoles, :name))
    missing = setdiff(wanted, known)
    isempty(missing) || error("Unknown target sextupoles: $(join(sort!(collect(missing)), ", "))")
    return [(index, entry) for (index, entry) in enumerate(sextupoles) if entry.name in wanted]
end

function derivative_value(value, k2_index, x_index, y_index, x_power, y_power)
    result = GTPSA.deriv(value, k2_index)
    for _ in 1:x_power
        result = GTPSA.deriv(result, x_index)
    end
    for _ in 1:y_power
        result = GTPSA.deriv(result, y_index)
    end
    return Float64(constant_term(result))
end

function parameter_only(map)
    return map ∘ zero(map)
end

function fixed_point_residual(one_turn, fixed_point)
    image = parameter_only(one_turn ∘ fixed_point)
    origin = parameter_only(fixed_point)
    nv = SciBmad.NNF.nvars(one_turn)
    return image.v[1:nv] .- origin.v[1:nv]
end

function residual_norm(residual)
    return sum(GTPSA.normTPS(value) for value in residual)
end

"""Solve z*(p) = M(z*(p),p) without a transverse normal-form gauge.

`normal(one_turn, 0)` performs only the analytic first-order parameter-dependent
fixed-point solve.  Each defect-correction iteration then removes the next
parameter order using the constant closed-orbit Jacobian.  This avoids the
linear eigenvector normalization used by full periodic Twiss, whose gauge can
be non-analytic at exactly zero coupling even when the closed orbit is regular.
"""
function parameter_dependent_fixed_point(one_turn, parameter_order)
    fixed_point = normal(one_turn, 0)
    linear = Matrix(SciBmad.NNF.jacobian(one_turn, SciBmad.NNF.HVARS))
    correction_operator = inv(I - linear)
    history = Float64[]
    for _ in 1:(parameter_order + 1)
        residual = fixed_point_residual(one_turn, fixed_point)
        push!(history, residual_norm(residual))
        correction = correction_operator * residual
        for coordinate in eachindex(correction)
            GTPSA.add!(
                fixed_point.v[coordinate],
                fixed_point.v[coordinate],
                correction[coordinate],
            )
        end
    end
    push!(history, residual_norm(fixed_point_residual(one_turn, fixed_point)))
    return fixed_point, history
end

function direct_fixed_point_orbits(ring, detectors, closed, descriptor, parameter_order)
    transport = calculate_direct_transport(
        ring,
        vcat(detectors, [ring.line[end]]),
        closed;
        descriptor,
    )
    transport_index = Dict(
        first(split(uppercase(String(name)), '!')) => index
        for (index, name) in enumerate(transport.names)
    )
    fixed_point, residual_history = parameter_dependent_fixed_point(
        transport.maps[end],
        parameter_order,
    )
    orbit_values = Matrix{Any}(undef, length(detectors), 2)
    for (detector_index, detector) in enumerate(detectors)
        name = uppercase(String(base_name(detector)))
        local_map = parameter_only(transport.maps[transport_index[name]] ∘ fixed_point)
        orbit_values[detector_index, 1] = local_map.v[1]
        orbit_values[detector_index, 2] = local_map.v[3]
    end
    return orbit_values, residual_history
end

function main(args=ARGS)
    defaults = Dict(
        "targets" => "all",
        "maximum-offset-order" => "3",
        "periodic-solver" => "fixed-point",
        "output-dir" => joinpath(@__DIR__, "results", "gtpsa_maps"),
        "overwrite" => "false",
    )
    options = parse_exact11_options(defaults, args)
    maximum_offset_order = parse(Int, options["maximum-offset-order"])
    1 <= maximum_offset_order <= 3 || error("--maximum-offset-order must be 1, 2, or 3")
    periodic_solver = lowercase(strip(options["periodic-solver"]))
    periodic_solver in ("fixed-point", "twiss-normal-form") ||
        error("--periodic-solver must be fixed-point or twiss-normal-form")
    output_dir = abspath(options["output-dir"])
    metadata_path = joinpath(output_dir, "map_metadata.toml")
    isfile(metadata_path) && lowercase(options["overwrite"]) != "true" &&
        error("Output exists; use --overwrite=true: $metadata_path")

    ring = cesr
    sextupoles = active_sextupole_inventory(ring)
    detectors = measurable_bpms(ring)
    selected = selected_target_entries(options["targets"], sextupoles)
    target_names = [entry.name for (_, entry) in selected]
    target_inventory_indices = [index for (index, _) in selected]
    bpm_names = String.(base_name.(detectors))
    nt, nd = length(selected), length(detectors)

    parameter_order = maximum_offset_order + 1
    descriptor = Descriptor(6, parameter_order, 3, parameter_order)
    parameters = params(descriptor)
    k2_parameter_index, x_parameter_index, y_parameter_index = 7, 8, 9
    monomials = [
        (x_power, total - x_power)
        for total in 0:maximum_offset_order
        for x_power in 0:total
    ]
    derivatives = fill(NaN, nt, length(monomials), 2 * nd)
    timing_rows = NamedTuple[]
    closed_timed = @timed solve_closed_orbit(ring)
    closed = closed_timed.value
    calculation_start = time()

    try
        for (target_counter, (_, target)) in enumerate(selected)
            element = ring.line[target.index]
            element.Kn2 = target.kn2_m3 + parameters[1]
            element.x_offset = target.x_offset_m + parameters[2]
            element.y_offset = target.y_offset_m + parameters[3]
            map_timed = if periodic_solver == "fixed-point"
                @timed direct_fixed_point_orbits(
                    ring,
                    detectors,
                    closed,
                    descriptor,
                    parameter_order,
                )
            else
                @timed begin
                    optics = calculate_twiss(
                        ring,
                        detectors,
                        closed;
                        descriptor=descriptor,
                    )
                    values = Matrix{Any}(undef, nd, 2)
                    values[:, 1] .= optics.table.orbit_x
                    values[:, 2] .= optics.table.orbit_y
                    (values, Float64[])
                end
            end
            orbit_values, residual_history = map_timed.value
            extraction_start = time()
            for detector_index in 1:nd
                for plane_index in 1:2
                    value = orbit_values[detector_index, plane_index]
                    channel = 2 * (detector_index - 1) + plane_index
                    for (monomial_index, (x_power, y_power)) in enumerate(monomials)
                        derivatives[target_counter, monomial_index, channel] = derivative_value(
                            value,
                            k2_parameter_index,
                            x_parameter_index,
                            y_parameter_index,
                            x_power,
                            y_power,
                        )
                    end
                end
            end
            all(isfinite, derivatives[target_counter, :, :]) ||
                error("Non-finite high-order map derivative for $(target.name)")
            extraction_seconds = time() - extraction_start
            push!(timing_rows, (;
                target_index=target_counter,
                sextupole_inventory_index=target_inventory_indices[target_counter],
                target=target.name,
                map_seconds=map_timed.time,
                extraction_seconds,
                initial_fixed_point_residual=isempty(residual_history) ? "" : first(residual_history),
                final_fixed_point_residual=isempty(residual_history) ? "" : last(residual_history),
            ))
            restore_sextupole!(ring, target)
            @printf(
                "%s GTPSA map %d/%d: %s %.3f s, extraction %.3f s\n",
                target.name,
                target_counter,
                nt,
                periodic_solver,
                map_timed.time,
                extraction_seconds,
            )
            flush(stdout)
            GC.gc()
        end
    finally
        for entry in sextupoles
            restore_sextupole!(ring, entry)
        end
    end

    mkpath(output_dir)
    write_npy(joinpath(output_dir, "k2_offset_derivatives.npy"), derivatives)
    write_lines(joinpath(output_dir, "target_names.txt"), target_names)
    write_lines(joinpath(output_dir, "bpm_names.txt"), bpm_names)
    write_rows(joinpath(output_dir, "offset_monomials.csv"), [
        (;
            monomial_index=index,
            x_offset_power=x_power,
            y_offset_power=y_power,
            factorial_divisor=factorial(x_power) * factorial(y_power),
        )
        for (index, (x_power, y_power)) in enumerate(monomials)
    ])
    write_rows(joinpath(output_dir, "map_timings.csv"), timing_rows)
    write_metadata(metadata_path, Dict(
        "format" => "cesr-latest-gtpsa-k2-offset-high-order-map-v1",
        "date" => string(Dates.today()),
        "engine" => periodic_solver == "fixed-point" ?
            "SciBmad/GTPSA direct one-turn map with order-by-order parameter-dependent fixed point" :
            "SciBmad/GTPSA parameterized periodic Twiss BPM orbit",
        "periodic_solver" => periodic_solver,
        "lattice" => LATEST_LATTICE,
        "descriptor" => "Descriptor(6, $parameter_order, 3, $parameter_order)",
        "parameters" => ["delta_Kn2_m3", "delta_x_offset_m", "delta_y_offset_m"],
        "maximum_offset_order" => maximum_offset_order,
        "maximum_total_parameter_order" => parameter_order,
        "derivative_semantics" => "saved values are true derivatives d^(1+a+b)O/dKn2/dx^a/dy^b; divide by a!b! for Taylor evaluation",
        "target_count" => nt,
        "target_inventory" => target_names,
        "bpm_count" => nd,
        "channel_order" => "BPM inventory order with x,y interleaved",
        "nominal_closed_orbit_seconds" => closed_timed.time,
        "calculation_wall_seconds" => time() - calculation_start,
    ))
    println("Wrote high-order GTPSA maps to $output_dir")
    return 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main())
end
