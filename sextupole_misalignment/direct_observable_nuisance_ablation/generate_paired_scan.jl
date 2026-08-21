#!/usr/bin/env julia

"""Generate a paired, nuisance-hidden sextupole-center scan with direct readbacks.

The stored transport maps are not reported as ideal Twiss observables.  They are
the raw ingredients used by `analyze_paired_scan.py` to emulate TBT BPM phase,
amplitude, spectral-line, tune, fixed-energy dispersion, and chromatic-tune
measurements.  Two actual CESR correctors are additionally finite-differenced
to produce a closed-orbit response measurement.
"""

include(joinpath(
    @__DIR__, "..", "quadrupole_affinity", "exact_11_triplet_validation", "common.jl",
))

const DEFAULT_BUMP_KNOBS = joinpath(
    @__DIR__, "..", "quadrupole_affinity", "exact_11_triplet_validation",
    "results", "bump_knobs", "local_bump_knobs.csv",
)

function scalar_warm(ring, v0)
    solution = find_closed_orbit(
        ring; v0=copy(v0), coasting_beam=false, batch=Val{false}(), warn=false,
    )
    all(solution.sol.retcode .== SciBmad.BatchSolve.RETCODE_SUCCESS) ||
        error("RF-on closed-orbit solve failed: $(solution.sol.retcode)")
    return solution
end

grid3(amplitude) = [
    (x, y) for x in (-amplitude, 0.0, amplitude)
           for y in (-amplitude, 0.0, amplitude)
]

function map_jacobian(map, output_coordinates, input_coordinates)
    return [
        Float64(constant_term(GTPSA.deriv(map.v[i], j)))
        for i in output_coordinates, j in input_coordinates
    ]
end

function map_hessian(map, output_coordinates, input_coordinates)
    return [
        Float64(constant_term(GTPSA.deriv(GTPSA.deriv(map.v[i], j), k)))
        for i in output_coordinates, j in input_coordinates, k in input_coordinates
    ]
end

function corrector_value(ring, control)
    element = ring.line[first(control.indices)]
    return Float64(constant_term(control.axis == :Kn0 ? element.Kn0 : element.Ks0))
end

function restore_correctors!(ring, controls)
    for control in controls
        for (index, original) in zip(control.indices, control.originals)
            if control.axis == :Kn0
                ring.line[index].Kn0 = original
            else
                ring.line[index].Ks0 = original
            end
        end
    end
end

function orbit_matrix(ring, closed, bpm_names)
    tracked = track_orbits_at_names(ring, closed, bpm_names)
    result = zeros(length(bpm_names), 2)
    for (i, name) in enumerate(bpm_names)
        result[i, 1] = tracked.horizontal[name][1]
        result[i, 2] = tracked.vertical[name][1]
    end
    return result
end

function choose_probe_correctors(controls, knob_by_key)
    ix = argmax([
        abs(knob_by_key[(control.name, String(control.axis))][1]) for control in controls
    ])
    iy = argmax([
        abs(knob_by_key[(control.name, String(control.axis))][2]) for control in controls
    ])
    return controls[[ix, iy]]
end

function main(args=ARGS)
    defaults = Dict(
        "target" => "SEX_09AW",
        "realizations" => "8",
        "seed" => "20260817",
        "target-halfwidth-m" => "3.5e-4",
        "other-sext-rms-m" => "3.0e-4",
        "quadrupole-fraction" => "0.01",
        "bump-amplitude-m" => "5.0e-4",
        "k2-step-m3" => "0.01",
        "corrector-probe-field" => "1.0e-6",
        "bump-knobs-csv" => DEFAULT_BUMP_KNOBS,
        "output-dir" => joinpath(@__DIR__, "results", "sex_09aw_paired_pilot"),
        "overwrite" => "false",
    )
    options = parse_exact11_options(defaults, args)
    target_name = uppercase(options["target"])
    realization_count = parse(Int, options["realizations"])
    seed = parse(Int, options["seed"])
    target_halfwidth = parse(Float64, options["target-halfwidth-m"])
    other_sext_rms = parse(Float64, options["other-sext-rms-m"])
    quad_fraction = parse(Float64, options["quadrupole-fraction"])
    bump_amplitude = parse(Float64, options["bump-amplitude-m"])
    k2_step = parse(Float64, options["k2-step-m3"])
    corrector_step = parse(Float64, options["corrector-probe-field"])
    output_dir = abspath(options["output-dir"])
    metadata_path = joinpath(output_dir, "scan_metadata.toml")
    isfile(metadata_path) && lowercase(options["overwrite"]) != "true" &&
        error("Output exists; use --overwrite=true: $metadata_path")

    ring = cesr
    sextupoles = active_sextupole_inventory(ring)
    quadrupoles = active_quadrupole_inventory(ring)
    detectors = measurable_bpms(ring)
    bpm_names = String.(base_name.(detectors))
    target_index = findfirst(entry -> entry.name == target_name, sextupoles)
    isnothing(target_index) && error("Unknown target sextupole: $target_name")
    target_entry = sextupoles[target_index]
    target_element = ring.line[target_entry.index]
    length(sextupoles) == 76 || error("Expected 76 active sextupoles")

    controls = independent_corrector_inventory(ring)
    knob_rows = read_bump_knobs(abspath(options["bump-knobs-csv"]), target_name)
    knob_by_key = Dict(
        (row["corrector"], row["field"]) => (
            parse(Float64, row["field_per_x_bump_m"]),
            parse(Float64, row["field_per_y_bump_m"]),
        ) for row in knob_rows
    )
    Set(keys(knob_by_key)) == Set((c.name, String(c.axis)) for c in controls) ||
        error("Bump knob and corrector inventory mismatch")
    probe_controls = choose_probe_correctors(controls, knob_by_key)

    bumps = grid3(bump_amplitude)
    k2_levels = [-2.0, -1.0, 0.0, 1.0, 2.0]
    nr, nb, nk, nd = realization_count, length(bumps), length(k2_levels), length(detectors)
    orbit = zeros(nr, nb, nk, nd, 2)
    target_orbit = zeros(nr, nb, nk, 2)
    bpm_jacobian = zeros(nr, nb, nk, nd, 2, 6)
    one_turn_jacobian = zeros(nr, nb, nk, 6, 6)
    one_turn_hessian = zeros(nr, nb, nk, 4, 6, 6)
    orm = zeros(nr, nb, nk, length(probe_controls), nd, 2)
    target_truth = zeros(nr, 2)
    sext_offsets = zeros(nr, length(sextupoles), 2)
    quadrupole_errors = zeros(nr, length(quadrupoles))
    timings = NamedTuple[]
    rng = MersenneTwister(seed)

    nominal_closed = solve_closed_orbit(ring)
    nominal_v0 = copy(nominal_closed.v0)
    transport_at = vcat(detectors, [ring.line[end]])

    try
        for realization in 1:nr
            target_truth[realization, :] .= (
                (2rand(rng) - 1) * target_halfwidth,
                (2rand(rng) - 1) * target_halfwidth,
            )
            for (sindex, entry) in enumerate(sextupoles)
                dx, dy = if sindex == target_index
                    target_truth[realization, 1], target_truth[realization, 2]
                else
                    other_sext_rms * randn(rng), other_sext_rms * randn(rng)
                end
                sext_offsets[realization, sindex, :] .= (dx, dy)
                element = ring.line[entry.index]
                element.x_offset = entry.x_offset_m + dx
                element.y_offset = entry.y_offset_m + dy
            end
            for (qindex, entry) in enumerate(quadrupoles)
                relative_error = (2rand(rng) - 1) * quad_fraction
                quadrupole_errors[realization, qindex] = relative_error
                set_quadrupole!(ring, entry, entry.kn1_m2 * (1 + relative_error))
            end

            current_v0 = copy(nominal_v0)
            realization_start = time()
            for (bump_index, (bump_x, bump_y)) in enumerate(bumps)
                for control in controls
                    cx, cy = knob_by_key[(control.name, String(control.axis))]
                    set_corrector_scalar!(
                        ring, control,
                        Float64(constant_term(first(control.originals))) + cx*bump_x + cy*bump_y,
                    )
                end
                for (k2_index, level) in enumerate(k2_levels)
                    state_start = time()
                    target_element.Kn2 = target_entry.kn2_m3 + level*k2_step
                    solve_timed = @timed scalar_warm(ring, current_v0)
                    closed = solve_timed.value
                    current_v0 .= closed.v0
                    requested = vcat(bpm_names, [target_name])
                    track_timed = @timed track_orbits_at_names(ring, closed, requested)
                    tracked = track_timed.value
                    for (i, name) in enumerate(bpm_names)
                        orbit[realization, bump_index, k2_index, i, 1] = tracked.horizontal[name][1]
                        orbit[realization, bump_index, k2_index, i, 2] = tracked.vertical[name][1]
                    end
                    target_orbit[realization, bump_index, k2_index, :] .= (
                        tracked.horizontal[target_name][1], tracked.vertical[target_name][1],
                    )

                    transport_timed = @timed calculate_direct_transport(
                        ring, transport_at, closed; descriptor=Descriptor(6, 2),
                    )
                    transport = transport_timed.value
                    index_by_name = Dict(
                        first(split(uppercase(String(name)), '!')) => i
                        for (i, name) in enumerate(transport.names)
                    )
                    for (i, name) in enumerate(bpm_names)
                        map = transport.maps[index_by_name[name]]
                        bpm_jacobian[realization, bump_index, k2_index, i, :, :] .=
                            map_jacobian(map, (1, 3), 1:6)
                    end
                    one_turn = transport.maps[end]
                    one_turn_jacobian[realization, bump_index, k2_index, :, :] .=
                        map_jacobian(one_turn, 1:6, 1:6)
                    one_turn_hessian[realization, bump_index, k2_index, :, :, :] .=
                        map_hessian(one_turn, 1:4, 1:6)

                    orm_seconds = 0.0
                    for (probe_index, control) in enumerate(probe_controls)
                        baseline = corrector_value(ring, control)
                        plus_timed = @timed begin
                            set_corrector_scalar!(ring, control, baseline + corrector_step)
                            plus_closed = scalar_warm(ring, closed.v0)
                            orbit_matrix(ring, plus_closed, bpm_names)
                        end
                        minus_timed = @timed begin
                            set_corrector_scalar!(ring, control, baseline - corrector_step)
                            minus_closed = scalar_warm(ring, closed.v0)
                            orbit_matrix(ring, minus_closed, bpm_names)
                        end
                        set_corrector_scalar!(ring, control, baseline)
                        orm[realization, bump_index, k2_index, probe_index, :, :] .=
                            (plus_timed.value .- minus_timed.value) ./ (2corrector_step)
                        orm_seconds += plus_timed.time + minus_timed.time
                    end
                    push!(timings, (;
                        realization, bump_index, k2_index,
                        closed_orbit_seconds=solve_timed.time,
                        tracking_seconds=track_timed.time,
                        direct_transport_seconds=transport_timed.time,
                        orm_probe_seconds=orm_seconds,
                        total_state_seconds=time() - state_start,
                    ))
                end
            end
            @printf(
                "%s realization %d/%d complete in %.1f s\n",
                target_name, realization, nr, time() - realization_start,
            )
            flush(stdout)
        end
    finally
        restore_correctors!(ring, controls)
        for entry in quadrupoles
            restore_quadrupole!(ring, entry)
        end
        for entry in sextupoles
            restore_sextupole!(ring, entry)
        end
    end

    mkpath(output_dir)
    write_npy(joinpath(output_dir, "bpm_orbits.npy"), orbit)
    write_npy(joinpath(output_dir, "target_orbits.npy"), target_orbit)
    write_npy(joinpath(output_dir, "bpm_jacobians.npy"), bpm_jacobian)
    write_npy(joinpath(output_dir, "one_turn_jacobians.npy"), one_turn_jacobian)
    write_npy(joinpath(output_dir, "one_turn_hessians.npy"), one_turn_hessian)
    write_npy(joinpath(output_dir, "corrector_orm.npy"), orm)
    write_npy(joinpath(output_dir, "target_truth.npy"), target_truth)
    write_npy(joinpath(output_dir, "sextupole_offsets.npy"), sext_offsets)
    write_npy(joinpath(output_dir, "quadrupole_relative_errors.npy"), quadrupole_errors)
    write_lines(joinpath(output_dir, "bpm_names.txt"), bpm_names)
    write_lines(joinpath(output_dir, "sextupole_names.txt"), getproperty.(sextupoles, :name))
    write_lines(joinpath(output_dir, "quadrupole_names.txt"), getproperty.(quadrupoles, :name))
    write_lines(joinpath(output_dir, "corrector_probe_names.txt"), getproperty.(probe_controls, :name))
    write_rows(joinpath(output_dir, "timings.csv"), timings)
    write_rows(joinpath(output_dir, "bump_points.csv"), [
        (; bump_index=i, bump_x_command_m=x, bump_y_command_m=y)
        for (i, (x, y)) in enumerate(bumps)
    ])
    write_metadata(metadata_path, Dict(
        "format" => "cesr-direct-observable-nuisance-ablation-v1",
        "date" => string(Dates.today()),
        "engine" => "SciBmad exact RF-on closed orbit/tracking plus direct order-2 maps",
        "lattice" => LATEST_LATTICE,
        "target_sextupole" => target_name,
        "realization_count" => nr,
        "random_seed" => seed,
        "other_sextupole_count" => 75,
        "other_sextupole_offset_rms_m" => other_sext_rms,
        "target_offset_distribution" => "independent uniform[-halfwidth,+halfwidth] in x and y",
        "target_offset_halfwidth_m" => target_halfwidth,
        "quadrupole_count" => length(quadrupoles),
        "quadrupole_relative_error_distribution" => "independent uniform[-fraction,+fraction], fixed within each scan tensor",
        "quadrupole_fraction" => quad_fraction,
        "bump_count" => nb,
        "bump_amplitude_m" => bump_amplitude,
        "k2_count" => nk,
        "k2_levels" => k2_levels,
        "k2_step_m3" => k2_step,
        "bpm_count" => nd,
        "corrector_probe_names" => getproperty.(probe_controls, :name),
        "corrector_probe_field" => corrector_step,
        "local_orbit_measurement_noise" => "none",
        "hidden_from_inverse" => "all non-target sextupole offsets and all quadrupole relative errors",
        "direct_observable_protocol" => "TBT readbacks synthesized from BPM/one-turn maps; fixed-delta energy readback from order-2 map; actual-corrector closed-orbit finite differences",
        "rf_frequency_qualification" => "uses fixed beam-energy delta rather than retuning harmon-master RF frequency",
    ))
    println("Wrote paired scan to $output_dir")
    return 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main())
end
