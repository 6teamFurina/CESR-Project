#!/usr/bin/env julia

"""Orbit-only paired K1 ablation on the same nuisance distribution as the pilot."""

include(joinpath(
    @__DIR__, "..", "..", "quadrupole_affinity", "exact_11_triplet_validation", "common.jl",
))

const BUMP_KNOBS = joinpath(
    @__DIR__, "..", "..", "quadrupole_affinity", "exact_11_triplet_validation",
    "results", "bump_knobs", "local_bump_knobs.csv",
)

function scalar_warm_k1(ring, v0)
    result = find_closed_orbit(
        ring; v0=copy(v0), coasting_beam=false, batch=Val{false}(), warn=false,
    )
    all(result.sol.retcode .== SciBmad.BatchSolve.RETCODE_SUCCESS) ||
        error("Closed-orbit solve failed")
    return result
end

function main(args=ARGS)
    defaults = Dict(
        "target" => "SEX_09AW",
        "realizations" => "8",
        "seed" => "20260817",
        "other-sext-rms-m" => "3e-4",
        "target-halfwidth-m" => "3.5e-4",
        "quadrupole-error-fraction" => "0.01",
        "k1-scan-fraction" => "0.01",
        "bump-amplitude-m" => "5e-4",
        "k2-step-m3" => "0.01",
        "selected-quadrupoles" => "QX4D,Q18W,Q24E",
        "output-dir" => joinpath(@__DIR__, "results", "sex_09aw_k1_orbit_ablation"),
        "overwrite" => "false",
    )
    options = parse_exact11_options(defaults, args)
    target_name = uppercase(options["target"])
    nr = parse(Int, options["realizations"])
    seed = parse(Int, options["seed"])
    other_rms = parse(Float64, options["other-sext-rms-m"])
    target_halfwidth = parse(Float64, options["target-halfwidth-m"])
    quad_error_fraction = parse(Float64, options["quadrupole-error-fraction"])
    k1_fraction = parse(Float64, options["k1-scan-fraction"])
    bump_amplitude = parse(Float64, options["bump-amplitude-m"])
    k2_step = parse(Float64, options["k2-step-m3"])
    selected_names = uppercase.(strip.(split(options["selected-quadrupoles"], ',')))
    output_dir = abspath(options["output-dir"])
    metadata_path = joinpath(output_dir, "scan_metadata.toml")
    isfile(metadata_path) && lowercase(options["overwrite"]) != "true" &&
        error("Output exists; use --overwrite=true")

    ring = cesr
    sextupoles = active_sextupole_inventory(ring)
    quadrupoles = active_quadrupole_inventory(ring)
    detectors = measurable_bpms(ring)
    bpm_names = String.(base_name.(detectors))
    target_index = findfirst(entry -> entry.name == target_name, sextupoles)
    isnothing(target_index) && error("Unknown target $target_name")
    target_entry = sextupoles[target_index]
    target_element = ring.line[target_entry.index]
    quadrupole_by_name = Dict(entry.name => entry for entry in quadrupoles)
    selected = [quadrupole_by_name[name] for name in selected_names]

    controls = independent_corrector_inventory(ring)
    knob_rows = read_bump_knobs(BUMP_KNOBS, target_name)
    knob_by_key = Dict(
        (row["corrector"], row["field"]) => (
            parse(Float64, row["field_per_x_bump_m"]),
            parse(Float64, row["field_per_y_bump_m"]),
        ) for row in knob_rows
    )
    bumps = [
        (-bump_amplitude, 0.0), (0.0, -bump_amplitude), (0.0, 0.0),
        (0.0, bump_amplitude), (bump_amplitude, 0.0),
    ]
    k2_levels = [-1.0, 0.0, 1.0]
    conditions = NamedTuple[(; name="nominal", quadrupole="nominal", sign=0)]
    for entry in selected, sign in (1, -1)
        push!(conditions, (;
            name="$(lowercase(entry.name))_$(sign > 0 ? "plus" : "minus")",
            quadrupole=entry.name, sign,
        ))
    end
    nc, nb, nk, nd = length(conditions), length(bumps), length(k2_levels), length(detectors)
    orbit = zeros(nr, nc, nb, nk, nd, 2)
    target_orbit = zeros(nr, nc, nb, nk, 2)
    truth = zeros(nr, 2)
    quad_errors = zeros(nr, length(quadrupoles))
    rng = MersenneTwister(seed)
    nominal_v0 = copy(solve_closed_orbit(ring).v0)
    start_all = time()

    try
        for realization in 1:nr
            truth[realization, :] .= (
                (2rand(rng)-1)*target_halfwidth, (2rand(rng)-1)*target_halfwidth,
            )
            for (sindex, entry) in enumerate(sextupoles)
                dx, dy = sindex == target_index ?
                    (truth[realization, 1], truth[realization, 2]) :
                    (other_rms*randn(rng), other_rms*randn(rng))
                ring.line[entry.index].x_offset = entry.x_offset_m + dx
                ring.line[entry.index].y_offset = entry.y_offset_m + dy
            end
            physical_k1 = Dict{String,Float64}()
            for (qindex, entry) in enumerate(quadrupoles)
                error_fraction = (2rand(rng)-1)*quad_error_fraction
                quad_errors[realization, qindex] = error_fraction
                physical_k1[entry.name] = entry.kn1_m2*(1+error_fraction)
                set_quadrupole!(ring, entry, physical_k1[entry.name])
            end
            current_v0 = copy(nominal_v0)
            for (condition_index, condition) in enumerate(conditions)
                for entry in selected
                    set_quadrupole!(ring, entry, physical_k1[entry.name])
                end
                if condition.sign != 0
                    entry = quadrupole_by_name[condition.quadrupole]
                    set_quadrupole!(
                        ring, entry,
                        physical_k1[entry.name]*(1 + condition.sign*k1_fraction),
                    )
                end
                for (bump_index, (bump_x, bump_y)) in enumerate(bumps)
                    for control in controls
                        cx, cy = knob_by_key[(control.name, String(control.axis))]
                        set_corrector_scalar!(
                            ring, control,
                            Float64(constant_term(first(control.originals))) + cx*bump_x + cy*bump_y,
                        )
                    end
                    for (k2_index, level) in enumerate(k2_levels)
                        target_element.Kn2 = target_entry.kn2_m3 + level*k2_step
                        closed = scalar_warm_k1(ring, current_v0)
                        current_v0 .= closed.v0
                        tracked = track_orbits_at_names(
                            ring, closed, vcat(bpm_names, [target_name]),
                        )
                        for (i, name) in enumerate(bpm_names)
                            orbit[realization, condition_index, bump_index, k2_index, i, :] .= (
                                tracked.horizontal[name][1], tracked.vertical[name][1],
                            )
                        end
                        target_orbit[realization, condition_index, bump_index, k2_index, :] .= (
                            tracked.horizontal[target_name][1], tracked.vertical[target_name][1],
                        )
                    end
                end
            end
            @printf("K1 orbit ablation realization %d/%d complete\n", realization, nr)
            flush(stdout)
        end
    finally
        for control in controls
            for (index, original) in zip(control.indices, control.originals)
                control.axis == :Kn0 ? (ring.line[index].Kn0 = original) :
                                       (ring.line[index].Ks0 = original)
            end
        end
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
    write_npy(joinpath(output_dir, "target_truth.npy"), truth)
    write_npy(joinpath(output_dir, "quadrupole_relative_errors.npy"), quad_errors)
    write_lines(joinpath(output_dir, "condition_names.txt"), getproperty.(conditions, :name))
    write_rows(joinpath(output_dir, "bump_points.csv"), [
        (; bump_index=i, bump_x_command_m=x, bump_y_command_m=y)
        for (i, (x, y)) in enumerate(bumps)
    ])
    write_metadata(metadata_path, Dict(
        "format" => "cesr-k1-orbit-nuisance-ablation-v1",
        "lattice" => LATEST_LATTICE,
        "target_sextupole" => target_name,
        "realization_count" => nr,
        "random_seed" => seed,
        "other_sextupole_count" => 75,
        "other_sextupole_offset_rms_m" => other_rms,
        "quadrupole_error_fraction" => quad_error_fraction,
        "selected_quadrupoles" => selected_names,
        "k1_scan_fraction" => k1_fraction,
        "k1_condition_count" => nc,
        "bump_count" => nb,
        "k2_count" => nk,
        "k2_levels" => k2_levels,
        "k2_step_m3" => k2_step,
        "wall_seconds" => time()-start_all,
    ))
    println("Wrote K1 orbit ablation to $output_dir")
    return 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main())
end
