#!/usr/bin/env julia

"""Generate exact local-orbit validation rays for bump/K2 Taylor maps.

The scan is deterministic and uses the validated latest SciBmad lattice.  It
varies the two model-based local-bump knobs and the target sextupole Kn2 around
their nominal settings.  The entrance and exit six-vectors of the target,
full-ring orbit extrema, and corrector demand are saved for every state.
"""

include(joinpath(
    @__DIR__, "..", "quadrupole_affinity", "exact_11_triplet_validation", "common.jl",
))

const DEFAULT_BUMP_KNOBS = joinpath(
    @__DIR__, "..", "quadrupole_affinity", "exact_11_triplet_validation",
    "results", "bump_knobs", "local_bump_knobs.csv",
)

function selected_entries(requested, sextupoles)
    text = lowercase(strip(requested))
    text == "all" && return collect(sextupoles)
    wanted = Set(uppercase(strip(name)) for name in split(requested, ',') if !isempty(strip(name)))
    known = Set(getproperty.(sextupoles, :name))
    missing = setdiff(wanted, known)
    isempty(missing) || error("Unknown targets: $(join(sort!(collect(missing)), ", "))")
    return [entry for entry in sextupoles if entry.name in wanted]
end

function scan_directions()
    rows = NamedTuple[]
    function add(
        label, family, qx, qy, qk;
        base_qx=0, base_qy=0, base_qk=0, maximum_radius=Inf,
    )
        push!(rows, (;
            label,
            family,
            base_qx=Float64(base_qx),
            base_qy=Float64(base_qy),
            base_qk=Float64(base_qk),
            qx=Float64(qx),
            qy=Float64(qy),
            qk=Float64(qk),
            maximum_radius=Float64(maximum_radius),
        ))
    end
    for sign in (-1, 1)
        suffix = sign < 0 ? "m" : "p"
        add("x_$suffix", "x_axis", sign, 0, 0)
        add("y_$suffix", "y_axis", 0, sign, 0)
        add("k_$suffix", "k_axis", 0, 0, sign)
    end
    for sx in (-1, 1), sy in (-1, 1)
        add("xy_$(sx < 0 ? "m" : "p")$(sy < 0 ? "m" : "p")", "xy_diagonal", sx, sy, 0)
    end
    # The 1.5:1 ratio passes through the maintained +/-1.5 mm, +/-0.10 m^-3
    # stochastic-inverse operating states when radius_scale == 1.
    for sb in (-1, 1), sk in (-1, 1)
        bsign = sb < 0 ? "m" : "p"
        ksign = sk < 0 ? "m" : "p"
        add("xk_$(bsign)$(ksign)", "xk_protocol", 1.5sb, 0, sk)
        add("yk_$(bsign)$(ksign)", "yk_protocol", 0, 1.5sb, sk)
    end
    # Conditional one-variable sweeps separate a bump limit at the maintained
    # |delta K2|=0.10 m^-3 from a K2 limit at the maintained 1.5 mm bump.
    for sb in (-1, 1), sk in (-1, 1)
        bsign = sb < 0 ? "m" : "p"
        ksign = sk < 0 ? "m" : "p"
        add(
            "x_at_k_$(bsign)$(ksign)", "x_bump_fixed_k", sb, 0, 0;
            base_qk=sk,
        )
        add(
            "y_at_k_$(bsign)$(ksign)", "y_bump_fixed_k", 0, sb, 0;
            base_qk=sk,
        )
        add(
            "k_at_x_$(bsign)$(ksign)", "k_fixed_x_bump", 0, 0, sk;
            base_qx=1.5sb,
        )
        add(
            "k_at_y_$(bsign)$(ksign)", "k_fixed_y_bump", 0, 0, sk;
            base_qy=1.5sb,
        )
    end
    # Inner three-parameter corners make the total-order-four empirical
    # fallback design full rank without paying for unnecessary outer x-y-K2
    # diagnostic rays.
    for sx in (-1, 1), sy in (-1, 1), sk in (-1, 1)
        add(
            "xyk_inner_$(sx < 0 ? "m" : "p")$(sy < 0 ? "m" : "p")$(sk < 0 ? "m" : "p")",
            "xyk_inner_fit",
            sx,
            sy,
            sk;
            maximum_radius=1.5,
        )
    end
    for sx in (-1, 1), sy in (-1, 1)
        add(
            "xy_asym_inner_$(sx < 0 ? "m" : "p")$(sy < 0 ? "m" : "p")",
            "xy_asymmetric_inner_fit",
            sx,
            2sy,
            0;
            maximum_radius=1.5,
        )
    end
    return rows
end

function batch_closed_orbit(ring, state_count, nominal_v0)
    solution = find_closed_orbit(
        ring;
        v0=repeat(reshape(vec(nominal_v0), 1, 6), state_count, 1),
        coasting_beam=false,
        batch=Val{true}(),
        warn=false,
    )
    converged = vec(Array(solution.sol.retcode .== SciBmad.BatchSolve.RETCODE_SUCCESS))
    return solution, converged
end

function local_and_ring_orbits(ring, closed, target_index, bpm_names)
    state_count = size(closed.v0, 1)
    bunch = Bunch(v=copy(closed.v0))
    SciBmad.BTBL.check_bl_bunch!(bunch, ring, false)
    entry = fill(NaN, state_count, 6)
    exit = fill(NaN, state_count, 6)
    maximum_abs_x = zeros(state_count)
    maximum_abs_y = zeros(state_count)
    maximum_bpm_abs_x = zeros(state_count)
    maximum_bpm_abs_y = zeros(state_count)
    for (index, element) in enumerate(ring.line)
        coordinates = Array(bunch.coords.v)
        maximum_abs_x .= max.(maximum_abs_x, abs.(coordinates[:, 1]))
        maximum_abs_y .= max.(maximum_abs_y, abs.(coordinates[:, 3]))
        if base_name(element) in bpm_names
            maximum_bpm_abs_x .= max.(maximum_bpm_abs_x, abs.(coordinates[:, 1]))
            maximum_bpm_abs_y .= max.(maximum_bpm_abs_y, abs.(coordinates[:, 3]))
        end
        index == target_index && (entry .= coordinates)
        track!(bunch, element)
        index == target_index && (exit .= Array(bunch.coords.v))
    end
    return (; entry, exit, maximum_abs_x, maximum_abs_y, maximum_bpm_abs_x, maximum_bpm_abs_y)
end

function main(args=ARGS)
    defaults = Dict(
        "targets" => "all",
        "bump-scale-m" => "1.0e-3",
        "k2-scale-m3" => "0.1",
        "radius-levels" => "0.25,0.5,0.75,1,1.25,1.5,2,2.5,3,4,5,6,8,10,12,15,20",
        "bump-knobs-csv" => DEFAULT_BUMP_KNOBS,
        "output-dir" => joinpath(@__DIR__, "results", "exact_validation"),
        "overwrite" => "false",
    )
    options = parse_exact11_options(defaults, args)
    bump_scale = parse(Float64, options["bump-scale-m"])
    k2_scale = parse(Float64, options["k2-scale-m3"])
    radii = parse.(Float64, split(options["radius-levels"], ','))
    all(isfinite, radii) && all(>(0), radii) || error("Radius levels must be positive and finite")
    issorted(radii) || error("Radius levels must be increasing")
    output_dir = abspath(options["output-dir"])
    metadata_path = joinpath(output_dir, "scan_metadata.toml")
    isfile(metadata_path) && lowercase(options["overwrite"]) != "true" &&
        error("Output exists; use --overwrite=true: $metadata_path")

    ring = cesr
    sextupoles = active_sextupole_inventory(ring)
    targets = selected_entries(options["targets"], sextupoles)
    controls = independent_corrector_inventory(ring)
    bpm_names = Set(String.(base_name.(measurable_bpms(ring))))
    directions = scan_directions()
    knob_rows = read_simple_csv(abspath(options["bump-knobs-csv"]))
    knobs_by_target = Dict{String,Dict{Tuple{String,String},Tuple{Float64,Float64}}}()
    for row in knob_rows
        target_knobs = get!(
            knobs_by_target,
            uppercase(row["target_sextupole"]),
            Dict{Tuple{String,String},Tuple{Float64,Float64}}(),
        )
        target_knobs[(row["corrector"], row["field"])] = (
            parse(Float64, row["field_per_x_bump_m"]),
            parse(Float64, row["field_per_y_bump_m"]),
        )
    end
    control_keys = Set((control.name, String(control.axis)) for control in controls)
    all(Set(keys(knobs_by_target[target.name])) == control_keys for target in targets) ||
        error("At least one target bump knob has a control mismatch")

    # The exported lattice retains some scalar strengths as DefExpr values.
    # Materialize the multipole parameter blocks that will be promoted to
    # BatchParam; mixing DefExpr{Any} and BatchParam is otherwise ambiguous.
    promoted_indices = Set(index for control in controls for index in control.indices)
    union!(promoted_indices, getproperty.(targets, :index))
    for element_index in promoted_indices
        element = ring.line[element_index]
        element.BMultipoleParams = Beamlines.deval(element.BMultipoleParams)
    end

    nominal_closed_timed = @timed solve_closed_orbit(ring)
    nominal_v0 = copy(nominal_closed_timed.value.v0)
    coordinate_names = ("x", "px", "y", "py", "z", "pz")
    result_rows = NamedTuple[]
    timing_rows = NamedTuple[]
    mkpath(joinpath(output_dir, "targets"))
    calculation_start = time()
    for (target_index, target) in enumerate(targets)
        target_start = time()
        first_target_row = length(result_rows) + 1
        target_knobs = knobs_by_target[target.name]
        states = [(; label="origin", family="origin", qx=0.0, qy=0.0, qk=0.0, radius=0.0)]
        append!(states, [
            (;
                label=direction.label,
                family=direction.family,
                qx=direction.base_qx + radius*direction.qx,
                qy=direction.base_qy + radius*direction.qy,
                qk=direction.base_qk + radius*direction.qk,
                radius,
            )
            for direction in directions for radius in radii
            if radius <= direction.maximum_radius
        ])
        state_count = length(states)
        maximum_abs_corrector_delta = zeros(state_count)
        maximum_abs_corrector_total = zeros(state_count)
        for control in controls
            cx, cy = target_knobs[(control.name, String(control.axis))]
            deltas = [bump_scale * (cx*state.qx + cy*state.qy) for state in states]
            maximum_abs_corrector_delta .= max.(maximum_abs_corrector_delta, abs.(deltas))
            totals = Float64[]
            for original in control.originals
                append!(totals, constant_term(original) .+ deltas)
            end
            maximum_abs_corrector_total .= max.(
                maximum_abs_corrector_total,
                maximum(abs.(reshape(totals, state_count, :)); dims=2)[:, 1],
            )
            set_corrector_values!(
                ring,
                control,
                constant_term(first(control.originals)) .+ deltas,
            )
        end
        delta_k2 = [k2_scale * state.qk for state in states]
        ring.line[target.index].Kn2 = BatchParam(target.kn2_m3 .+ delta_k2)
        solve_timed = @timed batch_closed_orbit(ring, state_count, nominal_v0)
        closed, converged = solve_timed.value
        orbit = local_and_ring_orbits(ring, closed, target.index, bpm_names)
        for (state_index, state) in enumerate(states)
            failure = converged[state_index] ? "" : string(closed.sol.retcode[state_index])
            entry = converged[state_index] ? orbit.entry[state_index, :] : fill(NaN, 6)
            exit = converged[state_index] ? orbit.exit[state_index, :] : fill(NaN, 6)
            entry_keys = Tuple(Symbol.("entry_" .* collect(coordinate_names)))
            exit_keys = Tuple(Symbol.("exit_" .* collect(coordinate_names)))
            entry_values = NamedTuple{entry_keys}(Tuple(entry))
            exit_values = NamedTuple{exit_keys}(Tuple(exit))
            push!(result_rows, merge((;
                target=target.name,
                target_s_m=target.s_m,
                target_length_m=target.length_m,
                target_nominal_kn2_m3=target.kn2_m3,
                state_index,
                direction=state.label,
                direction_family=state.family,
                radius_scale=state.radius,
                qx=state.qx,
                qy=state.qy,
                qk=state.qk,
                bump_x_command_m=bump_scale*state.qx,
                bump_y_command_m=bump_scale*state.qy,
                delta_k2_m3=delta_k2[state_index],
                converged=converged[state_index],
                failure,
                closed_orbit_seconds=solve_timed.time,
                maximum_abs_corrector_delta=maximum_abs_corrector_delta[state_index],
                maximum_abs_corrector_total=maximum_abs_corrector_total[state_index],
                maximum_ring_abs_x_m=converged[state_index] ? orbit.maximum_abs_x[state_index] : NaN,
                maximum_ring_abs_y_m=converged[state_index] ? orbit.maximum_abs_y[state_index] : NaN,
                maximum_bpm_abs_x_m=converged[state_index] ? orbit.maximum_bpm_abs_x[state_index] : NaN,
                maximum_bpm_abs_y_m=converged[state_index] ? orbit.maximum_bpm_abs_y[state_index] : NaN,
            ), entry_values, exit_values))
        end
        restore_sextupole!(ring, target)
        for control in controls
            for (element_index, original) in zip(control.indices, control.originals)
                if control.axis == :Kn0
                    ring.line[element_index].Kn0 = original
                else
                    ring.line[element_index].Ks0 = original
                end
            end
        end
        target_seconds = time() - target_start
        push!(timing_rows, (;
            target_index,
            target=target.name,
            state_count=length(states),
            converged_count=count(converged),
            seconds=target_seconds,
        ))
        write_rows(
            joinpath(output_dir, "targets", "$(lowercase(target.name)).csv"),
            result_rows[first_target_row:end],
        )
        @printf(
            "%s exact rays %d/%d: %d states in %.1f s (elapsed %.1f s)\n",
            target.name,
            target_index,
            length(targets),
            length(states),
            target_seconds,
            time() - calculation_start,
        )
        flush(stdout)
    end

    write_rows(joinpath(output_dir, "exact_local_orbit_states.csv"), result_rows)
    write_rows(joinpath(output_dir, "target_timings.csv"), timing_rows)
    write_rows(joinpath(output_dir, "scan_directions.csv"), directions)
    write_metadata(metadata_path, Dict(
        "format" => "cesr-sextupole-local-orbit-bump-k2-exact-rays-v1",
        "date" => string(Dates.today()),
        "engine" => "SciBmad exact scalar RF-on closed orbit and element-by-element tracking",
        "lattice" => LATEST_LATTICE,
        "target_count" => length(targets),
        "target_inventory" => getproperty.(targets, :name),
        "direction_count" => length(directions),
        "radius_levels" => radii,
        "state_count_per_target" => 1 + sum(
            count(radius -> radius <= direction.maximum_radius, radii)
            for direction in directions
        ),
        "bump_scale_m" => bump_scale,
        "k2_scale_m3" => k2_scale,
        "nominal_closed_orbit_seconds" => nominal_closed_timed.time,
        "calculation_wall_seconds" => time() - calculation_start,
        "machine_error_model" => "validated latest lattice nominal settings; no randomized latent errors",
        "interpretation_boundary" => "model-validity scan only; corrector fields are reported but no hardware or aperture limit is imposed",
    ))
    println("Wrote exact validation rays to $output_dir")
    return 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main())
end
