#!/usr/bin/env julia

"""
Compute the CESR control-to-closed-orbit response matrix with GTPSA and
compare it with the labeled Tao/Bmad reference matrix.

The 119 CESR HKICK/VKICK Overlay knobs are represented by first-order GTPSA
parameters.  If the one-turn map about the baseline closed orbit is

    F(z, k) = z0 + A * (z - z0) + B * k,

then the parameterized closed-orbit derivative is found without finite
differencing from

    X = dz_closed/dk = (I - A) \\ B.

For an RF-off coasting beam the fixed-z, fixed-pz four-dimensional form is
used.  The parameterized closed orbit is then tracked element by element and
the parameter coefficients of x and y are recorded at every DET_* marker.

Run from the CESR Project directory:

    julia --project=. bmad_comparison/test_control_response_scibmad.jl --mode=rf_on
    julia --project=. bmad_comparison/test_control_response_scibmad.jl --mode=rf_off
    julia --project=. bmad_comparison/test_control_response_scibmad.jl --mode=both
    julia --project=. bmad_comparison/test_control_response_scibmad.jl --mode=summary
"""

using LinearAlgebra
using Printf
using SciBmad
using GTPSA
using Beamlines

const HERE = @__DIR__
const PROJECT_ROOT = normpath(joinpath(HERE, ".."))

include(joinpath(PROJECT_ROOT, "cesr_model.jl"))

struct LabeledMatrix
    row_labels::Vector{String}
    column_labels::Vector{String}
    values::Matrix{Float64}
end

csv_quote(value) = '"' * replace(string(value), '"' => "\"\"") * '"'
function csv_unquote(value::AbstractString)
    text = String(value)
    if length(text) >= 2 && first(text) == '"' && last(text) == '"'
        return replace(text[2:(end - 1)], "\"\"" => "\"")
    end
    return text
end

function read_labeled_matrix(path::AbstractString)
    lines = readlines(path)
    isempty(lines) && error("Empty matrix CSV: $path")
    header = split(lines[1], ',')
    length(header) > 1 || error("Matrix CSV has no control columns: $path")
    column_labels = csv_unquote.(header[2:end])
    row_labels = String[]
    values = Matrix{Float64}(undef, length(lines) - 1, length(column_labels))
    for (row_index, line) in enumerate(lines[2:end])
        fields = split(line, ',')
        length(fields) == length(header) ||
            error("CSV row $row_index has $(length(fields)) fields; expected $(length(header))")
        push!(row_labels, csv_unquote(fields[1]))
        for column_index in eachindex(column_labels)
            values[row_index, column_index] = parse(
                Float64,
                csv_unquote(fields[column_index + 1]),
            )
        end
    end
    return LabeledMatrix(row_labels, column_labels, values)
end

function write_labeled_matrix(path, matrix::LabeledMatrix)
    open(path, "w") do io
        println(io, join(csv_quote.(vcat("observable", matrix.column_labels)), ','))
        for row_index in eachindex(matrix.row_labels)
            fields = Any[matrix.row_labels[row_index]]
            append!(fields, matrix.values[row_index, :])
            println(io, join(csv_quote.(fields), ','))
        end
    end
    return path
end

function write_entry_comparison(path, reference::LabeledMatrix, candidate::LabeledMatrix)
    reference.row_labels == candidate.row_labels || error("Response row labels differ")
    reference.column_labels == candidate.column_labels || error("Response column labels differ")
    open(path, "w") do io
        println(io, "observable,control,bmad,scibmad,difference,abs_difference")
        for column in axes(reference.values, 2), row in axes(reference.values, 1)
            bmad = reference.values[row, column]
            scibmad = candidate.values[row, column]
            difference = scibmad - bmad
            println(io, join(csv_quote.((
                reference.row_labels[row],
                reference.column_labels[column],
                bmad,
                scibmad,
                difference,
                abs(difference),
            )), ','))
        end
    end
    return path
end

function control_attribute_map()
    result = Dict{String,String}()
    for (name, attribute, _, _) in CESR_OVERLAY_TERMS
        if haskey(result, name) && result[name] != attribute
            error("Control $name has inconsistent attributes")
        end
        result[name] = attribute
    end
    return result
end

function detector_names(reference::LabeledMatrix)
    iseven(length(reference.row_labels)) || error("Response row count must be even")
    n = div(length(reference.row_labels), 2)
    x_labels = reference.row_labels[1:n]
    y_labels = reference.row_labels[(n + 1):end]
    all(endswith(":x"), x_labels) || error("First response block is not orbit.x")
    all(endswith(":y"), y_labels) || error("Second response block is not orbit.y")
    x_names = first.(split.(x_labels, ":"))
    y_names = first.(split.(y_labels, ":"))
    x_names == y_names || error("The orbit.x and orbit.y detector orders differ")
    return uppercase.(x_names)
end

function validate_reference(reference::LabeledMatrix)
    size(reference.values) == (198, 119) ||
        error("Expected a 198 x 119 Bmad matrix, got $(size(reference.values))")
    detectors = detector_names(reference)
    length(detectors) == 99 || error("Expected 99 detector positions")
    length(unique(detectors)) == length(detectors) || error("Detector labels are not unique")

    attributes = control_attribute_map()
    Set(reference.column_labels) == Set(keys(attributes)) ||
        error("Bmad control columns do not match the CESR Overlay controls")
    expected_attributes = vcat(fill("HKICK", 58), fill("VKICK", 61))
    actual_attributes = [attributes[name] for name in reference.column_labels]
    actual_attributes == expected_attributes ||
        error("Bmad control columns are not ordered as 58 HKICK then 61 VKICK")
    return detectors
end

function parameterized_model(rf_on::Bool, control_names::Vector{String})
    descriptor = Descriptor(6, 1, length(control_names), 1)
    variables = vars(descriptor)
    parameters = params(descriptor)
    model = load_cesr_model(zero_value=zero(parameters[1]), rf_on=rf_on)
    for (index, name) in enumerate(control_names)
        model.controls[name] = parameters[index]
    end
    return (; descriptor, variables, parameters, model)
end

function track_map(ring, coordinates)
    bunch = Bunch(v=reshape(coordinates, 1, 6))
    SciBmad.BTBL.check_bl_bunch!(bunch, ring, false)
    track!(bunch, ring)
    return vec(bunch.coords.v)
end

function coasting_closed_orbit(ring; tolerance=1e-13, max_iterations=20)
    # SciBmad's current coasting-beam solver differentiates through the CESR
    # implicit integrator with ForwardDiff, which faults on this lattice. Use a
    # Float64 finite-difference Jacobian only for the baseline 4D Newton solve;
    # the 119 requested control derivatives remain analytic GTPSA parameters.
    state = zeros(4)
    residual(state4) = track_map(ring, [state4; 0.0; 0.0])[1:4] - state4
    current = residual(state)

    for iteration in 0:max_iterations
        residual_norm = maximum(abs, current)
        residual_norm <= tolerance && return (
            orbit=[state; 0.0; 0.0],
            solution=(coasting_beam=true, iterations=iteration, residual=residual_norm),
        )
        iteration == max_iterations && break

        jacobian = zeros(4, 4)
        for coordinate in 1:4
            step = 1e-7 * max(1.0, abs(state[coordinate]))
            plus = copy(state)
            minus = copy(state)
            plus[coordinate] += step
            minus[coordinate] -= step
            jacobian[:, coordinate] .= (residual(plus) - residual(minus)) / (2step)
        end
        newton_step = -(jacobian \ current)

        accepted = false
        scale = 1.0
        for _ in 1:12
            trial_state = state + scale * newton_step
            trial_residual = residual(trial_state)
            if maximum(abs, trial_residual) < residual_norm
                state = trial_state
                current = trial_residual
                accepted = true
                break
            end
            scale *= 0.5
        end
        accepted || error(
            "Float64 coasting closed-orbit Newton line search failed at iteration $iteration",
        )
    end
    error(
        "Float64 coasting closed-orbit solve did not converge; " *
        "residual=$(maximum(abs, current))",
    )
end

function closed_orbit_vector(ring, rf_on::Bool)
    if !rf_on
        result = coasting_closed_orbit(ring)
        return result.orbit, result.solution
    end

    solution = find_closed_orbit(
        ring;
        coasting_beam=false,
        batch=Val{false}(),
        warn=true,
    )
    bad = findall(solution.sol.retcode .!= SciBmad.BatchSolve.RETCODE_SUCCESS)
    isempty(bad) || error("Closed-orbit solve did not converge: retcode=$(solution.sol.retcode)")
    return vec(Float64.(solution.v0[1, :])), solution
end

function one_turn_response(ring, closed_orbit, variables, n_controls, rf_on)
    input_map = [closed_orbit[index] + copy(variables[index]) for index in 1:6]
    output_map = track_map(ring, input_map)
    full_jacobian = Matrix(GTPSA.jacobian(output_map; include_params=true))
    size(full_jacobian) == (6, 6 + n_controls) ||
        error("Unexpected one-turn GTPSA Jacobian size: $(size(full_jacobian))")
    A = full_jacobian[:, 1:6]
    B = full_jacobian[:, 7:end]
    response = zeros(6, n_controls)
    if rf_on
        response .= (I - A) \ B
        closure_residual = (I - A) * response - B
    else
        response[1:4, :] .= (I - A[1:4, 1:4]) \ B[1:4, :]
        closure_residual = (I - A[1:4, 1:4]) * response[1:4, :] - B[1:4, :]
    end
    return (; A, B, response, closure_residual, output_map)
end

function detector_response(
    ring,
    closed_orbit,
    closed_orbit_response,
    parameters,
    detectors::Vector{String},
)
    n_controls = length(parameters)
    initial_map = [
        closed_orbit[coordinate] + sum(
            closed_orbit_response[coordinate, control] * parameters[control]
            for control in 1:n_controls
        )
        for coordinate in 1:6
    ]
    bunch = Bunch(v=reshape(initial_map, 1, 6))
    SciBmad.BTBL.check_bl_bunch!(bunch, ring, false)

    detector_index = Dict(name => index for (index, name) in enumerate(detectors))
    found = falses(length(detectors))
    horizontal = zeros(length(detectors), n_controls)
    vertical = zeros(length(detectors), n_controls)

    for element in ring.line
        track!(bunch, element)
        name = uppercase(String(element.name))
        haskey(detector_index, name) || continue
        index = detector_index[name]
        found[index] && error("Detector $name occurs more than once in the tracking line")
        jac = Matrix(GTPSA.jacobian(vec(bunch.coords.v); include_params=true))
        horizontal[index, :] .= jac[1, 7:end]
        vertical[index, :] .= jac[3, 7:end]
        found[index] = true
    end

    all(found) || error("Missing SciBmad detectors: $(detectors[.!found])")
    return vcat(horizontal, vertical), vec(bunch.coords.v)
end

function normalized_difference(difference, reference)
    scale = norm(reference)
    return iszero(scale) ? (iszero(norm(difference)) ? 0.0 : Inf) : norm(difference) / scale
end

function matrix_correlation(a, b)
    av = vec(a)
    bv = vec(b)
    denom = norm(av) * norm(bv)
    return iszero(denom) ? NaN : dot(av, bv) / denom
end

function block_metrics(bmad, scibmad, row_range, column_range)
    reference = bmad[row_range, column_range]
    candidate = scibmad[row_range, column_range]
    difference = candidate - reference
    reference_max = maximum(abs, reference)
    difference_max = maximum(abs, difference)
    return (
        bmad_fro=norm(reference),
        scibmad_fro=norm(candidate),
        difference_fro=norm(difference),
        relative_fro=normalized_difference(difference, reference),
        max_abs_difference=difference_max,
        bmad_max_abs=reference_max,
        max_normalized=iszero(reference_max) ? Inf : difference_max / reference_max,
        correlation=matrix_correlation(reference, candidate),
    )
end

function comparison_metrics(reference::LabeledMatrix, candidate::LabeledMatrix)
    bmad = reference.values
    scibmad = candidate.values
    difference = scibmad - bmad
    worst = argmax(abs.(difference))
    column_relative = [
        normalized_difference(difference[:, column], bmad[:, column])
        for column in axes(bmad, 2)
    ]
    worst_column = argmax(column_relative)
    singular_bmad = svdvals(bmad)
    singular_scibmad = svdvals(scibmad)
    singular_difference = singular_scibmad - singular_bmad
    blocks = (
        xH=block_metrics(bmad, scibmad, 1:99, 1:58),
        xV=block_metrics(bmad, scibmad, 1:99, 59:119),
        yH=block_metrics(bmad, scibmad, 100:198, 1:58),
        yV=block_metrics(bmad, scibmad, 100:198, 59:119),
    )
    return (
        difference=difference,
        relative_fro=normalized_difference(difference, bmad),
        max_abs_difference=maximum(abs, difference),
        bmad_max_abs=maximum(abs, bmad),
        max_normalized=maximum(abs, difference) / maximum(abs, bmad),
        correlation=matrix_correlation(bmad, scibmad),
        worst_row=worst[1],
        worst_column=worst[2],
        column_relative=column_relative,
        worst_column_index=worst_column,
        singular_bmad=singular_bmad,
        singular_scibmad=singular_scibmad,
        singular_relative=normalized_difference(singular_difference, singular_bmad),
        blocks=blocks,
    )
end

function write_column_metrics(path, reference, metrics)
    open(path, "w") do io
        println(io, "column,control,relative_2norm")
        for column in eachindex(reference.column_labels)
            println(io, join(csv_quote.((
                column,
                reference.column_labels[column],
                metrics.column_relative[column],
            )), ','))
        end
    end
end

function format_vector(values; count=6)
    n = min(count, length(values))
    return join((@sprintf("%.9e", values[index]) for index in 1:n), ", ")
end

function write_mode_summary(
    path,
    mode,
    reference,
    candidate,
    metrics,
    closed_orbit,
    coasting_beam,
    closure_residual,
    runtime_seconds,
)
    worst_row_label = reference.row_labels[metrics.worst_row]
    worst_control = reference.column_labels[metrics.worst_column]
    worst_column_control = reference.column_labels[metrics.worst_column_index]
    open(path, "w") do io
        println(io, "# Bmad–SciBmad CESR control-response comparison ($(uppercase(replace(mode, '_' => '-'))))")
        println(io)
        println(io, "The response matrix maps 119 Bmad-compatible CESR control Overlay knobs to horizontal and vertical closed orbit at 99 `DET_*` markers. SciBmad computes all control derivatives simultaneously with a first-order GTPSA parameter map; no control finite differences are used.")
        println(io)
        println(io, "- Matrix shape: `$(size(candidate.values, 1)) x $(size(candidate.values, 2))`")
        println(io, "- Units: `m/rad`")
        println(io, "- Closed-orbit model: `$(coasting_beam ? "4D coasting beam with fixed z and pz" : "6D RF-confined beam")`")
        println(io, "- SciBmad baseline closed orbit: `[$(format_vector(closed_orbit))]`")
        println(io, @sprintf("- GTPSA closure-equation residual: `%.6e`", maximum(abs, closure_residual)))
        println(io, @sprintf("- Runtime: `%.3f s`", runtime_seconds))
        println(io)
        println(io, "## Overall agreement")
        println(io)
        println(io, @sprintf("- Relative Frobenius difference: `%.9e` (`%.6f%%`)", metrics.relative_fro, 100metrics.relative_fro))
        println(io, @sprintf("- Maximum absolute entry difference: `%.9e m/rad`", metrics.max_abs_difference))
        println(io, @sprintf("- Maximum difference normalized by the Bmad matrix maximum: `%.9e` (`%.6f%%`)", metrics.max_normalized, 100metrics.max_normalized))
        println(io, @sprintf("- Full-matrix cosine correlation: `%.12f`", metrics.correlation))
        println(io, @sprintf("- Singular-value relative 2-norm difference: `%.9e`", metrics.singular_relative))
        println(io, "- Worst entry: `$(worst_row_label)` versus `$(worst_control)`")
        println(io, @sprintf("- Worst column relative 2-norm: `%.9e` for `%s`", metrics.column_relative[metrics.worst_column_index], worst_column_control))
        println(io)
        println(io, "## Plane blocks")
        println(io)
        println(io, "| Block | Relative Frobenius | Max absolute difference (m/rad) | Max-normalized difference | Correlation |")
        println(io, "|---|---:|---:|---:|---:|")
        for name in (:xH, :xV, :yH, :yV)
            block = getproperty(metrics.blocks, name)
            println(io, @sprintf(
                "| %s | %.9e | %.9e | %.9e | %.12f |",
                String(name), block.relative_fro, block.max_abs_difference,
                block.max_normalized, block.correlation,
            ))
        end
        println(io)
        println(io, "## Leading singular values")
        println(io)
        println(io, "| Index | Bmad | SciBmad | Difference |")
        println(io, "|---:|---:|---:|---:|")
        for index in 1:min(12, length(metrics.singular_bmad))
            println(io, @sprintf(
                "| %d | %.9e | %.9e | %.9e |",
                index,
                metrics.singular_bmad[index],
                metrics.singular_scibmad[index],
                metrics.singular_scibmad[index] - metrics.singular_bmad[index],
            ))
        end
    end
    return path
end

function mode_paths(mode::String)
    directory = joinpath(HERE, "bmad_control_response_$mode")
    return (
        directory=directory,
        bmad=joinpath(directory, "bmad_control_response_$mode.csv"),
        scibmad=joinpath(directory, "scibmad_control_response_$mode.csv"),
        entries=joinpath(directory, "bmad_scibmad_control_response_entries_$mode.csv"),
        columns=joinpath(directory, "bmad_scibmad_control_response_columns_$mode.csv"),
        summary=joinpath(directory, "bmad_scibmad_control_response_summary_$mode.md"),
    )
end

function compare_mode(mode::String)
    mode in ("rf_on", "rf_off") || error("Unknown mode: $mode")
    rf_on = mode == "rf_on"
    paths = mode_paths(mode)
    isfile(paths.bmad) || error("Bmad response CSV is missing: $(paths.bmad)")
    reference = read_labeled_matrix(paths.bmad)
    detectors = validate_reference(reference)

    @printf("\n[%s] Loading Float64 CESR model and finding baseline closed orbit...\n", mode)
    float_model = load_cesr_model(rf_on=rf_on)
    closed_orbit, closed_orbit_solution = closed_orbit_vector(float_model.ring, rf_on)
    @printf("[%s] Baseline closed orbit: [%s]\n", mode, format_vector(closed_orbit))

    elapsed = @elapsed begin
        @printf("[%s] Attaching %d GTPSA control parameters...\n", mode, length(reference.column_labels))
        parameterized = parameterized_model(rf_on, reference.column_labels)
        @printf("[%s] Tracking the parameterized one-turn map...\n", mode)
        one_turn = one_turn_response(
            parameterized.model.ring,
            closed_orbit,
            parameterized.variables,
            length(reference.column_labels),
            rf_on,
        )
        @printf("[%s] Closure residual: %.3e\n", mode, maximum(abs, one_turn.closure_residual))
        @printf("[%s] Tracking the parameterized closed orbit to %d detectors...\n", mode, length(detectors))
        response_values, final_map = detector_response(
            parameterized.model.ring,
            closed_orbit,
            one_turn.response,
            parameterized.parameters,
            detectors,
        )
    end

    candidate = LabeledMatrix(reference.row_labels, reference.column_labels, response_values)
    metrics = comparison_metrics(reference, candidate)
    write_labeled_matrix(paths.scibmad, candidate)
    write_entry_comparison(paths.entries, reference, candidate)
    write_column_metrics(paths.columns, reference, metrics)
    write_mode_summary(
        paths.summary,
        mode,
        reference,
        candidate,
        metrics,
        closed_orbit,
        closed_orbit_solution.coasting_beam,
        one_turn.closure_residual,
        elapsed,
    )

    @printf("[%s] Relative Frobenius difference: %.6e (%.6f%%)\n", mode, metrics.relative_fro, 100metrics.relative_fro)
    @printf("[%s] Maximum absolute difference:  %.6e m/rad\n", mode, metrics.max_abs_difference)
    @printf("[%s] Full-matrix correlation:      %.12f\n", mode, metrics.correlation)
    println("[$mode] SciBmad CSV: $(paths.scibmad)")
    println("[$mode] Summary:     $(paths.summary)")
    return (; mode, status="complete", paths, metrics, elapsed)
end

function write_overall_summary(path, results, failures)
    open(path, "w") do io
        println(io, "# CESR Bmad–SciBmad control-response summary")
        println(io)
        println(io, "SciBmad uses 119 first-order GTPSA parameters to differentiate the closed orbit at all 99 detector markers with respect to the Bmad-compatible HKICK/VKICK Overlay knobs.")
        println(io)
        println(io, "| Mode | Status | Relative Frobenius | Max absolute difference (m/rad) | Correlation | Runtime (s) |")
        println(io, "|---|---|---:|---:|---:|---:|")
        for result in results
            println(io, @sprintf(
                "| %s | complete | %.9e | %.9e | %.12f | %.3f |",
                result.mode,
                result.metrics.relative_fro,
                result.metrics.max_abs_difference,
                result.metrics.correlation,
                result.elapsed,
            ))
        end
        for failure in failures
            println(io, "| $(failure.mode) | failed | — | — | — | — |")
        end
        if !isempty(results)
            println(io)
            println(io, "## Interpretation")
            println(io)
            worst_relative = maximum(result.metrics.relative_fro for result in results)
            worst_correlation = minimum(result.metrics.correlation for result in results)
            println(io, @sprintf(
                "All completed matrices agree with Bmad to at most `%.6f%%` in relative Frobenius norm; the lowest full-matrix correlation is `%.12f`.",
                100worst_relative,
                worst_correlation,
            ))
            println(io)
            println(io, "- Bmad `HKICK`/`VKICK` are laboratory-frame kicks. For element alignment tilt `t`, SciBmad uses `HKICK -> (Kn0L, Ks0L) = (-cos(t), -sin(t)) HKICK` and `VKICK -> (Kn0L, Ks0L) = (-sin(t), cos(t)) VKICK`.")
            println(io, "- RF-on uses the 6D RF-confined closed orbit. RF-off uses a 4D coasting closed orbit with fixed `z=pz=0`.")
            println(io, "- The RF-off baseline orbit is found with a Float64 finite-difference Newton Jacobian to avoid the current ForwardDiff/implicit-integrator fault. All 119 control derivatives in both modes are still computed simultaneously with first-order GTPSA parameters.")
        end
        if !isempty(failures)
            println(io)
            println(io, "## Failures")
            println(io)
            for failure in failures
                println(io, "- `$(failure.mode)`: $(failure.error)")
            end
        end
        println(io)
        println(io, "Detailed per-mode summaries and labeled matrices are stored in the corresponding `bmad_control_response_rf_on` and `bmad_control_response_rf_off` directories.")
    end
    return path
end

function existing_mode_result(mode)
    paths = mode_paths(mode)
    all(isfile, (paths.bmad, paths.scibmad, paths.summary)) || return nothing
    reference = read_labeled_matrix(paths.bmad)
    candidate = read_labeled_matrix(paths.scibmad)
    reference.row_labels == candidate.row_labels || return nothing
    reference.column_labels == candidate.column_labels || return nothing
    runtime_match = match(r"- Runtime: `([0-9.eE+\-]+) s`", read(paths.summary, String))
    elapsed = isnothing(runtime_match) ? NaN : parse(Float64, runtime_match.captures[1])
    return (;
        mode,
        status="complete",
        paths,
        metrics=comparison_metrics(reference, candidate),
        elapsed,
    )
end

function option_value(args, prefix, default)
    hit = findfirst(argument -> startswith(argument, prefix), args)
    return isnothing(hit) ? default : split(args[hit], '='; limit=2)[2]
end

function main(args=ARGS)
    requested = lowercase(option_value(args, "--mode=", "both"))
    modes = requested == "both" ? ["rf_on", "rf_off"] :
            requested == "summary" ? String[] : [requested]
    all(mode -> mode in ("rf_on", "rf_off"), modes) ||
        error("--mode must be rf_on, rf_off, both, or summary")

    results = NamedTuple[]
    failures = NamedTuple[]
    for mode in modes
        try
            push!(results, compare_mode(mode))
        catch exception
            message = sprint(showerror, exception, catch_backtrace())
            push!(failures, (; mode, error=message))
            @error "Control-response comparison failed" mode exception=(exception, catch_backtrace())
        end
    end


    # Keep the root summary complete when only one mode is rerun and the other
    # mode already has valid labeled Bmad/SciBmad matrices on disk.
    handled_modes = Set(vcat([result.mode for result in results],
                             [failure.mode for failure in failures]))
    for mode in ("rf_on", "rf_off")
        mode in handled_modes && continue
        existing = existing_mode_result(mode)
        !isnothing(existing) && push!(results, existing)
    end
    sort!(results; by=result -> result.mode == "rf_on" ? 1 : 2)

    overall = joinpath(HERE, "bmad_scibmad_control_response_summary.md")
    write_overall_summary(overall, results, failures)
    println("\nOverall summary: $overall")
    isempty(failures) || return 1
    return 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main())
end
