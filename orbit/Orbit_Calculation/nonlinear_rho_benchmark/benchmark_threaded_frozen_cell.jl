#!/usr/bin/env julia

const HERE_THREAD_TEST = @__DIR__
const CALCULATION_DIR_THREAD_TEST = normpath(joinpath(HERE_THREAD_TEST, ".."))
include(joinpath(CALCULATION_DIR_THREAD_TEST, "benchmark_scibmad.jl"))

const INPUT_THREAD_TEST = joinpath(
    HERE_THREAD_TEST,
    "shared_input",
    "latest_cesr",
    "nonlinear_rho_correctors.csv",
)
const MANIFEST_THREAD_TEST = joinpath(
    HERE_THREAD_TEST,
    "shared_input",
    "latest_cesr",
    "sample_manifest.csv",
)
const RESPONSE_THREAD_TEST = joinpath(
    ORBIT_ROOT,
    "reference",
    "latest_cesr",
    "gtpsa",
    "closed_orbit_response.csv",
)
const OUTPUT_DIR_THREAD_TEST = joinpath(
    HERE_THREAD_TEST,
    "results",
    "latest_cesr",
    "thread_scaling",
)
const MODEL_FACTORY_THREAD_TEST = (; kwargs...) ->
    load_ring_model(; ring=:latest, kwargs...)

function selected_rows_thread_test()
    rows = Int[]
    for (row, line) in enumerate(readlines(MANIFEST_THREAD_TEST)[2:end])
        fields = split(line, ',')
        scenario = fields[2]
        rho = parse(Float64, fields[3])
        scenario == "all" && isapprox(rho, 9.05; atol=1.0e-12, rtol=0.0) &&
            push!(rows, row)
    end
    length(rows) == 600 || error("Expected 600 all/rho=9.05 samples")
    return rows
end

function run_cell_thread_test(model, guess, use_threads)
    timed = @timed frozen_solve_and_track(
        model,
        size(guess.v0, 1),
        guess.nominal_jacobian;
        initial_v0=guess.v0,
        reltol=1.0e-8,
        abstol=1.0e-10,
        maxiter=100,
        use_cpu_multithreading=use_threads,
    )
    return timed.value, timed.time
end

function main_thread_test()
    Threads.nthreads() > 1 || error("Launch this check with more than one Julia thread")
    BLAS.set_num_threads(1)
    samples = read_samples(INPUT_THREAD_TEST)
    selected = selected_rows_thread_test()
    values = Matrix(samples.values[selected, :])
    guess = prepare_initial_guess(
        samples.names,
        values,
        "response-linear";
        response_matrix_cache=RESPONSE_THREAD_TEST,
        recompute_response=false,
        response_method="gtpsa",
        reltol=1.0e-8,
        abstol=1.0e-10,
        maxiter=100,
        model_factory=MODEL_FACTORY_THREAD_TEST,
    )
    model = prepare_batch_model(
        samples.names,
        values;
        model_factory=MODEL_FACTORY_THREAD_TEST,
    )

    run_cell_thread_test(model, guess, false)
    run_cell_thread_test(model, guess, true)
    serial, serial_seconds = run_cell_thread_test(model, guess, false)
    threaded, threaded_seconds = run_cell_thread_test(model, guess, true)

    max_observable_difference = maximum(abs, serial.observables .- threaded.observables)
    max_orbit_difference = maximum(abs, serial.final_v0 .- threaded.final_v0)
    max_closure_difference = maximum(abs, serial.closure_norms .- threaded.closure_norms)
    identical_iterations = serial.iterations == threaded.iterations
    identical_convergence = serial.converged == threaded.converged
    all(serial.converged) || error("Serial reference did not converge every sample")
    all(threaded.converged) || error("Threaded calculation did not converge every sample")
    maximum(threaded.closure_norms) <= 1.0e-10 ||
        error("Threaded calculation exceeded the closure threshold")

    mkpath(OUTPUT_DIR_THREAD_TEST)
    output_path = joinpath(
        OUTPUT_DIR_THREAD_TEST,
        "all_rho9p05_threads$(Threads.nthreads()).toml",
    )
    open(output_path, "w") do io
        TOML.print(io, Dict(
            "ring" => "latest_cesr",
            "scenario" => "all",
            "rho" => 9.05,
            "sample_count" => length(selected),
            "julia_threads" => Threads.nthreads(),
            "blas_threads" => BLAS.get_num_threads(),
            "serial_seconds" => serial_seconds,
            "threaded_seconds" => threaded_seconds,
            "speedup" => serial_seconds / threaded_seconds,
            "serial_solve_seconds" => serial.solve_seconds,
            "threaded_solve_seconds" => threaded.solve_seconds,
            "serial_track_seconds" => serial.track_seconds,
            "threaded_track_seconds" => threaded.track_seconds,
            "max_observable_difference" => max_observable_difference,
            "max_orbit_difference" => max_orbit_difference,
            "max_closure_difference" => max_closure_difference,
            "identical_iterations" => identical_iterations,
            "identical_convergence" => identical_convergence,
            "threaded_max_closure" => maximum(threaded.closure_norms),
        ); sorted=true)
    end
    @printf(
        "Thread scaling all/rho=9.05: serial %.6f s, threaded %.6f s, speedup %.3fx\n",
        serial_seconds,
        threaded_seconds,
        serial_seconds / threaded_seconds,
    )
    @printf(
        "Differences: observables %.3e, orbit %.3e, closure %.3e; output %s\n",
        max_observable_difference,
        max_orbit_difference,
        max_closure_difference,
        output_path,
    )
end

main_thread_test()
