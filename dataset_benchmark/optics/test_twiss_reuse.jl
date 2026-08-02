include(joinpath(@__DIR__, "benchmark_scibmad_chromatic_optics.jl"))

samples = read_samples(joinpath(ORBIT_DIR, "inputs", "cesr_corrector_samples_1000.csv"))
values = Matrix(samples.values[1:2, :])
response_cache = joinpath(ORBIT_DIR, "reference", "closed_orbit_response_6x119.csv")

rf_on = solve_input_closed_orbits(
    samples.names,
    values;
    response_matrix_cache=response_cache,
    reltol=1e-8,
    abstol=1e-10,
    maxiter=100,
)
setup = prepare_chromatic_optics_model()
coasting = solve_coasting_closed_orbits!(
    setup.model,
    samples.names,
    values,
    Matrix(rf_on.result.final_v0);
    reltol=1e-8,
    abstol=1e-10,
    maxiter=100,
)

apply_sample!(setup.model, samples.names, view(values, 1, :))
reference = twiss(
    setup.model.ring;
    GTPSA_descriptor=setup.descriptor,
    at=setup.detectors,
    v0_and_coast=(reshape(copy(view(coasting.orbits, 1, :)), 1, 6), true),
    spin=false,
    RDTs=false,
    normalizing_map=false,
)
workspace = prepare_reusable_twiss_workspace(
    setup.model.ring,
    setup.detectors,
    setup.descriptor,
    view(coasting.orbits, 1, :),
)

for row in 1:2
    apply_sample!(setup.model, samples.names, view(values, row, :))
    reused = twiss!(workspace, view(coasting.orbits, row, :))
    fresh = twiss(
        setup.model.ring;
        GTPSA_descriptor=setup.descriptor,
        at=setup.detectors,
        v0_and_coast=(reshape(copy(view(coasting.orbits, row, :)), 1, 6), true),
        spin=false,
        RDTs=false,
        normalizing_map=false,
    )
    differences = Float64[]
    for column in (CHROMATIC_TWISS_COLUMNS..., CHROMATIC_ORBIT_COLUMNS...)
        fresh_column = getproperty(fresh.table, column)
        reused_column = getproperty(reused.table, column)
        for index in eachindex(fresh_column)
            push!(differences, abs(constant_term(fresh_column[index]) - constant_term(reused_column[index])))
            push!(differences, abs(delta_derivative(fresh_column[index]) - delta_derivative(reused_column[index])))
        end
    end
    for index in eachindex(fresh.tunes)
        push!(differences, abs(constant_term(fresh.tunes[index]) - constant_term(reused.tunes[index])))
        push!(differences, abs(delta_derivative(fresh.tunes[index]) - delta_derivative(reused.tunes[index])))
    end
    println("row=$row max_abs_difference=$(maximum(differences))")
end
