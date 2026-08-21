#!/usr/bin/env julia

"""Focused regression for the maintained latest-ring GTPSA response path.

This intentionally tracks one full latest lattice map and its detector map;
it does not generate the 1000-sample benchmark dataset.  The test also checks
that the 21 non-steering controls in the complete registry remain primitive
Float64 zeros in the parameterized model.
"""

using Test
using LinearAlgebra
using SciBmad

const HERE = @__DIR__
include(joinpath(HERE, "benchmark_scibmad.jl"))

const LATEST_RING = :latest

@testset "latest CESR selected-control GTPSA response" begin
    model = load_ring_model(; ring=LATEST_RING, zero_value=0.0, rf_on=true)
    specs = read_control_specs()
    steering_names = [
        spec.name for spec in specs if spec.plane in (:horizontal, :vertical)
    ]

    # Inventory assertions are regression targets for the maintained export;
    # production code obtains these dimensions from model metadata.
    @test length(model.ring.line) == 1177
    @test length(model.metadata.detector_names) == 144
    @test length(steering_names) == 103
    @test length(model.metadata.all_control_names) == 124

    coordinate_count = orbit_coordinate_count(model)
    nominal_solution = find_closed_orbit(
        model.ring;
        v0=zeros(1, coordinate_count),
        coasting_beam=false,
        batch=Val{false}(),
        reltol=1.0e-8,
        abstol=1.0e-10,
        maxiter=100,
        warn=false,
    )
    @test nominal_solution.sol.retcode == SciBmad.BatchSolve.RETCODE_SUCCESS
    nominal_orbit = vec(Float64.(copy(nominal_solution.v0)))
    @test length(nominal_orbit) == 6
    @test all(isfinite, nominal_orbit)

    model_factory = (; kwargs...) -> load_ring_model(; ring=LATEST_RING, kwargs...)
    response = gtpsa_first_order_responses(
        steering_names;
        nominal_orbit,
        rf_on=true,
        reltol=1.0e-8,
        abstol=1.0e-10,
        maxiter=100,
        model_factory,
    )

    @test size(response.A) == (6, 6)
    @test size(response.B) == (6, 103)
    @test size(response.closed_orbit_response) == (6, 103)
    @test size(response.detector_response) == (288, 103)
    @test length(response.detectors) == 144
    @test length(response.observable_labels) == 288
    @test all(isfinite, response.A)
    @test all(isfinite, response.B)
    @test all(isfinite, response.closed_orbit_response)
    @test all(isfinite, response.detector_response)
    @test maximum(abs, response.closure_residual) < 1.0e-10

    # The inactive registry entries are deliberately left as primitive values;
    # this is the condition that avoids promoting unused combined-multipole
    # controls into the latest-ring GTPSA map.
    @test response.model.controls["SK_Q14W"] isa Float64
    @test response.model.controls["RAW_XQUNEING_1"] isa Float64
    @test response.model.controls["RAW_XQUNEING_2"] isa Float64
    @test !(response.model.controls[steering_names[1]] isa Float64)

    incompatible_error = try
        gtpsa_first_order_responses(
            ["SK_Q14W"];
            nominal_orbit,
            rf_on=true,
            reltol=1.0e-8,
            abstol=1.0e-10,
            maxiter=100,
            model_factory,
        )
        nothing
    catch exception
        exception
    end
    @test incompatible_error isa Exception
    @test occursin("SK_Q14W", sprint(showerror, incompatible_error))
    @test occursin("central-difference", sprint(showerror, incompatible_error))

    println((
        ring=:latest_cesr,
        line_elements=length(model.ring.line),
        controls=length(steering_names),
        detectors=length(response.detectors),
        closed_orbit_response_shape=size(response.closed_orbit_response),
        detector_response_shape=size(response.detector_response),
        closure_residual_max=maximum(abs, response.closure_residual),
    ))
end
