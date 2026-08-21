#!/usr/bin/env julia

using LinearAlgebra
using Printf
using SciBmad
using Test

include(joinpath(@__DIR__, "..", "cesr.jl"))
include(joinpath(@__DIR__, "..", "scibmad_coasting_forwarddiff_patch.jl"))
using .SciBmadCoastingForwardDiffPatch

function finite_difference_closed_orbit(
    ring;
    v0=zeros(6),
    z=v0[5],
    pz=v0[6],
    abstol=1e-13,
    maxiter=20,
    max_backtracks=12,
)
    state = Float64.(v0[1:4])

    for iteration in 0:maxiter
        residual = coasting_residual(ring, state; z, pz)
        residual_norm = norm(residual, Inf)
        residual_norm <= abstol && return (;
            orbit=[state; Float64(z); Float64(pz)],
            iterations=iteration,
            residual=residual_norm,
        )
        iteration == maxiter && break

        jacobian = zeros(4, 4)
        for coordinate in 1:4
            step_size = 1e-7 * max(1.0, abs(state[coordinate]))
            plus = copy(state)
            minus = copy(state)
            plus[coordinate] += step_size
            minus[coordinate] -= step_size
            jacobian[:, coordinate] .= (
                coasting_residual(ring, plus; z, pz) -
                coasting_residual(ring, minus; z, pz)
            ) / (2step_size)
        end

        step = -(jacobian \ residual)
        accepted = false
        scale = 1.0
        for _ in 0:max_backtracks
            trial = state + scale * step
            if norm(coasting_residual(ring, trial; z, pz), Inf) < residual_norm
                state = trial
                accepted = true
                break
            end
            scale *= 0.5
        end
        accepted || error(
            "finite-difference coasting closed-orbit line search failed at iteration $iteration",
        )
    end

    error("finite-difference coasting closed-orbit solve did not converge")
end

ring = load_cesr()
set_cesr_rf!(ring; on=false)

println("Computing RF-off orbit with patched 4D ForwardDiff...")
forwarddiff = find_closed_orbit_coasting_forwarddiff(ring; coasting_beam=true)

println("Computing RF-off orbit with Float64 central differences...")
finite_difference = finite_difference_closed_orbit(ring)

fd_orbit = forwarddiff.orbit
finite_orbit = finite_difference.orbit
orbit_difference = norm(fd_orbit - finite_orbit, Inf)
fd_closure = norm(coasting_residual(ring, fd_orbit[1:4]), Inf)
finite_closure = norm(coasting_residual(ring, finite_orbit[1:4]), Inf)

_, forwarddiff_jacobian =
    SciBmadCoastingForwardDiffPatch.residual_and_jacobian(ring, fd_orbit[1:4])
finite_jacobian = zeros(4, 4)
for coordinate in 1:4
    step_size = 1e-7 * max(1.0, abs(fd_orbit[coordinate]))
    plus = copy(fd_orbit[1:4])
    minus = copy(fd_orbit[1:4])
    plus[coordinate] += step_size
    minus[coordinate] -= step_size
    finite_jacobian[:, coordinate] .= (
        coasting_residual(ring, plus) - coasting_residual(ring, minus)
    ) / (2step_size)
end
jacobian_difference = norm(forwarddiff_jacobian - finite_jacobian, Inf)

@printf("Patched ForwardDiff orbit:  [%s]\n", join((@sprintf("%.15e", x) for x in fd_orbit), ", "))
@printf("Finite-difference orbit:    [%s]\n", join((@sprintf("%.15e", x) for x in finite_orbit), ", "))
@printf("Maximum orbit difference:  %.6e\n", orbit_difference)
@printf("ForwardDiff closure:        %.6e\n", fd_closure)
@printf("Finite-difference closure:  %.6e\n", finite_closure)
@printf("Maximum Jacobian row-sum difference: %.6e\n", jacobian_difference)
@printf("Iterations (ForwardDiff / finite): %d / %d\n", forwarddiff.iterations, finite_difference.iterations)

@testset "RF-off patched ForwardDiff closed orbit" begin
    @test forwarddiff.coasting_beam
    @test forwarddiff.sol.retcode == SciBmad.BatchSolve.RETCODE_SUCCESS
    @test fd_orbit[5:6] == [0.0, 0.0]
    @test fd_closure <= 1e-12
    @test finite_closure <= 1e-12
    @test orbit_difference <= 1e-10
    @test jacobian_difference <= 1e-6
end

println("RF-off patched ForwardDiff comparison passed.")
