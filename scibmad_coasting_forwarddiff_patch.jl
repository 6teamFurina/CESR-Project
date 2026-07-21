module SciBmadCoastingForwardDiffPatch

using LinearAlgebra
using SciBmad

export find_closed_orbit_coasting_forwarddiff, coasting_residual

const FD = SciBmad.ForwardDiff

"""Track one six-dimensional phase-space vector through one turn."""
function track_vector(ring, coordinates::AbstractVector)
    length(coordinates) == 6 || throw(ArgumentError("expected six coordinates"))
    bunch = Bunch(v=reshape(copy(coordinates), 1, 6))
    SciBmad.BTBL.check_bl_bunch!(bunch, ring, false)
    track!(bunch, ring)
    return vec(bunch.coords.v)
end

"""
Return the four-dimensional coasting-beam residual with fixed `z` and `pz`.

Only `(x, px, y, py)` participate in the closed-orbit equation.
"""
function coasting_residual(ring, state4::AbstractVector; z=0.0, pz=0.0)
    length(state4) == 4 || throw(ArgumentError("expected four transverse coordinates"))
    input = [state4; z; pz]
    return track_vector(ring, input)[1:4] - state4
end

"""
Construct the 4D residual Jacobian using ForwardDiff while seeding six directions.

BeamTracking's current implicit-integrator ForwardDiff path assumes that every
`Dual` carries six partial derivatives. SciBmad's coasting solver seeds only
four and therefore faults on the CESR lattice. This adapter differentiates the
six-dimensional one-turn residual, then selects its transverse 4-by-4 block.
The Newton solve itself remains the fixed-z, fixed-pz 4D coasting problem.
"""
function residual_and_jacobian(ring, state4::AbstractVector; z=0.0, pz=0.0)
    state6 = [state4; z; pz]
    residual6(v) = track_vector(ring, v) - v
    jacobian6 = FD.jacobian(residual6, state6)
    residual = residual6(state6)[1:4]
    return residual, jacobian6[1:4, 1:4]
end

"""
    find_closed_orbit_coasting_forwarddiff(ring; kwargs...)

Find an RF-off closed orbit with a 4D Newton solve at fixed `z` and `pz`.
The adapter above supplies the Jacobian, while SciBmad's native
`BatchSolve.newton!` performs the Newton iteration.
"""
function find_closed_orbit_coasting_forwarddiff(
    ring;
    v0=zeros(6),
    coasting_beam::Bool=true,
    z=v0[5],
    pz=v0[6],
    reltol=1e-13,
    abstol=1e-13,
    maxiter=100,
)
    coasting_beam || throw(ArgumentError(
        "this patch is only for coasting_beam=true; use SciBmad.find_closed_orbit for RF-on",
    ))
    length(v0) == 6 || throw(ArgumentError("v0 must contain six coordinates"))
    state = Float64.(v0[1:4])
    residual = zeros(4)
    jacobian = zeros(4, 4)

    function value_and_jacobian!(residual_out, jacobian_out, state_in)
        residual_value, jacobian_value =
            residual_and_jacobian(ring, state_in; z, pz)
        residual_out .= residual_value
        jacobian_out .= jacobian_value
        return nothing
    end

    solution = SciBmad.BatchSolve.newton!(
        value_and_jacobian!,
        residual,
        jacobian,
        state;
        reltol,
        abstol,
        maxiter,
        batchdim=nothing,
    )
    solution.retcode == SciBmad.BatchSolve.RETCODE_SUCCESS || error(
        "BatchSolve coasting closed-orbit solve did not converge; " *
        "retcode=$(solution.retcode), residual=$(norm(solution.f, Inf))",
    )

    orbit = [solution.u; Float64(z); Float64(pz)]
    return (;
        v0=reshape(orbit, 1, 6),
        orbit,
        coasting_beam=true,
        iterations=solution.iters,
        residual=norm(solution.f, Inf),
        jacobian=solution.jac,
        sol=solution,
    )
end

end # module
