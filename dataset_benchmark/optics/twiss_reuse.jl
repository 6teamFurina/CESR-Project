"""
Prototype an allocation-reducing `twiss!` interface using SciBmad's current
internal Twiss stages.

SciBmad 0.4.x does not expose a public `twiss!` method.  This workspace caches
the detector/step lookup, the identity DAMap, and the segment DAMaps/TPS
storage.  The normal-form and returned Twiss table are still constructed for
each call, so this is a conservative measurement of the benefit available
from a future fully in-place implementation.
"""

using LinearAlgebra

struct ReusableTwissWorkspace{B,E,M,S,N,I,ZL,ZP,ZO,ZH}
    beamline::B
    eye::E
    maps::M
    s::S
    names::N
    indices::I
    step_save::Vector{Int}
    zero_lattice_function::ZL
    zero_phase::ZP
    zero_orbit::ZO
    zero_h::ZH
    symplectic_tol::Float64
    in_body_coordinates::Bool
end

function prepare_reusable_twiss_workspace(
    beamline,
    at,
    descriptor,
    closed_orbit;
    symplectic_tol=1e-8,
)
    orbit_matrix = reshape(copy(closed_orbit), 1, 6)
    s, names, indices, step_save = SciBmad._twiss_1(beamline, at)
    concat, eye, zero_lf, zero_phase, zero_orbit, zero_h = SciBmad._twiss_2(
        step_save,
        (orbit_matrix, true),
        descriptor,
        Val{false}(),
        Val{false}(),
    )
    concat || error(
        "Reusable twiss prototype currently requires the concatenating Twiss path",
    )
    maps = SciBmad._twiss_concat_preallocate(step_save, eye)
    return ReusableTwissWorkspace(
        beamline,
        eye,
        maps,
        s,
        names,
        indices,
        step_save,
        zero_lf,
        zero_phase,
        zero_orbit,
        zero_h,
        Float64(symplectic_tol),
        false,
    )
end

function reset_identity_damap!(map, closed_orbit; absolute_scalar=false)
    nv = SciBmad.NNF.nvars(map)
    length(closed_orbit) >= nv || error("Closed orbit is shorter than DAMap v0")
    map.v0 .= view(closed_orbit, 1:nv)
    for component in map.v
        SciBmad.TI.clear!(component)
    end
    SciBmad.NNF.setray!(map.v; v_matrix=I)
    if absolute_scalar
        SciBmad.NNF.setscalar!(map, map.v0)
    end
    return map
end

function twiss!(workspace::ReusableTwissWorkspace, closed_orbit)
    reset_identity_damap!(workspace.eye, closed_orbit)
    if !isempty(workspace.step_save) && first(workspace.step_save) == 0
        reset_identity_damap!(workspace.maps[1], closed_orbit; absolute_scalar=true)
    end

    callback = SciBmad._twiss_3(
        workspace.step_save,
        workspace.maps,
        workspace.in_body_coordinates,
    )
    tracked = SciBmad._twiss_4(workspace.eye, callback, workspace.beamline)
    one_turn_map = SciBmad._twiss_5!(workspace.eye, tracked, workspace.maps)
    tunes, normalizing_map = SciBmad._twiss_6(one_turn_map)
    table = SciBmad._twiss_7(
        normalizing_map,
        one_turn_map,
        workspace.maps,
        workspace.s,
        workspace.names,
        workspace.indices,
        workspace.symplectic_tol,
        workspace.zero_lattice_function,
        workspace.zero_phase,
        workspace.zero_orbit,
        workspace.zero_h,
        Val{false}(),
        Val{false}(),
    )
    return SciBmad.Twiss(true, tunes, table)
end
