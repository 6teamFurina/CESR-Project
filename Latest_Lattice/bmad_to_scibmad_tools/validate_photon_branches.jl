#!/usr/bin/env julia

"""Validate the explicit photon-branch registry emitted by the repair."""

using Printf
using SciBmad

const LATTICE_DIR = normpath(joinpath(@__DIR__, ".."))

include(joinpath(LATTICE_DIR, "latest_cesr_scibmad_repaired.jl"))

const BRANCH_CSV = joinpath(LATTICE_DIR, "bmad_reference", "inventory", "bmad_branches.csv")
const INVENTORY_CSV = joinpath(LATTICE_DIR, "bmad_reference", "inventory", "bmad_element_inventory.csv")
const SUMMARY = joinpath(LATTICE_DIR, "scibmad_validation", "PHOTON_BRANCH_VALIDATION.md")

function csv_rows(path)
    lines = readlines(path)
    header = Symbol.(split(first(lines), ','))
    return [NamedTuple{Tuple(header)}(Tuple(split(line, ','))) for line in @view lines[2:end]]
end

function main()
    branch_rows = filter(row -> parse(Int, row.branch) > 0, csv_rows(BRANCH_CSV))
    inventory = filter(row -> parse(Int, row.branch) > 0, csv_rows(INVENTORY_CSV))
    length(branch_rows) == 11 || error("Expected eleven Bmad photon branches")
    length(latest_photon_branches) == 11 || error("Expected eleven registered photon branches")
    length(latest_photon_fork_targets) == 11 || error("Expected eleven registered photon forks")

    rows = NamedTuple[]
    for branch in branch_rows
        index = parse(Int, branch.branch)
        name = branch.branch_name
        line = latest_photon_branch(name)
        bmad_elements = filter(row -> parse(Int, row.branch) == index, inventory)
        bmad_length = maximum(parse(Float64, row.s) for row in bmad_elements)
        expected_count = parse(Int, branch.n_ele_track)
        length_error = abs(Float64(line.line[end].s_downstream) - bmad_length)

        reference_output = track_latest_photon_branch(name)
        probe = LatestPhotonRay(x=1e-3, xprime=2e-4, y=-8e-4, yprime=1.5e-4)
        probe_output = track_latest_photon_branch(name; ray=probe)
        push!(rows, (
            branch=index,
            name,
            expected_count,
            scibmad_count=length(line.line),
            bmad_length,
            scibmad_length=Float64(line.line[end].s_downstream),
            length_error,
            reference_output,
            probe_output,
        ))
    end

    all(row.expected_count == row.scibmad_count for row in rows) || error("Photon branch element-count mismatch")
    maximum(row.length_error for row in rows) <= 1e-12 || error("Photon branch length mismatch")
    all(row.reference_output.x == 0.0 && row.reference_output.y == 0.0 for row in rows) ||
        error("Reference photon left the local branch axis")
    maximum(abs(row.reference_output.path_length - row.bmad_length) for row in rows) <= 1e-12 ||
        error("Reference photon path-length mismatch")
    all(isfinite(row.probe_output.x) && isfinite(row.probe_output.y) for row in rows) ||
        error("Non-finite paraxial photon output")

    open(SUMMARY, "w") do io
        println(io, "# Latest CESR photon branch registry validation")
        println(io)
        println(io, "All eleven Bmad photon fork targets are registered as independently " *
                    "queryable Beamlines branches and runnable branch-local paraxial ray lines.")
        println(io)
        println(io, "| Bmad branch | Name | Elements | Length (m) | Length error (m) |")
        println(io, "|---:|---|---:|---:|---:|")
        for row in rows
            @printf(
                io, "| %d | `%s` | %d | %.12g | %.3e |\n",
                row.branch, row.name, row.scibmad_count,
                row.scibmad_length, row.length_error,
            )
        end
        println(io)
        println(io, "`latest_photon_branch(name)` and `latest_photon_branch_for_fork(name)` " *
                    "provide lookup, while `track_latest_photon_branch(name; ray=...)` " *
                    "propagates a paraxial ray through the branch drifts. The reference ray " *
                    "remains on axis and accumulates the exact Bmad branch length.")
        println(io)
        println(io, "The archived flat mirrors have `REF_TILT=-pi/2`, " *
                    "`GRAZE_ANGLE=0.004 rad`, and 10 keV reference energy. They remain " *
                    "identity markers in branch-local coordinates because the Bmad branch " *
                    "already defines the reflected reference frame. This interface does not " *
                    "model reflectivity, apertures, curvature, or off-axis mirror scattering; " *
                    "SciBmad/Beamlines currently has no photon `Fork`/`Mirror` tracker.")
    end
    println("Validated $(length(rows)) photon branches")
    println("Fork registry entries: $(length(latest_photon_fork_targets))")
    println("Mirror placeholders: $(latest_photon_mirror_placeholders)")
    println("Wrote $SUMMARY")
end

main()
