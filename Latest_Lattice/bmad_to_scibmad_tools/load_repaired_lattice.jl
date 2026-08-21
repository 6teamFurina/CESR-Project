"""Load and inventory the first repaired SciBmad representation."""

include(joinpath(@__DIR__, "..", "latest_cesr_scibmad_repaired.jl"))

wigglers = filter(element -> uppercase(String(element.kind)) == "WIGGLER", cesr.line)
println("repaired_export_loaded=true")
println("cesr_type=", typeof(cesr))
println("cesr_elements=", length(cesr.line))
println("cesr_length=", sum(element.L for element in cesr.line))
println("pc_ref=", cesr.pc_ref)
println("species_ref=", cesr.species_ref)
println("wiggler_segments=", length(wigglers))
for element in wigglers
    println(
        "wiggler name=", element.name,
        " L=", element.L,
        " phase=", element.four_potential_params[3],
        " steps=", element.tracking_method.n_steps,
    )
end
