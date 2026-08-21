#!/usr/bin/env julia

"""Run resumable truth-only exact-11 scans for one shard of target sextupoles."""

include(joinpath(@__DIR__, "run_exact_11_scan.jl"))

function main(args=ARGS)
    defaults = Dict(
        "targets" => "all",
        "shard-index" => "1",
        "shard-count" => "1",
        "quadrupole-fraction" => "0.001",
        "candidate-csv" => joinpath(
            AFFINITY_HERE,
            "results",
            "scibmad_latest",
            "selection",
            "quadrupole_sets_by_sextupole.csv",
        ),
        "output-dir" => joinpath(@__DIR__, "results", "scans"),
    )
    options = parse_exact11_options(defaults, args)
    shard_index = parse(Int, options["shard-index"])
    shard_count = parse(Int, options["shard-count"])
    1 <= shard_index <= shard_count || error("shard-index must lie in 1:shard-count")
    inventory_names = String.(getproperty.(active_sextupole_inventory(cesr), :name))
    requested = lowercase(options["targets"]) == "all" ? inventory_names :
        uppercase.(strip.(split(options["targets"], ',')))
    unknown = setdiff(requested, inventory_names)
    isempty(unknown) || error("Unknown targets: $(join(unknown, ", "))")
    selected = [
        name for (index, name) in enumerate(inventory_names)
        if name in requested && mod(index - 1, shard_count) + 1 == shard_index
    ]
    output_root = abspath(options["output-dir"])
    completed = 0
    skipped = 0
    for (index, target) in enumerate(selected)
        metadata_path = joinpath(output_root, safe_name(target), "scan_metadata.toml")
        if isfile(metadata_path)
            skipped += 1
            println("Skip existing $target ($index/$(length(selected)))")
            continue
        end
        println("Start $target ($index/$(length(selected))) shard $shard_index/$shard_count")
        run_scan([
            "--target=$target",
            "--scenario-mode=truth-only",
            "--quadrupole-fraction=$(options["quadrupole-fraction"])",
            "--candidate-csv=$(options["candidate-csv"])",
            "--output-dir=$output_root",
        ])
        completed += 1
    end
    println("Shard $shard_index/$shard_count complete: generated=$completed skipped=$skipped")
    return 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main())
end
