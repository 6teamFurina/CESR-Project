"""Bmad-linearized coherent alignment controls for the twelve CESR girders."""

const LATEST_GIRDER_COEFFICIENT_PATH = joinpath(
    @__DIR__, "..", "bmad_reference", "girder", "alignment_coefficients.csv",
)
const LATEST_GIRDER_PARAMETER_NAMES = (
    :x_offset, :y_offset, :z_offset, :x_pitch, :y_pitch, :tilt,
)

function load_latest_girder_coefficients(path=LATEST_GIRDER_COEFFICIENT_PATH)
    lines = readlines(path)
    header = split(first(lines), ',')
    column = Dict(name => index for (index, name) in enumerate(header))
    result = NamedTuple[]
    for line in @view lines[2:end]
        values = split(line, ',')
        coefficient = parse(Float64, values[column["coefficient"]])
        abs(coefficient) > 1e-15 || continue
        push!(result, (
            girder=values[column["girder"]],
            member_index=parse(Int, values[column["member_index"]]),
            member_name=values[column["member_name"]],
            girder_parameter=Symbol(values[column["girder_parameter"]]),
            scibmad_property=Symbol(values[column["scibmad_property"]]),
            coefficient,
        ))
    end
    return result
end

const LATEST_GIRDER_COEFFICIENTS = load_latest_girder_coefficients()
const LATEST_GIRDER_NAMES = sort!(unique(row.girder for row in LATEST_GIRDER_COEFFICIENTS))

"""
    set_latest_girder!(ring, name; x_offset=0, y_offset=0, z_offset=0,
                       x_pitch=0, y_pitch=0, tilt=0)

Apply a coherent small-displacement girder transform to the corresponding
SciBmad member elements. The mapping is the central-difference Jacobian of
Bmad's own girder geometry at the nominal lattice, including curved-floor
coordinate rotations and longitudinally varying offsets. Values are absolute
girder settings; calling with all zeros resets the girder.
"""
function set_latest_girder!(
    ring,
    name::AbstractString;
    x_offset=0.0,
    y_offset=0.0,
    z_offset=0.0,
    x_pitch=0.0,
    y_pitch=0.0,
    tilt=0.0,
)
    girder = uppercase(String(name))
    girder in LATEST_GIRDER_NAMES || throw(ArgumentError(
        "Unknown latest CESR girder $name; expected one of $(LATEST_GIRDER_NAMES)",
    ))
    settings = Dict(
        :x_offset => x_offset,
        :y_offset => y_offset,
        :z_offset => z_offset,
        :x_pitch => x_pitch,
        :y_pitch => y_pitch,
        :tilt => tilt,
    )
    totals = Dict{Tuple{Int,Symbol},Float64}()
    member_names = Dict{Int,String}()
    for row in LATEST_GIRDER_COEFFICIENTS
        row.girder == girder || continue
        key = (row.member_index, row.scibmad_property)
        totals[key] = get(totals, key, 0.0) + row.coefficient * settings[row.girder_parameter]
        member_names[row.member_index] = row.member_name
    end
    isempty(totals) && error("No alignment coefficients found for $girder")

    member_indices = sort!(collect(keys(member_names)))
    for index in member_indices
        element = ring.line[index]
        expected = lowercase(replace(member_names[index], "#" => "!s"))
        lowercase(String(element.name)) == expected || error(
            "$girder member index $index is $(element.name), expected $(member_names[index])",
        )
        for property in (:x_offset, :y_offset, :z_offset, :x_rot, :y_rot, :tilt)
            setproperty!(element, property, get(totals, (index, property), 0.0))
        end
    end
    return member_indices
end
