#!/usr/bin/env julia

using Printf
using SciBmad
using GTPSA

const SCRIPT_DIR = @__DIR__
const PROJECT_ROOT = normpath(joinpath(SCRIPT_DIR, "..", ".."))

include(joinpath(PROJECT_ROOT, "cesr.jl"))
include(joinpath(PROJECT_ROOT, "scibmad_coasting_forwarddiff_patch.jl"))
using .SciBmadCoastingForwardDiffPatch

const DELTA_INDEX = 6

"Return the constant coefficient of a scalar or TPS value as Float64."
constant_term(value) = Float64(GTPSA.scalar(value))

"Return the first derivative with respect to relative momentum delta."
delta_derivative(value) = Float64(value[DELTA_INDEX])

function csv_field(value::AbstractString)
    return '"' * replace(value, '"' => "\"\"") * '"'
end

csv_field(value) = string(value)

function write_results_csv(path, table, tune_1, tune_2, xi_1, xi_2)
    columns = (
        "sample_index",
        "beamline_index",
        "s_m",
        "element_name",
        "tune_1",
        "tune_2",
        "chromaticity_1",
        "chromaticity_2",
        "phi_1_turn",
        "phi_2_turn",
        "accumulated_chromaticity_1",
        "accumulated_chromaticity_2",
        "beta_1_m",
        "beta_2_m",
        "dbeta_1_ddelta_m",
        "dbeta_2_ddelta_m",
        "alpha_1",
        "alpha_2",
        "dalpha_1_ddelta",
        "dalpha_2_ddelta",
    )

    open(path, "w") do io
        println(io, join(columns, ','))
        for i in eachindex(table.s)
            row = (
                i,
                table.beamline_index[i],
                Float64(table.s[i]),
                csv_field(table.name[i]),
                tune_1,
                tune_2,
                xi_1,
                xi_2,
                constant_term(table.phi_1[i]),
                constant_term(table.phi_2[i]),
                delta_derivative(table.phi_1[i]),
                delta_derivative(table.phi_2[i]),
                constant_term(table.beta_1[i]),
                constant_term(table.beta_2[i]),
                delta_derivative(table.beta_1[i]),
                delta_derivative(table.beta_2[i]),
                constant_term(table.alpha_1[i]),
                constant_term(table.alpha_2[i]),
                delta_derivative(table.alpha_1[i]),
                delta_derivative(table.alpha_2[i]),
            )
            println(io, join(csv_field.(row), ','))
        end
    end
end

function nice_ticks(lo, hi, count=6)
    lo == hi && return [lo]
    raw_step = (hi - lo) / max(count - 1, 1)
    magnitude = 10.0^floor(log10(raw_step))
    normalized = raw_step / magnitude
    step = (normalized <= 1 ? 1.0 : normalized <= 2 ? 2.0 : normalized <= 5 ? 5.0 : 10.0) * magnitude
    first_tick = ceil(lo / step) * step
    last_tick = floor(hi / step) * step
    return collect(first_tick:step:last_tick)
end

function svg_escape(text)
    return replace(string(text), '&' => "&amp;", '<' => "&lt;", '>' => "&gt;", '"' => "&quot;")
end

function write_svg_plot(path, s, series; title, ylabel)
    width, height = 1400, 720
    left, right, top, bottom = 112, 40, 62, 88
    plot_width = width - left - right
    plot_height = height - top - bottom

    xlo, xhi = extrema(s)
    values = reduce(vcat, (entry.values for entry in series))
    finite_values = filter(isfinite, values)
    isempty(finite_values) && error("No finite values available for $title")
    ylo, yhi = extrema(finite_values)
    if ylo == yhi
        padding = max(abs(ylo), 1.0) * 0.05
        ylo -= padding
        yhi += padding
    else
        padding = 0.06 * (yhi - ylo)
        ylo -= padding
        yhi += padding
    end

    xpixel(x) = left + (x - xlo) / (xhi - xlo) * plot_width
    ypixel(y) = top + (yhi - y) / (yhi - ylo) * plot_height
    colors = ("#1769aa", "#d1495b", "#2a9d8f", "#7b2cbf")
    x_ticks = nice_ticks(xlo, xhi, 8)
    y_ticks = nice_ticks(ylo, yhi, 7)

    open(path, "w") do io
        println(io, "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"$width\" height=\"$height\" viewBox=\"0 0 $width $height\">")
        println(io, "<rect width=\"100%\" height=\"100%\" fill=\"white\"/>")
        println(io, "<style>text{font-family:Arial,sans-serif;fill:#202124}.tick{font-size:15px}.label{font-size:18px}.title{font-size:24px;font-weight:600}.legend{font-size:16px}</style>")
        println(io, "<text class=\"title\" x=\"$(width / 2)\" y=\"32\" text-anchor=\"middle\">$(svg_escape(title))</text>")

        for value in y_ticks
            y = ypixel(value)
            println(io, "<line x1=\"$left\" y1=\"$y\" x2=\"$(left + plot_width)\" y2=\"$y\" stroke=\"#d9dce1\" stroke-width=\"1\"/>")
            println(io, "<text class=\"tick\" x=\"$(left - 12)\" y=\"$(y + 5)\" text-anchor=\"end\">$(@sprintf("%.5g", value))</text>")
        end
        for value in x_ticks
            x = xpixel(value)
            println(io, "<line x1=\"$x\" y1=\"$top\" x2=\"$x\" y2=\"$(top + plot_height)\" stroke=\"#eef0f2\" stroke-width=\"1\"/>")
            println(io, "<text class=\"tick\" x=\"$x\" y=\"$(top + plot_height + 28)\" text-anchor=\"middle\">$(@sprintf("%.0f", value))</text>")
        end

        println(io, "<rect x=\"$left\" y=\"$top\" width=\"$plot_width\" height=\"$plot_height\" fill=\"none\" stroke=\"#5f6368\" stroke-width=\"1.2\"/>")
        if ylo <= 0 <= yhi
            yzero = ypixel(0.0)
            println(io, "<line x1=\"$left\" y1=\"$yzero\" x2=\"$(left + plot_width)\" y2=\"$yzero\" stroke=\"#777\" stroke-width=\"1.2\"/>")
        end

        for (index, entry) in enumerate(series)
            points = join((@sprintf("%.3f,%.3f", xpixel(s[i]), ypixel(entry.values[i])) for i in eachindex(s)), ' ')
            println(io, "<polyline points=\"$points\" fill=\"none\" stroke=\"$(colors[index])\" stroke-width=\"2.2\" stroke-linejoin=\"round\"/>")
        end

        legend_x = left + 20
        legend_y = top + 26
        for (index, entry) in enumerate(series)
            x = legend_x + (index - 1) * 260
            println(io, "<line x1=\"$x\" y1=\"$legend_y\" x2=\"$(x + 34)\" y2=\"$legend_y\" stroke=\"$(colors[index])\" stroke-width=\"3\"/>")
            println(io, "<text class=\"legend\" x=\"$(x + 43)\" y=\"$(legend_y + 5)\">$(svg_escape(entry.label))</text>")
        end

        println(io, "<text class=\"label\" x=\"$(left + plot_width / 2)\" y=\"$(height - 24)\" text-anchor=\"middle\">s [m]</text>")
        println(io, "<text class=\"label\" transform=\"translate(28 $(top + plot_height / 2)) rotate(-90)\" text-anchor=\"middle\">$(svg_escape(ylabel))</text>")
        println(io, "</svg>")
    end
end

function compute_rf_off()
    mode_name = "RF_off"
    output_dir = joinpath(SCRIPT_DIR, mode_name)
    mkpath(output_dir)

    ring = load_cesr()
    set_cesr_rf!(ring; on=false)

    println("[$mode_name] Finding closed orbit...")
    closed_orbit = find_closed_orbit_coasting_forwarddiff(ring; coasting_beam=true)

    println("[$mode_name] Computing second-order GTPSA Twiss table...")
    optics = twiss(
        ring;
        GTPSA_descriptor=Descriptor(6, 2),
        at=:,
        v0=closed_orbit.v0,
        v0_and_coast=(closed_orbit.v0, closed_orbit.coasting_beam),
    )

    table = optics.table
    xi_1 = delta_derivative(optics.tunes[1])
    xi_2 = delta_derivative(optics.tunes[2])
    tune_1 = constant_term(optics.tunes[1])
    tune_2 = constant_term(optics.tunes[2])

    csv_path = joinpath(output_dir, "nonlinear_twiss.csv")
    write_results_csv(csv_path, table, tune_1, tune_2, xi_1, xi_2)

    s = Float64.(table.s)
    accumulated_xi_1 = delta_derivative.(table.phi_1)
    accumulated_xi_2 = delta_derivative.(table.phi_2)
    dbeta_1 = delta_derivative.(table.beta_1)
    dbeta_2 = delta_derivative.(table.beta_2)
    dalpha_1 = delta_derivative.(table.alpha_1)
    dalpha_2 = delta_derivative.(table.alpha_2)

    write_svg_plot(
        joinpath(output_dir, "chromaticity_along_ring.svg"),
        s,
        [(label="mode 1", values=accumulated_xi_1), (label="mode 2", values=accumulated_xi_2)];
        title="$mode_name accumulated chromaticity",
        ylabel="dphi / ddelta [turn]",
    )
    write_svg_plot(
        joinpath(output_dir, "beta_derivative_along_ring.svg"),
        s,
        [(label="d beta_1 / d delta", values=dbeta_1), (label="d beta_2 / d delta", values=dbeta_2)];
        title="$mode_name beta derivative",
        ylabel="d beta / d delta [m]",
    )
    write_svg_plot(
        joinpath(output_dir, "alpha_derivative_along_ring.svg"),
        s,
        [(label="d alpha_1 / d delta", values=dalpha_1), (label="d alpha_2 / d delta", values=dalpha_2)];
        title="$mode_name alpha derivative",
        ylabel="d alpha / d delta",
    )

    endpoint_xi = (accumulated_xi_1[end], accumulated_xi_2[end])
    @printf("[%s] tunes: Q1=%.12g, Q2=%.12g\n", mode_name, tune_1, tune_2)
    @printf("[%s] chromaticities: xi1=%.12g, xi2=%.12g\n", mode_name, xi_1, xi_2)
    @printf(
        "[%s] end accumulated chromaticities: xi1=%.12g, xi2=%.12g\n",
        mode_name,
        endpoint_xi...,
    )
    println("[$mode_name] Wrote $(length(s)) longitudinal samples to $output_dir")

    return (; mode_name, tune_1, tune_2, xi_1, xi_2, endpoint_xi, samples=length(s))
end

function main(args=ARGS)
    valid_args = ("--mode=rf_off", "--mode=off")
    all(arg -> lowercase(arg) in valid_args, args) || error(
        "Nonlinear Twiss derivatives are computed with RF cavities off; " *
        "run without arguments or use --mode=rf_off.",
    )
    result = compute_rf_off()
    println("Completed nonlinear Twiss calculation for ", result.mode_name)
end

main()
