#!/usr/bin/env julia

"""Exercise representative Overlay and Group DefExpr controls dynamically."""

using Printf

include(joinpath(@__DIR__, "latest_cesr_scibmad_repaired.jl"))

const SUMMARY = joinpath(@__DIR__, "CONTROL_VALIDATION.md")

function check_slope(variable, element, attribute, delta, expected_slope)
    Core.eval(Main, Expr(:(=), variable, 0.0))
    baseline = Float64(getproperty(element, attribute))
    Core.eval(Main, Expr(:(=), variable, delta))
    changed = Float64(getproperty(element, attribute))
    Core.eval(Main, Expr(:(=), variable, 0.0))
    measured_slope = (changed - baseline) / delta
    error = abs(measured_slope - expected_slope)
    return (;
        variable=String(variable),
        element=String(element.name),
        attribute=String(attribute),
        baseline,
        delta,
        measured_slope,
        expected_slope,
        error,
    )
end

function main()
    checks = [
        check_slope(:H11W_hkick, b11w, :Kn0, 1e-5, -1.52108327900530887e-1),
        check_slope(:H12W_hkick, b12w, :Kn0, 1e-5, -1.52108327900530887e-1 / 2),
        check_slope(:H12W_hkick, b13w, :Kn0, 1e-5, -1.52084639360173829e-1 / 2),
        check_slope(:V09AW_vkick, sex_09aw, :Ks0, 1e-5, 3.67647058823529393),
        check_slope(:SK_Q14W_k1, sex_14w, :Ks1, 1e-5, -1.0),
        check_slope(:H_CANT_S3_hkick, hs3a, :Kn0, 1e-5, -9.09090909090909172 * 0.6506),
        check_slope(:H_CANT_S3_hkick, hs3c, :Kn0, 1e-5, -9.09090909090909172 * 0.6506),
        check_slope(:H_CANT_S3_hkick, hs3b!s1, :Kn0, 1e-5, 20.0),
        check_slope(:H_CANT_S3_hkick, hs3b!s2, :Kn0, 1e-5, 20.0),
        check_slope(:H48W_hkick, b48w!s1, :Kn0, 1e-5, -3.39522373505846908e-1),
        check_slope(:H48W_hkick, b48w!s2, :Kn0, 1e-5, -3.39522373505846908e-1),
        check_slope(:RAW_XQUNEING_1_command, sex_12w, :Kn2, 1e-3, -0.008718),
        check_slope(:RAW_XQUNEING_1_command, sex_27w, :Kn2, 1e-3, -0.517422),
        check_slope(:RAW_XQUNEING_2_command, sex_47e, :Kn2, 1e-3, 0.299675),
    ]
    maximum_error = maximum(check.error for check in checks)

    open(SUMMARY, "w") do io
        println(io, "# Latest CESR control validation")
        println(io)
        println(io, "Representative Overlay and Group variables were changed after lattice " *
                    "construction to verify that `DefExpr` relationships remain live.")
        println(io)
        println(io, "| Control variable | Element | Attribute | Expected slope | Measured slope | Abs. error |")
        println(io, "|---|---|---|---:|---:|---:|")
        for check in checks
            @printf(
                io, "| `%s` | `%s` | `%s` | %.12e | %.12e | %.3e |\n",
                check.variable, check.element, check.attribute,
                check.expected_slope, check.measured_slope, check.error,
            )
        end
        println(io)
        @printf(io, "Maximum tested slope error: `%.3e`.\n", maximum_error)
        println(io)
        println(io, "The tests cover a one-slave bend overlay, a two-slave overlay, vertical " *
                    "and skew correctors, the repaired split-superlord cases, and both " *
                    "sextupole Group commands.")
    end

    for check in checks
        @printf(
            "%s -> %s.%s expected=% .12e measured=% .12e error=%.3e\n",
            check.variable, check.element, check.attribute,
            check.expected_slope, check.measured_slope, check.error,
        )
    end
    @printf("Maximum slope error: %.3e\n", maximum_error)
    maximum_error <= 1e-10 || error("A repaired control relationship failed")
    println("Wrote $SUMMARY")
end

main()
