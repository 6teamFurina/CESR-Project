using Test
using Beamlines

const PROJECT_ROOT = normpath(joinpath(@__DIR__, ".."))
include(joinpath(PROJECT_ROOT, "cesr_model.jl"))

element(ring, name) = only(filter(ele -> uppercase(String(ele.name)) == name, ring.line))
elements_with_prefix(ring, prefix) =
    filter(ele -> startswith(uppercase(String(ele.name)), prefix), ring.line)

@testset "CESR DefExpr controls" begin
    ring = load_cesr()

    hv01w = element(ring, "HV01W")
    b04w = element(ring, "B04W")
    b05w = element(ring, "B05W")
    sex08w = element(ring, "SEX_08W")
    q00w_slices = elements_with_prefix(ring, "Q00W!")
    q00e_slices = elements_with_prefix(ring, "Q00E!")

    baselines = (
        hv01w_kn0l = hv01w.Kn0L,
        b04w_kn0l = b04w.Kn0L,
        b05w_kn0l = b05w.Kn0L,
        sex08w_kn2 = sex08w.Kn2,
        q00w_tilt = [ele.tilt for ele in q00w_slices],
        q00e_tilt = [ele.tilt for ele in q00e_slices],
    )

    controls = attach_cesr_controls!(ring)

    @test length(controls.overlays) == 119
    @test length(controls.groups) == 26
    @test hv01w.Kn0L == baselines.hv01w_kn0l
    @test b04w.Kn0L == baselines.b04w_kn0l
    @test b05w.Kn0L == baselines.b05w_kn0l
    @test sex08w.Kn2 == baselines.sex08w_kn2
    @test [ele.tilt for ele in q00w_slices] == baselines.q00w_tilt
    @test [ele.tilt for ele in q00e_slices] == baselines.q00e_tilt

    controls["H01W"] = 1.0e-4
    @test hv01w.Kn0L ≈ baselines.hv01w_kn0l + 1.0e-4

    controls["HB01"] = 2.0e-4
    controls["RAW_PRETZING_4"] = 1.0
    expected_hb01_delta = 0.5 * (2.0e-4 - 5.5886e-4)
    @test b04w.Kn0L ≈ baselines.b04w_kn0l + expected_hb01_delta
    @test b05w.Kn0L ≈ baselines.b05w_kn0l + expected_hb01_delta

    controls["RAW_XQUNE_1"] = 1.0
    @test sex08w.Kn2 ≈ baselines.sex08w_kn2 - 0.6415e-3

    controls["ASYM_Q0"] = 0.01
    @test [ele.tilt for ele in q00w_slices] ≈ baselines.q00w_tilt .+ 0.01
    @test [ele.tilt for ele in q00e_slices] ≈ baselines.q00e_tilt .- 0.01

    # Attaching controls to one copied ring must not mutate the static template.
    untouched = load_cesr()
    @test element(untouched, "HV01W").Kn0L == baselines.hv01w_kn0l
    @test element(untouched, "B04W").Kn0L == baselines.b04w_kn0l

    model = load_cesr_model()
    @test length(model.controls.overlays) == 119
    @test length(model.controls.groups) == 26
    @test element(model.ring, "HV01W").Kn0L == baselines.hv01w_kn0l
end
