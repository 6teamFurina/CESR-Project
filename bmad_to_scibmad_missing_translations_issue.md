<!--
GitHub issue title:
bmad_to_scibmad cannot translate groups, kicker controls, some expressions, and wigglers
-->

`bmad_to_scibmad` currently cannot produce a complete, loadable SciBmad
lattice for lattices using several common Bmad features.

Tested with:

- Bmad: `20260806-0`
- SciBmad.jl: `0.4.1`
- Beamlines.jl: `0.9.2`

The following problems were encountered while converting a CESR lattice.

## 1. Group definitions are not translated

Bmad `group` definitions are skipped during conversion.

The converter prints:

```text
GROUP ELEMENTS CANNOT YET BE TRANSLATED!
```

However, conversion continues and an output file is still produced. The group
and its control relationships are absent from the generated SciBmad lattice.

For example, a Bmad definition such as:

```text
g: group = {...}
```

has no equivalent definition in the generated Julia file.

## 2. `HKICK` and `VKICK` controls are not translated

Overlay or group controls targeting the `HKICK` and `VKICK` attributes are not
handled by the SciBmad attribute-name translation.

For `HKICK`, the converter prints:

```text
Attribute not yet coded for translation: HKICK
Please report this.
```

For `VKICK`, it prints:

```text
Attribute not yet coded for translation: VKICK
Please report this.
```

The converter still continues instead of returning an error. Because no valid
SciBmad attribute name is returned, the corresponding generated assignment can
be missing or invalid.

This affects corrector-magnet controls in the CESR lattice. In this case, 138
overlay-to-element `HKICK`/`VKICK` target terms cannot be translated directly.

Tilted kickers also require conversion between Bmad kick attributes and SciBmad
multipole components. The expected component mappings are:

```text
HKICK:
    Kn0L = -cos(tilt) * HKICK
    Ks0L = -sin(tilt) * HKICK

VKICK:
    Kn0L = -sin(tilt) * VKICK
    Ks0L =  cos(tilt) * VKICK
```

Therefore, translating only the attribute name would not be sufficient for a
tilted kicker; the control must be distributed between both normal and skew
dipole components.

## 3. Some Bmad expressions are not translated into valid Julia expressions

Control expressions are reconstructed using Bmad function and constant names.
Some emitted names are undefined in Julia, while others have different
semantics.

| Bmad expression | Generated Julia | Result |
|---|---|---|
| `ATAN2(y, x)` | `atan2(y, x)` | `UndefVarError`; Julia uses `atan(y, x)` |
| `MODULO(a, b)` | `modulo(a, b)` | `UndefVarError`; Julia uses `mod(a, b)` |
| `INT(x)` | `int(x)` | `UndefVarError` |
| `NINT(x)` | `nint(x)` | `UndefVarError` |
| `RAN()` | `ran()` | `UndefVarError` |
| `RAN_GAUSS()` | `ran_gauss()` | `UndefVarError` |
| `MASS_OF(...)` | `mass_of(...)` | `UndefVarError` |
| `CHARGE_OF(...)` | `charge_of(...)` | `UndefVarError` |
| `ANTIPARTICLE(...)` | `antiparticle(...)` | `UndefVarError` |
| `PI` | `PI` | `UndefVarError`; Julia provides `pi` |
| `C_LIGHT` | `C_LIGHT` | `UndefVarError` |
| `SINC(x)` | `sinc(x)` | No exception, but different numerical semantics |

For example, generated code containing:

```julia
value = atan2(y, x)
```

fails while loading:

```text
ERROR: LoadError: UndefVarError: `atan2` not defined in `Main`
```

Similarly:

```julia
value = C_LIGHT
```

fails with:

```text
ERROR: LoadError: UndefVarError: `C_LIGHT` not defined in `Main`
```

`SINC` is more problematic because it does not raise an exception:

- Bmad: `SINC(x) = sin(x) / x`
- Julia: `sinc(x) = sin(pi*x) / (pi*x)`

For example:

```text
Bmad SINC(1.0)  = 0.8414709848...
Julia sinc(1.0) = 0.0
```

Thus, the generated lattice can load successfully while silently producing a
different control value.

## 4. Wiggler elements produce undefined Julia constructors

A Bmad wiggler is currently emitted using a `Wiggler(...)` constructor:

```julia
wiggler_name = Wiggler(...)
```

However, `Wiggler` is not defined or exported by the tested Beamlines.jl
version.

The converter does not report an error when writing this element. Instead, the
error occurs later when loading the generated file:

```julia
using Beamlines
include("converted_lattice.jl")
```

which produces:

```text
ERROR: LoadError: UndefVarError: `Wiggler` not defined in `Main`
```

Therefore, `bmad_to_scibmad` reports that the conversion completed, but the
resulting Julia lattice cannot be loaded.
