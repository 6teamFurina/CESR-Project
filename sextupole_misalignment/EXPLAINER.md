# Inferring CESR Sextupole Magnetic Centers from BPM Orbit Responses

> **Core idea.** A static closed orbit cannot uniquely identify a sextupole
> offset. Instead, we deliberately vary the target sextupole strength $K_2$
> while steering the beam through signed local orbit bumps. Taking the
> $K_2$-odd and bump-odd response leaves a two-parameter inverse problem for
> the horizontal and vertical magnetic-center coordinates.

> **Current result.** In the synthetic latest-lattice SciBmad benchmark,
> finite BPM information predicts the local bump at sub-micrometer accuracy.
> A fixed physics template, covariance-matched estimation, and time-series
> drift filtering give a two-dimensional center RMSE of
> $20.297~\mu\mathrm{m}$ with the assumed BPM white noise and random-walk
> drift. A separate noise-free, all-SciBmad-forward observation-derivative
> benchmark reaches $5.128~\mu\mathrm{m}$, but has not yet been validated
> with measurement or calibration errors. Quadrupole misalignment remains the
> largest unresolved problem because it can move the closed orbit outside the
> present bump and inverse ranges before the sextupole scan begins.

---

## 1. What quantity are we estimating?

Let the two-dimensional magnetic center of target sextupole $s$ be

$$
\mathbf c_s=(c_{x,s},c_{y,s}),
$$

and let the local closed orbit at zero bump and nominal $K_2$ be

$$
\mathbf z_{s0}=(x_{s0},y_{s0}).
$$

The quantity directly identifiable by beam-based alignment is the
**beam-relative center**

$$
\mathbf d_s=\mathbf c_s-\mathbf z_{s0}.
$$

It is the local bump displacement required to place the beam on the magnetic
center. Recovering the absolute mechanical offset $\mathbf c_s$ additionally
requires an estimate of $\mathbf z_{s0}$ and an alignment reference. Without
that reference, a common translation of the orbit and magnetic center is
non-identifiable.

A static orbit is also insufficient because offsets and errors in other
sextupoles, quadrupoles, correctors, and BPMs can create similar observations.
The identifiable signal is therefore the response to a known intervention,
not a single machine snapshot.

---

## 2. How is the dataset acquired?

```text
Select one target sextupole
        |
        v
Apply four signed local bumps: +x, -x, +y, -y
        |
        v
Scan the target K2 at negative, zero, and positive values
        |
        v
Record x/y closed orbit at all 111 BPMs: 222 channels
        |
        v
Take the K2-odd contrast, followed by the bump-odd contrast
        |
        v
Infer the two-dimensional beam-relative magnetic center
```

The current high-signal protocol is:

| Item | Current setting |
|---|---|
| Lattice | `Latest_Lattice/latest_cesr_scibmad_repaired.jl`, RF on |
| Targets | All 76 active normal sextupoles, scanned one at a time |
| Hidden machines | Four realizations per target; target offset uniform over $[-350,+350]~\mu\mathrm{m}$, other 75 sextupoles at $300~\mu\mathrm{m}$ RMS per plane |
| Local bumps | $(0,0),(\pm1.5~\mathrm{mm},0),(0,\pm1.5~\mathrm{mm})$ |
| $K_2$ levels | $-0.10,0,+0.10~\mathrm{m}^{-3}$ |
| Observation per state | 111 BPMs times two planes |
| Atomic inverse example | One complete bump-by-$K_2$ scan tensor for one target in one hidden machine, not one BPM row |

One baseline tensor contains

$$
76\times4\times5\times3=4560
$$

exact SciBmad closed-orbit states. A paired drift-response secant scan doubles
the exact-state total to 9,120. All 608 signed maintained-protocol states pass
the fourth-order local-model validity gate of
$1~\mu\mathrm{m}/1~\mu\mathrm{rad}$. This is a **model-validity** result,
not approval of corrector or power-supply ranges, aperture, lifetime,
interlocks, or machine operation.

### Repeated time-series acquisition

The four signed bumps combined with $K_2=\pm0.10~\mathrm{m}^{-3}$ form the
eight signal states that carry the center information. The selected
time-series protocol uses 3,072 reads per signal state:

- $8\times3072=24{,}576$ core BPM acquisitions per target;
- every 256 signal cycles, and at the final cycle, a fixed-bump
  $K_2:0,+,0,-,0$ reference block;
- 32 calibration reads for each of the four $K_2=0$ reference baselines,
  with their finite uncertainty retained in the filter;
- 24,860 acquisitions per target, or 1,889,360 for all 76 targets before
  switching and settling overhead.

These are simulated BPM acquisitions, not 1,889,360 independent SciBmad
closed-orbit solves. Exact SciBmad states define the deterministic signal and
drift secant; repeated white noise and random-walk histories are propagated
through analytic covariance and time-series models.

---

## 3. From neighboring BPMs to the orbit at the sextupole

A real machine measures position only at BPMs, not inside the target
sextupole. We infer the relative local bump using the nearest upstream and
downstream BPMs together with nominal SciBmad transport.

For each nonzero bump, first subtract the zero-bump state:

$$
\Delta\mathbf y_{\mathrm{BPM}}(b)
=\mathbf y_{\mathrm{BPM}}(b,K_{2,0})
-\mathbf y_{\mathrm{BPM}}(0,K_{2,0}).
$$

Known corrector commands provide the nominal prediction at the BPMs and target.
Let $\mathbf r_u$ and $\mathbf r_d$ be measured-minus-model position
residuals at the upstream and downstream BPMs. Partition the transverse
upstream-to-downstream transport so that

$$
\mathbf r_d=A_{ud}\mathbf r_u+B_{ud}\mathbf p_u,
\qquad
\mathbf p_u=B_{ud}^{+}(\mathbf r_d-A_{ud}\mathbf r_u),
$$

where $B_{ud}^{+}$ is the pseudoinverse and
$\mathbf p_u=(\Delta p_x,\Delta p_y)$ contains the two unmeasured angles.
Transport the reconstructed residual to the sextupole:

$$
\Delta\mathbf z_{s,\mathrm{res}}
=A_{us}\mathbf r_u+B_{us}\mathbf p_u,
$$

and add it to the command prediction:

$$
\Delta\mathbf z_{s,\mathrm{pred}}
=\Delta\mathbf z_{s,\mathrm{command}}
+\Delta\mathbf z_{s,\mathrm{res}}.
$$

### Local-orbit prediction error

The following noise-free comparison covers all 76 targets, eight hidden
machines per target, and four nonzero bumps per machine:

| Predictor | x RMSE | y RMSE | 2D RMSE | P90 | Maximum |
|---|---:|---:|---:|---:|---:|
| Corrector command only | $10.302~\mu\mathrm{m}$ | $20.451~\mu\mathrm{m}$ | $22.899~\mu\mathrm{m}$ | $37.964~\mu\mathrm{m}$ | $91.689~\mu\mathrm{m}$ |
| Nearest upstream/downstream BPM transport | $0.237~\mu\mathrm{m}$ | $0.185~\mu\mathrm{m}$ | **$0.301~\mu\mathrm{m}$** | **$0.250~\mu\mathrm{m}$** | $4.385~\mu\mathrm{m}$ |
| Full-ring effective-corrector MAP | $7.048~\mu\mathrm{m}$ | $13.665~\mu\mathrm{m}$ | $15.375~\mu\mathrm{m}$ | $24.095~\mu\mathrm{m}$ | $96.678~\mu\mathrm{m}$ |

The two-sided BPM method improves all 76 targets. Its largest transverse
momentum-block condition number is 4.290. The global MAP result is worse not
because more BPMs are intrinsically harmful, but because a local bump is
designed to leak very little orbit to distant BPMs. A full-ring effective-
corrector basis can fit the visible residual while losing the unobserved local
null-space component at the sextupole.

When the two-sided coordinates are passed to the earlier all-111-BPM center
inverse, the center RMSE is $5.864~\mu\mathrm{m}$, essentially identical to
the $5.870~\mu\mathrm{m}$ internal-orbit oracle baseline. Their center-error
vectors differ by only $0.192~\mu\mathrm{m}$ RMS. In that noise-free study,
local-orbit reconstruction is no longer the limiting error.

Stable BPM offsets cancel in the zero-bump difference. BPM gain, roll, white
noise, missing channels, corrector calibration, and transport-model errors do
not cancel automatically. The sub-micrometer result is therefore not a
model-free machine-accuracy claim.

---

## 4. Basic magnetic-center inverse

Changing $K_2$ in a normal sextupole converts the beam-center displacement
into normal- and skew-quadrupole feed-down. SciBmad/GTPSA transports the two
local sources

$$
q_n=\frac{x^2-y^2}{2},
\qquad
q_s=xy
$$

to every BPM, producing fixed full-ring templates $P_n$ and $P_s$. Let
$\mathbf O$ be the full BPM orbit vector. First form the symmetric
$K_2$ slope

$$
\mathbf S(b)=
\frac{\mathbf O(b,K_+)-\mathbf O(b,K_-)}{K_+-K_-},
$$

then take the bump-odd gradients

$$
\mathbf G_x=\frac{\mathbf S(+b_x)-\mathbf S(-b_x)}{2b},
\qquad
\mathbf G_y=\frac{\mathbf S(+b_y)-\mathbf S(-b_y)}{2b}.
$$

With the current SciBmad sign convention, the beam-relative center
$\mathbf d=(d_x,d_y)$ enters linearly:

$$
\mathbf G_x=-P_n d_x-P_s d_y,
\qquad
\mathbf G_y=-P_s d_x+P_n d_y.
$$

Stacking both gradients and all BPM x/y channels gives

$$
\mathbf g=A\mathbf d+\boldsymbol\epsilon.
$$

Only $d_x$ and $d_y$ are unknown. Given measurement covariance $C$, the
generalized least-squares estimator is

$$
\widehat{\mathbf d}
=(A^{T}C^{-1}A)^{-1}A^{T}C^{-1}\mathbf g.
$$

Three projections are essential:

1. The $K_2$-odd contrast removes $K_2$-independent closed orbit and stable
   BPM offsets.
2. The bump-odd contrast removes bump-even and bump-independent terms.
3. The design matrix $A$ remains fixed by the physics model. It is not
   re-normalized from each noisy scan, which would allow noise to define the
   apparent scale of individual BPM channels.

---

## 5. Noise-free SciBmad observation-map method

This is a separate structural benchmark from the stochastic fixed-template
inverse discussed below. Every forward state is an exact RF-on SciBmad
closed-orbit solution on the repaired latest lattice. The only machine errors
are fixed x/y offsets on all 76 active normal sextupoles. The benchmark has:

- no BPM noise, offset, gain, roll, or missing channels;
- no corrector or target-$K_2$ calibration error;
- no time drift;
- no quadrupole strength, roll, or alignment error;
- one hidden machine realization per target.

The target offset is uniform over $[-350,+350]~\mu\mathrm{m}$, while the
other 75 sextupoles have independent $300~\mu\mathrm{m}$ RMS Gaussian
offsets per plane. For each target, the scan uses a $5\times5$ two-plane bump
grid over $\pm0.5~\mathrm{mm}$ and target-$K_2$ levels
$-0.02,-0.01,0,+0.01,+0.02~\mathrm{m}^{-3}$. This gives 125 states per
target and 9,500 exact SciBmad states in total. The deployable inputs are the
111-BPM orbit tensor and the known bump and $K_2$ commands. Exact
target-local orbit and alignment are loaded only after fitting and are used
solely for scoring.

### Observation-derivative construction

For BPM channel $j$, first fit the measured $K_2$ derivative as a
quadratic function of the BPM-predicted local orbit $\mathbf z=(x,y)$:

$$
S_j(\mathbf z)
=a_j+\mathbf g_j^T\mathbf z
+\frac{1}{2}\mathbf z^T H_j\mathbf z.
$$

The common magnetic center satisfies

$$
H_j\mathbf d=-\mathbf g_j
$$

for all retained BPM channels. The scaled stacked system is solved for the
same two center coordinates. Unlike the fixed local-source model, this method
learns the finite-amplitude observable derivative directly from the exact
SciBmad scan.

A related empirical observation-map method fits the raw surface
$O_j(x,y,K_2)$ with all identifiable monomials through total order $N$,
differentiates that polynomial analytically with respect to $K_2$, and finds
the common two-dimensional root.

### Noise-free result

| Machine-facing method | Beam-relative 2D RMSE | Median | P90 |
|---|---:|---:|---:|
| Physical two-source inverse, linear $K_2$ derivative | $6.394706~\mu\mathrm{m}$ | $5.000718~\mu\mathrm{m}$ | $8.956349~\mu\mathrm{m}$ |
| **Direct quadratic observable derivative** | **$5.128185~\mu\mathrm{m}$** | **$3.738924~\mu\mathrm{m}$** | **$7.549159~\mu\mathrm{m}$** |
| Command-space derivative transformed by the local-orbit chain rule | $5.465936~\mu\mathrm{m}$ | $3.981950~\mu\mathrm{m}$ | $7.771729~\mu\mathrm{m}$ |
| Empirical observation Taylor map, order 3 | $33.008657~\mu\mathrm{m}$ | $27.946649~\mu\mathrm{m}$ | $47.192103~\mu\mathrm{m}$ |
| Empirical observation Taylor map, order 4 | $6.362654~\mu\mathrm{m}$ | $4.656776~\mu\mathrm{m}$ | $9.719398~\mu\mathrm{m}$ |
| Empirical observation Taylor map, order 5 | $5.724636~\mu\mathrm{m}$ | $4.490103~\mu\mathrm{m}$ | $8.573020~\mu\mathrm{m}$ |

The BPM-predicted local orbit itself has $0.039371~\mu\mathrm{m}$ 2D RMSE
over all 125 states, confirming that local-orbit inference is not the limiting
error in this idealized test. Total order 3 is inadequate; order 4 removes most
of its truncation bias, and order 5 improves further. The direct quadratic
observable-derivative fit remains the best full-inventory method in this
particular scan.

This $5.128~\mu\mathrm{m}$ result is **only a noise-free structural
validation** with one hidden realization per target. It has not been tested
against BPM noise or calibration, drift, missing channels, actuator errors, or
quadrupole errors and must not be presented as expected machine accuracy or
compared directly with the $20.297~\mu\mathrm{m}$ stochastic result.

Direct high-order SciBmad/GTPSA $K_2$-offset maps remain a subset diagnostic.
The lower-level GTPSA path still terminates on `invalid domain sqrt(0)` for at
least `SEX_14W` and `SEX_44E`. Therefore the full 76-target high-order
comparison uses empirical observation Taylor maps fitted to exact SciBmad scan
states, not a silently incomplete set of direct GTPSA coefficients.

---

## 6. Input-amplitude and model-validity envelope

The maintained $\pm1.5~\mathrm{mm}$ bump and
$\Delta K_2=\pm0.10~\mathrm{m}^{-3}$ are protocol design points, not input
limits. A separate scan expands the two model-based bump knobs and target
$\Delta K_2$ along axes, diagonals, fixed-input rays, and joint-scaling rays.
It evaluates all six coordinates at the target entrance and exit.

The working validity gate requires, for every signed state at a sampled radius:

- exact RF-on closed-orbit convergence;
- no more than $1~\mu\mathrm{m}$ maximum transverse position error at the
  target entrance or exit;
- no more than $1~\mu\mathrm{rad}$ maximum transverse slope error.

The scan contains 651 exact states per target, or 49,476 in total. Of those,
49,362 converge. The 114 failures occur only on large outer rays and are
counted as validity failures. All maintained-protocol states converge.

### Conservative fourth-order all-ring limits

The table reports the **minimum last sampled value that passes across all 76
sextupoles**. The next column is the next sampled value that fails. These are
discrete brackets, not continuously located physical limits.

| Varied input family | All-ring last pass | Next sampled fail | Median last pass | Limiting target(s) |
|---|---:|---:|---:|---|
| x bump at $\lvert\Delta K_2\rvert=0.10~\mathrm{m}^{-3}$ | $3.0~\mathrm{mm}$ | $4.0~\mathrm{mm}$ | $10.0~\mathrm{mm}$ | `SEX_39W` |
| y bump at $\lvert\Delta K_2\rvert=0.10~\mathrm{m}^{-3}$ | $4.0~\mathrm{mm}$ | $5.0~\mathrm{mm}$ | $10.0~\mathrm{mm}$ | `SEX_18W` |
| $\lvert\Delta K_2\rvert$ at $\lvert x_{\mathrm{bump}}\rvert=1.5~\mathrm{mm}$ | $0.60~\mathrm{m}^{-3}$ | $0.80~\mathrm{m}^{-3}$ | $2.0~\mathrm{m}^{-3}$ | `SEX_18W` |
| $\lvert\Delta K_2\rvert$ at $\lvert y_{\mathrm{bump}}\rvert=1.5~\mathrm{mm}$ | $0.80~\mathrm{m}^{-3}$ | $1.00~\mathrm{m}^{-3}$ | $2.0~\mathrm{m}^{-3}$ | `SEX_18W`, `SEX_24E` |
| Joint scaling from $(1.5~\mathrm{mm},0.10~\mathrm{m}^{-3})$ | $2.5\times=(3.75~\mathrm{mm},0.25~\mathrm{m}^{-3})$ | $3.0\times=(4.5~\mathrm{mm},0.30~\mathrm{m}^{-3})$ | $4.5\times$ | `SEX_18W`, `SEX_39W` |

All 608 signed maintained-protocol states pass the fourth-order gate. On the
60-target direct-GTPSA subset, the worst maintained-protocol discrepancies
against exact SciBmad are $0.00753~\mu\mathrm{m}$ in position and
$0.00199~\mu\mathrm{rad}$ in slope. Sixteen targets require a full-rank
fourth-order fit to the inner exact SciBmad scan because direct high-order
GTPSA propagation is unavailable there.

Second order is not uniformly valid at the maintained excitation: it passes
only 586 of 608 signed states. Its conservative all-ring joint scale is 0.5,
equivalent to $0.75~\mathrm{mm}$ and
$0.05~\mathrm{m}^{-3}$. Fourth order is therefore required for the common
maintained protocol unless target-specific smaller second-order bounds are
used.

A guarded synthetic common envelope at scale 2.0--a
$3.0~\mathrm{mm}$ bump and $0.20~\mathrm{m}^{-3}$ along the maintained
joint ray--remains one sampled step inside the worst fourth-order boundary.
This is still only a Taylor-truncation and closed-orbit-convergence result.
Even at the maintained scale, exact full-ring orbit extrema reach
$5.67~\mathrm{mm}$ horizontally and $8.56~\mathrm{mm}$ vertically.
Consequently, none of these values is a CESR corrector, power-supply, aperture,
lifetime, interlock, settling, or operator-safety limit.

---

## 7. Mitigating the two dominant stochastic errors

### A. BPM white noise: fixed templates and covariance-matched GLS

The current sensitivity setting assigns independent
$5~\mu\mathrm{m}$ RMS white noise to every BPM plane and read. With $R$
repeats per signed state, $K_2$ span $\Delta K_2$, and bump span $2b$,
the variance of each gradient channel is

$$
\mathrm{Var}(G)=
\left[
\frac{2\sigma}{\sqrt{R}\,\Delta K_2\,(2b)}
\right]^2.
$$

White-noise suppression therefore combines:

- a larger model-validated signal product $\Delta K_2 b$;
- repeated acquisition, giving the expected $1/\sqrt{R}$ reduction;
- whitening and matched filtering of the full BPM response with $C^{-1}$.

For the present equal, independent BPM-noise model, GLS is numerically the same
as the OLS matched filter. The implementation can instead accept a nonuniform
diagonal or full measured BPM covariance.

### B. Scan-time drift: balanced order, references, and a state-space inverse

Drift is modeled as a Gaussian random walk along one fixed two-dimensional
local-orbit mode, with $10~\mu\mathrm{m}$ RMS endpoint change over the core
scan. The mitigation has three layers:

1. Acquire the four signed states for each plane in a time-balanced
   `+,-,-,+` contrast order. Common constant and linear additive drift have
   zero total contrast weight.
2. Insert fixed-bump `0,+,0,-,0` $K_2$ reference blocks every 256 signal
   cycles and at the final cycle, making the random walk observable during the
   scan.
3. Retain finite reference-calibration errors as nuisance states and apply a
   random-walk state-space functional filter. The equivalent covariance is
   also propagated exactly with a reverse-cumulative closed form, without
   storing the full raw BPM tensor.

Acquisition order or per-point averaging alone is insufficient. The earlier
estimator normalized each BPM channel using the noisy scan itself, while the
typical $K_2$-dependent orbit signal was only at the nanometer scale. Fixed
physics templates and explicit covariance are the key changes.

---

## 8. Results under different noise conditions

The current eight-state time-series benchmark covers 76 targets, four hidden
all-sextupole-offset machines per target, 3,072 reads per signal state, and 512
stochastic measurement draws. Quadrupole strength, roll, and misalignment are
disabled in this table to isolate white-noise and drift mitigation.

| Case | Treatment | 2D RMSE | Median | P90 | P99 |
|---|---|---:|---:|---:|---:|
| No stochastic noise | Fixed-template inverse | $12.761~\mu\mathrm{m}$ | $8.323~\mu\mathrm{m}$ | $20.627~\mu\mathrm{m}$ | $38.946~\mu\mathrm{m}$ |
| BPM white noise only | Covariance-matched full-BPM filter | $20.385~\mu\mathrm{m}$ | $15.846~\mu\mathrm{m}$ | $31.464~\mu\mathrm{m}$ | $49.146~\mu\mathrm{m}$ |
| Random-walk drift only | Periodic-reference state-space filter | $12.764~\mu\mathrm{m}$ | $8.270~\mu\mathrm{m}$ | $20.729~\mu\mathrm{m}$ | $39.448~\mu\mathrm{m}$ |
| White noise + drift | Balanced eight-state sequence, no drift inverse | $21.108~\mu\mathrm{m}$ | $16.565~\mu\mathrm{m}$ | $32.503~\mu\mathrm{m}$ | $50.502~\mu\mathrm{m}$ |
| White noise + drift | **Periodic-reference state-space inverse** | **$20.297~\mu\mathrm{m}$** | **$15.817~\mu\mathrm{m}$** | **$31.289~\mu\mathrm{m}$** | **$48.681~\mu\mathrm{m}$** |

The stochastic drift contribution falls from $5.757$ to
$0.319~\mu\mathrm{m}$. The filtered result has a worst-target RMSE of
$44.561~\mu\mathrm{m}$ and passes the statistical gate requiring aggregate
RMSE, P99, and every target-level RMSE to remain below
$50~\mu\mathrm{m}$. Individual Monte Carlo draws can still exceed the
threshold because Gaussian noise is unbounded. The slightly smaller combined
RMSE than white-noise-only RMSE is finite Monte Carlo variation, not evidence
that adding drift improves the estimator.

The clean RMSE remains $12.761~\mu\mathrm{m}$, so finite-amplitude and
fixed-template model error remain even without stochastic measurement errors.

### Compound static nuisances

Simultaneously enabling 1% RMS BPM gain, 1% RMS corrector gain, 1% RMS
$K_2$ gain, independent $\pm1\%$ quadrupole-strength errors,
$1~\mathrm{mrad}$ RMS quadrupole roll, white noise, and filtered drift--while
still **excluding quadrupole misalignment**--gives:

- aggregate 2D RMSE: $30.334~\mu\mathrm{m}$;
- P99: $73.738~\mu\mathrm{m}$;
- worst-target RMSE: $66.560~\mu\mathrm{m}$;
- targets failing the current target-level gate: `SEX_09AW` and
  `SEX_38E`.

This compound case fails the current acceptance gates. Quadrupole strength is
the dominant static component, with a paired error-vector increment of
$19.384~\mu\mathrm{m}$ RMS. The combined white-noise and filtered-drift
stochastic component remains approximately $15.817~\mu\mathrm{m}$, showing
that stochastic mitigation itself does not become unstable. Quadrupole-
strength conditioning or machine-specific response relinearization is still
required.

### Why does another clean result report $5.864~\mu\mathrm{m}$?

The $5.864~\mu\mathrm{m}$ result uses the earlier
$0.5~\mathrm{mm}$ bump, $\Delta K_2=\pm0.02~\mathrm{m}^{-3}$, two-sided
BPM local-orbit estimate, and shared thin-source fit. The
$12.761~\mu\mathrm{m}$ result uses the present
$1.5~\mathrm{mm}/\pm0.10~\mathrm{m}^{-3}$ fixed-derivative template
benchmark. The excitation, hidden-machine ensemble, and response model differ.
The two values answer different validation questions and are not an
algorithm ranking on one common dataset.

---

## 9. Largest unresolved issue: quadrupole misalignment

A transversely misaligned quadrupole produces a dipole feed-down kick and can
create a large closed orbit. The primary problem is not merely a small
transport-matrix mismatch: the orbit can move the beam-relative sextupole
center beyond the available bump and inverse search ranges.

In the earlier one-at-a-time, uncorrected-orbit stress test, all 113
quadrupoles received independent $100~\mu\mathrm{m}$ RMS x/y offsets:

| Case | 2D center RMSE | Median | P90 | Maximum |
|---|---:|---:|---:|---:|
| No added nuisance | $6.051~\mu\mathrm{m}$ | $4.056~\mu\mathrm{m}$ | $8.945~\mu\mathrm{m}$ | $24.161~\mu\mathrm{m}$ |
| Quadrupole misalignment | **$1181.981~\mu\mathrm{m}$** | $224.169~\mu\mathrm{m}$ | $1917.196~\mu\mathrm{m}$ | $6204.899~\mu\mathrm{m}$ |

Of the 304 paired inversions:

- 290 beam-relative truths lie outside the $0.5~\mathrm{mm}$ bump radius;
- 176 have at least one plane outside the current
  $\pm1.5~\mathrm{mm}$ fit box.

The $1181.981~\mu\mathrm{m}$ RMSE is therefore an **uncorrected-orbit and
dynamic-range failure**, not a pure measurement of transfer-matrix sensitivity
to quadrupole offset. The current compound benchmark deliberately sets
quadrupole misalignment to zero, so this problem has not been solved or bypassed.

### Required next experiment

1. **Correct the orbit before the sextupole scan.** Quadrupole offsets must
   remain hidden from both correction and center inference. Correction may use
   only machine-available BPM observations, corrector actions, and a measured
   or model response matrix.
2. **Choose a defensible correction reference.** Candidates include BPM
   electrical zero, an operational golden orbit, or a beam-based-alignment-
   anchored reference. The simulator's latent ideal mechanical-center orbit
   must not be used as an unlabeled oracle.
3. **Restore dynamic range before attributing residual error.** Once the orbit
   is returned to the bump and fit ranges, residual response-matrix mismatch,
   local-orbit error, and center-inversion bias can be separated.
4. **Run one integrated end-to-end validation.** Combine two-sided BPM local-
   orbit reconstruction, full BPM covariance, drift states, and static lattice
   nuisances in a test that includes quadrupole misalignment.

Orbit correction restores a usable and repeatable trajectory. It does not by
itself identify quadrupole magnetic centers or provide an absolute mechanical
reference for the sextupoles.

---

## 10. Interpretation boundary

### Supported conclusions

- The $K_2$-odd and bump-odd intervention protocol makes the two-dimensional
  beam-relative sextupole center identifiable in the nominal response model.
- Nearest upstream/downstream BPM transport reduces local-bump prediction
  error to $0.301~\mu\mathrm{m}$ 2D RMSE in the current noise-free ensemble.
- The noise-free all-SciBmad-forward observation-derivative benchmark reaches
  $5.128~\mu\mathrm{m}$ 2D RMSE, but this is structural validation without
  measurement, calibration, drift, or quadrupole errors.
- The maintained $1.5~\mathrm{mm}/0.10~\mathrm{m}^{-3}$ excitation is
  inside the fourth-order all-ring model-validity envelope; the discrete
  envelope is not a machine-operating limit.
- Fixed templates, covariance-matched GLS, and periodic-reference state-space
  filtering control the assumed $5~\mu\mathrm{m}$-per-read white noise and
  $10~\mu\mathrm{m}$ endpoint drift at approximately
  $20.3~\mu\mathrm{m}$ combined RMSE.
- Quadrupole strength is the largest static response-model error in the
  current compound test. Quadrupole misalignment is the largest excluded and
  unresolved dynamic-range problem.

### Conclusions not yet supported

- These results are not demonstrated CESR machine precision. Noise amplitudes,
  drift mode, and cadence are synthetic sensitivity settings rather than
  measured priors.
- The non-misalignment nuisances are not collectively solved; the compound
  case fails in P99 and worst-target performance.
- Orbit correction is not equivalent to locating an absolute mechanical
  center.
- The current excitation has not been approved against power-supply, aperture,
  lifetime, interlock, or operator constraints.
- Four hidden machines per target do not establish the real machine
  distribution, and the present benchmark does not yet cover missing or
  outlier BPMs, multidirectional colored drift, hysteresis, or polarity
  asymmetry.

---

## Results and implementation

- [Finite-BPM local-orbit and center inversion](finite_bpm_inversion/README.md)
- [Noise-free all-SciBmad BPM/Taylor-map benchmark](sextupole_misalignment_only_bpm_taylor_map/README.md)
- [Fixed GTPSA derivative templates and stochastic inverse](gtpsa_derivative_stochastic_inverse/README.md)
- [Eight-state time-series result](gtpsa_derivative_stochastic_inverse/results/time_series_analysis/SUMMARY.md)
- [Compound nuisance result excluding quadrupole misalignment](gtpsa_derivative_stochastic_inverse/results/compound_nuisance_analysis/SUMMARY.md)
- [Quadrupole-misalignment one-at-a-time stress test](real_machine_nuisance_ablation/README.md)
- [SciBmad model-validity envelope for the excitation](sextupole_excitation_validity_envelope/README.md)
