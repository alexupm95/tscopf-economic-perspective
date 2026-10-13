#=
================================================================================
 numerics.jl  —  numerical constants of the dynamic models
================================================================================
 Every smoothing width and formulation threshold used by the transient-stability
 builders, in one place. None of these is configuration:

 - user-facing limits, tolerances and model options live in the config structs
   (`TsBoundLimitsConfig`, `DynModelConfig`, `TsBuilderConfig`, the solver configs)
   and are documented in `docs/parameter_reference.md`;
 - structural metadata (equality-family registries, dual-export specs, required
   CSV columns) stays next to the code that defines those structures;
 - zero guards (`isapprox(x, 0; atol = 1e-12)`, `Tf <= 1e-9`) stay inline.

 Each value below has one correct order of magnitude, explained in its docstring.
 Changing one changes the formulation. These are `const`: Julia compiles them into
 the builders, and Revise does not replace them, so **restart Julia after editing
 this file**. A run in a live session keeps the old value.

 Included from `src/TSCOPF.jl` before any `_transient_stability` builder.
================================================================================
=#

# ===================================================================================
# Governor (TGOV1)
# ===================================================================================

"""Smoothing of the GOV_SMOOTH sqrt min/max on the valve (keeps the gradient finite)."""
const _GOV_SMOOTH_RHO = 1.0e-4

# ===================================================================================
# AVR (exciter)
# ===================================================================================

"""Smoothing of the exciter sqrt min/max (matches `_GOV_SMOOTH_RHO`)."""
const _AVR_SMOOTH_RHO = 1.0e-4

# ===================================================================================
# Grid-forming converters (REGFM_A1)
# ===================================================================================

# --- Formulation switch ------------------------------------------------------------
"""
Current-limiter bypass threshold [pu, system base].

`Imax` at or above this value cannot bind (`scale ≡ 1`), so `_add_gfm_current_limiter!`
emits unsaturated stator algebra `Xl·Id = ΔVd`, `Xl·Iq = ΔVq` instead of the nonconvex
smooth-max clamp, dropping one sqrt expression per (gen, t).

This is a **formulation switch driven by input data**, not a physical limit: the model
changes shape when a `gfm_dynamic_data.csv` row crosses the threshold. Units that take
the bypass are recorded on `meta[:gfm_limiter_bypassed]` and logged once per window.
"""
const _GFM_IMAX_NO_LIMIT = 20.0

# --- Smoothing: Q–V PI and current limiter (fixed values) --------------------------
# Each has one correct order of magnitude and no study varies it.

"""Smoothing of the E-clip sqrt composition (anti-windup on the Q–V PI)."""
const _GFM_EPS_E = 1.0e-4

"""√ regularisation inside the limiter norms — keeps `∂|I|/∂I` finite at `I = 0`."""
const _GFM_EPS_NORM_V = 1.0e-8

"""Smooth-max scale of the current limiter, in current units (multiplied by `Xl` at use)."""
const _GFM_EPS_LIM_I = 1.0e-3

# --- Smoothing: P/Q overload PIs (resolved per converter and axis) -----------------
# The overload branch has two smoothed clips, and they get different ε because they
# fail in different ways (resolved by `_gfm_pq_eps` / `_gfm_pq_eps_out`, recorded in
# `run_manifest.toml` under `[resolved] gfm_pq_eps_*`).

"""
Overload **state** clip ε, as a fraction `c` of one integrator step.

The clipped state `x = smooth_max(x_raw, 0)` is carried from step to step. Each
branch gets ε = c · s, where s = Δt · Ki · E_ref is the step the integrator takes
for a 1 pu error on the converter's own base. Two errors trade against each other,
and both scale as ε²/s:

- **Bias.** Inside the band the smoothed clip never returns exactly 0. The clipped
  state then settles at x* ≈ ε²/(4s) instead of 0, and every in-band branch feeds a
  false Δω (P axis) or voltage offset (Q axis) into the droop. With a fixed
  ε = 1e-2 and Kipmax = 0.1, that was ≈ 0.1 pu Δω on the P axis, and an all-GFM
  fleet's frequency ramped for seconds.
- **Conditioning.** With the raw state one step below the rail, ∂x_clip/∂x_raw ≈
  ε²/(4s²). In principle a very small slope stops the clip from tying the raw column
  down. In practice the raw column is pinned by its own integrator equality and this
  never bit (see below).

With ε = c·s the bias is c²·s/4 for a 1 pu error. Near the rail the error |e| is
small and the bias grows as (s/|e|)·c²/4. That growth is what sets c: the false state
has to stay below the proportional term Kp·|e|, or the branch holds the unit at a
spurious equilibrium inside the band.

c = 0.01, measured on the all-GFM 9-bus load trip (Δt = 5 ms, REGFM_A1 gains):
- At c = 0.2, a converter with Qmin = −5 MVAr settled at −0.63 MVAr, 4.4 MVAr inside
  its rail, and a unit pinned at Pmin = 0 sat at 0.5 MW.
- At c = 0.01 they settle at −4.95 MVAr and 0.001 MW, with no more Ipopt iterations
  (25 and 31 against 26 and 32).
- At c = 1 the P and Q corners were both visibly soft (10 MW, 24 MVAr).

A single shared ε cannot serve both axes: REGFM_A1 gains put Kiq about 200× above Kip.
"""
const _GFM_PQLIM_EPS_STEP_FRAC = 0.01

"""
Overload **output** clip width, as an error tolerance [pu, converter base].

The branch output `u = smooth_max(Kp·e + x, 0)` carries nothing to the next step, so
it has neither the drift nor the raw-column conditioning problem of the state clip.
Its only cost is a soft corner: the √ smoothing starts acting about ε/Kp before the
rail, measured in error units. ε_out = Kp · tol keeps that zone at `tol` of the
converter's own rating on both axes. The value 1e-3 is 0.3 MVAr on a 300 MVA unit.
With the state-clip ε used here too, the Q axis (Kpq = 3) acted about 2 MVAr early at
c = 0.2 and about 10 MVAr early at c = 1 (state-clip ε, before it moved to c = 0.01).
"""
const _GFM_PQLIM_OUT_ERR_TOL = 1.0e-3

"""Floor on every overload ε, for branches with a zero gain (no integrator / no P term)."""
const _GFM_PQLIM_EPS_FLOOR = 1.0e-6
