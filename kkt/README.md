# KKT stationarity decomposition for `P_g`

Tier 1 of the dual analysis: split `∂L/∂P_g` into its constituent dual terms for
every model path (Kron / FULL_BUS / DQ_4TH × AVR × governor × GFM) and check that
the parts add up. Julia computes the coefficients, Python groups and plots.

## The rule

Every row of the model that touches `P_g` reduces to the same statement. With each
constraint stored in JuMP canonical form (`expr == 0` or `expr ≤ 0`):

```
contribution(family) = −dual × ∂expr/∂P_g

Σ contribution + dC/dP_g = residual ≈ 0
```

The cost gradient is `c1 + 2·c2·P_MW` €/MW for `cost_type = "quadratic"` and `c1`
for `"linear"` (`src/_common/functions_4_obj_func.jl:8-67`). No rescaling is needed:
the objective carries a `1/base_MVA` factor and the constraints are stated in p.u.,
and the two cancel, so every term below already lands in €/MW.

## Every row containing `P_g`

Verified across all four builders. A family only contributes when its condition
holds; `Save_KKT_Pg_Terms!` skips absent containers and warns about unknown ones.

| family | container key | `∂expr/∂P_g` | condition |
|---|---|---|---|
| `p_balance` | `opf_dict[:eq_const][:eq_const_p_balance][bus(g)]` | `+1` | always |
| `pg_lower` / `pg_upper` | `:ineq_const_pg_lower` / `_upper` | `−1` / `+1` | `bound_P_g`, CONSTRAINT encoding |
| `sg_upper` | `:ineq_const_sg_upper` | `2·P_g*` | `ineq_sg_upper` |
| `gfm_Imax` | `:ineq_const_gfm_Imax` | `2·P_g*` | `allow_gfm` |
| `gfm_Emax` | `:ineq_const_gfm_Emax` | `2·P_g*·(Xl/V*)²` | `allow_gfm` |
| `Pe_init` | `dyn_model_dict[:eq_const][:eq_const_P_init][g]` | `+1` | all TS paths |
| `Pm_init` | `:eq_const_Pm_init[g]` | `−1` | `mech_power_mode = USE_PM` |
| `gfm_Pm_init` | `:eq_const_gfm_Pm_init[g]` | `−1` | `allow_gfm` |
| `gfm_Pmeas_init` | `:eq_const_gfm_Pmeas_init[g]` | `−1` | `allow_gfm` |
| `swing_first_step` | `:eq_const_Δω_tf[g][1]` | `+Δt/(4H_g)` | FULL_BUS **and** `ode_first_step = :trapezoidal`; zero under `:backward_euler` |

Source anchors for the `Pe_init` row: Kron `functions_4_TS_kron_eqconst.jl:144`,
Kron-linear `:358`, dq `functions_4_TS_dq_eqconst.jl:66`, GFM
`functions_4_TS_gfm.jl:269`. `Pm_init` is `functions_4_TS_kron_eqconst.jl:20`; the
GFM twins are `functions_4_TS_gfm.jl:274` and `:283`. The trapezoidal first step is
`functions_4_TS_fullbus_eqconst.jl:265` and `:394`.

Coefficients are read back from the model with `JuMP.normalized_coefficient`
wherever the row is affine, so a builder rewrite that renormalises a constraint
cannot silently flip a sign. Nonlinear rows (`sg_upper`, the GFM limiters) have no
such method and fall back to the analytic entry above — which the residual
assertion then polices.

## Trap: `Pm_init` is a capacity rent, not dynamics

`eq_const_Pm_init` is `P_m = P_g`. Stationarity in `P_m` (contribution
`= −dual × ∂expr/∂P_m`) gives

```
−ν − η_u + η_l = 0     ⇒     ν = −η_u + η_l
```

so ν is the rent of whatever box sits on `P_m` — with `bound_P_m = true` and
`P_m_source = :dgen_pg_limits`, that is the **same DGEN Pmax** the dispatch already
imposes on `P_g` through `bound_P_g`. The unit's capacity rent therefore enters the
`P_g` row twice: once as `η⁺`, once through the pin.

The feasible set is unchanged and λ is correct. What is not determined is the
*attribution*: any split of the rent between `η⁺` and ν satisfies stationarity, and
an interior-point solver typically splits it evenly, so the two columns can come out
numerically identical. **Only `scarcity + mech_pin` is meaningful.**

Verified on a 39-bus DQ + AVR + governor case: `max |ν − (P_m box rent)| = 3.8e-7` across all ten machines,
i.e. the pin carries the box and nothing else, governor notwithstanding. In that run
the δ corridor never binds (`max |dual| ≈ 1e-11`), so the genuine dynamic column is
zero — reporting the pin as "dynamic" would have claimed a 12 €/MW stability rent
that does not exist.

For a case study needing a clean attribution, set `bound_P_m = false` or widen the
`P_m` box beyond the dispatch limits, and re-solve.

## Where the corridor duals are

In this repository δ and Δω are **explicit variables** with their own defining
equalities, so the corridor duals Γ (`dual_delta_COI_*`, `dual_delta_ref_*`,
`dual_Delta_Omega_*`) do **not** appear in the `P_g` row at all. That differs from
the substituted formulation used in the older linear-test notebook, where `h_P`
injected `c₁Γ` and `c₃Γ` straight into it. Here their influence reaches `P_g`
through the adjoint chain, into `ν` (the `Pe_init` dual) and `φ¹` (the first swing
step).

Tier 1 therefore gives the exact split of λ into **cost / own-limit scarcity /
dynamic coupling**. Attributing the dynamic part back to individual corridor duals
needs the implicit-function-theorem sensitivity `S = ∂h/∂P` and the `F_x` assembly —
that is Tier 2, and it is out of scope here.

## Tier 2 — `S = ∂h/∂P` and the causer-pays identity

With `h_g^t := δ_g^t − δ_COI^t` (what the corridor bounds) and `F(x, P) = 0` the
dynamic equality block, the implicit function theorem gives

```
dx/dP = −F_x⁻¹ F_P        S_P = h_P + h_x·(dx/dP)
```

and stationarity in the state, `F_xᵀα + h_xᵀγ̃ = 0` with `γ̃ = Γ⁺ − Γ⁻`, collapses the
weighted sensitivity to

```
γ̃ᵀS_P = γ̃ᵀh_P + αᵀF_P = −dyn
```

**`h_P = 0` here.** δ is an explicit variable, so the corridor row contains no `P_g`
at all — unlike the substituted formulation of the reference notebook, where `h_P`
carried `c₁`/`c₃`. The identity is therefore just `γ̃ᵀS_P = αᵀF_P`, i.e. minus the sum
of the `Pe_init`, `Pm_init` and `swing_first_step` contributions in
`kkt_Pg_terms.csv`. Tier 2 is the numerical confirmation of that.

### Every derivative comes from the model

`Save_KKT_Sensitivity!` never differentiates anything by hand. `F_x`, `F_P`, `h_x`,
`h_P` are slices of the **exact** constraint Jacobian at the solution, obtained from
`build_constraint_evaluator` (MOI automatic differentiation of the actual JuMP
model). Consequences: no calculus to re-derive when a builder changes, and extending
to FULL_BUS / DQ / GFM means adding a `DynStateSpec` — a list of family names — not
new maths.

The work is therefore *indexing*: `jacobian_rows.csv` and `jacobian_cols.csv` label
every Jacobian row and column with the model's own family names and stamp a `block`
column (`F` / `h` / `h_lower` / `other`, and `x` / `P` / `Q` / `other`), so the
notebook slices by following the model rather than by re-declaring the partition.
Anything the labeller cannot place is marked `UNMAPPED` and counted in the meta file
— an unmapped `F` row would silently corrupt `S`.

### Squareness

`F_x` must be square over the state. On the Kron classical path both sides count
`3n + 3nT + T`. The exporter asserts it and prints the per-family counts on failure
rather than producing a plausible wrong answer.

### Boxes on state variables

The clean identity needs the corridor to be the **only** inequality touching the
state. `bound_E`, `bound_δ` and `bound_P_m` each add one, and their duals then appear
in that state's column of the adjoint residual. The exporter reports the residual per
state family, so such a run is explainable rather than merely wrong; the proof run
turns those boxes off.

Related: this is the same mechanism as the Tier-1 `mech_pin` finding above — the
`P_m` box reaching the `P_g` row through `Pm_init`.

### The joint identity residual is a solver-tolerance quantity, not an error

`identity_residual_joint` is **not** the number to read. Substituting the three
sensitivity definitions — `S_P = h_P + h_x·dxdP`, `S_ω = hω_P + hω_x·dxdP`,
`S_b,P = b_P + b_x·dxdP`, all built from the single `dxdP` — into the joint identity and
using `dxdP = −F_x⁻¹F_P` collapses it exactly to

```
joint − (γ̃ᵀh_P + γ̃ωᵀhω_P + αᵀF_P + Σμ_bᵀb_P)  ≡  rᵀ·dxdP,
    r = F_xᵀα + h_xᵀγ̃ + hω_xᵀγ̃ω + Σb_xᵀμ_b
```

`r` is the solver's own stationarity residual on the state columns — `adjoint_residual_full_max`
reproduces Ipopt's reported *Dual infeasibility* digit for digit. So the joint residual is
Ipopt's dual tolerance amplified by `‖dx/dP‖`, and it **grows with the state count**:
~2e-11 at `n_x ≈ 9e3`, and 1.8e-5 … 3.1e-3 across the FULL_BUS DQ runs at `n_x` 23k–32k
over a 5 s horizon. Bounding it by a fixed tolerance fails on large *correct* runs.

The exporter therefore splits it:

| key | meaning |
|---|---|
| `identity_residual_explained` | `rᵀ·dxdP` with `r` assembled from the **full Jacobian against every row's dual**, with no reference to the F / h / block partition |
| `identity_residual_unexplained` | `joint − explained` = the charge of rows the partition missed. **This is the number that must be ~0.** |
| `stationarity_residual_max` | `max|r|`; cross-check against `solver_log.txt` |

The independence of `r_true` is what gives the split teeth. Measured on the five
FULL_BUS DQ runs of 2026-08-15 (plain, +AVR, +governor, +both, Δω-corridor):
`explained` reproduces `joint` to 8–10 significant digits and `unexplained` sits at
1e-13 … 1e-11 throughout. Negative control on the 179-bus `P_m`-box run: dropping
`ineq_const_P_m_lower` from block discovery leaves `explained` unchanged at 5.0e-11 and
moves the entire 44.3786 into `unexplained`.

Per generator, `charges_by_channel.csv` carries the same split as `residual_explained`
and `residual_unexplained`.

### The corridor must actually bind — and on 9-bus Kron it never does

If `γ̃ ≈ 0` the identity is `0 = 0` and proves nothing.

**Measured on this path, and it is a property of the case rather than of the code:**
the δ–COI corridor on `9bus` Kron is not merely slack, it is *unpriceable*. It goes
from slack straight to `LOCALLY_INFEASIBLE` with no binding window between — slack at
56°, infeasible at 52° (0.3 s clearing); γ̃ ≤ 4e-6 up to `load_factor` 1.8, infeasible
at 2.0.

The mechanism: at `δ_tol = 32°` against an unconstrained peak of 40.85°, the
trajectory *does* move (peak → 27.3°) at **Δobj = 0**. The swing is reshaped through
**Q and V, which cost nothing in the objective**, so the constraint is met for free
and never acquires a multiplier.

Consequence for Tier 2: on this case every check below passes, but the identity is
verified at ~1e-7 €/MW — real, and not discriminating. The finite-difference check is
unaffected, since it validates `S` without reference to any dual. A discriminating
run needs a case where the corridor is priced; the known-binding lever in this
repository is **GFM `Imax` on `9bus_gfm`** (≈157 €/pu), which needs a `DynStateSpec`
for the DQ FULL_BUS path — a list of family names, not new calculus. The notebook
prints `*** NOT A DISCRIMINATING TEST ***` rather than passing quietly on noise.

## Files

The title of this file says `P_g`, but the directory now holds three chains: Tier 1 on the
generator side, Tier 2 (`S = ∂h/∂P`, documented above), and the demand side.

**Tier 1 — `∂L/∂P_g`**

| file | role |
|---|---|
| `kkt_pg_terms.jl` | `Save_KKT_Pg_Terms!` — writes `Dispatch/CSV_duals/kkt_Pg_terms.csv` |
| `runtests_kkt_terms.jl` | residual assertion, run by hand |
| `KKT_Pg_Duals.ipynb` | reads `run_manifest.toml` + the CSVs, plots the decomposition |

**Tier 2 — `S = ∂h/∂P` and the causer-pays identity**

| file | role |
|---|---|
| `kkt_sensitivity.jl` | `Save_KKT_Sensitivity!` — the `DynStateSpec` registry and `Sensitivity/` |
| `S_Sensitivity_Pg.ipynb` | channel charges vs `−dyn`, the explained / unexplained residual split |
| `runtests_kkt_sensitivity_blocks.jl` | replays N exported runs from their Jacobian CSVs, no solve |

**Tier 1-D — `∂L/∂P_d`, the demand side**

| file | role |
|---|---|
| `param_jacobian.jl` | `Save_Param_Residuals!` + `param_column`: `∂c/∂P_d` by model rebuild, since demand is a parameter and MOI cannot differentiate it. A *validation* tool — the notebooks do not need it |
| `pd_coeffs.py` | `coeff_arrays` / `f_pd_matrix` — the single implementation of `∂expr/∂P_d`, shared by Tier 1-D and Tier 2-D so the ZIP and Schur formulas cannot drift apart |
| `kkt_pd_terms.jl` | `Save_KKT_Pd_Terms!` — the energy / stability ledger, and the `\|energy − λ\|` reduction check that pins the sign convention |
| `KKT_Pd_Duals.ipynb` | the report: per-bus `π_k` split, from one run's existing exports, no solve |

**Tier 2-D — `S_Pd` and the demand-side causer-pays identity**

| file | role |
|---|---|
| `S_Sensitivity_Pd.ipynb` | splits each bus's stability price across the machines whose corridor its load moves: `γ̃ᵀS_Pd + Σ_b μ_bᵀS_b,Pd = αᵀF_Pd = −stability_k` |

Tier 2-D needs no `DynStateSpec` and no new calculus. Every row carrying `P_d` is already
inside the `F` block (`eq_const_Pbalance_tf/_tpf` on FULL_BUS, `eq_const_Pe_tf/_tpf` on Kron),
and `h_{P_d} = b_{P_d} = 0` because corridors bound δ/Δω and the state boxes bound `E`/`δ`/`P_m`
— none of which contain load. The notebook asserts that rather than assuming it. `F_x` and its
LU are the ones Tier 2 already builds, so a bus costs one back-substitution.

Each side has two notebooks: `KKT_Pg_Duals` then `S_Sensitivity_Pg`, and `KKT_Pd_Duals` then
`S_Sensitivity_Pd`. The second one of each pair needs the `Sensitivity/` folder of the run, which
exists only when the driver attached the Tier-2 hook.

## Output schema

`kkt_Pg_terms.csv`, `;`-delimited, one row per (generator, family, time):

```
gen_id;bus_id;family;store;t_index;dual;coeff;contribution
```

`store ∈ {opf, dyn, objective, check}` travels with every row so a consumer never has
to re-derive which families are dynamic coupling and which are steady-state scarcity —
`KKT_Settlement.ipynb` discovers the dynamic family list from it rather than hardcoding
names, which is how it stays correct when a family is added.

`t_index` is empty for time-independent families. Each generator also gets one
`cost_gradient` row (`store = objective`, `dual` and `coeff` empty,
`contribution = dC/dP_g`) and one `RESIDUAL` row (`store = check`) holding the sum of
everything above it.

## Settlement — putting the two sides against the actual dispatch

Tier 1 and Tier 1-D each produce a price. The settlement report runs them against the dispatch
and reports the money: `KKT_Settlement.ipynb`.

A generator maximising its own profit at a posted price knows only its own limits — it has
never heard of a rotor-angle corridor:

```
max  pi * P_g - C_g(P_g)   s.t.  P_min <= P_g <= P_max
  =>  pi = C'_g + mu_upper - mu_lower
```

Subtract that from the market's stationarity and the residual is exactly the dynamic term:

```
lambda_g - C'_g - mu_upper + mu_lower  =  dyn_g
```

so **the uniform nodal price does not support the dispatch** whenever a corridor reaches the
active-power columns. The price that restores support is

```
pi_hat_g = lambda_g - dyn_g = lambda_g + charge_g
```

with `charge_g` the Tier 2 causer charge, which the notebook cross-checks against
`Sensitivity/charges_by_channel.csv` (agreement 8e-11 EUR/MW on the 9-bus Kron run).

The sharp version: in a classical ACOPF an interior generator always satisfies
`lambda_g = C'_g`. In a TSC-OPF it need not — `mu = 0` and yet the gap is `dyn_g`, a corridor
repricing a unit nowhere near its own capacity.

### What it reports

Four pricing rules, active power only, as a 2x2 of `R = what load pays - what generation is
paid` (`+` TSO surplus, `-` TSO tops up), decomposed into

```
R(pi_d, pi_hat) = R(lambda, lambda) + sum_k stability_k * P_d,k - sum_g charge_g * P_g
```

plus per-generator profit and any make-whole owed. Measured on `Results - 2026-08-16 130202`
(9-bus, Kron, classical, `:coi_box` delta, SC contingency 2, `max|gamma|` 4.50), EUR/h:

|  | gen @ lambda | gen @ pi_hat |
|---|---|---|
| **load @ lambda** | -534.82 | +6595.93 |
| **load @ pi_d** | -81.83 | +7048.93 |

`= -534.82 + 452.998 + 7130.756`. Self-test: `cost_type = "linear"` with every unit interior
makes `pi_hat = C'` exactly, so `sum pi_hat * P_g` must reproduce `[status] objective_MVA` —
13286.0360 both ways.

### Three things it is there to make visible

- **`pi_hat` as a linear tariff is confiscatory.** It strips the whole inframarginal rent
  (7130.76 EUR/h here) and leaves every unit at exactly zero profit. Restoring incentive
  compatibility with a linear price and preserving that rent are the same instrument pulled in
  opposite directions; a lump-sum side payment fixes cost recovery without touching the
  marginal condition. Report the four cells as a menu, not a ranking.
- **A unit can lose money strictly inside its own limits.** G1 runs 161.9 MW at `C' = 50` while
  `lambda_1 = 42.29`, losing 1248.39 EUR/h — the "penalised at `P_min`" case with the corridor,
  not `P_min`, doing the holding. A stability-blind settlement has no vocabulary for it.
- **The residual is not all stability.** `R(lambda, lambda)` is already -534.82 from losses,
  congestion and voltage rents, plus an unsettled reactive side. Stability moves the residual;
  it does not create it.

Linear cost makes the support claim weak — `pi_hat = C'` leaves the unit indifferent across its
whole range. A quadratic curve pins `P_g*` uniquely and is the sharper case study.

It depends on `kkt_Pd_prices.csv`, which `KKT_Pd_Duals.ipynb` writes: run that notebook on the
same run first. Verified on `KRON_REDUCED` + `CLASSICAL_2ND` and on `FULL_BUS` +
`DQ_4TH` + AVR + governor with no edit — the dynamic family list is discovered from the
`store` column, never hardcoded.

## Tier 1-P — the demand pin, and why it needs no census

Tier 1-D obtains `∂c_r/∂P_d` by rebuilding the model at `p_d ± ε` and differencing
(`param_jacobian.jl`). That exists because demand was never a JuMP object: `p_d` was an
inlined `Float64`, so `build_constraint_evaluator` had no column to differentiate.

`DispatchConfig.pd_as_variable = true` removes the obstacle instead of working around it.
The ACOPF then builds

```
P_d[k]  a variable,      eq_const_pd_pin[k] :  P_d[k] - p_d⁰[k] == 0   : σ_k
```

and stationarity of the Lagrangian in `P_d[k]` reads

```
∂f/∂P_d,k - Σ_r λ_r ∂h_r/∂P_d,k - Σ_s μ_s ∂g_s/∂P_d,k - σ_k = 0     ⇒     σ_k = π_k
```

so **the pin's own dual is the price**, summed by the solver over whatever rows happen to
contain the demand. That is the property the FD route had to reconstruct by census: here it
is structural, and nothing can go stale when a builder changes.

| | Tier 1-D, rebuild | Tier 1-P, pin |
|---|---|---|
| generic over `h`, `g` | yes, by census | yes, by construction |
| accuracy | `O(ε²)` truncation + dual convergence | solver tolerance |
| cost | 3 solves **per bus** | 1 solve, **all** buses |
| traps | `ε` tuning, the `V²` row-rescale branch, unconverged legs | none |
| Kron TSC | works | **cannot** — `Y_red` is assembled outside JuMP |

The pin does not replace `param_jacobian.jl`. Where the demand is absorbed into a numeric
matrix before JuMP sees it — the Kron constant-impedance shunt propagating through
`Reduce_Matrix` into every entry of `Y_red` — there is no coefficient to attach a variable
to, and the rebuild remains the only route. Where the demand survives as an inlined
coefficient (the ACOPF balance today, the FULL_BUS ZIP right-hand side if extended) the pin
is exact and covers every bus at once.

### Sign, measured not assumed

`σ_k = +λ^P_k` on a dispatch-only ACOPF, where the balance and the pin are the only rows
containing the demand. Measured on both `use_matrix` forms and both fixtures:

| case | λ range (€/MW) | max \|π − σ\| |
|---|---|---|
| 9bus, `load_factor` 1.5 | 50.0 flat (lossless, uncongested, `c1 = 50` on all units) | 0.0 |
| 39bus, `load_factor` 1.0 | 6.72 … 35.96 | 7.1e-15 |

The 9bus fixture cannot discriminate a locational bug — every branch has `r = 0`, `rateA`
is 9999 and all three units bid 50, so λ is 50 everywhere by construction. 39bus has 42
resistive branches and a real merit order, and is the case that shows σ tracking λ bus by
bus. Note 39bus is `LOCALLY_INFEASIBLE` at the `RunConfig` default `load_factor = 1.5`,
with or without the pin.

### The ledger — `kkt_pd_pin.jl`

`dual_Pd_pin.csv` carries `π_k` on every run with the flag on and needs nothing else.
`Save_KKT_Pd_Pin!` adds the per-row breakdown: it reads the `P_d[k]` **column** of the exact
Jacobian, labels each row through `label_jacobian`, and writes
`Dispatch/CSV_duals/kkt_Pd_pin.csv` in the Tier 1-D schema, so `KKT_Pd_Duals.ipynb` reads it
unchanged. `label_jacobian` needed no edit — it walks the model's own containers, so storing
the pin under `opf_dict[:eq_const]` and the variable under `opf_dict[:vars]` was enough to
have both labelled rather than `UNMAPPED`.

**The pin row is excluded from the ledger, deliberately.** It is the instrument, not part of
the model being priced: folding `σ_k` back into the sum returns `π_k − σ_k`, i.e. the
stationarity condition, which is ≈ 0 by construction and says nothing about the price. It is
written out anyway as the `STATIONARITY` row, where it is a free check on the solver's dual
convergence. (This was a real bug in the first draft of the test, and it fails in the most
misleading way possible — a clean zero.)

Include `kkt_sensitivity.jl` first, as `param_jacobian.jl` requires, then attach:

```julia
hook(; path_names, opf_dict, dyn_model_dict, model, kwargs...) =
    Save_KKT_Pd_Pin!(path_names, opf_dict, dyn_model_dict, model)
```

## The ACOPF baseline — `ACOPF_Settlement.ipynb`

Every chain above prices a transient-stability run, and the settlement report's headline
is that `R` does **not** close. That reads as a finding only against a case where it
provably does. `ACOPF_Settlement.ipynb` is that case: the same accounting on
`trans_stab = false`, where the residual is an algebraic identity and any gap is an
export bug or a loose solver tolerance, never economics.

| file | role |
|---|---|
| `ACOPF_Settlement.ipynb` | the report: two identities, four CSVs, four figures |

The run it reads is an ACOPF with `pd_as_variable = true`, `ineq_sbranch_upper = true` and
`post_solve_hook = Save_KKT_Pg_Terms!`.

### The two identities, kept apart

The point of the notebook is that these are different statements and carry different duals.

**Network identity — the operator's money.** Multiply the balance row by `λ_k` and sum:

```
MS_P ≡ Σ_k λ_k (P_d,k − P_g,k)
     = − Σ_k λ_k P_sh,k − Σ_branches (λ_i P_ik + λ_k P_ki)
```

and each branch term splits exactly, two ways:

```
λ_i P_ik + λ_k P_ki = ½(λ_i+λ_k)(P_ik+P_ki) + ½(λ_i−λ_k)(P_ik−P_ki)     symmetric
                    = λ_k (P_ik+P_ki)       + (λ_i−λ_k) P_ik            directed
```

The symmetric split is orientation-independent; the directed one reproduces the textbook
congestion-rent expression but flips sign if the branch is entered from the other end.
Both are written to `acopf_settlement_branch.csv`; the symmetric one is what the closure
asserts. `μ_V`, `μ_θ` and `μ_Pg` appear nowhere in this identity — a voltage limit changes
the dispatch and therefore `λ`, but collects no rent of its own. Charging for one on top
of LMP settlement would collect the same money twice, the steady-state twin of the
`γ̃ᵀS_P` warning in `S_Sensitivity_Pg.ipynb`.

**Stationarity identity — the generator's money.** With the bound senses of
`bound_builders.jl:160,169` (`lb − P_g ≤ 0`, `P_g − ub ≤ 0`):

```
λ_bus(g) = C'_g(P_g) + μ_lower − μ_upper
```

so profit `(λ − C')·P_g` **is** the scarcity rent. Under a linear cost an interior unit
earns exactly zero, so a nonzero make-whole says some unit is pinned against the merit
order. That number is the baseline the TSC uplift is read against — the difference between
the two is what stability costs, not the TSC figure alone.

### What it asserts

Six checks, and a failed one fails the `nbconvert --execute` render. That is the test; no
`runtests_*.jl` accompanies this notebook.

1. `buses_report` satisfies its own balance row (`PG − PD − Psh = P_MW`), 1e-6 MW.
2. `λ_k = σ_k`. Trivially true here, and that is the point — proving the gap is zero
   without transient rows is what licenses reading `σ − λ` as stability content when they
   are present.
3. `MS_P` closes on the branch decomposition, 1e-6 €/h.
4. `MS_Q` closes likewise, at `μ^Q`.
5. Both branch splits reproduce the same total, 1e-9 €/h.
6. `kkt_Pg_terms.csv` `RESIDUAL` is zero and its `λ` and `C'` match the dual export and
   the manifest cost curve, 1e-8. Skipped with a message when the driver set no hook.

An independent cell rebuilds `Y` from `Inputs/line_data.csv` and recomputes the polar
injection, so a bad flow column cannot hide inside a self-consistent identity. It is
advisory, and it skips itself on a case with off-nominal taps or phase shifters, which
that build does not model.

### Run resolution

`ACOPF_RUN_DIR` names a run folder or a `RESULTS` root, newest wins — same shape as the
other notebooks. The notebook refuses a run whose manifest says `trans_stab = true`.

Nothing has to be rendered first, because there is no demand-side price to import.

### Writes

Into the run folder, all `;`-delimited: `acopf_settlement_bus.csv` (prices, payments,
`σ`, `μ_V` per bus), `acopf_settlement_branch.csv` (both rent splits, `μ_S`, loading,
binding flag), `acopf_settlement_gen.csv` (cost curve, `λ`, `C'`, `μ_Pg`, profit,
stationarity residual), `acopf_settlement_summary.csv` (the scalars), and
`ACOPF_Settlement.html`.
