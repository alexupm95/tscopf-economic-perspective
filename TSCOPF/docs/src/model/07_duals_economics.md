# 7. Duals, KKT, and the economics

Dual variables are the language TSCOPF uses to connect optimisation algebra to market intuition. A multiplier on a power-balance row is an energy price. A multiplier on a binding generator limit is a capacity rent. A multiplier on a rotor-angle corridor is a stability scarcity signal. This page states the sign convention the code enforces, the stationarity conditions that define economic prices, and how exported CSV/XLSX columns map back to constraint families.

The coupled model is on [6. The TSC-OPF, assembled](06_tsc_opf_assembled.md). Runnable export flags and folder layout are in [Running a case — outputs](../running_a_case.md).

---

## JuMP Lagrangian and inequality sense

TSCOPF writes every constraint in **≤ form** for inequalities and **= 0** for equalities. JuMP's Lagrangian (simplified) is:

```math
\mathcal{L}(x, \lambda, \mu) = f(x) - \sum_i \lambda_i\, h_i(x) - \sum_j \mu_j\, g_j(x),
\qquad g_j(x) \leq 0 .
```

!!! warning "Sign convention — read this before opening a CSV"
    For an **active upper-bound** inequality coded as $(\mathrm{LHS} - \mathrm{RHS}) \leq 0$, `JuMP.dual()` returns a value **$\leq 0$**. A negative dual on $P_g \leq P_g^{\max}$ at a binding limit is correct, not a bug.

    **Balance rows are the trap, and there are two different $\lambda$.** On the **primal JuMP export** the nodal LMP is $\pi_k = +\lambda_k$: the balance is coded $P_g - P_d - \sum(\cdot) = 0$ and the Lagrangian above carries $-\lambda$, so the two minus signs cancel. On the **explicit dual LP** (`Dispatch_Dual/`) the multiplier is the negative of that, and there $\pi_k = -\lambda_k$. Both routes give the same price from opposite-signed multipliers. Cross-check on case9, contingency 2, before plotting rents.

---

## Steady-state stationarity and LMPs

Consider a linear DCOPF (the cleanest LP case). Primal:

```math
\min \sum_g c_g P_g
\quad\text{s.t.}\quad
P_{g(k)} - \sum_m B_{km}(\theta_k - \theta_m) = P_d(k),
\quad P_g^{\min} \leq P_g \leq P_g^{\max},
\quad \ldots
```

Let $\lambda_k$ be the multiplier on balance at bus $k$, and $\underline{\eta}_g, \overline{\eta}_g$ on generator limits. Stationarity w.r.t. $P_g$ at bus $k(g)$:

!!! note "Model 7.1 (dispatch stationarity)"
    ```math
    \begin{equation}
    \label{eq:stationarity-pg}
    c_g - \lambda_{k(g)} - \underline{\eta}_g + \overline{\eta}_g = 0
    \end{equation}
    ```

!!! note "Interpretation (LMP decomposition)"
    Define the **economic energy price** $\pi_{k(g)} = +\lambda_{k(g)}$. Rearranging $\eqref{eq:stationarity-pg}$:

    ```math
    \pi_{k(g)} = c_g - \underline{\eta}_g + \overline{\eta}_g .
    ```

    At bus $k$, $\pi_k$ is the marginal cost of serving one more MW of load: fuel cost $c_g$ at the marginal unit, plus any capacity rent $\overline{\eta}_g$ if the upper bound binds, minus any must-run subsidy $\underline{\eta}_g$ if the lower bound binds. Congestion and angle limits enter through $\lambda_k$ indirectly when they force a re-dispatch.

ACOPF uses the same $\pi_k = +\lambda_k$ rule on the active balance rows from [2. Steady-state OPF](02_steady_state_opf.md). Reactive balance multipliers $\beta_k$ play the analogous role for reactive support scarcity.

### Two ways to obtain steady-state duals

| Mechanism | Flag | Output | When |
|---|---|---|---|
| Primal JuMP duals | `RunConfig.save_duals = true` | `Dispatch/CSV_duals/`, `Dispatch_Duals.xlsx` | Default; all dispatch types |
| Explicit dual LP | `dispatch.solve_explicit_dual = true` | `Dispatch_Dual/` | Linear ED or DCOPF only; strong-duality check |

The explicit LP treats prices as decision variables and recovers the same **price** when strong duality holds, but from a multiplier of the opposite sign: it codes stationarity as $c_g + \lambda_k = 0$ (`eq_const_dcopf_dual_stationarity_Pg!` in `src/_dcopf/functions_build_DC_OPF_Dual.jl`), so its $\lambda_k = -c_g$ and there $\pi_k = -\lambda_k$. In short, $\lambda^{\text{explicit}}_k = -\lambda^{\text{primal}}_k$ and both give $\pi_k = c_g$ at the marginal unit.

---

## Transient duals are local KKT multipliers

When `trans_stab = true` and `save_duals = true`, TSCOPF exports JuMP duals from the **primal NLP** into `Transient_Stability/CSV_duals/` and `OPF_Duals_Results.xlsx`. These are multipliers on the transient constraints from [6. The TSC-OPF, assembled](06_tsc_opf_assembled.md): initial links, trapezoidal swing rows, $P_e$ definitions, COI equalities, and $\delta$ corridor inequalities.

!!! info "Assumption"
    Nonlinear TSC-ACOPF duals characterise the **specific** KKT point Ipopt converged to. They are not guaranteed to be unique global shadow prices the way an LP dual is. Use them for local sensitivity and economic interpretation at the solved trajectory, not as proof of convex market equilibrium.

### Which multipliers matter economically

| Family | Registry name | Constraint | Typical reading when binding |
|---|---|---|---|
| Initial P–δ link | `dual_Pe_init` | [Model 6.1](06_tsc_opf_assembled.md) | Cost of enforcing consistency between dispatch and initial swing state |
| COI definition | `dual_δCOI` | [Model 6.6](06_tsc_opf_assembled.md) | Gauge on the inertia-weighted reference angle |
| Swing discretisation | `dual_δ`, `dual_Δω` | [Models 6.3–6.4](06_tsc_opf_assembled.md) | Shadow on integrating the swing one more step |
| Electrical power | `dual_Pe` | [Model 6.5](06_tsc_opf_assembled.md) | Network torque on the rotor at $(g,t)$ |
| **Stability corridor (COI)** | `dual_δ_COI_upper`, `dual_δ_COI_lower` | [Model 6.7](06_tsc_opf_assembled.md) | **Stability scarcity**: margin on rotor angle w.r.t. COI |
| **Stability corridor (reference machine)** | `dual_δ_ref_upper`, `dual_δ_ref_lower` | `bound_style_δ ∈ (:highest_H, :ref_gen)` | Same scarcity, priced against one live machine instead of the COI |
| **Stability corridor (pairwise)** | `dual_δ_pair_upper`, `dual_δ_pair_lower` | `bound_style_δ = :pairwise` | Same scarcity with no reference at all; columns are `Pair_G<i>_G<j>`, so a non-zero column names the separating pair |
| Speed COI box (optional) | `dual_Δω_COI_upper`, `dual_Δω_COI_lower` | `constrain_Δω` with `:coi_box` | Frequency-coherence scarcity: spread between machines |
| Absolute speed box (optional) | `dual_Δω_abs_upper`, `dual_Δω_abs_lower` | `constrain_Δω` with `:abs` | Also prices the fleet's common-mode frequency excursion |

When the upper $\delta$ bound binds for generator $g$ at time $t$, `dual_δ_COI_upper` is the incremental cost (in the Lagrangian sense) of tightening the corridor by one unit. The reference-machine family `dual_δ_ref_upper` answers the same question against a different yardstick, so the two are not interchangeable and never appear in the same run: the COI is an inertia-weighted average that moves with the whole fleet, while a reference machine is one physical unit whose own swing enters every row. Comparing the price surfaces over the same contingency is the point of running more than one. The pairwise family `dual_δ_pair_upper` answers the question a third way, and is the only one of the three whose value does not depend on the frame: shifting every angle by a constant leaves each $\delta_i - \delta_j$ untouched, so the price is attached to a *pair of units* rather than to a unit and a convention. Its columns are mostly zero by construction — at any step only the widest-spread pair can be active — and that sparsity is the readout, not a defect. That is the object you would trace when asking how stability constraints re-shape effective prices, not the nodal $\lambda_k$ from dispatch alone.

---

## Steady-state export map

Authoritative source: `src/_manage_outputs/DispatchDualRegistry.jl` (`STEADY_STATE_DUAL_SPECS`). The save layer writes only families present in the solved model.

| XLSX sheet | Constraint family | Dual / price | Id column |
|---|---|---|---|
| `P_Balance` | Active nodal balance | $\lambda_k$; LMP $= +\lambda_k$ | `Bus_ID` |
| `Q_Balance` | Reactive nodal balance | $\beta_k$ | `Bus_ID` |
| `Pg_Lo` / `Pg_Up` | Generator P limits | $\underline{\eta}_g$ / $\overline{\eta}_g$ | `Gen_ID` |
| `Qg_Lo` / `Qg_Up` | Generator Q limits | reactive bound multipliers | `Gen_ID` |
| `V_Mag_Lo` / `V_Mag_Up` | Voltage magnitude | $\alpha_k^-$ / $\alpha_k^+$ | `Bus_ID` |
| `V_Ang_Lo` / `V_Ang_Up` | Voltage angle box | $\mu_k^-$ / $\mu_k^+$ | `Bus_ID` |
| `Ang_Diff_Lo` / `Ang_Diff_Up` | Branch angle difference | $\rho_{km}^-$ / $\rho_{km}^+$ | `Branch_ID` |
| `Sik_Upper` / `Ski_Upper` | Branch thermal limits | congestion multipliers | `Branch_ID` |
| `Sg_Upper` | Gen capability curve | capability rent | `Gen_ID` |

CSV mirrors: one file per family under `Dispatch/CSV_duals/` (e.g. `dual_P_balance.csv`).

---

## Transient export map (classical Kron)

Authoritative source: `src/_transient_stability/DynDualRegistry.jl`, which holds **three** spec lists selected by the active path — `CLASSICAL_KRON_DUAL_SPECS` (the table below), `CLASSICAL_FULL_BUS_DUAL_SPECS` (adds nodal balance and voltage-bound families), and `DQ_4TH_FULL_BUS_DUAL_SPECS` (drops the classical-only entries and adds the two-axis machine, AVR, governor, and GFM pre-fault families). The control and converter duals — `dual_Vref_init`, `dual_avr_E_fd`, `dual_avr_E_fd_sat`, `dual_Pref_init`, `dual_gov_valve`, `dual_gov_mech`, and the limiter-specific ones — appear only on that third list. Read them alongside [8. Machine controls](08_controls_avr_governor.md) and [9. Grid-forming inverters](09_grid_forming.md), which explain what each multiplier prices.

A sheet is named after its registry export with `dual_` dropped and `_Duals` appended; the CSV files use ASCII names.

| XLSX sheet | Registry export | CSV |
|---|---|---|
| `Pe_init_Duals` | `dual_Pe_init` | `dual_Pe_init.csv` |
| `δCOI_Duals` | `dual_δCOI` | `dual_delta_COI.csv` |
| `Pe_Duals` | `dual_Pe` | `dual_Pe.csv` |
| `δ_Duals` | `dual_δ` | `dual_delta.csv` |
| `Δω_Duals` | `dual_Δω` | `dual_Delta_Omega.csv` |
| `δ_COI_lower_Duals` | `dual_δ_COI_lower` | `dual_delta_COI_lower.csv` |
| `δ_COI_upper_Duals` | `dual_δ_COI_upper` | `dual_delta_COI_upper.csv` |
| `δ_ref_lower_Duals` | `dual_δ_ref_lower` | `dual_delta_ref_lower.csv` |
| `δ_ref_upper_Duals` | `dual_δ_ref_upper` | `dual_delta_ref_upper.csv` |
| `δ_pair_lower_Duals` | `dual_δ_pair_lower` | `dual_delta_pair_lower.csv` |
| `δ_pair_upper_Duals` | `dual_δ_pair_upper` | `dual_delta_pair_upper.csv` |
| `Δω_COI_lower_Duals` | `dual_Δω_COI_lower` | `dual_Delta_Omega_COI_lower.csv` |
| `Δω_COI_upper_Duals` | `dual_Δω_COI_upper` | `dual_Delta_Omega_COI_upper.csv` |
| `Δω_abs_lower_Duals` | `dual_Δω_abs_lower` | `dual_Delta_Omega_abs_lower.csv` |
| `Δω_abs_upper_Duals` | `dual_Δω_abs_upper` | `dual_Delta_Omega_abs_upper.csv` |

Governor and AVR paths add `dual_Pref_init`, `dual_gov_valve`, `dual_gov_mech`, `dual_gov_valve_sat` (or `dual_LB_gov_valve` / `dual_UB_gov_valve`, depending on `governor_limiter`), `dual_Vref_init`, `dual_avr_E_fd`, etc., when those constraints are built. What each of those multipliers prices is in [8. Machine controls](08_controls_avr_governor.md); the file-by-file output inventory is in [Dynamic controls (AVR / governor)](../dynamic_controls_avr_governor.md).

**CSV export is registry-driven.** Every entry in the catalog that declares a `csv_file` and whose constraint family exists in the run is written to `Transient_Stability/CSV_duals/`. Adding a `DualRegistryEntry` is the only step needed to export a new family — the save layer needs no edit. That includes the optional variable-bound duals (`dual_LB_Ed`, `dual_UB_V_ref`, `dual_LB_Pv`, …), which appear as soon as the matching `TsBuilderConfig.bound_*` toggle builds the ≤-row.

---

## Worked reading (9-bus, contingency 2)

!!! tip "Example 7.1 · tracing a binding stability dual"
    Run the Kron TSC-ACOPF example from [6. The TSC-OPF, assembled](06_tsc_opf_assembled.md) with `δ_tol_deg = 120` and `save_duals = true`. After a successful solve:

    1. Open `RESULTS/Results - <timestamp>/Dispatch/Dispatch_Duals.xlsx` and note $\lambda_k$ on the marginal buses (LMP $= +\lambda_k$, so the value should read positive and near the marginal unit's `g_cost_1`).
    2. Open `Transient_Stability/OPF_Duals_Results.xlsx` → `δ_COI_upper_Duals`. Find generator 2 in the fault window; a non-zero entry means the upper $\delta$ corridor bound was active at that time step.
    3. Compare magnitude and sign with `Transient_Stability/CSV_duals/dual_delta_COI_upper.csv` (same registry row).

    If you enable `dispatch.solve_explicit_dual = true` on a **dispatch-only** linear DCOPF of the same case, $\pi_k$ from `Dispatch_Dual/` should match $+\lambda_k$ from the primal dispatch duals on balance rows (up to solver tolerance) — the two multipliers differ by a sign, the two prices do not.

---

## Where this lives in the code

| Task | Module | Config |
|---|---|---|
| Steady-state dual registry | `DispatchDualRegistry.jl` | `save_duals`, `DispatchConfig` toggles |
| Transient dual registry | `DynDualRegistry.jl` | `transient.builder`, `dyn_model` |
| Primal dual export | `functions_2_save_dispatch_duals.jl`, `Save_Duals_Dynamic_Model_*` | `save_duals = true` |
| Explicit ED / DC dual LP | `_ed/`, `_dcopf/` dual builders | `solve_explicit_dual = true` |

Parameter-level export flags: [Parameter reference — `RunConfig`](../parameter_reference.md).
