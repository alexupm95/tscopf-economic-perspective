# 5. Swing dynamics and discretization

Transient stability in TSCOPF is not a separate time-domain simulation pasted onto an OPF result. The swing equations are **equality constraints** inside the same JuMP model as dispatch. A feasible solution is a dispatch **and** a discrete trajectory of rotor angles and speed deviations that satisfy the trapezoidal discretisation of the classical swing equation, together with the electrical-power map from [4. The network: Kron vs `FULL_BUS`](04_network_kron_fullbus.md).

This page states the continuous-time model, derives the trapezoidal rows the builders install, defines the centre-of-inertia (COI) reference, and explains how rotor-angle corridors are enforced. The full coupled optimisation problem is assembled on [6. The TSC-OPF, assembled](06_tsc_opf_assembled.md).

| Symbol | Meaning |
|---|---|
| $H_g$, $D_g$ | Inertia constant [s] and damping on machine $g$ (`H`, `D` in `gen_dynamic_data.csv`) |
| $\omega_s$ | Synchronous speed [rad/s]; $2\pi f_{\mathrm{syn}}$ from `f_syn` |
| $\Delta t$ | Trapezoidal step [s]; `t_step` |
| $\delta_g^t$, $\Delta\omega_g^t$ | Rotor angle and per-unit speed deviation at step $t$ |
| $P_{m,g}$, $P_g^{\mathrm{ele},t}$ | Mechanical and electrical power [p.u.] |
| $\delta^{\mathrm{COI},t}$ | Inertia-weighted COI angle |
| $\delta^{\max}$ | Half-width of the stability corridor; from `δ_tol_deg` |

Machine parameters and timing live on `TransientConfig` — see [Parameter reference — `TransientConfig`](../parameter_reference.md).

---

## Continuous swing equation

Each classical machine is a second-order rotor:

!!! note "Model 5.1 (continuous swing)"
    For generator $g$,

    ```math
    \begin{equation}
    \label{eq:swing-cont}
    \frac{2H_g}{\omega_s}\frac{d\Delta\omega_g}{dt}
    = P_{m,g} - P_g^{\mathrm{ele},t} - D_g \Delta\omega_g,
    \qquad
    \frac{d\delta_g}{dt} = \omega_s \Delta\omega_g .
    \end{equation}
    ```

$P_g^{\mathrm{ele},t}$ couples to the network (Kron sum in Model 6.5 on [page 6](06_tsc_opf_assembled.md), or nodal classical maps on `FULL_BUS`). During a short-circuit contingency the same ODE holds in the fault-on and post-fault windows; only $Y_{\mathrm{bus}}$ (and hence electrical power) changes at clearing.

!!! info "Assumption"
    Governor and AVR dynamics are not modelled on the default classical path. $P_{m,g}$ is either fixed at dispatch $P_g$ (`mech_power_mode = USE_PG`, Kron default) or an explicit variable (`USE_PM`, required on `FULL_BUS`). See [8. Machine controls: AVR and turbine governor](08_controls_avr_governor.md) for the extensions that make $P_{m,g}$ and $E_{fd,g}$ time-varying.

---

## Trapezoidal discretisation

TSCOPF replaces $\eqref{eq:swing-cont}$ with the **implicit trapezoidal rule** on each time window. The angle update is explicit in $\Delta\omega$; the swing equation is implicit in $\Delta\omega^t$ and $P_g^{\mathrm{ele},t}$.

### Angle update

Integrating $d\delta/dt = \omega_s \Delta\omega$ with trapezoids over $[t-\Delta t, t]$ gives:

!!! note "Model 5.2 (trapezoidal angle update)"
    For $t \geq 1$ inside a window,

    ```math
    \begin{equation}
    \label{eq:trap-delta-5}
    \delta_g^t - \delta_g^{t-1}
    = \omega_s \frac{\Delta t}{2}\left(\Delta\omega_g^t + \Delta\omega_g^{t-1}\right) .
    \end{equation}
    ```

    The **first** step in a window uses the pre-fault (or antecedent-window) state $(\delta_g^0, \Delta\omega_g^0)$ instead of $\delta_g^{t-1}$. Post-fault steps chain from the last fault-on values.

    Implemented in `eq_const_kron_δ_swingeq_generic!` (`functions_4_TS_kron_eqconst.jl`).

### Swing equation (implicit trapezoid)

Apply the trapezoidal rule to $\frac{2H_g}{\omega_s}\,\frac{d\Delta\omega}{dt} = P_{m,g} - P_{\mathrm{ele}} - D_g\Delta\omega$. Rearranging the implicit form

$$\left(1 + \frac{D_g \Delta t}{4H_g}\right)\Delta\omega_g^t
- \left(1 - \frac{D_g \Delta t}{4H_g}\right)\Delta\omega_g^{t-1}
= \frac{\Delta t}{4H_g}\left(2P_{m,g} - P_g^{\mathrm{ele},t} - P_g^{\mathrm{ele},t-1}\right)$$

yields the row stamped at each $(g,t)$:

!!! note "Model 5.3 (trapezoidal swing equation)"
    ```math
    \begin{equation}
    \label{eq:trap-swing-5}
    \left(1 + \frac{D_g \Delta t}{4H_g}\right)\Delta\omega_g^t
    - \left(1 - \frac{D_g \Delta t}{4H_g}\right)\Delta\omega_g^{t-1}
    - \frac{\Delta t}{4H_g}\left(2P_{m,g} - P_g^{\mathrm{ele},t} - P_g^{\mathrm{ele},t-1}\right) = 0 .
    \end{equation}
    ```

    `eq_const_tsred_Δω_swingeq_generic!` implements $\eqref{eq:trap-swing-5}$ for Kron TSC-ACOPF. TSC-DCOPF uses a linearised variant (`eq_const_tsredlinear_Δω_swingeq_generic_modified!`) with the same trapezoidal structure on affine electrical power.

!!! warning "Time-step sensitivity"
    The trapezoidal rule is second-order accurate in $\Delta t$ for smooth trajectories. In practice, too large a `t_step` can smooth peak swings and weaken the stability constraint. If margins look optimistic, halve `t_step` and compare `Transient_Stability/CSV/` angle traces before trusting duals near the corridor.

---

## Centre of inertia

Absolute rotor angles are gauge-dependent. TSCOPF measures stability relative to an inertia-weighted mean:

!!! note "Model 5.4 (COI angle)"
    ```math
    \begin{equation}
    \label{eq:coi-5}
    \delta^{\mathrm{COI},t} = \frac{\sum_{g \in \mathcal{G}_{\mathrm{act}}} H_g \delta_g^t}{\sum_{g \in \mathcal{G}_{\mathrm{act}}} H_g} .
    \end{equation}
    ```

    `eq_const_kron_COI_generic!` enforces $\eqref{eq:coi-5}$ at every step. $\mathcal{G}_{\mathrm{act}}$ is the **active** generator set for the disturbance window (`active_gen` passed to the builder). For SC faults it equals all machines; for **GL generator trip** it excludes disconnected units so COI is defined on the synchronized survivors only.

    The weights come from `_coi_weights`, which looks `H_g` up **by generator id** rather than by row position — the two coincide only while `DGEN_DYN.id == 1:n`, which `validate_dgen_dyn_row_ids!` guarantees on the SG-only paths but not on a mixed fleet. On a mixed SG + GFM fleet $\mathcal{G}_{\mathrm{act}}$ is narrowed to the machines unless `DynModelConfig.gfm_in_coi` is set, in which case each converter joins with the equivalent inertia of its droop loop, $H^{\mathrm{eff}}_c = T^P_{f,c}/(2 m_{p,c})$ — see chapter 9.

Synchronism is preserved when each machine stays within a corridor around COI, not around zero.

---

## Rotor-angle stability bounds

Set a symmetric half-width with `δ_tol_deg`, or independent below/above limits with `δ_tol_deg_lower` and `δ_tol_deg_upper` (degrees, each positive; unset fields fall back to `δ_tol_deg`). Internally `common_ts_parameters` builds `δ_tol = (−lower, +upper)` in radians.

Two COI-referenced styles exist. They bound $\delta_g^t - \delta^{\mathrm{COI},t}$ inside $[\delta^{\mathrm{low}}, \delta^{\mathrm{hi}}]$ (symmetric when only `δ_tol_deg` is set) but produce different JuMP rows (and therefore different dual exports).

### Direct COI box (`:coi_box`)

!!! note "Model 5.5 (COI-referenced angle corridor)"
    ```math
    \begin{equation}
    \label{eq:delta-coi-box-5}
    \delta^{\mathrm{low}} \leq \delta_g^t - \delta^{\mathrm{COI},t} \leq \delta^{\mathrm{hi}} .
    \end{equation}
    ```

    `ineq_const_kron_δ_COI_generic!` stamps $\eqref{eq:delta-coi-box-5}$ directly on $(\delta_g^t, \delta^{\mathrm{COI},t})$ using `δ_tol[1]` / `δ_tol[2]` from `common_ts_parameters`. **Required** (or another box-form style) when `mech_power_mode = USE_PM` (explicit $P_{m,g}$), and therefore on every `FULL_BUS` run. It is not the default: set `bound_style_δ = :coi_box` explicitly.

### Swing-propagated (`:swing_propagated`)

For $t > 1$, substitute the trapezoidal updates into $\eqref{eq:delta-coi-box-5}$ so the inequality involves $(\delta_g^{t-1}, \Delta\omega_g^{t-1}, P_g^{\mathrm{ele},t}, P_g^{\mathrm{ele},t-1}, P_{m,g})$ instead of $\delta_g^t$ alone. The propagated rows are built by `ineq_const_kron_δ_COI_generic_modified!`.

| `bound_style_δ` | Reference | Constraint function |
|---|---|---|
| `:swing_propagated` | $\delta^{\mathrm{COI}}$ | `ineq_const_kron_δ_COI_generic_modified!` |
| `:coi_box` | $\delta^{\mathrm{COI}}$ | `ineq_const_kron_δ_COI_generic!` |
| `:highest_H` | largest-$H$ surviving machine | `ineq_const_kron_δ_COI_generic!` with `δ[ref]` as the reference series |
| `:ref_gen` | `δ_ref_gen_id` | same |
| `:pairwise` | none — every unordered pair | `ineq_const_kron_δ_pairwise!` |

Every network form dispatches on `bound_style_δ` through the same pair of entry points,
`_add_δ_bounds_fault!` / `_add_δ_bounds_postf!` in `functions_4_TS_kron_ineqconst.jl` — Kron,
Kron-linear, classical `FULL_BUS` and DQ `FULL_BUS` alike. The default is `:swing_propagated`
everywhere; `FULL_BUS` requires `USE_PM`, which rejects it, so a `FULL_BUS` run must name a
box-form style. The choice is honoured on every network form.

### Machine-referenced (`:highest_H`, `:ref_gen`)

Instead of the inertia-weighted average, the band is measured against one live machine:

```math
\delta^{\mathrm{tol}}_{\mathrm{lower}} \;\le\; \delta_g^t - \delta_{\mathrm{ref}}^t \;\le\; \delta^{\mathrm{tol}}_{\mathrm{upper}},
\qquad g \neq \mathrm{ref}.
```

The reference is excluded from its own corridor, so an $n$-unit fleet yields $n-1$ rows
per time step. `:highest_H` resolves the reference over the machines that survive the
disturbance, so a unit the run trips cannot be selected; `:ref_gen` takes the id from
`DynModelConfig.δ_ref_gen_id`, validated against the system data before the warm-start
solve.

On a mixed fleet these two styles span the grid-forming converters as well: measuring
$\delta_g - \delta_{\mathrm{ref}}$ needs no inertia from $g$, only an angle in the same
frame, which the converter angle is. The *reference* under `:highest_H` is still always a
machine, since that style ranks by $H$ and a converter carries none — a converter becomes
the reference only through an explicit `:ref_gen`. See
[9. Grid-forming converters](09_grid_forming.md) for what remains machine-only.

Because nothing references $\delta^{\mathrm{COI}}$ under these styles, it is built as a JuMP
expression rather than a variable with a defining equality — the solver never sees it, but
the COI trajectory is still written to `angle_rel_COI.csv` for comparison against a
COI-referenced run. The duals are exported separately as `dual_δ_ref_lower` /
`dual_δ_ref_upper`.

### Pairwise separation (`:pairwise`)

The classical loss-of-synchronism criterion, stated without a reference at all:

```math
-\Delta^{\max} \;\le\; \delta_i^t - \delta_j^t \;\le\; +\Delta^{\max},
\qquad i < j,\ i,j \in \mathcal{G}_{\mathrm{SG}} \cup \mathcal{G}_{\mathrm{GFM}} .
```

$\Delta^{\max}$ comes from its own knob, `TsSimulationConfig.δ_pair_tol_deg` (default
180°), not from `δ_tol_deg*`: the two live on different scales, and a separation limit is
symmetric by construction because the pair $(i,j)$ and the pair $(j,i)$ are the same
constraint read from opposite sides. The rows are built by `ineq_const_kron_δ_pairwise!` and
stored under the `δ_pair_*` keys.

An $n$-unit fleet yields $n(n-1)/2$ rows per time step instead of $n-1$, and every unit
carries rows — converters included, for the same reason the machine-referenced styles bound
them. Only the widest-spread pair can be active at a given step, so most multipliers are
zero; that is the point, because the non-zero column of `dual_delta_pair_upper.csv` names
the two units the corridor is holding apart. When two pairs bind simultaneously with the
same spread the split between their multipliers is not unique, and only their sum is
meaningful.

!!! warning "Thresholds are not comparable across styles"
    $|\delta_i - \delta_{\mathrm{ref}}| \le \Delta$ for all $i$ implies only
    $|\delta_i - \delta_j| \le 2\Delta$. A pairwise run at $\Delta$ is a *different*
    feasible set from a reference-framed run at the same number; comparisons must be made on
    the realised spread.

Like the machine-referenced styles, `:pairwise` leaves $\delta^{\mathrm{COI}}$ as a JuMP
expression, so no `dual_δCOI` exists and `bound_δCOI_*` is rejected at validation. Duals
export as `dual_δ_pair_lower` / `dual_δ_pair_upper`, with columns labelled
`Pair_G<i>_G<j>` — the only two-index family in the model.

!!! note "Interpretation"
    When the upper bound in $\eqref{eq:delta-coi-box-5}$ binds, the dispatch is marginal with respect to machine $g$'s swing margin at time $t$. The corresponding dual (`dual_δ_COI_upper` in CSV export) prices **stability scarcity** at $(g,t)$ — not a nodal LMP. Sign conventions are on [7. Duals, KKT, and the economics](07_duals_economics.md).

### Optional speed corridor

Set `DynModelConfig.constrain_Δω = true` to bound the speed deviation (`ineq_const_kron_Δω_COI_generic!`). This is off by default, and `bound_style_Δω` picks the reference: `:coi_box` bounds $\Delta\omega_g^t - \Delta\omega^{\mathrm{COI},t}$, `:abs` bounds the raw $\Delta\omega_g^t$ and builds no COI at all.

The two are not the same constraint with the COI pinned at zero. $\Delta\omega^{\mathrm{COI}}$ is a free variable, so the COI box prices only the spread between machines; the absolute box additionally pins the fleet's common-mode drift from synchronous speed. That drift is set mostly by the accelerating power during the fault, which redispatch cannot remove, so an absolute band tighter than it renders the problem infeasible rather than expensive.

The corridor need not be symmetric. `Δω_tol_pu` sets a symmetric half-width; `Δω_tol_pu_lower` and `Δω_tol_pu_upper` override each side independently as positive magnitudes, so

```math
-\Delta\omega^{\mathrm{tol}}_{\mathrm{lower}} \;\le\; \Delta\omega_g^t - \Delta\omega^{\mathrm{COI},t} \;\le\; \Delta\omega^{\mathrm{tol}}_{\mathrm{upper}} .
```

This mirrors `δ_tol_deg_lower` / `δ_tol_deg_upper` for the angle corridor, and lets under- and over-frequency excursions be priced against different limits (`dual_Δω_COI_lower` / `dual_Δω_COI_upper`).

---

## Time windows

Discrete steps are built from `TransientConfig.simulation`:

| Field | Role |
|---|---|
| `t_start_fault` | Start of disturbance |
| `clearing_time` | Fault duration (SC only); clearing at $t_{\mathrm{start}} + t_{\mathrm{clear}}$. The fault-on network still holds at that sample, so the fault effectively lasts `clearing_time + t_step` |
| `t_end_sim` | End of post-fault horizon |
| `t_step` | $\Delta t$ in Models 5.2–5.3 |
| `f_syn`, `Δω_0` | $\omega_s = 2\pi f_{\mathrm{syn}}$; initial speed deviation |

For a **short-circuit** (`FaultConfig` type SC), `time_windows_sc` builds:

- **Fault-on** `*_tf`: $t \in [t_{\mathrm{start}}, t_{\mathrm{start}} + t_{\mathrm{clear}}]$ stepped by $\Delta t$
- **Post-fault** `*_tpf`: $t \in [t_{\mathrm{start}} + t_{\mathrm{clear}} + \Delta t, t_{\mathrm{end}}]$ stepped by $\Delta t$

The fault-on network therefore applies through the clearing sample, and the post-fault
network takes over one step later: a 100 ms fault on a 5 ms grid is `clearing_time = 0.095`.

Generator/load trips (GL) use a **single** window from fault start to `t_end_sim`. Network topology and $Y_{\mathrm{bus}}$ change between windows as described on [page 4](04_network_kron_fullbus.md).

---

## Worked example

!!! tip "Example 5.1 · reading the discretisation on IEEE 9-bus"
    Contingency 2 (fault at bus 7, trip branch 5–7) with a generous corridor for debugging:

    ```julia
    cfg = RunConfig(
        trans_stab = true, case = "9bus", solver_name = "Ipopt",
        dispatch = DispatchConfig(type_model = "ACOPF"),
        transient = TransientConfig(
            simulation = TsSimulationConfig(
                δ_tol_deg = 100.0,
                t_step    = 0.01,
                clearing_time = 0.1,
            ),
            dyn_model = DynModelConfig(
                fault = FaultConfig(contingency_id = 2),
            ),
        ),
    )
    ```

    After `run_case!(cfg, sys, path_main, path_results)`:

    1. **`Transient_Stability/CSV/angle_rel_COI.csv`** — check $\delta_g^t - \delta^{\mathrm{COI},t}$ stays inside $\pm\delta^{\max}$.
    2. **`Transient_Stability/CSV/speed_dev.csv`** against **`angle_abs.csv`** — verify trapezoidal consistency of the speed and angle traces.
    3. Halve `t_step` and re-run; if the optimal dispatch changes, the coarse grid was masking swing peaks.

    `FULL_BUS` runs need `network_form = FULL_BUS`, `mech_power_mode = USE_PM`, and `bound_style_δ = :coi_box` (Example 4.1 on [page 4](04_network_kron_fullbus.md)).

---

## Where this lives in the code

| Piece | Function | File |
|---|---|---|
| $\eqref{eq:trap-delta-5}$ | `eq_const_kron_δ_swingeq_generic!` | `functions_4_TS_kron_eqconst.jl` |
| $\eqref{eq:trap-swing-5}$ | `eq_const_tsred_Δω_swingeq_generic!` | `functions_4_TS_kron_eqconst.jl` |
| $\eqref{eq:coi-5}$ | `eq_const_kron_COI_generic!` | `functions_4_TS_kron_eqconst.jl` |
| $\eqref{eq:delta-coi-box-5}$ direct | `ineq_const_kron_δ_COI_generic!` | `functions_4_TS_kron_ineqconst.jl` |
| Propagated corridor | `ineq_const_kron_δ_COI_generic_modified!` | `functions_4_TS_kron_ineqconst.jl` |
| Kron window assembly | `Define_Fault_Dynamic_Model_tsred!`, `Define_PostFault_Dynamic_Model_tsred!` | `functions_2_build_TS_model_w_Kron.jl` |
| `FULL_BUS` bound dispatch | `_add_fullbus_δ_COI_bounds_fault!` | `functions_4_TS_fullbus_ineqconst.jl` |
| Timing / `δ_tol` tuple | `common_ts_parameters`, `time_windows_sc` | `TsConfig.jl` |

Electrical power during the swing and the steady-state initial links are Models 6.1–6.5 on [6. The TSC-OPF, assembled](06_tsc_opf_assembled.md).
