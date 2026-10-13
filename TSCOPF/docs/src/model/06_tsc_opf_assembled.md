# 6. The TSC-OPF, assembled

A TSC-OPF run is one optimisation problem: minimise fuel cost at the steady state, subject to power flow **and** a trajectory of rotor angles and speeds that survives a disturbance. The dispatch you would have chosen from a plain OPF is still there. What changes is the feasible set. Any dispatch that would lose synchronism after a fault is ruled out before you see a solution.

This page assembles the pieces from [2. Steady-state OPF](02_steady_state_opf.md) with the transient constraints implemented in TSCOPF today: Kron-reduced TSC-ACOPF (classical 2nd-order machines), the FULL_BUS nodal path, and TSC-DCOPF (linearised electrical power). Swing discretisation and Kron reduction get fuller treatment on [5. Swing dynamics and discretization](05_swing_dynamics.md) and [4. The network: Kron vs FULL_BUS](04_network_kron_fullbus.md); here they appear where they belong in the coupled model.

| Symbol | Meaning |
|---|---|
| $\mathcal{T}$, $t$ | Discrete time steps in a fault or post-fault window |
| $\delta_g^t$, $\Delta\omega_g^t$ | Rotor angle and speed deviation at step $t$ |
| $P_g^{\mathrm{ele},t}$ | Electrical power delivered to the network at $t$ |
| $\delta^{\mathrm{COI},t}$ | Inertia-weighted centre-of-inertia angle |
| $\delta^{\max}$ | Half-width of the stability corridor, from `δ_tol_deg` |
| $G^{\mathrm{red}}, B^{\mathrm{red}}$ | Kron-reduced admittance among generator buses |
| $Y_{\mathrm{bus}}$ | Full nodal admittance (FULL_BUS path) |
| $\omega_s$, $\Delta t$ | Synchronous speed [rad/s] and trapezoidal step [s] |

Set `RunConfig.trans_stab = true` and provide `RunConfig.transient` to activate this block on top of the steady-state formulation.

---

## Problem structure

At a high level, TSCOPF builds a single JuMP model with two layers:

1. **Dispatch (pre-fault)** — ACOPF, DCOPF, or ED as in [2. Steady-state OPF](02_steady_state_opf.md).
2. **Transient (fault and post-fault windows)** — classical swing dynamics, network-dependent $P_e$, COI definition, and optional rotor-angle limits relative to COI.

The two layers meet at **initial conditions**: the dispatch $P_g$, $Q_g$, $V_k$, $\theta_k$ must be consistent with initial rotor states $\delta_g^0$, $E_g$, and (when used) $P_{m,g}$. Data flows: steady-state OPF → initial $P$–$\delta$–$E$ link → trapezoidal swing and $P_e$ → COI definition and $\delta$ corridor limits.

---

## TSC-ACOPF with Kron reduction

The default path (`DynModelConfig.network_form = KRON_REDUCED`) keeps only generator-internal buses in the transient network. Loads and passive buses are eliminated into $Y^{\mathrm{red}}$; electrical power during the swing is a function of generator angles only.

### Steady-state block

Identical to Model 2.1–2.2 on the [steady-state page](02_steady_state_opf.md): minimise fuel cost subject to AC nodal balance and limits. TSC-ACOPF requires `dispatch.type_model = "ACOPF"` and a nonlinear solver (Ipopt or MadNLP).

### Initial conditions (OPF → transient)

At $t=0$, dispatch and classical-machine states must agree:

!!! note "Model 6.1 (initial active power link)"
    For each generator $g$ at terminal bus $k(g)$,

    ```math
    \begin{equation}
    \label{eq:pe-init-p}
    P_{g(k)} - \frac{E_g V_{k(g)} \sin(\delta_g^0 - \theta_{k(g)})}{x'_g} = 0
    \end{equation}
    ```

!!! note "Model 6.2 (initial reactive power link)"
    ```math
    \begin{equation}
    \label{eq:pe-init-q}
    Q_{g(k)} - \frac{E_g V_{k(g)} \cos(\delta_g^0 - \theta_{k(g)})}{x'_g} + \frac{V_{k(g)}^2}{x'_g} = 0
    \end{equation}
    ```

$E_g$ and $\delta_g^0$ are decision variables at the pre-fault point. They are the bridge between the OPF voltage solution and the swing that follows.

### Swing dynamics (trapezoidal)

Continuous-time classical swing:

```math
\frac{2H_g}{\omega_s}\frac{d\Delta\omega_g}{dt} = P_{m,g} - P_g^{\mathrm{ele},t} - D_g \Delta\omega_g, \qquad
\frac{d\delta_g}{dt} = \omega_s \Delta\omega_g .
```

TSCOPF discretises with the trapezoidal rule on each time window (fault `*_tf`, post-fault `*_tpf` for short-circuit contingencies):

!!! note "Model 6.3 (trapezoidal angle update)"
    For $t \geq 1$ in a window,

    ```math
    \begin{equation}
    \label{eq:trap-delta}
    \delta_g^t - \delta_g^{t-1} = \omega_s \frac{\Delta t}{2}\left(\Delta\omega_g^t + \Delta\omega_g^{t-1}\right)
    \end{equation}
    ```

    The first step in a window couples to the pre-fault state $(\delta_g^0, \Delta\omega_g^0)$ instead of $\delta_g^{t-1}$.

!!! note "Model 6.4 (trapezoidal swing equation)"
    ```math
    \begin{equation}
    \label{eq:trap-swing}
    \left(1 + \frac{D_g \Delta t}{4H_g}\right)\Delta\omega_g^t
    - \left(1 - \frac{D_g \Delta t}{4H_g}\right)\Delta\omega_g^{t-1}
    - \frac{\Delta t}{4H_g}\left(2P_{m,g} - P_g^{\mathrm{ele},t} - P_g^{\mathrm{ele},t-1}\right) = 0
    \end{equation}
    ```

!!! info "Assumption"
    Mechanical power $P_{m,g}$ is fixed at the dispatch $P_g$ when `mech_power_mode = USE_PG` (Kron default). With `USE_PM`, $P_m$ is an explicit variable (requires `bound_style_δ = :coi_box` on FULL_BUS paths).

### Electrical power during the transient (nonlinear Kron)

!!! note "Model 6.5 (Kron electrical power)"
    At each time $t$ and for each active generator $g$,

    ```math
    \begin{equation}
    \label{eq:pe-kron}
    P_g^{\mathrm{ele},t} = E_g \sum_{i \in \mathcal{G}} E_i \left(
      G^{\mathrm{red}}_{gi} \cos(\delta_g^t - \delta_i^t)
      + B^{\mathrm{red}}_{gi} \sin(\delta_g^t - \delta_i^t)
    \right)
    \end{equation}
    ```

$G^{\mathrm{red}}$ and $B^{\mathrm{red}}$ are rebuilt when the network topology changes: augmented $Y_{\mathrm{bus}}$ during a bus fault, cleared branch after fault removal, or modified admittance after a generator/load trip.

### Centre of inertia and stability bounds

!!! note "Model 6.6 (COI angle)"
    ```math
    \begin{equation}
    \label{eq:coi}
    \delta^{\mathrm{COI},t} = \frac{\sum_{g \in \mathcal{G}_{\mathrm{act}}} H_g \delta_g^t}{\sum_{g \in \mathcal{G}_{\mathrm{act}}} H_g}
    \end{equation}
    ```

    $\mathcal{G}_{\mathrm{act}}$ is the active generator set for the window (all machines on SC; survivors only after GL gen-trip). See [5. Swing dynamics](05_swing_dynamics.md#centre-of-inertia).

Rotor angles are limited relative to COI, not to an absolute reference. With `bound_style_δ = :coi_box`, the direct inequalities are:

!!! note "Model 6.7 (COI-referenced angle corridor)"
    ```math
    \begin{equation}
    \label{eq:delta-coi-box}
    \delta^{\mathrm{low}} \leq \delta_g^t - \delta^{\mathrm{COI},t} \leq \delta^{\mathrm{hi}}
    \end{equation}
    ```

    Defaults: $\delta^{\mathrm{low}} = -\texttt{δ\_tol\_deg}$, $\delta^{\mathrm{hi}} = +\texttt{δ\_tol\_deg}$. Override with `δ_tol_deg_lower` / `δ_tol_deg_upper` on `TsSimulationConfig`.

    $\delta^{\max}$ comes from `transient.simulation.δ_tol_deg` converted to radians. With `:swing_propagated`, the same corridor is enforced through propagated swing constraints (algebraically equivalent intent, different constraint rows for dual extraction).

!!! note "Interpretation (binding stability bound)"
    When the upper bound in $\eqref{eq:delta-coi-box}$ is active for generator $g$ at time $t$, the dispatch is marginal with respect to that machine's swing margin: a small increase in $P_g$ or a slightly weaker network would violate the corridor. The corresponding dual entry in `Transient_Stability/CSV_duals/` (`dual_δ_COI_upper`) is the local KKT multiplier on that inequality. It is not a nodal LMP, but it prices **stability scarcity** at $(g,t)$. Full dual sign conventions are on [7. Duals, KKT, and the economics](07_duals_economics.md).

### Time windows and disturbances

| Disturbance | `FaultConfig` | Network change | Time structure |
|---|---|---|---|
| **SC (bus fault)** | `fault_type = SC`, `contingency_id` from `contingencies.csv` | Augmented $Y$ at faulted bus; branch out after clearing | Fault window $[t_f, t_c]$ + post-fault $(t_c, t_{\mathrm{end}}]$ |
| **GL (gen trip)** | `gl_gen_ids` | Generator removed; admittance rebuilt on a copied case | Single window $[t_f, t_{\mathrm{end}}]$ |
| **GL (load trip)** | `gl_load_bus_ids`, `gl_percent_power` | Scaled $P_d, Q_d$ on copied `DBUS` | Single window $[t_f, t_{\mathrm{end}}]$ |

Timing (`t_start_fault`, `clearing_time`, `t_end_sim`, `t_step`) lives on `TransientConfig.simulation`, not in `contingencies.csv`.

---

## TSC-ACOPF with FULL_BUS

`network_form = FULL_BUS` keeps the full sparse $Y_{\mathrm{bus}}$ and nodal states $V_k^t$, $\theta_k^t$ at every transient step. Generator buses still carry $\delta_g$, $E_g$, and optional $P_m$; $P_e$ and $Q_e$ couple to the terminal bus through the same classical maps as $\eqref{eq:pe-init-p}$–$\eqref{eq:pe-init-q}$.

What changes relative to Kron:

- **KCL at every bus** — active and reactive nodal balance at each $t$ (not a reduced-gen-only network).
- **ZIP loads** — `zip_load_p` / `zip_load_q` split constant-Z/I/P components for active and reactive demand independently; default `(1,0,0)` on both is constant impedance.
- **Mandatory warm start** — TSCOPF solves ACOPF first and passes `SteadyStateHints` into the transient model; `FULL_BUS` requires `mech_power_mode = USE_PM`.
- **Post-fault window** — SC with branch trip uses `*_tpf` steps; GL uses a single window only.

Disturbances modify **copies** of `DBUS`/`DGEN`/`DCIR` inside the builder; CSV files on disk are untouched.

---

## TSC-DCOPF (linearised transient)

`dispatch.type_model = "DCOPF"` with `trans_stab = true` selects the linear benchmark path (Gurobi or HiGHS). Steady state is the DC model from $\eqref{eq:dcopf-balance}$ on the [steady-state page](02_steady_state_opf.md) (cross-page reference: eq. for DC balance).

Transient simplifications:

- $V = 1$ p.u. everywhere; $P_{m,g} = P_g$ during the swing.
- A reference angle $\delta_g^{\mathrm{ref}}$ from the solved DCOPF linearises electrical power:

!!! note "Model 6.8 (linearised electrical power, schematic)"
    ```math
    \begin{equation}
    \label{eq:pe-linear}
    P_g^{\mathrm{ele},t} \approx P_{g}^{\mathrm{ele,ref}} +
    \sum_{i \in \mathcal{G}} K_{gi}\left(\delta_g^t - \delta_i^t - \delta_g^{\mathrm{ref}} + \delta_i^{\mathrm{ref}}\right)
    \end{equation}
    ```

    Coefficients $K_{gi}$ come from $G^{\mathrm{red}}, B^{\mathrm{red}}$ at the reference point (`Make_Dynamic_Model_tsredlinear!`).

Initial link: $P_{g(k)} - (\delta_g^0 - \theta_{k(g)})/x'_g = 0$. COI and $\eqref{eq:delta-coi-box}$ match the Kron TSC-ACOPF structure. `FULL_BUS` and `DQ_4TH` are not implemented on this path.

---

## Worked example

!!! tip "Example 6.1 · IEEE 9-bus, fault at bus 7"
    Contingency 2 in `INPUT_FILES/9bus/contingencies.csv` applies a fault at bus 7 and trips branch 5–7 after clearing. With `δ_tol_deg = 120`, a binding upper bound on machine 2 often appears in the fault window; its multiplier is exported in `Transient_Stability/CSV_duals/` (row family `dual_δ_COI_upper`).

    ```julia
    cfg = RunConfig(
        trans_stab = true,
        case = "9bus",
        load_factor = 1.5,
        solver_name = "Ipopt",
        dispatch = DispatchConfig(type_model = "ACOPF"),
        transient = TransientConfig(
            simulation = TsSimulationConfig(δ_tol_deg = 120.0),
            dyn_model = DynModelConfig(
                fault = FaultConfig(contingency_id = 2),
            ),
        ),
    )
    sys = load_system(cfg, @__DIR__)
    res = run_case!(cfg, sys, @__DIR__, joinpath(@__DIR__, "RESULTS"))
    ```

    More presets (FULL_BUS, GL trips, δ_tol sweeps): [Running a case](../running_a_case.md).

---

## Where this lives in the code

| Model piece | Builder / file | Config knob |
|---|---|---|
| Steady-state OPF | `_acopf/`, `_dcopf/` | `RunConfig.dispatch` |
| Route Kron / FULL_BUS / linear | `functions_2_build_TS_model_common.jl` | `dyn_model.network_form`, `dispatch.type_model` |
| $\eqref{eq:pe-init-p}$–$\eqref{eq:pe-init-q}$ | `functions_4_TS_kron_eqconst.jl` | always built on AC paths; `Q` omitted on Kron linear |
| $\eqref{eq:trap-delta}$–$\eqref{eq:trap-swing}$ | `functions_4_TS_kron_eqconst.jl` | always built in tf/tpf windows |
| $\eqref{eq:pe-kron}$ | `eq_const_tsred_Pe_generic!` | always built (Kron) |
| $\eqref{eq:coi}$ | `eq_const_kron_COI_generic!` | always built |
| $\eqref{eq:delta-coi-box}$ | `functions_4_TS_kron_ineqconst.jl` | `ineq_δ_COI_tf_*`, `ineq_δ_COI_tpf_*`; `bound_style_δ` |
| FULL_BUS nodal KCL | `functions_2_build_TS_model_w_FullBus.jl` | `network_form = FULL_BUS` |
| $\eqref{eq:pe-linear}$ | `functions_2_build_TS_model_w_Kron_Linear.jl` | `dispatch.type_model = "DCOPF"` |
| SC / GL topology | `FaultConfig.jl`, fault builders | `dyn_model.fault` |
| Dual export | `DynDualRegistry.jl` | `RunConfig.save_duals = true` |

Next in the reading order: [7. Duals, KKT, and the economics](07_duals_economics.md) for sign conventions and sheet-by-sheet dual interpretation.
