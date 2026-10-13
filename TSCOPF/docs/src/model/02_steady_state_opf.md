# 2. Steady-state OPF

Before you attach swing dynamics or fault windows, every TSC run solves a **steady-state dispatch**: how much each generator produces, at what terminal voltage and angle, while respecting network physics and equipment limits. That dispatch is the economic problem market operators already recognise. TSC-OPF keeps it and adds transient feasibility on top.

This page states the steady-state models TSCOPF implements today: ACOPF (the nonlinear workhorse for TSC-ACOPF), DCOPF (linear benchmark and the steady-state block inside TSC-DCOPF), and ED (scalar-balance dispatch with explicit dual LP). The notation table below is shared with the transient pages.

| Symbol | Meaning |
|---|---|
| $\mathcal{N}$, $k,m$ | Buses |
| $\mathcal{G}$, $g$ | Generators |
| $P_g, Q_g$ | Active / reactive generator dispatch [p.u.] |
| $V_k, \theta_k$ | Bus voltage magnitude / angle [p.u., rad] |
| $G_{km}, B_{km}$ | Entries of the nodal admittance matrix |
| $P_d(k), Q_d(k)$ | Active / reactive demand at bus $k$ |
| $c_{2,g}, c_{1,g}, c_{0,g}$ | Quadratic fuel-cost coefficients |

---

## AC optimal power flow

ACOPF is the steady-state block inside TSC-ACOPF (Kron and FULL_BUS). You pick the fuel-cost curve, then minimise total generation cost subject to nonlinear power flow.

!!! note "Model 2.1 (ACOPF objective)"
    Minimise quadratic (or linear) generator fuel cost:

    ```math
    \begin{equation}
    \label{eq:acopf-objective}
    \min \sum_{g \in \mathcal{G}} \left(c_{2,g} P_g^2 + c_{1,g} P_g + c_{0,g}\right)
    \end{equation}
    ```

    Set `dispatch.cost_type = "linear"` for a linear objective; the default quadratic form matches MATPOWER `gencost` model 2.

Active and reactive power must balance at every bus. In matrix form (the default when `dispatch.use_matrix = true`):

!!! note "Model 2.2 (AC nodal balance)"
    For each bus $k \in \mathcal{N}$,

    ```math
    \begin{align}
    \label{eq:acopf-p-balance}
    P_{g(k)} - V_k \sum_{m \in \mathcal{N}} V_m \left(G_{km}\cos(\theta_k - \theta_m) + B_{km}\sin(\theta_k - \theta_m)\right) &= P_d(k), \\
    \label{eq:acopf-q-balance}
    Q_{g(k)} - V_k \sum_{m \in \mathcal{N}} V_m \left(G_{km}\sin(\theta_k - \theta_m) - B_{km}\cos(\theta_k - \theta_m)\right) &= Q_d(k).
    \end{align}
    ```

    $P_{g(k)}$ and $Q_{g(k)}$ are zero when no generator sits at bus $k$.

The slack bus fixes the voltage-angle gauge. Generator active/reactive limits, bus voltage magnitude limits, and optional branch thermal and angle-difference limits close the feasible set. TSCOPF codes upper and lower bounds as explicit inequalities in $(\mathrm{LHS} - \mathrm{RHS}) \leq 0$ form so JuMP duals are well defined.

| Dual | Constraint (schematic) | Role when binding |
|---|---|---|
| $\lambda_k$ | Active power balance, $\eqref{eq:acopf-p-balance}$ | Nodal energy price (LMP) |
| $\beta_k$ | Reactive power balance, $\eqref{eq:acopf-q-balance}$ | Reactive support scarcity |
| $\zeta$ | $\theta^{\mathrm{ref}} = 0$ | Angle reference (gauge) |
| $\eta_g^-,\eta_g^+$ | $P_g^{\min} \leq P_g \leq P_g^{\max}$ | Must-run subsidy / capacity rent |
| $\kappa_g^-,\kappa_g^+$ | $Q_g^{\min} \leq Q_g \leq Q_g^{\max}$ | Reactive capability rent |
| $\alpha_k^-,\alpha_k^+$ | $V_k^{\min} \leq V_k \leq V_k^{\max}$ | Voltage support scarcity |
| branch duals | Thermal / angle-difference limits | Congestion rent |

!!! info "Assumption"
    Branch flow limits can be enforced either through explicit branch-flow variables or through matrix balance plus angle-difference inequalities, depending on `DispatchConfig` toggles. The economic reading is the same: a binding limit raises the marginal cost of serving load downstream.

!!! note "Interpretation (LMP from AC balance)"
    Stationarity with respect to $P_g$ at bus $k(g)$ ties the marginal fuel cost to the balance multiplier. TSCOPF exports $\lambda_k$ from the primal, and the **economic price is $\pi_k = +\lambda_k$** — no flip. The balance is coded $P_{g(k)} - P_d(k) - \sum(\cdot) = 0$ and JuMP's Lagrangian carries $-\lambda$, so the two signs cancel. This holds for both steady-state and transient primal exports. The **explicit dual LP** (`Dispatch_Dual/`) is the one object with the opposite convention; see [7. Duals, KKT, and the economics](07_duals_economics.md).

!!! tip "Example 2.1 · IEEE 9-bus, plain ACOPF"
    A dispatch-only run with `trans_stab = false` and `dispatch.type_model = "ACOPF"` solves Model 2.1–2.2 on the 9-bus case at `load_factor = 1.5`. Ipopt typically converges in a few seconds; primal dispatch lands in `RESULTS/.../Dispatch/`.

    ```julia
    cfg = RunConfig(
        trans_stab = false,
        case = "9bus",
        load_factor = 1.5,
        solver_name = "Ipopt",
        dispatch = DispatchConfig(type_model = "ACOPF"),
    )
    sys = load_system(cfg, @__DIR__)
    res = run_case!(cfg, sys, @__DIR__, joinpath(@__DIR__, "RESULTS"))
    ```

---

## DC optimal power flow

DCOPF drops reactive power and voltage magnitude, linearises the network around small angles, and solves an LP. It is the steady-state block for TSC-DCOPF and the cleanest place to cross-check dual signs against an explicit dual LP.

!!! note "Model 2.3 (DCOPF balance)"
    With $B_{km}$ the branch susceptance (see below),

    ```math
    \begin{equation}
    \label{eq:dcopf-balance}
    P_{g(k)} - \sum_{m \in \mathcal{N}} B_{km}(\theta_k - \theta_m) = P_d(k), \qquad \forall k \in \mathcal{N}.
    \end{equation}
    ```

!!! info "Assumption"
    Two susceptance conventions exist, selected by `dispatch.susceptance_model` (DC-OPF only):

    - `SIMPLE`: $B_{km} = 1/x_{km}$ (textbook DC, ignores series $r$)
    - `POWERMODELS`: $B_{km} = \operatorname{Im}\!\big(1/(r_{km}+jx_{km})\big) = x_{km}/(r_{km}^2+x_{km}^2)$

    They coincide when $r_{km}=0$. With `POWERMODELS`, in-house DC-OPF matches PowerModels `DCPPowerModel` on primal quantities (see `test/runtests_powermodels_crosscheck.jl`).

Angle limits and optional branch angle-difference limits replace thermal constraints. The same LMP convention applies: $\pi_k = +\lambda_k$ on $\eqref{eq:dcopf-balance}$ for the primal duals.

---

## Economic dispatch

ED collapses the network to a single active-power balance (no angles, no branches). It is useful for sanity-checking the dual LP machinery before you open the DC or AC models.

!!! note "Model 2.4 (ED)"
    ```math
    \begin{align}
    \label{eq:ed-objective}
    \min \;& \sum_g c_g P_g \\
    \text{s.t.}\;& \sum_g P_g = \sum_k P_d(k), \\
    & P_g^{\min} \leq P_g \leq P_g^{\max}.
    \end{align}
    ```

    Requires `cost_type = "linear"`. A scalar balance multiplier $\lambda$ gives the system marginal price $\pi_{\mathrm{SMP}} = +\lambda$ from the primal duals. The explicit ED dual LP reports the same price from its own multiplier, which carries the opposite sign — see [7. Duals, KKT, and the economics](07_duals_economics.md).

---

## Where this lives in the code

| Equation / object | Builder / module | Config knob |
|---|---|---|
| $\eqref{eq:acopf-objective}$–$\eqref{eq:acopf-q-balance}$ | `_acopf/` primal builder | `dispatch.type_model = "ACOPF"` |
| $\eqref{eq:dcopf-balance}$ | `_dcopf/` | `dispatch.type_model = "DCOPF"`, `susceptance_model` |
| $\eqref{eq:ed-objective}$ | `_ed/` | `dispatch.type_model = "ED"`, `cost_type = "linear"` |
| Explicit dual LP (ED / DC) | `_ed/`, `_dcopf/` dual builders | `dispatch.solve_explicit_dual = true` |
| JuMP primal dual export | `DispatchDualRegistry.jl` | `RunConfig.save_duals = true` |

Field-level descriptions: [Parameter reference](../parameter_reference.md). The coupled steady-state + transient model is assembled on [6. The TSC-OPF, assembled](06_tsc_opf_assembled.md).
