# 4. The network: Kron vs FULL_BUS

The transient network model answers one question: how does electrical torque on each machine depend on rotor angles (and, on the nodal path, bus voltages) during and after a fault? TSCOPF offers two answers. **Kron reduction** collapses loads and passive buses into a small $n_{\mathrm{gen}} \times n_{\mathrm{gen}}$ admittance among generator internal nodes. **FULL_BUS** keeps every network bus, enforces KCL at each time step, and couples generators to their terminal buses explicitly.

This page gives the only full network derivation in Part I (Schur complement / Kron). The assembled optimisation problem is on [6. The TSC-OPF, assembled](06_tsc_opf_assembled.md). Electrical power during the swing uses Model 6.5 there for Kron and nodal balances for FULL_BUS.

| Symbol | Meaning |
|---|---|
| $Y_{\mathrm{bus}}$ | $n_{\mathrm{bus}} \times n_{\mathrm{bus}}$ sparse nodal admittance |
| $Y^{\mathrm{aug}}$ | Augmented matrix with generator internal nodes |
| $Y^{\mathrm{red}}$ | Kron-reduced admittance among generator internals |
| $G^{\mathrm{red}}, B^{\mathrm{red}}$ | $\operatorname{Re}(Y^{\mathrm{red}})$, $\operatorname{Im}(Y^{\mathrm{red}})$ |
| $x'_g$ | Transient reactance (`Xd_tr` in `gen_dynamic_data.csv`) |

Select the path with `DynModelConfig.network_form`: `KRON_REDUCED` (default) or `FULL_BUS`.

---

## Augmented network before reduction

Classical machines are modelled as internal voltage $E_g\angle\delta_g$ behind transient reactance $x'_g$. TSCOPF stamps each generator as a terminal bus link plus an **internal node** (the $\delta_g$ bus in the reduced picture).

Loads enter as constant admittances on the network diagonal (the same $P_d + jQ_d$ scaling used when building fault-on and post-fault $Y$). For a short-circuit, a large shunt is added at the faulted bus; after clearing, the tripped branch is removed and $Y_{\mathrm{bus}}$ is rebuilt.

!!! note "Model 4.1 (generator internal stamp, schematic)"
    For each in-service generator $g$, augment $Y_{\mathrm{bus}}$ with one internal row/column and a $2 \times 2$ admittance block coupling terminal bus $k(g)$ to the internal node:

    ```math
    \begin{equation}
    \label{eq:aug-stamp}
    y_{g} = \frac{1}{\mathrm{j}\, x'_{g}},
    \qquad
    Y^{\mathrm{aug}} =
    \begin{bmatrix}
      Y_{\mathrm{bus}} & Y_{12} \\
      Y_{21} & Y_{22}
    \end{bmatrix}
    \in \mathbb{C}^{(n_{\mathrm{bus}}+n_{\mathrm{gen}})\times(n_{\mathrm{bus}}+n_{\mathrm{gen}})} .
    \end{equation}
    ```

    In code, `augment_with_gen_internals_Kron` builds this sparse augmented matrix without densifying $Y_{\mathrm{bus}}$ (`functions_4_admittance_matrices.jl`).

Partition nodes into **network** indices $\{1,\ldots,n_{\mathrm{bus}}\}$ and **generator internal** indices $\{n_{\mathrm{bus}}+1,\ldots,n_{\mathrm{bus}}+n_{\mathrm{gen}}\}$.

---

## Kron reduction (Schur complement)

Linear network theory eliminates passive nodes algebraically. Solve $I_{\mathrm{net}} = Y_{11} V_{\mathrm{net}} + Y_{12} V_{\mathrm{gen}}$ for $V_{\mathrm{net}}$ given generator internal voltages, substitute into the generator-block equations, and you obtain an equivalent $n_{\mathrm{gen}} \times n_{\mathrm{gen}}$ admittance.

!!! note "Model 4.2 (Kron reduction)"
    Partition $Y^{\mathrm{aug}}$ into blocks $Y_{11}, Y_{12}, Y_{21}, Y_{22}$ on network vs generator nodes. The Kron-reduced matrix is the Schur complement

    ```math
    \begin{equation}
    \label{eq:kron-schur}
    Y^{\mathrm{red}} = Y_{22} - Y_{21}\, Y_{11}^{-1}\, Y_{12} .
    \end{equation}
    ```

    `Reduce_Matrix` implements $\eqref{eq:kron-schur}$ with sparse LU on $Y_{11}$ and a thin solve for $Y_{11}^{-1} Y_{12}$, so memory stays $O(\mathrm{nnz} + n_{\mathrm{bus}}\, n_{\mathrm{gen}})$ rather than $O(n_{\mathrm{bus}}^2)$.

During the transient, electrical power on the Kron path is Model 6.5 on [page 6](06_tsc_opf_assembled.md): a sum over $G^{\mathrm{red}}_{gi}$, $B^{\mathrm{red}}_{gi}$ and pairwise angle differences $\delta_g^t - \delta_i^t$.

### When is Kron appropriate?

Kron is a good fit when:

- Every dynamic machine sits behind a classical $E'\angle\delta$ model and loads are representable as fixed admittances on the pre-fault network.
- You want the smallest NLP: only generator $\delta_g^t$, $\Delta\omega_g^t$, and $P_e$ trajectories, not $V_k^t$, $\theta_k^t$ at every load bus.
- Faults change $Y_{\mathrm{bus}}$ but the **same reduction procedure** applies to fault-on and post-fault topologies (SC augments the fault bus; post-fault rebuilds `DCIR` without the cleared line).

Kron is **not** appropriate when you need ZIP load voltage dependence, explicit reactive balance at each bus during the swing, or dq machine models that require terminal $V$, $\theta$ at every step. Those paths use FULL_BUS.

!!! info "Assumption"
    On the Kron path, reactive power during the swing does not enter the swing equation. $Q_e$ can be reconstructed post-solve as a diagnostic (`expressions[:Qe_tf]` on Kron TSC-ACOPF) but is not a primal constraint.

---

## FULL_BUS nodal path

`network_form = FULL_BUS` skips $\eqref{eq:kron-schur}$. The builder keeps the full sparse $Y_{\mathrm{bus}}$ at each time step and introduces $V_k^t$, $\theta_k^t$ for **every** bus in the fault and post-fault windows.

!!! note "Model 4.3 (nodal active-power balance, schematic)"
    At each time $t$ and bus $k$,

    ```math
    \begin{equation}
    \label{eq:fullbus-p-balance}
    P^{\mathrm{inj}}_k(t) - V_k^t \sum_{m} V_m^t \left( G_{km}\cos(\theta_k^t - \theta_m^t) + B_{km}\sin(\theta_k^t - \theta_m^t) \right) = 0
    \end{equation}
    ```

    with a matching reactive balance row. Generator injections $P^{\mathrm{inj}}, Q^{\mathrm{inj}}$ couple through the classical terminal maps (Models 6.1–6.2 on [page 6](06_tsc_opf_assembled.md)). Loads use the ZIP splits `zip_load_p = (a,b,c)` and `zip_load_q = (a,b,c)` on `DynModelConfig`: constant-Z, constant-I, and constant-P fractions that must each sum to 1. The active and reactive splits are **independent** — e.g. the Spanish TSO models active demand as constant current (`(0,1,0)`) and reactive demand as constant admittance (`(1,0,0)`).

Mandatory differences from Kron:

| Topic | Kron (`KRON_REDUCED`) | FULL_BUS |
|---|---|---|
| Transient states | $\delta_g$, $\Delta\omega_g$, $P_e$ (and optional $P_m$) | Above **plus** $V_k$, $\theta_k$ at all buses |
| Network matrix | $Y^{\mathrm{red}}$ among generators | Full $Y_{\mathrm{bus}}$ per fault topology |
| Mechanical power | `USE_PG` default ($P_m = P_g$) | `USE_PM` required (explicit $P_m$) |
| Warm start | Optional hints from dispatch | **Mandatory** ACOPF solve → `SteadyStateHints` |
| GL disturbance | Rebuild $Y^{\mathrm{red}}$ on copied case | Rebuild full $Y_{\mathrm{bus}}$; scale loads on copied `DBUS` |
| TSC-DCOPF | Supported (linearised $P_e$) | **Not implemented** |
| DQ 4th-order | Not on Kron path | `DQ_4TH` + `FULL_BUS` only |

Disturbances always modify **copies** of `DBUS`/`DGEN`/`DCIR` inside the builder; files under `INPUT_FILES/` are unchanged.

---

## Fault-dependent topologies

Both paths rebuild admittance when the network changes.

| Phase | SC (bus fault) | GL (gen / load trip) |
|---|---|---|
| Pre-fault | Base $Y_{\mathrm{bus}}$ + gen internals | Same |
| Fault-on | Augmented $Y$ with fault shunt at `faulted_bus` | Trip gens (`g_status`) or scale $P_d,Q_d$ on `DBUS` copy |
| Post-fault | Rebuild without cleared branch (`circuit_id` in `contingencies.csv`) | *(no separate post-fault window)* |

Kron: `Calculate_Ybus_fault_SC_Kron` → `Reduce_Matrix` → $G^{\mathrm{red}}, B^{\mathrm{red}}$ for the fault window; `Calculate_Ybus_postf_ClearFault_Kron` for post-fault. FULL_BUS: parallel routines in `functions_4_admittance_matrices.jl` / full-bus builders without Schur complement.

Exported matrices land in `RESULTS/.../Bus_Matrices/` when `trans_stab = true`.

---

## Side-by-side on IEEE 9-bus

!!! tip "Example 4.1 · same contingency, two network forms"
    Contingency 2 (fault at bus 7, trip branch 5–7) on the 9-bus case:

    **Kron (default)** — classical preset, no extra flags:

    ```julia
    cfg_kron = RunConfig(
        trans_stab = true, case = "9bus", solver_name = "Ipopt",
        dispatch = DispatchConfig(type_model = "ACOPF"),
        transient = TransientConfig(
            simulation = TsSimulationConfig(δ_tol_deg = 100.0),
            dyn_model = DynModelConfig(fault = FaultConfig(contingency_id = 2)),
        ),
    )
    ```

    **FULL_BUS** — requires `USE_PM`, `:coi_box`, and accepts ZIP loads:

    ```julia
    cfg_full = RunConfig(
        trans_stab = true, case = "9bus", solver_name = "Ipopt",
        dispatch = DispatchConfig(type_model = "ACOPF"),
        transient = TransientConfig(
            simulation = TsSimulationConfig(δ_tol_deg = 100.0),
            dyn_model = DynModelConfig(
                network_form = FULL_BUS,
                mech_power_mode = USE_PM,
                bound_style_δ = :coi_box,
                zip_load_p = (1.0, 0.0, 0.0),  # (Z, I, P) — active demand, constant impedance
                zip_load_q = (1.0, 0.0, 0.0),  # (Z, I, P) — reactive demand, constant impedance
                fault = FaultConfig(contingency_id = 2),
            ),
        ),
    )
    ```

    Compare `Transient_Stability/CSV/` trajectories: Kron writes generator-centred angles; FULL_BUS adds `bus_*.csv` and `generator_*.csv` with nodal $V$, $\theta$. Primal objectives and dispatch can differ when ZIP load matters; on constant-Z defaults they should be close but not identical.

---

## Where this lives in the code

| Step | Function / module | Notes |
|---|---|---|
| Augment $Y$ with gen internals | `augment_with_gen_internals_Kron` | `functions_4_admittance_matrices.jl` |
| $\eqref{eq:kron-schur}$ | `Reduce_Matrix` | Sparse LU on $Y_{11}$ |
| SC / post-fault Kron $Y$ | `Calculate_Ybus_fault_SC_Kron`, `Calculate_Ybus_postf_ClearFault_Kron` | Fault shunt + branch trip |
| GL Kron $Y$ | `Calculate_Ybus_fault_LG_Kron` | After gen/load scaling on copy |
| Kron TS builder | `Make_Dynamic_Model_tsred!` | `functions_2_build_TS_model_w_Kron.jl` |
| FULL_BUS TS builder | `Make_Dynamic_Model_fullbus!` | `functions_2_build_TS_model_w_FullBus.jl` |
| Network form knob | `DynModelConfig.network_form` | `KRON_REDUCED` \| `FULL_BUS` |

Swing equations and time discretisation: [5. Swing dynamics and discretization](05_swing_dynamics.md). Generator internals: [3. Generator models](03_generator_models.md).
