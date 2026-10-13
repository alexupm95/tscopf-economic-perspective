# TSCOPF.jl

[![Docs](https://img.shields.io/badge/docs-dev-blue.svg)](https://alexupm95.github.io/tscopf-economic-perspective/dev/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://github.com/alexupm95/tscopf-economic-perspective/blob/main/TSCOPF/LICENSE)
![Julia >= 1.10](https://img.shields.io/badge/julia-%E2%89%A5%201.10-9558B2.svg)

**Transient-Stability-Constrained Optimal Power Flow** (TSC-OPF) in Julia/JuMP.

TSCOPF builds and solves steady-state OPF (ED, DCOPF, ACOPF) and TSC formulations with Kron-reduced or **full-bus admittance** dynamics, short-circuit and generator/load trip contingencies, and structured dual export for economic interpretation.

## What it is for

Two capabilities that most OPF codes keep separate:

1. **Model fidelity as a dial.** The same contingency runs through a Kron-reduced classical machine model, a full nodal (`FULL_BUS`) network, or fourth-order dq machines with optional AVR and TGOV1 governor, so the effect of dynamic-model detail on the dispatch is measurable rather than assumed. A linearised TSC-DCOPF is kept as a benchmark, not as the target formulation.
2. **Duals as first-class output.** Every constraint family — steady-state and transient — exports its dual through a registry, so shadow prices on the rotor-angle and speed bounds arrive alongside the usual nodal prices and the cost of a binding stability limit can be read off directly.

## At a glance

| | |
|---|---|
| Steady-state formulations | ED, DC-OPF, AC-OPF (quadratic or linear cost) |
| Transient network models | Kron-reduced classical 2nd-order; full sparse-`Ybus` nodal (`FULL_BUS`) |
| Disturbances | Three-phase bus short-circuit + branch trip; generator or load disconnection |
| Dual export | Registry-driven, for both steady-state and transient constraint families |
| Solvers | Ipopt or MadNLP (NLP); Gurobi or HiGHS (LP/MILP) |

## Where to start

**Read the model first** if you know power systems but haven't seen this formulation:

| Page | Content |
|------|---------|
| [2. Steady-state OPF](model/02_steady_state_opf.md) | AC / DC / ED objectives and balance equations |
| [6. The TSC-OPF, assembled](model/06_tsc_opf_assembled.md) | Coupled dispatch + transient model |
| [7. Duals, KKT, and the economics](model/07_duals_economics.md) | Sign convention and dual export |
| [8. Machine controls](model/08_controls_avr_governor.md) | AVR and TGOV1 governor: equations, discretisation, limiters |
| [9. Grid-forming inverters](model/09_grid_forming.md) | Droop, Q–V PI, current limiter, and what stays SG-only |

**Run the code:**

| Page | Content |
|------|---------|
| [Install](install.md) | Julia environment and dependencies |
| [Quick start](quickstart.md) | The runs of the paper, and a single run via `run_case!` |
| [Running a case](running_a_case.md) | Example configs, outputs, validation |
| [Parameter reference](parameter_reference.md) | Every `RunConfig` knob explained |
| [Complete example](complete_example.md) | One run with **every** field of every config struct set explicitly |
| [Configuration map](configuration_map.md) | Struct field reference for developers |
| [Developer architecture](developer_architecture.md) | Source layout for transient-stability builders |
| [API](api.md) | Docstrings from the `TSCOPF` module |

!!! tip "New here?"
    Skim [2. Steady-state OPF](model/02_steady_state_opf.md), then [Install](install.md) and [Quick start](quickstart.md). Search [Parameter reference](parameter_reference.md) when you need a specific field.

## Conventions

!!! note "Cases and scripts named in the manuals"
    The examples in these pages use the demonstration cases of the package's
    development tree (`9bus`, `9bus_gfm`, …) and refer to its entry scripts
    (`main.jl`, `main_loop.jl`), examples and test suite. This repository does not
    ship them. It ships the three cases of the paper under `INPUT_FILES/` and their
    run scripts under `runs/`: see [Quick start](quickstart.md), and put one of
    those case names where an example says `case = "9bus"`.

## Repository

[github.com/alexupm95/tscopf-economic-perspective](https://github.com/alexupm95/tscopf-economic-perspective). The package is the `TSCOPF/` folder; the
repository also holds the cases, runs and post-processing of the paper
"Transient Stability Constrained Optimal Power Flow: An Economic Perspective".

## Citing

Cite the repository. Machine-readable metadata lives in
[`CITATION.cff`](https://github.com/alexupm95/tscopf-economic-perspective/blob/main/TSCOPF/CITATION.cff)
in `TSCOPF/`, which GitHub renders as a *Cite this repository* button.

**Author:** Alex Junior da Cunha Coelho — Technical University of Madrid (UPM)
