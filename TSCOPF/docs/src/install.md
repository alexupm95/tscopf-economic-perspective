# Install

Julia **≥ 1.10** is required (`Project.toml` `[compat]`). The committed `Manifest.toml` was
resolved with Julia 1.11 and pins the versions used for the runs of the paper.

!!! note "Research code, not a registered package"
    TSCOPF.jl is not on the General registry: `Pkg.add("TSCOPF")` will not work. The
    package is the `TSCOPF/` folder of the repository.

## From a clone of the repository

```bash
git clone https://github.com/alexupm95/tscopf-economic-perspective
cd tscopf-economic-perspective
julia --project=TSCOPF -e "using Pkg; Pkg.instantiate()"
```

Then `julia --project=TSCOPF` starts a session where `using TSCOPF` works.

## Solvers

| Problem class | Typical solver | Notes |
|---------------|----------------|-------|
| ACOPF, TSC-ACOPF | Ipopt, MadNLP | NLP |
| DCOPF, TSC-DCOPF, ED | HiGHS, Gurobi | LP |
| UC toy | Gurobi | MILP; optional |

The runs of the paper use **Ipopt**. No Gurobi license is required.

!!! warning "Solver ↔ formulation pairing is enforced, not a suggestion"
    `run_case!` calls `Check_Coherence_Input_Data` before building anything: an
    ACOPF or TSC-ACOPF run with `solver_name = "Gurobi"` or `"HiGHS"` (LP/MILP
    solvers) throws immediately, and likewise a DCOPF/ED run with Ipopt/MadNLP.
    See the solver-compatibility table in the [Parameter reference](parameter_reference.md).

## Build these docs locally

From the root of the repository:

```bash
julia TSCOPF/docs/setup_docs_env.jl          # first time only
julia --project=TSCOPF/docs --color=yes TSCOPF/docs/make.jl
```

Open `TSCOPF/docs/build/index.html` in a browser.

!!! tip
    `docs/parameter_reference.md`, `docs/running_a_case.md`, and `docs/configuration_map.md` are the
    single source of truth for those pages — `make.jl` copies them into
    `docs/src/` at build time. Edit the top-level `docs/*.md` files, not the
    copies under `docs/src/` (those are gitignored and get overwritten).
