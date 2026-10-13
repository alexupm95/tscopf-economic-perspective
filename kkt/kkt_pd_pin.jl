# =====================================================================================
#  Tier 1-P — ∂L/∂P_d(k) read straight off a pinned demand variable
# =====================================================================================
#
#  The demand-side price is
#
#      π_k = ∂f/∂P_d(k) − Σ_r λ_r ∂h_r/∂P_d(k) − Σ_s μ_s ∂g_s/∂P_d(k)
#
#  and the difficulty has never been the formula, it is obtaining ∂c_r/∂P_d for EVERY
#  row without first knowing which rows contain the demand. A hand-maintained list of
#  "families that carry P_d" is exactly the thing that goes stale when a builder changes.
#
#  `DispatchConfig.pd_as_variable` removes the problem at the source. With the demand
#  promoted to a variable pinned by `P_d[k] − p_d⁰[k] == 0 : σ_k`, stationarity of the
#  Lagrangian in `P_d[k]` reads `π_k − σ_k = 0`, so
#
#      σ_k = π_k,
#
#  summed by the solver over whatever rows happen to exist. `dual_Pd_pin.csv` already
#  carries that number on every run with the flag on, and needs nothing from this file.
#
#  What this file adds is the LEDGER: which rows the demand actually reaches, and what
#  each contributes. `P_d` is now a real Jacobian column, so those coefficients come
#  from MOI's exact derivative of the real model — analytic, one solve, no ε, and no
#  family list anywhere. Rows are labelled by `label_jacobian`, which walks the model's
#  own containers, so a family nobody anticipated still arrives named rather than
#  UNMAPPED.
#
#  RELATION TO param_jacobian.jl. That file answers the same question by rebuilding the
#  model at `p_d ± ε` and differencing, at three solves per bus. It remains the ONLY
#  route where the demand is absorbed into a numeric matrix before JuMP sees it — the
#  Kron `Y_red` constant-impedance shunt. Where the demand survives as an inlined
#  coefficient (the ACOPF balance; the FULL_BUS ZIP right-hand side) the pin is exact
#  and covers every bus in a single solve.
#
#  THE PIN ROW IS EXCLUDED from the ledger, deliberately. It is the instrument, not part
#  of the model being priced. Folding σ_k back in returns the stationarity condition
#  itself — a number that is ≈ 0 by construction and says nothing about the price. It is
#  still worth reporting, as a free check on the solver's dual convergence, and this file
#  writes it as `STATIONARITY`.
#
#  Runs inside `RunConfig.post_solve_hook`: `release_solver_backend!` empties the model
#  immediately afterwards, so no ConstraintRef survives the return.
#
#  Requires `kkt_sensitivity.jl` to be included first (for `label_jacobian`), the same
#  way `param_jacobian.jl` does.
# =====================================================================================

using JuMP
using CSV, DataFrames, DataStructures, SparseArrays
using TSCOPF
import MathOptInterface as MOI
import TSCOPF: build_constraint_evaluator

"""
    Save_KKT_Pd_Pin!(path_names, opf_dict, dyn_model_dict, model;
                     tol = 1e-12, rtol_check = 1e-6) -> String

Write the per-row ledger of `π_k = ∂L/∂P_d(k)` for every pinned bus to
`Dispatch/CSV_duals/kkt_Pd_pin.csv`, and check each bus total against the pin's own dual.

Requires the model to have been built with `DispatchConfig.pd_as_variable = true`;
returns `""` with a warning otherwise, so a driver can attach the hook unconditionally.

Columns (the Tier 1-D schema, so `KKT_Pd_Duals.ipynb` reads this file unchanged):

    bus_id;param;source;family;store;window;t_index;dual;coeff;contribution

`coeff` is `∂c_r/∂P_d(k)` from the exact Jacobian; `contribution` is `−dual · coeff`.
Each bus also gets `ENERGY` (dispatch rows), `STABILITY` (dynamic rows), `OTHER`,
`TOTAL`, `SIGMA_PIN` (the solver's own answer) and `STATIONARITY` summary rows.
"""
function Save_KKT_Pd_Pin!(path_names, opf_dict, dyn_model_dict, model;
                          tol::Float64 = 1e-12, rtol_check::Float64 = 1e-6)::String

    if !haskey(opf_dict[:vars], :P_d)
        @warn "Save_KKT_Pd_Pin!: no P_d variable in the model — set " *
              "DispatchConfig.pd_as_variable = true to enable the demand pin. Skipping."
        return ""
    end

    evaluator, rows, xvars = build_constraint_evaluator(model; features = [:Jac])
    row_of = Dict(rows[i] => i for i in eachindex(rows))
    col_of = Dict(xvars[j] => j for j in eachindex(xvars))
    rows_df, _cols_df = label_jacobian(opf_dict, dyn_model_dict, rows, xvars, row_of, col_of)

    # One pass over the COO triplets, bucketed BY COLUMN. Never assemble a
    # SparseMatrixCSC and slice rows from it: CSC stores columns, so a per-row view
    # rescans the whole structure and the loop degenerates to O(n_rows · nnz).
    sparsity = MOI.jacobian_structure(evaluator)
    vals = zeros(length(sparsity))
    MOI.eval_constraint_jacobian(evaluator, vals, JuMP.value.(xvars))

    by_col = Dict{Int, Vector{Tuple{Int, Float64}}}()
    for (t, (r, c)) in enumerate(sparsity)
        abs(vals[t]) > tol || continue
        push!(get!(by_col, c, Tuple{Int, Float64}[]), (r, vals[t]))
    end

    pins = opf_dict[:eq_const][:eq_const_pd_pin]
    out  = NamedTuple[]

    println("\n[KKT] demand pin — π_k = ∂L/∂P_d(k), from the exact Jacobian column")

    for bus in sort(collect(keys(opf_dict[:vars][:P_d])))
        var     = opf_dict[:vars][:P_d][bus]
        j       = col_of[var]
        pin_row = row_of[pins[bus]]
        σ       = JuMP.dual(pins[bus])

        bus_rows = NamedTuple[]
        stationarity = 0.0

        for (r, coeff) in sort(get(by_col, j, Tuple{Int, Float64}[]); by = first)
            dual = rows_df.dual[r]
            stationarity -= dual * coeff
            r == pin_row && continue

            src    = rows_df.source[r]
            family = rows_df.family[r]
            store  = src == "dispatch" ? "opf" : src == "dynamic" ? "dyn" : src
            window = endswith(family, "_tpf") ? "tpf" :
                     endswith(family, "_tf")  ? "tf"  : ""
            push!(bus_rows, (bus_id = bus, param = "P_d", source = src, family = family,
                             store = store, window = window,
                             t_index = rows_df.key2[r] < 0 ? missing : rows_df.key2[r],
                             dual = dual, coeff = coeff, contribution = -dual * coeff))
        end

        energy    = sum(x.contribution for x in bus_rows if x.store == "opf"; init = 0.0)
        stability = sum(x.contribution for x in bus_rows if x.store == "dyn"; init = 0.0)
        other     = sum(x.contribution for x in bus_rows
                        if !(x.store in ("opf", "dyn")); init = 0.0)
        total     = energy + stability + other

        append!(out, bus_rows)
        for (fam, val) in (("ENERGY", energy), ("STABILITY", stability),
                           ("OTHER", other), ("TOTAL", total),
                           ("SIGMA_PIN", σ), ("STATIONARITY", stationarity))
            push!(out, (bus_id = bus, param = "P_d", source = "check", family = fam,
                        store = "check", window = "", t_index = missing,
                        dual = missing, coeff = missing, contribution = val))
        end

        # The ledger and the pin's dual are the same quantity by two routes. A mismatch
        # means the Jacobian column and the solver disagree about what P_d touches —
        # which is the one failure this whole construction is meant to make impossible.
        if !isapprox(total, σ; rtol = rtol_check, atol = rtol_check)
            @warn "bus $bus: ledger $total does not match the pin dual $σ"
        end

        n_dyn = count(x -> x.store == "dyn", bus_rows)
        println("  bus $bus  π = $(round(total, digits = 6)) €/MW  = energy " *
                "$(round(energy, digits = 6)) + stability $(round(stability, digits = 6))" *
                "   [$(length(bus_rows)) rows, $(n_dyn) dynamic]  " *
                "|π − σ| = $(round(abs(total - σ), sigdigits = 3))")
    end

    out_dir = path_names[:pf_dispatch_CSV_duals]
    mkpath(out_dir)
    file_path = joinpath(out_dir, "kkt_Pd_pin.csv")
    CSV.write(file_path, DataFrame(out); delim = ';')
    println("[KKT] demand-pin ledger written to $file_path\n")
    return file_path
end
