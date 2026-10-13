# =====================================================================================
#  Tier 1-D — the demand-side decomposition of ∂L/∂P_d
# =====================================================================================
#
#  The mirror of `kkt_pg_terms.jl`, with one structural difference that has to be stated
#  up front because it changes what the file can prove.
#
#  On the generator side `P_g` is a decision variable, so stationarity holds and
#  `Σ contribution + dC/dP_g = 0` is a free correctness check: the residual must vanish.
#
#  **`P_d` is a parameter, not a variable. There is no ∂L/∂P_d = 0 to satisfy.** This is a
#  VALUE statement, not a stationarity statement:
#
#      π_k = ∂L/∂P_d,k = Σ_rows ( −dual_row × ∂expr_row/∂P_d,k )
#
#  and by the envelope theorem (Fiacco), evaluated at the optimal primal-dual point with
#  x* and the multipliers HELD FIXED, that equals dJ*/dP_d,k — the nodal price. Nothing
#  here is required to sum to zero, so correctness has to come from outside: the
#  reduction check against the exported LMP (below), and a finite difference on the
#  objective.
#
#  The whole point of the exercise: in a conventional OPF the balance row is the only
#  place load appears, so π_k collapses to the balance dual. In a TSC-OPF it does not,
#  and the extra rows are the stability content of the price.
#
#  Runs entirely offline from the CSVs `param_jacobian.jl` produced — it never touches a
#  JuMP model, so the dead-model trap does not apply.
# =====================================================================================

using CSV, DataFrames

"""
    Save_KKT_Pd_Terms!(out_dir, col, bus; lambda_csv = nothing) -> String

Turn one parametric column (the output of `param_column`) into the Tier 1-D ledger for
bus `bus`, and write `kkt_Pd_terms_bus<bus>.csv`.

Columns:

    bus_id;param;source;family;store;window;t_index;dual;coeff;contribution

`store` partitions the price the way the generator-side file does — off a column the
exporter wrote, never off a list of family names, because that is precisely how a stale
name silently survived a rename on the `P_g` side:

    "opf"  → the steady-state balance  → the incumbent LMP
    "dyn"  → every transient row       → the stability content

Two synthetic rows close the file: `ENERGY`/`STABILITY` subtotals and a `TOTAL` holding
`π_k`. There is deliberately no `RESIDUAL` row — see the header.

`lambda_csv` optionally points at `Dispatch/CSV_duals/dual_P_balance.csv`; when given, the
reduction check `π_k − λ_k = stability` is printed, which is the cheapest available
confirmation that the sign convention is right.
"""
function Save_KKT_Pd_Terms!(out_dir::String, col::DataFrame, bus::Int;
                            param::String = "P_d", lambda_csv = nothing)::String

    rows = NamedTuple[]
    for r in eachrow(col)
        store  = r.source == "dispatch" ? "opf" : "dyn"
        window = endswith(r.family, "_tpf") ? "tpf" :
                 endswith(r.family, "_tf")  ? "tf"  : ""
        t_index = r.key2 < 0 ? missing : r.key2
        # The sign convention, in one line. JuMP's Lagrangian is
        #   L = f − Σ dual·(LHS − RHS),  every row stored as expr {==,≤} 0,
        # so a row contributes −dual × ∂expr/∂p. Identical to the P_g rule; the only
        # thing that changed is which parameter the derivative is taken with respect to.
        push!(rows, (bus_id = bus, param = param, source = r.source, family = r.family,
                     store = store, window = window, t_index = t_index,
                     dual = r.dual, coeff = r.dcdp, contribution = -r.dual * r.dcdp))
    end

    energy    = sum(x.contribution for x in rows if x.store == "opf"; init = 0.0)
    stability = sum(x.contribution for x in rows if x.store == "dyn"; init = 0.0)
    total     = energy + stability

    for (fam, val) in (("ENERGY", energy), ("STABILITY", stability), ("TOTAL", total))
        push!(rows, (bus_id = bus, param = param, source = "check", family = fam,
                     store = "check", window = "", t_index = missing,
                     dual = missing, coeff = missing, contribution = val))
    end

    file_path = joinpath(out_dir, "kkt_Pd_terms_bus$(bus).csv")
    CSV.write(file_path, DataFrame(rows); delim = ';')

    n_dyn = count(x -> x.store == "dyn", rows)
    println("  bus $bus  π = $(round(total, digits = 6)) €/MW  " *
            "= energy $(round(energy, digits = 6)) + stability $(round(stability, digits = 6))" *
            "   [$(n_dyn) dynamic rows]")
    if stability != 0.0
        println("        stability share = $(round(100 * stability / total, digits = 3)) %")
    end

    # Reduction check: the energy channel must BE the exported LMP, to solver precision.
    # If this fails, the sign convention is wrong and every number above is wrong with it.
    if lambda_csv !== nothing && isfile(lambda_csv)
        lam = CSV.read(lambda_csv, DataFrame; delim = ';')
        idx = findfirst(==(bus), lam[!, 1])
        if idx !== nothing
            λ = lam[idx, 2]
            println("        exported λ_$bus = $(round(λ, digits = 6)) €/MW   " *
                    "|energy − λ| = $(round(abs(energy - λ), sigdigits = 3))")
        end
    end
    return file_path
end
