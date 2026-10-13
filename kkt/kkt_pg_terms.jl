#=
================================================================================
 kkt_pg_terms.jl  —  KKT stationarity decomposition of ∂L/∂P_g
================================================================================
 Not part of the TSCOPF module: this file is loaded by the scripts in `runs/` and
 reaches the model through `RunConfig.post_solve_hook`.

 With every constraint in JuMP canonical form (`expr == 0` or `expr ≤ 0`):

     contribution(family) = −dual × ∂expr/∂P_g
     Σ contribution + dC/dP_g = residual ≈ 0

 The family table below is the complete list of rows containing `P_g`, with the
 derivation and source anchors in kkt/README.md.
================================================================================
=#

using JuMP
using CSV, DataFrames, DataStructures
using TSCOPF
import TSCOPF: extract_bound_manifest_duals, bound_manifest_entry_present

# --- the family table -------------------------------------------------------------
#
# `store`   :opf (steady-state dict) | :dyn (dynamic dict)
# `source`  :eq_const | :ineq_const
# `index`   :gen        container[gen]        → ConstraintRef
#           :bus        container[bus(gen)]   → ConstraintRef
#           :first_step container[gen][1]     → ConstraintRef
# `affine`  true  → ∂expr/∂P_g is read back from JuMP (immune to a builder
#                   renormalising the row)
#           false → the analytic `coeff` is used. Required for the squared
#                   capability limits: `normalized_coefficient` returns the AFFINE
#                   part of a QuadExpr, which is 0 for `P_g²`, so trusting JuMP
#                   there would silently drop the term.
# `coeff`   Float64, or (ctx, gen) -> Float64 for the state-dependent ones.

const _PG_FAMILIES = [
    (family = "p_balance", store = :opf, source = :eq_const,
     key = :eq_const_p_balance, index = :bus, affine = true, coeff = 1.0),

    (family = "pg_lower", store = :opf, source = :ineq_const,
     key = :ineq_const_pg_lower, index = :gen, affine = true, coeff = -1.0),
    (family = "pg_upper", store = :opf, source = :ineq_const,
     key = :ineq_const_pg_upper, index = :gen, affine = true, coeff = 1.0),

    # |S_g|² ≤ S_max² — quadratic in P_g.
    (family = "sg_upper", store = :opf, source = :ineq_const,
     key = :ineq_const_sg_upper, index = :gen, affine = false,
     coeff = (ctx, gen) -> 2 * ctx.pg_value[gen]),

    # GFM converter capability (allow_gfm; attach_gfm_acopf_limits!).
    (family = "gfm_Imax", store = :opf, source = :ineq_const,
     key = :ineq_const_gfm_Imax, index = :gen, affine = false,
     coeff = (ctx, gen) -> 2 * ctx.pg_value[gen]),
    (family = "gfm_Emax", store = :opf, source = :ineq_const,
     key = :ineq_const_gfm_Emax, index = :gen, affine = false,
     coeff = (ctx, gen) -> 2 * ctx.pg_value[gen] * _gfm_Xl_over_V(ctx, gen)^2),

    # Pre-fault coupling: P_g pinned to the machine's terminal injection. Present on
    # every TS path (Kron, Kron-linear, dq, GFM), always with P_g on the +1 side.
    (family = "Pe_init", store = :dyn, source = :eq_const,
     key = :eq_const_P_init, index = :gen, affine = false, coeff = 1.0),

    # P_m − P_g == 0. One family for the whole fleet: the GFM builder writes its
    # converters into the SAME :eq_const_Pm_init container as the machines
    # (functions_4_TS_gfm.jl:225), so there is no separate gfm_Pm_init key.
    (family = "Pm_init", store = :dyn, source = :eq_const,
     key = :eq_const_Pm_init, index = :gen, affine = true, coeff = -1.0),

    # GFM droop power measurement: P_meas − P_g == 0.
    (family = "gfm_Pmeas_init", store = :dyn, source = :eq_const,
     key = :eq_const_gfm_Pmeas_init, index = :gen, affine = true, coeff = -1.0),

    # The swing rows. Scanned at EVERY step rather than only the first, because how
    # many of them carry P_g depends on `mech_power_mode`:
    #
    #   USE_PM — P_m is its own variable, so P_g enters only the trapezoidal first
    #            step, as the pre-fault electrical power:
    #              … − (Δt/4H)·(2·P_m − Pe[1] − P_g) == 0
    #   USE_PG — P_mech IS P_g, so P_g sits in every row of both windows.
    #
    # These rows are affine, so the coefficient is read back from the model per row
    # and steps that do not contain P_g return exactly 0 and emit nothing. One entry
    # therefore covers both modes without branching on configuration — and covers
    # `ode_first_step = :backward_euler`, whose first row has no P_g either.
    (family = "swing_tf", store = :dyn, source = :eq_const,
     key = :eq_const_Δω_tf, index = :all_steps, affine = true, coeff = 0.0),
    (family = "swing_tpf", store = :dyn, source = :eq_const,
     key = :eq_const_Δω_tpf, index = :all_steps, affine = true, coeff = 0.0),

    # The stability corridor itself, which contains P_g on ONE of the two styles:
    #
    #   :coi_box            — a plain box on the δ variable, |δ_g − δ_COI| ≤ δ_tol.
    #                         No P_g in the row, so these coefficients are exactly 0
    #                         and nothing is emitted. This is the h_P = 0 case.
    #   :swing_propagated   — the SUBSTITUTED form: δ is replaced by its swing
    #                         expansion, so P_g appears in the corridor row and the
    #                         corridor prices active power directly (h_P ≠ 0).
    #
    # Both sides of the corridor are listed; the sign of each is read back from the
    # model. As with the swing rows, one entry covers both styles because the
    # coefficient decides.
    (family = "corridor_tf_lower", store = :dyn, source = :ineq_const,
     key = :ineq_const_δ_COI_tf_lower, index = :all_steps, affine = true, coeff = 0.0),
    (family = "corridor_tf_upper", store = :dyn, source = :ineq_const,
     key = :ineq_const_δ_COI_tf_upper, index = :all_steps, affine = true, coeff = 0.0),
    (family = "corridor_tpf_lower", store = :dyn, source = :ineq_const,
     key = :ineq_const_δ_COI_tpf_lower, index = :all_steps, affine = true, coeff = 0.0),
    (family = "corridor_tpf_upper", store = :dyn, source = :ineq_const,
     key = :ineq_const_δ_COI_tpf_upper, index = :all_steps, affine = true, coeff = 0.0),

    # The machine-referenced δ corridor (:highest_H / :ref_gen) and the two speed
    # corridors. These are built as direct boxes on the δ / Δω variables, so today
    # they carry no P_g and contribute nothing. They are listed anyway: it costs one
    # zero-coefficient probe per row, and it means a builder that ever propagates
    # dispatch into them shows up as a contribution instead of as a silent gap in
    # the decomposition.
    (family = "corridor_ref_tf_lower", store = :dyn, source = :ineq_const,
     key = :ineq_const_δ_ref_tf_lower, index = :all_steps, affine = true, coeff = 0.0),
    (family = "corridor_ref_tf_upper", store = :dyn, source = :ineq_const,
     key = :ineq_const_δ_ref_tf_upper, index = :all_steps, affine = true, coeff = 0.0),
    (family = "corridor_ref_tpf_lower", store = :dyn, source = :ineq_const,
     key = :ineq_const_δ_ref_tpf_lower, index = :all_steps, affine = true, coeff = 0.0),
    (family = "corridor_ref_tpf_upper", store = :dyn, source = :ineq_const,
     key = :ineq_const_δ_ref_tpf_upper, index = :all_steps, affine = true, coeff = 0.0),
    (family = "omega_COI_tf_lower", store = :dyn, source = :ineq_const,
     key = :ineq_const_Δω_COI_tf_lower, index = :all_steps, affine = true, coeff = 0.0),
    (family = "omega_COI_tf_upper", store = :dyn, source = :ineq_const,
     key = :ineq_const_Δω_COI_tf_upper, index = :all_steps, affine = true, coeff = 0.0),
    (family = "omega_COI_tpf_lower", store = :dyn, source = :ineq_const,
     key = :ineq_const_Δω_COI_tpf_lower, index = :all_steps, affine = true, coeff = 0.0),
    (family = "omega_COI_tpf_upper", store = :dyn, source = :ineq_const,
     key = :ineq_const_Δω_COI_tpf_upper, index = :all_steps, affine = true, coeff = 0.0),
    (family = "omega_abs_tf_lower", store = :dyn, source = :ineq_const,
     key = :ineq_const_Δω_abs_tf_lower, index = :all_steps, affine = true, coeff = 0.0),
    (family = "omega_abs_tf_upper", store = :dyn, source = :ineq_const,
     key = :ineq_const_Δω_abs_tf_upper, index = :all_steps, affine = true, coeff = 0.0),
    (family = "omega_abs_tpf_lower", store = :dyn, source = :ineq_const,
     key = :ineq_const_Δω_abs_tpf_lower, index = :all_steps, affine = true, coeff = 0.0),
    (family = "omega_abs_tpf_upper", store = :dyn, source = :ineq_const,
     key = :ineq_const_Δω_abs_tpf_upper, index = :all_steps, affine = true, coeff = 0.0),
]

const _PG_FAMILY_KEYS = Set(f.key for f in _PG_FAMILIES)

# --- helpers ----------------------------------------------------------------------

"""`Xl / V*` for GFM unit `gen` — the factor the Emax limit scales P_g by."""
function _gfm_Xl_over_V(ctx, gen::Int)::Float64
    ctx.DGFM === nothing && return 0.0
    row = findfirst(==(gen), Int.(ctx.DGFM.id))
    row === nothing && return 0.0
    return Float64(ctx.DGFM.Xl[row]) / ctx.v_value[ctx.gen_bus[gen]]
end

"""
Rows of `family` that belong to generator `gen`, as `(t_index, ConstraintRef)` pairs.

Empty when this run did not build the family. Time-independent families return a
single pair with `missing` for the step.
"""
function _family_entries(ctx, family, gen::Int)
    empty_result = Tuple{Union{Missing, Int}, JuMP.ConstraintRef}[]
    dict = family.store == :opf ? ctx.opf_dict : ctx.dyn_model_dict
    dict === nothing && return empty_result
    haskey(dict, family.source) || return empty_result
    store = dict[family.source]
    haskey(store, family.key) || return empty_result
    container = store[family.key]

    if family.index == :gen
        con = get(container, gen, nothing)
        return con === nothing ? empty_result : [(missing, con)]
    elseif family.index == :bus
        con = get(container, ctx.gen_bus[gen], nothing)
        return con === nothing ? empty_result : [(missing, con)]
    elseif family.index == :first_step
        inner = get(container, gen, nothing)
        inner === nothing && return empty_result
        con = get(inner, 1, nothing)
        return con === nothing ? empty_result : [(1, con)]
    elseif family.index == :all_steps
        inner = get(container, gen, nothing)
        inner === nothing && return empty_result
        return Tuple{Union{Missing, Int}, JuMP.ConstraintRef}[(t, c) for (t, c) in inner]
    end
    return empty_result
end

"""
`∂expr/∂P_g` for one row.

Affine families are asked of the model itself, so a builder that rewrites a
constraint in a different but equivalent form cannot put this export out of step
with it. Everything else uses the analytic entry, which the residual row polices.
"""
function _pg_coefficient(ctx, family, gen::Int, con)::Float64
    if family.affine
        try
            return JuMP.normalized_coefficient(con, ctx.P_g[gen])
        catch
            # Nonlinear row: no such method. Fall through to the analytic value.
        end
    end
    return family.coeff isa Function ? family.coeff(ctx, gen) : Float64(family.coeff)
end

"""
Duals of a `gen`-indexed bound family that was built with `bound_encoding=VARIABLE`.

Those runs carry no constraint rows — the bound lives on the variable and its dual
comes through the bound manifest — so without this the `pg_lower` / `pg_upper` terms
would vanish and every residual would blow up by exactly the scarcity rent.
"""
function _manifest_duals(dict, key::Symbol)
    dict === nothing && return nothing
    bound_manifest_entry_present(dict, key) || return nothing
    ids, vals = extract_bound_manifest_duals(dict, key)
    return OrderedDict{Int, Float64}(Int(i) => Float64(v) for (i, v) in zip(ids, vals))
end

"""Marginal cost `dC/dP_g` in €/MW (the objective's p.u. scaling cancels — see README)."""
function _cost_gradient(ctx, gen::Int)::Float64
    row = findfirst(==(gen), Int.(ctx.DGEN.id))
    row === nothing && return 0.0
    c1 = Float64(ctx.DGEN.g_cost_1[row])
    ctx.cost_type == "quadratic" || return c1
    c2 = Float64(ctx.DGEN.g_cost_2[row])
    return c1 + 2 * c2 * ctx.pg_value[gen] * ctx.base_MVA
end

"""
Warn about any *other* constraint that couples `P_g`, so a new builder row cannot
silently skew the decomposition.

Only affine couplings are detectable this way: `normalized_coefficient` returns the
affine part of a QuadExpr (0 for a `P_g²` term) and throws on a nonlinear row. The
residual column is the net that catches the rest — this pass just names the culprit
when it can.
"""
function _warn_unlisted_pg_couplings(ctx)
    for dict in (ctx.opf_dict, ctx.dyn_model_dict)
        dict === nothing && continue
        for source in (:eq_const, :ineq_const)
            haskey(dict, source) || continue
            for (key, container) in dict[source]
                key in _PG_FAMILY_KEYS && continue
                for (gen, P_g_var) in ctx.P_g
                    coupled = false
                    for con in _leaf_constraints(container)
                        coefficient = try
                            JuMP.normalized_coefficient(con, P_g_var)
                        catch
                            0.0
                        end
                        if !iszero(coefficient)
                            coupled = true
                            break
                        end
                    end
                    if coupled
                        @warn "Constraint family couples P_g but is absent from the " *
                              "KKT family table; the decomposition is incomplete." key gen
                        break
                    end
                end
            end
        end
    end
    return nothing
end

"""Flatten a constraint container (scalar or time-nested) to its ConstraintRefs."""
function _leaf_constraints(container)
    container isa JuMP.ConstraintRef && return (container,)
    container isa AbstractDict || return ()
    refs = JuMP.ConstraintRef[]
    for value in values(container)
        if value isa JuMP.ConstraintRef
            push!(refs, value)
        elseif value isa AbstractDict
            append!(refs, _leaf_constraints(value))
        end
    end
    return refs
end

# --- entry point -------------------------------------------------------------------

"""
    Save_KKT_Pg_Terms!(path_names, cfg, sys, opf_dict, dyn_model_dict) -> String

Write `Dispatch/CSV_duals/kkt_Pg_terms.csv` and return its path.

Long format, one row per (generator, family, time):

    gen_id;bus_id;family;t_index;dual;coeff;contribution

plus, per generator, a `cost_gradient` row and a `RESIDUAL` row holding the sum of
everything above it. Must be called while the model is alive — after
`release_solver_backend!` every ConstraintRef here is dead.
"""
function Save_KKT_Pg_Terms!(path_names, cfg, sys, opf_dict, dyn_model_dict)::String
    P_g = opf_dict[:vars][:P_g]
    V = get(opf_dict[:vars], :V, nothing)

    gen_bus = OrderedDict{Int, Int}(
        Int(sys.DGEN.id[r]) => Int(sys.DGEN.bus[r]) for r in 1:nrow(sys.DGEN))

    ctx = (
        opf_dict = opf_dict, dyn_model_dict = dyn_model_dict,
        P_g = P_g, gen_bus = gen_bus,
        pg_value = OrderedDict{Int, Float64}(g => JuMP.value(v) for (g, v) in P_g),
        v_value = V === nothing ? OrderedDict{Int, Float64}() :
                  OrderedDict{Int, Float64}(b => JuMP.value(v) for (b, v) in V),
        DGEN = sys.DGEN, DGFM = sys.DGFM,
        cost_type = cfg.dispatch.cost_type, base_MVA = cfg.base_MVA,
    )

    # VARIABLE bound encoding keeps the P_g box off the constraint list.
    manifest_duals = Dict(
        :ineq_const_pg_lower => _manifest_duals(opf_dict, :ineq_const_pg_lower),
        :ineq_const_pg_upper => _manifest_duals(opf_dict, :ineq_const_pg_upper),
    )

    # `store` travels with every row so consumers never have to re-derive which families are
    # dynamic coupling and which are steady-state scarcity. `_PG_FAMILIES` already carries
    # the answer; a reader that reconstructs it from a list of names goes stale the moment a
    # family is added, which is how the notebook ended up carrying a `swing_first_step` entry
    # matching nothing. The two synthetic rows get their own labels so the partition over
    # this column is total.
    rows = DataFrame(gen_id = Int[], bus_id = Int[], family = String[],
                     store = String[],
                     t_index = Union{Missing, Int}[],
                     dual = Union{Missing, Float64}[],
                     coeff = Union{Missing, Float64}[],
                     contribution = Float64[])

    for gen in keys(P_g)
        bus = gen_bus[gen]
        total = 0.0

        for family in _PG_FAMILIES
            entries = _family_entries(ctx, family, gen)

            # No constraint row at all: either the family is off for this run, or the
            # bound was encoded on the variable and its dual lives in the manifest.
            if isempty(entries)
                fallback = get(manifest_duals, family.key, nothing)
                (fallback === nothing || !haskey(fallback, gen)) && continue
                coefficient = family.coeff isa Function ?
                    family.coeff(ctx, gen) : Float64(family.coeff)
                iszero(coefficient) && continue
                contribution = -fallback[gen] * coefficient
                total += contribution
                push!(rows, (gen, bus, family.family, String(family.store), missing,
                             fallback[gen], coefficient, contribution))
                continue
            end

            for (t_index, con) in entries
                coefficient = _pg_coefficient(ctx, family, gen, con)
                iszero(coefficient) && continue      # this step does not contain P_g
                dual = JuMP.dual(con)
                contribution = -dual * coefficient
                total += contribution
                push!(rows, (gen, bus, family.family, String(family.store), t_index,
                             dual, coefficient, contribution))
            end
        end

        gradient = _cost_gradient(ctx, gen)
        total += gradient
        push!(rows, (gen, bus, "cost_gradient", "objective",
                     missing, missing, missing, gradient))
        push!(rows, (gen, bus, "RESIDUAL", "check", missing, missing, missing, total))
    end

    _warn_unlisted_pg_couplings(ctx)

    file_path = joinpath(path_names[:pf_dispatch_CSV_duals], "kkt_Pg_terms.csv")
    mkpath(dirname(file_path))
    CSV.write(file_path, rows; delim = ';')
    return file_path
end
