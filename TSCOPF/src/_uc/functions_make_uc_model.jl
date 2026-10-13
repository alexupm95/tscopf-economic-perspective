# ==============================================================================
#                    UNIT COMMITMENT (UC) — single-period MILP + restricted LP
# ==============================================================================
#
# Mathematical form (single period, active generators only):
#
#   min  Σ_g c_g P_g
#   s.t. Σ_g P_g = D                          (system active-power balance)
#        P_g^min u_g ≤ P_g ≤ P_g^max u_g     (on/off generation limits)
#        u_g ∈ {0, 1}
#
# Bounds are written as explicit ≤ constraints (ED convention) so that
# JuMP.dual() on the restricted LP exposes α_lower / α_upper cleanly.
#
# Pricing / duals:
#   MILP optimal (P*, u*) → fix u = u* → solve LP with ED-style bounds
#   → JuMP.dual(balance) gives SMP (€/MW); JuMP.dual(Pg bounds) gives scarcity rents.
#   MILP KKT multipliers are NOT used (mixed-integer problem).
#
# Future hooks:
#   opf_dict[:meta][:n_periods] = 1  →  time-index u[g,t], P[g,t], ramp, startup later.
#   Reuses ed_const_power_balance! and obj_minimize_marginal_cost! from ED.
#
# Solver: Gurobi only (MILP + LP restricted pricing).
# ==============================================================================

# Variable / equality / inequality builders live in dedicated files for
# modularity (future ramp, min up/down, multi-period constraints):
#   _uc/uc_variables.jl         var_uc_dispatch!, var_uc_commitment!
#   _uc/uc_eq_constraints.jl    eq_const_uc_power_balance!
#   _uc/uc_ineq_constraints.jl  ineq_const_uc_pg_limits!, ineq_const_uc_pg_limits_fixed!

"""
Build the single-period UC MILP and export model structure before solve.

Returns `(model, obj_function, obj_function_MVA, opf_dict)` — same tuple pattern as ED.
`obj_function` is in €/MVA (JuMP objective); `obj_function_MVA` is the MVA-scaled form
used in model_details.txt exports.
"""
function Make_UC_Model!(
    model::Model,
    path_names::OrderedDict{Symbol, String},
    DBUS::DataFrame,
    DGEN::DataFrame,
    bus_gen_circ_dict_ON::OrderedDict,
    base_MVA::Float64,
    nBUS::Int64,
    nGEN::Int64,
    opf_input_param::OrderedDict{Symbol, Any},
)
    # Only in-service generators participate (g_status == 1).
    active_gen = findall(x -> x == 1, DGEN.g_status)
    gen_ids = DGEN.id[active_gen]

    # Standard dispatch container (vars / eq / ineq / meta).
    opf_dict = OrderedDict{Symbol, Any}()
    opf_dict[:vars] = OrderedDict{Symbol, Any}()
    opf_dict[:eq_const] = OrderedDict{Symbol, OrderedDict{Int, JuMP.ConstraintRef}}()
    opf_dict[:ineq_const] = OrderedDict{Symbol, OrderedDict{Int, JuMP.ConstraintRef}}()
    opf_dict[:meta] = OrderedDict{Symbol, Any}(:n_periods => 1)  # multi-period hook

    # Continuous dispatch P_g (pu); commitment u_g added below.
    P_g = var_uc_dispatch!(model, gen_ids)
    opf_dict[:vars][:P_g] = P_g

    u = var_uc_commitment!(model, gen_ids)
    opf_dict[:vars][:u_commit] = u

    # Limits in per-unit (same scaling as ED).
    pg_min_pu = DGEN.pg_min[active_gen] ./ base_MVA
    pg_max_pu = DGEN.pg_max[active_gen] ./ base_MVA
    ineq_lo, ineq_up = ineq_const_uc_pg_limits!(model, P_g, u, gen_ids, pg_min_pu, pg_max_pu)
    opf_dict[:ineq_const][:ineq_const_pg_lower] = ineq_lo
    opf_dict[:ineq_const][:ineq_const_pg_upper] = ineq_up

    # UC is MILP → linear cost only (quadratic would need MIQP / piecewise linearisation).
    if opf_input_param[:obj_function][:type] != "linear"
        throw(ArgumentError("UC requires cost_type=\"linear\" (MILP)."))
    end
    obj_function, obj_function_MVA = obj_minimize_marginal_cost!(
        model, P_g, gen_ids, DGEN.g_cost_1[active_gen], base_MVA)

    # Single-bus or multi-bus: UC balance wrapper (delegates to shared ED helper).
    eq_const_p_balance = eq_const_uc_power_balance!(model, P_g, bus_gen_circ_dict_ON, gen_ids, base_MVA)
    opf_dict[:eq_const][:eq_const_p_balance] = eq_const_p_balance

    # Ensure Dispatch/ exists before writing model_summary.txt at build time.
    mkpath(path_names[:pf_dispatch])

    # Archive MILP structure before optimize! (model_summary + model_details).
    println("--------------------------------------------------------------------------------------------------------------------------------------")
    Export_UC_Model(model, path_names, obj_function_MVA, opf_dict)
    println("--------------------------------------------------------------------------------------------------------------------------------------")

    return model, obj_function, obj_function_MVA, opf_dict
end

"""
Restricted pricing: fix u* from the MILP, rebuild as a pure LP, solve for shadow prices.

For u*_g = 1:  P_min ≤ P_g ≤ P_max  (standard ED bounds).
For u*_g = 0:  0 ≤ P_g ≤ 0           (generator forced off).

The returned `(lp_model, lp_dict)` is consumed by Save_Duals_UC_Model and
Export_UC_Restricted_LP_Model! after optimize!.
"""
function Solve_UC_Restricted_Pricing!(
    DGEN::DataFrame,
    bus_gen_circ_dict_ON::OrderedDict,
    base_MVA::Float64,
    u_star::OrderedDict{Int, Float64};
    silent_solver::Bool = true,
    gurobi::GurobiSolverConfig = GurobiSolverConfig(),
)
    active_gen = findall(x -> x == 1, DGEN.g_status)
    gen_ids = DGEN.id[active_gen]
    pg_min_pu = DGEN.pg_min[active_gen] ./ base_MVA
    pg_max_pu = DGEN.pg_max[active_gen] ./ base_MVA

    lp_model = Setup_Optim_Model("Gurobi"; gurobi=gurobi, silent=silent_solver)
    lp_dict = OrderedDict{Symbol, Any}()
    lp_dict[:vars] = OrderedDict{Symbol, Any}()
    lp_dict[:eq_const] = OrderedDict{Symbol, OrderedDict{Int, JuMP.ConstraintRef}}()
    lp_dict[:ineq_const] = OrderedDict{Symbol, OrderedDict{Int, JuMP.ConstraintRef}}()

    P_g = var_uc_dispatch!(lp_model, gen_ids)
    lp_dict[:vars][:P_g] = P_g

    # Substitute u* into the bilinear limits → fixed ED-style box constraints.
    ineq_lo, ineq_up = ineq_const_uc_pg_limits_fixed!(
        lp_model, P_g, gen_ids, pg_min_pu, pg_max_pu, u_star)
    lp_dict[:ineq_const][:ineq_const_pg_lower] = ineq_lo
    lp_dict[:ineq_const][:ineq_const_pg_upper] = ineq_up

    obj_minimize_marginal_cost!(lp_model, P_g, gen_ids, DGEN.g_cost_1[active_gen], base_MVA)
    lp_dict[:eq_const][:eq_const_p_balance] =
        eq_const_uc_power_balance!(lp_model, P_g, bus_gen_circ_dict_ON, gen_ids, base_MVA)

    JuMP.optimize!(lp_model)
    termination_status(lp_model) == MOI.OPTIMAL ||
        @warn "UC restricted LP did not terminate optimally: $(termination_status(lp_model))"

    return lp_model, lp_dict
end

"""Extract solved commitment values u*_g from the MILP (post optimize!)."""
function uc_commitment_values(u::OrderedDict{Int, JuMP.VariableRef})
    return OrderedDict(g => JuMP.value(v) for (g, v) in u)
end

# Dump MILP structure at build time (before solve). Layout mirrors Export_OPF_Model (ED subset).
function Export_UC_Model(
    model::Model,
    path_names::OrderedDict{Symbol, String},
    obj_function::JuMP.AbstractJuMPScalar,
    opf_dict::OrderedDict{Symbol, Any},
)
    open(joinpath(path_names[:pf_dispatch], "model_summary.txt"), "w") do io
        println(io, "=== UC MILP Model Summary ===")
        show(io, model)
    end

    open(joinpath(path_names[:pf_dispatch], "model_details.txt"), "w") do io
        println(io, "type_model: UC (single-period MILP; binary commitment u_g)")
        println(io, "n_periods: ", get(opf_dict[:meta], :n_periods, 1))
        println(io, "pricing_note: shadow prices come from restricted LP (u* fixed), not MILP KKT")
        println(io)

        println(io, "=========")
        println(io, "Objective ")
        println(io, "=========")
        println(io, obj_function)
        println(io, "\n")

        println(io, "=========")
        println(io, "Variables")
        println(io, "=========")
        for (g, v) in opf_dict[:vars][:P_g]
            println(io, "$g: ", v)
        end
        for (g, v) in opf_dict[:vars][:u_commit]
            println(io, "u[$g]: ", v)
        end
        println(io, "\n")

        if haskey(opf_dict[:eq_const], :eq_const_p_balance)
            println(io, "===================================================")
            println(io, "Equality Constraints Active Power Balance for Buses ")
            println(io, "===================================================")
            for (i, info) in opf_dict[:eq_const][:eq_const_p_balance]
                println(io, "$i: ", info)
            end
            println(io, "\n")
        end

        println(io, "=========================================================")
        println(io, "Inequality Constraints Inferior Limits Decision Variables ")
        println(io, "=========================================================")
        if haskey(opf_dict[:ineq_const], :ineq_const_pg_lower)
            println(io, "====================================")
            println(io, "Active Power Generated - Lower Bound (P_min u - P_g ≤ 0)")
            println(io, "====================================")
            for (i, info) in opf_dict[:ineq_const][:ineq_const_pg_lower]
                println(io, "$i: ", info)
            end
            println(io, "\n")
        end

        if haskey(opf_dict[:ineq_const], :ineq_const_pg_upper)
            println(io, "====================================")
            println(io, "Active Power Generated - Upper Bound (P_g - P_max u ≤ 0)")
            println(io, "====================================")
            for (i, info) in opf_dict[:ineq_const][:ineq_const_pg_upper]
                println(io, "$i: ", info)
            end
            println(io, "\n")
        end
    end
    println("UC MILP model saved as TXT in: ", path_names[:pf_dispatch])
end

# After restricted LP solve: archive fixed u* and LP structure for audit / debugging.
function Export_UC_Restricted_LP_Model!(
    path_names::OrderedDict{Symbol, String},
    lp_model::Model,
    lp_dict::OrderedDict{Symbol, Any},
    u_star::OrderedDict{Int, Float64},
)
    mkpath(path_names[:pf_dispatch])

    open(joinpath(path_names[:pf_dispatch], "restricted_lp_model_summary.txt"), "w") do io
        println(io, "=== UC Restricted-Pricing LP (u* fixed from MILP) ===")
        show(io, lp_model)
    end

    open(joinpath(path_names[:pf_dispatch], "restricted_lp_model_details.txt"), "w") do io
        println(io, "u_star (fixed commitment):")
        for (g, u) in sort(collect(u_star))
            println(io, "  gen $g: u = $u")
        end
        println(io)
        println(io, "========= Objective =========")
        println(io, objective_function(lp_model))
        println(io)
        if haskey(lp_dict[:eq_const], :eq_const_p_balance)
            println(io, "========= Equality: power balance =========")
            for (i, c) in lp_dict[:eq_const][:eq_const_p_balance]
                println(io, "$i: ", c)
            end
            println(io)
        end
        if haskey(lp_dict[:ineq_const], :ineq_const_pg_lower)
            println(io, "========= Inequality: P_g lower =========")
            for (i, c) in lp_dict[:ineq_const][:ineq_const_pg_lower]
                println(io, "$i: ", c)
            end
            println(io)
        end
        if haskey(lp_dict[:ineq_const], :ineq_const_pg_upper)
            println(io, "========= Inequality: P_g upper =========")
            for (i, c) in lp_dict[:ineq_const][:ineq_const_pg_upper]
                println(io, "$i: ", c)
            end
        end
    end
    println("UC restricted-pricing LP saved in: ", path_names[:pf_dispatch])
end
