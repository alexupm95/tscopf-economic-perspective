# ==============================================================================
#                    ED DUAL PROBLEM — Builder and Utilities
# ==============================================================================
#
# SIGN CONVENTION (read before modifying):
#   λ  is the KKT multiplier of the system active-power balance
#      ∑_g P_g = ∑_k P_{d,k}  (single equality, no network).
#   Stationarity w.r.t. each P_g (linear cost c_g):
#      c_g + λ - α̲_g + ᾱ_g = 0
#   System marginal price (SMP):  π = -λ  (€/MWh in normal merit-order dispatch).
#
#   Bound multipliers α̲, ᾱ are ≥ 0 (standard form for  g(x) ≤ 0 ).
#
#   NOTE: this explicit dual is a maximization LP.  Its λ is NOT identical to
#   JuMP.dual(eq_const_p_balance) on the primal (which follows L = f - λh - μg).
#   Economic price:  π_SMP = -λ_dual  ≈  JuMP.dual(balance)  on the primal.
#
# PAIRING WITH PRIMAL:
#   This explicit dual LP corresponds to the ED primal built by `Make_ED_Model!`
#   with `cost_type = "linear"` and explicit generator bound inequalities.
#   Quadratic ED primal → use `save_duals` on the NLP primal instead.
# ==============================================================================


# ==============================================================================
#                       SECTION 1 — DUAL VARIABLE CREATORS
# ==============================================================================

"""Scalar multiplier for the system-wide active-power balance."""
function var_ed_dual_lambda!(model::Model)
    return @variable(model, base_name = "λ_balance")
end

function var_ed_dual_alpha!(
    model::Model,
    gen_id::Vector,
)
    α_lower = OrderedDict{Int, JuMP.VariableRef}()
    α_upper = OrderedDict{Int, JuMP.VariableRef}()
    for gen in gen_id
        α_lower[gen] = @variable(model, base_name = "α_lower[$gen]")
        @constraint(model, α_lower[gen] ≥ 0.0)
        α_upper[gen] = @variable(model, base_name = "α_upper[$gen]")
        @constraint(model, α_upper[gen] ≥ 0.0)
    end
    return α_lower, α_upper
end


# ==============================================================================
#                   SECTION 2 — DUAL EQUALITY CONSTRAINTS
# ==============================================================================

function eq_const_ed_dual_stationarity_Pg!(
    model::Model,
    λ::JuMP.VariableRef,
    α_lower::Union{OrderedDict{Int, JuMP.VariableRef}, Nothing},
    α_upper::Union{OrderedDict{Int, JuMP.VariableRef}, Nothing},
    DGEN::DataFrame,
    active_gen::Vector,
    bound_P_g::Bool,
)
    eq_const_stat_Pg = OrderedDict{Int, JuMP.ConstraintRef}()

    for (idx, gen) in enumerate(DGEN.id[active_gen])
        c_g = DGEN.g_cost_1[active_gen[idx]]
        if bound_P_g
            eq_const_stat_Pg[gen] = @constraint(model,
                c_g + λ - α_lower[gen] + α_upper[gen] == 0.0)
        else
            eq_const_stat_Pg[gen] = @constraint(model, c_g + λ == 0.0)
        end
    end
    return eq_const_stat_Pg
end


# ==============================================================================
#                    SECTION 3 — TOP-LEVEL DUAL MODEL BUILDER
# ==============================================================================

function _ed_total_demand_pu(
    bus_gen_circ_dict_ON::OrderedDict,
    base_MVA::Float64,
)::Float64
    return sum(bus_gen_circ_dict_ON[bus][:pd_tot] for bus in keys(bus_gen_circ_dict_ON)) / base_MVA
end

function Make_ED_Dual_Model!(
    model::Model,
    path_names::OrderedDict{Symbol, String},
    DGEN::DataFrame,
    bus_gen_circ_dict_ON::OrderedDict,
    base_MVA::Float64,
    nGEN::Int64,
    opf_input_param::OrderedDict{Symbol, Any},
)
    active_gen = findall(x -> x == 1, DGEN.g_status)
    bound_P_g = opf_input_param[:var_bounds][:P_g]
    p_d_total = _ed_total_demand_pu(bus_gen_circ_dict_ON, base_MVA)

    dual_dict = OrderedDict{Symbol, Any}()
    dual_dict[:vars] = OrderedDict{Symbol, Any}()
    dual_dict[:eq_const] = OrderedDict{Symbol, OrderedDict{Int, JuMP.ConstraintRef}}()
    dual_dict[:ineq_const] = OrderedDict{Symbol, OrderedDict{Int, JuMP.ConstraintRef}}()
    dual_dict[:p_d_total_pu] = p_d_total

    λ = var_ed_dual_lambda!(model)
    dual_dict[:vars][:λ] = λ

    α_lower = α_upper = nothing
    if bound_P_g
        α_lower, α_upper = var_ed_dual_alpha!(model, DGEN.id[active_gen])
        dual_dict[:vars][:α_lower] = α_lower
        dual_dict[:vars][:α_upper] = α_upper
    end

    d_obj = JuMP.AffExpr(0.0)
    JuMP.add_to_expression!(d_obj, -p_d_total, λ)

    if bound_P_g
        for (idx, gen) in enumerate(DGEN.id[active_gen])
            pg_min_pu = DGEN.pg_min[active_gen[idx]] / base_MVA
            pg_max_pu = DGEN.pg_max[active_gen[idx]] / base_MVA
            JuMP.add_to_expression!(d_obj,  pg_min_pu, α_lower[gen])
            JuMP.add_to_expression!(d_obj, -pg_max_pu, α_upper[gen])
        end
    end

    @objective(model, Max, d_obj)

    eq_const_stat_Pg = eq_const_ed_dual_stationarity_Pg!(
        model, λ, α_lower, α_upper, DGEN, active_gen, bound_P_g)
    dual_dict[:eq_const][:eq_const_stat_Pg] = eq_const_stat_Pg

    println("--------------------------------------------------------------------------------------------------------------------------------------")
    println("  ED DUAL MODEL SUMMARY")
    println("  Generators (ON): $(length(active_gen))  |  System demand (p.u.): $(round(p_d_total, digits=6))")
    println("  Dual variables : $(num_variables(model))")
    println("  Equality const : $(length(eq_const_stat_Pg))")
    println("--------------------------------------------------------------------------------------------------------------------------------------")

    Export_ED_Dual_Model(model, path_names, d_obj, dual_dict)

    return model, dual_dict
end


# ==============================================================================
#               SECTION 4 — STRONG DUALITY VERIFICATION
# ==============================================================================

function Verify_Strong_Duality_ED!(
    primal_model::Model,
    primal_dict::OrderedDict{Symbol, Any},
    dual_model::Model,
    dual_dict::OrderedDict{Symbol, Any},
    DGEN::DataFrame;
    primal_f_star::Union{Nothing, Float64} = nothing,
    cross_check_primal_duals::Bool = true,
)
    primal_status = termination_status(primal_model)
    dual_status   = termination_status(dual_model)

    primal_status != MOI.OPTIMAL && @warn "Primal did not terminate optimally: $primal_status"
    dual_status != MOI.OPTIMAL && @warn "Dual did not terminate optimally: $dual_status"

    f_star = primal_f_star === nothing ? objective_value(primal_model) : primal_f_star
    d_star = objective_value(dual_model)

    gap_abs = abs(f_star - d_star)
    gap_rel = gap_abs / max(abs(f_star), 1.0)

    λ_val = JuMP.value(dual_dict[:vars][:λ])
    smp = -λ_val

    α_lo_vals = haskey(dual_dict[:vars], :α_lower) ?
        OrderedDict(g => JuMP.value(v) for (g, v) in dual_dict[:vars][:α_lower]) :
        OrderedDict{Int, Float64}()
    α_up_vals = haskey(dual_dict[:vars], :α_upper) ?
        OrderedDict(g => JuMP.value(v) for (g, v) in dual_dict[:vars][:α_upper]) :
        OrderedDict{Int, Float64}()

    primal_π = NaN
    if cross_check_primal_duals && haskey(primal_dict[:eq_const], :eq_const_p_balance)
        primal_π = JuMP.dual(primal_dict[:eq_const][:eq_const_p_balance][1])
    end

    println()
    println("══════════════════════════════════════════════════════════════════════")
    println("              ED STRONG DUALITY VERIFICATION REPORT")
    println("══════════════════════════════════════════════════════════════════════")
    @printf("  Primal objective  f* = %14.6f  (p.u.)\n", f_star)
    @printf("  Dual   objective  d* = %14.6f  (p.u.)\n", d_star)
    @printf("  Duality gap  |f*-d*| = %14.2e  (rel: %.2e)\n", gap_abs, gap_rel)
    if gap_rel < 1e-6
        println("  Strong duality verified  (gap < 1e-6)")
    else
        println("  WARNING: gap exceeds 1e-6 — check model formulation")
    end

    println()
    println("  ── System marginal price  π_SMP = -λ ──")
    @printf("    λ (dual LP)     = %12.6f\n", λ_val)
    @printf("    π_SMP           = %12.6f  €/MWh\n", smp)
    if cross_check_primal_duals && isfinite(primal_π)
        @printf("    π (JuMP primal) = %12.6f  (|Δπ| = %.2e)\n", primal_π, abs(smp - primal_π))
    end

    println()
    println("  ── Marginal / bounded generators ──")
    active_gen = findall(x -> x == 1, DGEN.g_status)
    for (idx, gen) in enumerate(DGEN.id[active_gen])
        α_lo = get(α_lo_vals, gen, 0.0)
        α_up = get(α_up_vals, gen, 0.0)
        c_g  = DGEN.g_cost_1[active_gen[idx]]
        if α_lo < 1e-4 && α_up < 1e-4
            @printf("    Gen %3d (c = %8.2f €/MWh) : MARGINAL\n", gen, c_g)
        elseif α_up > 1e-4
            @printf("    Gen %3d (c = %8.2f €/MWh) : AT MAX    (ᾱ = %.4f)\n", gen, c_g, α_up)
        elseif α_lo > 1e-4
            @printf("    Gen %3d (c = %8.2f €/MWh) : AT MIN    (α̲ = %.4f)\n", gen, c_g, α_lo)
        end
    end
    println("══════════════════════════════════════════════════════════════════════")
    println()

    return (;
        f_star, d_star, gap_abs, gap_rel, smp, λ = λ_val,
        α_lower = α_lo_vals, α_upper = α_up_vals,
    )
end


# ==============================================================================
#               SECTION 5 — EXPORT / SAVE DUAL RESULTS
# ==============================================================================

function Export_ED_Dual_Model(
    model::Model,
    path_names::OrderedDict{Symbol, String},
    obj_function::JuMP.AbstractJuMPScalar,
    dual_dict::OrderedDict{Symbol, Any},
)
    open(joinpath(path_names[:pf_dispatch], "model_summary.txt"), "w") do io
        show(io, model)
    end

    open(joinpath(path_names[:pf_dispatch], "model_details.txt"), "w") do io
        println(io, "=========")
        println(io, "Objective (max dual)")
        println(io, "=========")
        println(io, obj_function)
        println(io, "\n")
        println(io, "System demand (p.u.): ", dual_dict[:p_d_total_pu])
        println(io, "\n")
        println(io, "=========")
        println(io, "Dual variables")
        println(io, "=========")
        println(io, "λ: ", dual_dict[:vars][:λ])
        for sym in (:α_lower, :α_upper)
            haskey(dual_dict[:vars], sym) || continue
            for (j, info) in dual_dict[:vars][sym]
                println(io, "$sym[$j]: ", info)
            end
        end
        println(io, "\n")
        println(io, "=========")
        println(io, "D1 stationarity P_g")
        println(io, "=========")
        for (j, cref) in dual_dict[:eq_const][:eq_const_stat_Pg]
            println(io, "$j: ", cref)
        end
        println(io, "\n")
    end
    println("ED Dual model saved as TXT in: ", path_names[:pf_dispatch])
end

function Save_Duals_ED_Dual_Model(
    path_names::OrderedDict{Symbol, String},
    dual_model::Model,
    dual_dict::OrderedDict{Symbol, Any},
    verification::NamedTuple,
)
    λ_val = JuMP.value(dual_dict[:vars][:λ])
    smp = -λ_val

    α_lo_vals = haskey(dual_dict[:vars], :α_lower) ?
        OrderedDict(g => JuMP.value(v) for (g, v) in dual_dict[:vars][:α_lower]) :
        OrderedDict{Int, Float64}()
    α_up_vals = haskey(dual_dict[:vars], :α_upper) ?
        OrderedDict(g => JuMP.value(v) for (g, v) in dual_dict[:vars][:α_upper]) :
        OrderedDict{Int, Float64}()

    open(joinpath(path_names[:pf_dispatch], "dual_ED_dual_model.txt"), "w") do io
        println(io, "ED EXPLICIT DUAL — optimal multipliers")
        println(io, "  d* = $(objective_value(dual_model))  |  f* = $(verification.f_star)  |  gap = $(verification.gap_rel)")
        @printf(io, "  λ = %.8f   π_SMP = %.8f €/MWh\n", λ_val, smp)
    end

    open(joinpath(path_names[:pf_dispatch], "duality_verification.txt"), "w") do io
        @printf(io, "f* (primal) = %.8f\n", verification.f_star)
        @printf(io, "d* (dual)   = %.8f\n", verification.d_star)
        @printf(io, "gap_abs     = %.3e\n", verification.gap_abs)
        @printf(io, "gap_rel     = %.3e\n", verification.gap_rel)
    end

    df_smp = DataFrame(λ = [λ_val], SMP = [smp])
    CSV.write(joinpath(path_names[:pf_dispatch_CSV_duals], "dual_ed_SMP.csv"), df_smp; delim = ';')

    df_alpha = DataFrame()
    if !isempty(α_lo_vals)
        df_alpha = DataFrame(
            Gen_ID = collect(keys(α_lo_vals)),
            α_lower = collect(values(α_lo_vals)),
            α_upper = [α_up_vals[g] for g in keys(α_lo_vals)],
        )
        CSV.write(joinpath(path_names[:pf_dispatch_CSV_duals], "dual_ed_alpha.csv"), df_alpha; delim = ';')
    end

    df_opt = DataFrame(
        metric = ["primal_objective_pu", "dual_objective_pu", "gap_abs", "gap_rel", "SMP", "lambda"],
        value = [verification.f_star, verification.d_star, verification.gap_abs, verification.gap_rel, smp, λ_val],
    )
    CSV.write(joinpath(path_names[:pf_dispatch_CSV], "optimization_report.csv"), df_opt; delim = ';')
    open(joinpath(path_names[:pf_dispatch], "optimization_report.txt"), "w") do io
        println(io, "ED explicit dual optimisation report")
        @printf(io, "Primal f* (p.u.) : %.8f\n", verification.f_star)
        @printf(io, "Dual   d* (p.u.) : %.8f\n", verification.d_star)
        @printf(io, "Gap (relative)   : %.3e\n", verification.gap_rel)
        @printf(io, "SMP (π = -λ)     : %.8f €/MWh\n", smp)
    end

    sheets = Pair{String, DataFrame}[]
    push!(sheets, "SMP" => df_smp)
    !isempty(α_lo_vals) && push!(sheets, "Alpha_Gen" => df_alpha)
    push!(sheets, "Duality" => df_opt)

    excel_file = joinpath(path_names[:pf_dispatch], "ED_Dual_Results.xlsx")
    XLSX.writetable(excel_file, sheets...; overwrite = true)
    println("ED explicit dual results saved in: ", path_names[:pf_dispatch])
end


# ==============================================================================
#               SECTION 6 — ORCHESTRATION (build → solve → verify → save)
# ==============================================================================

function lp_solver_for_explicit_dual(primary_solver::String)::String
    if primary_solver in ("Gurobi", "HiGHS")
        return primary_solver
    end
    return gurobi_available() ? "Gurobi" : "HiGHS"
end

function Run_Explicit_ED_Dual!(
    path_names::OrderedDict{Symbol, String},
    dispatch::DispatchConfig,
    primal_model::Model,
    primal_dict::OrderedDict{Symbol, Any},
    DGEN::DataFrame,
    bus_gen_circ_dict_ON::OrderedDict,
    base_MVA::Float64,
    nGEN::Int64;
    primal_f_star::Union{Nothing, Float64} = nothing,
    solver_name::String = "Gurobi",
    silent_solver::Bool = false,
    highs::HiGHSSolverConfig = HiGHSSolverConfig(),
    gurobi::GurobiSolverConfig = GurobiSolverConfig(),
    cross_check_primal_duals::Bool = true,
)
    dual_paths = dispatch_dual_path_names(path_names)
    opf_input_param = build_opf_input_param(dispatch)

    lp_solver = lp_solver_for_explicit_dual(solver_name)
    dual_model = if lp_solver == "Gurobi"
        Setup_Optim_Model(lp_solver; gurobi=gurobi, silent=silent_solver)
    else
        Setup_Optim_Model(lp_solver; highs=highs, silent=silent_solver)
    end
    dual_model, dual_dict = Make_ED_Dual_Model!(
        dual_model, dual_paths, DGEN, bus_gen_circ_dict_ON,
        base_MVA, nGEN, opf_input_param)

    dual_log = joinpath(dual_paths[:pf_dispatch], "solver_log.txt")
    set_solver_log_path!(dual_model, lp_solver, dual_log)
    JuMP.optimize!(dual_model)

    verification = Verify_Strong_Duality_ED!(
        primal_model, primal_dict, dual_model, dual_dict, DGEN;
        primal_f_star = primal_f_star,
        cross_check_primal_duals = cross_check_primal_duals,
    )

    Save_Duals_ED_Dual_Model(dual_paths, dual_model, dual_dict, verification)

    release_solver_backend!(dual_model)

    return dual_model, dual_dict, verification
end
