# ==============================================================================
#                    DC-OPF DUAL PROBLEM — Builder and Utilities
# ==============================================================================
#
# SIGN CONVENTION (read before modifying):
#   λ_k  is the KKT multiplier of the nodal power balance (P1), defined so that
#   the stationarity condition w.r.t. P_g gives:  c_g + λ_{k(g)} - α̲_g + ᾱ_g = 0
#   The Locational Marginal Price is therefore  π_k = -λ_k  (positive in normal
#   operation because c_g > 0 and the unconstrained LMP equals the marginal cost).
#
#   All multipliers for inequality constraints (α̲, ᾱ, β̲, β̄, μ⁺, μ⁻, ν⁺, ν⁻)
#   are non-negative (KKT standard form for  g(x) ≤ 0  constraints).
#
#   All primal quantities are in per-unit (divided by base_MVA), matching the
#   convention used in Make_DCOPF_Model_w_Bbus!.  Angle bounds are in radians.
#
# PAIRING WITH PRIMAL:
#   This explicit dual LP corresponds to the Bbus (matrix) DC-OPF primal built by
#   Make_DCOPF_Model_w_Bbus!.  Constraint toggles in `opf_input_param` mirror the
#   primal builder (bound_P_g, bound_θ, ineq_sbranch_upper, ineq_ang_diff_branch).
# ==============================================================================


# ==============================================================================
#                       SECTION 1 — DUAL VARIABLE CREATORS
# ==============================================================================

function var_dcopf_dual_lambda!(
    model::Model,
    bus_id::Vector,
)
    λ = OrderedDict{Int, JuMP.VariableRef}()
    for bus in bus_id
        λ[bus] = @variable(model, base_name = "λ[$bus]")
    end
    return λ
end

function var_dcopf_dual_mu_ref!(model::Model)
    return @variable(model, base_name = "μ_ref")
end

function var_dcopf_dual_alpha!(
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

function var_dcopf_dual_beta!(
    model::Model,
    bus_id::Vector,
)
    β_lower = OrderedDict{Int, JuMP.VariableRef}()
    β_upper = OrderedDict{Int, JuMP.VariableRef}()
    for bus in bus_id
        β_lower[bus] = @variable(model, base_name = "β_lower[$bus]")
        @constraint(model, β_lower[bus] ≥ 0.0)
        β_upper[bus] = @variable(model, base_name = "β_upper[$bus]")
        @constraint(model, β_upper[bus] ≥ 0.0)
    end
    return β_lower, β_upper
end

function var_dcopf_dual_mu_branch!(
    model::Model,
    branch_id::Vector,
    from_bus::Vector,
    to_bus::Vector,
)
    μ_plus  = OrderedDict{Int, JuMP.VariableRef}()
    μ_minus = OrderedDict{Int, JuMP.VariableRef}()
    for (idx, branch) in enumerate(branch_id)
        i = from_bus[idx]
        k = to_bus[idx]
        μ_plus[branch]  = @variable(model, base_name = "μ_plus[$branch, ($i, $k)]")
        @constraint(model, μ_plus[branch]  ≥ 0.0)
        μ_minus[branch] = @variable(model, base_name = "μ_minus[$branch, ($i, $k)]")
        @constraint(model, μ_minus[branch] ≥ 0.0)
    end
    return μ_plus, μ_minus
end

function var_dcopf_dual_nu_anglediff!(
    model::Model,
    branch_id::Vector,
    from_bus::Vector,
    to_bus::Vector,
)
    ν_plus  = OrderedDict{Int, JuMP.VariableRef}()
    ν_minus = OrderedDict{Int, JuMP.VariableRef}()
    for (idx, branch) in enumerate(branch_id)
        i = from_bus[idx]
        k = to_bus[idx]
        ν_plus[branch]  = @variable(model, base_name = "ν_plus[$branch, ($i, $k)]")
        @constraint(model, ν_plus[branch]  ≥ 0.0)
        ν_minus[branch] = @variable(model, base_name = "ν_minus[$branch, ($i, $k)]")
        @constraint(model, ν_minus[branch] ≥ 0.0)
    end
    return ν_plus, ν_minus
end


# ==============================================================================
#                   SECTION 2 — DUAL EQUALITY CONSTRAINTS
# ==============================================================================

function eq_const_dcopf_dual_stationarity_Pg!(
    model::Model,
    λ::OrderedDict{Int, JuMP.VariableRef},
    α_lower::Union{OrderedDict{Int, JuMP.VariableRef}, Nothing},
    α_upper::Union{OrderedDict{Int, JuMP.VariableRef}, Nothing},
    DGEN::DataFrame,
    active_gen::Vector,
    bus_gen_circ_dict_ON::OrderedDict,
    bound_P_g::Bool,
)
    eq_const_stat_Pg = OrderedDict{Int, JuMP.ConstraintRef}()
    gen_bus_map = OrderedDict{Int, Int}()
    for bus in keys(bus_gen_circ_dict_ON)
        for gen in bus_gen_circ_dict_ON[bus][:gen_ids]
            gen_bus_map[gen] = bus
        end
    end

    for (idx, gen) in enumerate(DGEN.id[active_gen])
        bus_k = gen_bus_map[gen]
        c_g   = DGEN.g_cost_1[active_gen[idx]]
        if bound_P_g
            eq_const_stat_Pg[gen] = @constraint(model,
                c_g + λ[bus_k] - α_lower[gen] + α_upper[gen] == 0.0)
        else
            eq_const_stat_Pg[gen] = @constraint(model, c_g + λ[bus_k] == 0.0)
        end
    end
    return eq_const_stat_Pg
end

function eq_const_dcopf_dual_stationarity_theta!(
    model::Model,
    λ::OrderedDict{Int, JuMP.VariableRef},
    μ_ref::JuMP.VariableRef,
    β_lower::Union{OrderedDict{Int, JuMP.VariableRef}, Nothing},
    β_upper::Union{OrderedDict{Int, JuMP.VariableRef}, Nothing},
    μ_plus::Union{OrderedDict{Int, JuMP.VariableRef}, Nothing},
    μ_minus::Union{OrderedDict{Int, JuMP.VariableRef}, Nothing},
    ν_plus::Union{OrderedDict{Int, JuMP.VariableRef}, Nothing},
    ν_minus::Union{OrderedDict{Int, JuMP.VariableRef}, Nothing},
    DBUS::DataFrame,
    DCIR::DataFrame,
    active_branches::Vector,
    slack_bus::Int,
    susc_model::SusceptanceModel,
    bound_θ::Bool,
    include_branch_capacity::Bool,
    include_ang_diff::Bool,
)
    eq_const_stat_θ = OrderedDict{Int, JuMP.ConstraintRef}()
    bus_expr = OrderedDict{Int, JuMP.AffExpr}()
    for bus in DBUS.bus
        bus_expr[bus] = JuMP.AffExpr(0.0)
    end

    for branch in active_branches
        i = DCIR.from_bus[branch]
        k = DCIR.to_bus[branch]
        b = dc_branch_susceptance(DCIR.l_res[branch], DCIR.l_reac[branch], susc_model)
        JuMP.add_to_expression!(bus_expr[i],  b, λ[k])
        JuMP.add_to_expression!(bus_expr[i], -b, λ[i])
        JuMP.add_to_expression!(bus_expr[k],  b, λ[i])
        JuMP.add_to_expression!(bus_expr[k], -b, λ[k])
    end

    if include_branch_capacity && μ_plus !== nothing
        for branch in keys(μ_plus)
            i = DCIR.from_bus[branch]
            k = DCIR.to_bus[branch]
            b = dc_branch_susceptance(DCIR.l_res[branch], DCIR.l_reac[branch], susc_model)
            JuMP.add_to_expression!(bus_expr[i],  b,  μ_plus[branch])
            JuMP.add_to_expression!(bus_expr[i], -b,  μ_minus[branch])
            JuMP.add_to_expression!(bus_expr[k], -b,  μ_plus[branch])
            JuMP.add_to_expression!(bus_expr[k],  b,  μ_minus[branch])
        end
    end

    if include_ang_diff && ν_plus !== nothing
        for branch in active_branches
            i = DCIR.from_bus[branch]
            k = DCIR.to_bus[branch]
            JuMP.add_to_expression!(bus_expr[i],  1.0, ν_plus[branch])
            JuMP.add_to_expression!(bus_expr[i], -1.0, ν_minus[branch])
            JuMP.add_to_expression!(bus_expr[k], -1.0, ν_plus[branch])
            JuMP.add_to_expression!(bus_expr[k],  1.0, ν_minus[branch])
        end
    end

    for bus in DBUS.bus
        if bound_θ && β_lower !== nothing
            JuMP.add_to_expression!(bus_expr[bus], -1.0, β_lower[bus])
            JuMP.add_to_expression!(bus_expr[bus],  1.0, β_upper[bus])
        end
        if bus == slack_bus
            JuMP.add_to_expression!(bus_expr[bus], 1.0, μ_ref)
        end
        eq_const_stat_θ[bus] = @constraint(model, bus_expr[bus] == 0.0)
    end
    return eq_const_stat_θ
end


# ==============================================================================
#                    SECTION 3 — TOP-LEVEL DUAL MODEL BUILDER
# ==============================================================================

function Make_DC_OPF_Dual_Model!(
    model::Model,
    path_names::OrderedDict{Symbol, String},
    DBUS::DataFrame,
    DGEN::DataFrame,
    DCIR::DataFrame,
    bus_gen_circ_dict_ON::OrderedDict,
    base_MVA::Float64,
    nBUS::Int64,
    nGEN::Int64,
    nCIR::Int64,
    opf_input_param::OrderedDict{Symbol, Any},
)
    active_gen      = findall(x -> x == 1, DGEN.g_status)
    active_branches = findall(x -> x == 1, DCIR.l_status)
    pair_info, pair_circ_map = get_angle_limits(DCIR, active_branches)
    susc_model = opf_input_param[:susceptance_model]
    # Same limit values as the primal (θ box, DC angle-diff clamp). Branch capacity
    # duals exist only for rated branches (l_cap_1 != 0) when include_branch_capacity.
    limits = opf_input_param[:limits]

    bound_P_g = opf_input_param[:var_bounds][:P_g]
    bound_θ   = opf_input_param[:var_bounds][:θ]
    include_branch_capacity = opf_input_param[:ineq_cons][:sbranch_upper]
    include_ang_diff        = opf_input_param[:ineq_cons][:ang_diff_branch]

    SW = findall(x -> x == 3, DBUS.type)
    if isempty(SW)
        throw(ArgumentError("You must define one bus as the SLACK BUS (type 3)."))
    elseif length(SW) > 1
        throw(ArgumentError("This code still does not support more than one SLACK BUS (type 3)."))
    end
    slack_bus = SW[1]

    _, P_max_pu = branch_flow_limit_vectors(DCIR, base_MVA)
    rated_active_branches = [b for b in active_branches if isfinite(P_max_pu[b])]

    pair_circ_map_dict = OrderedDict(pair_circ_map)
    θ_diff_max_pu = OrderedDict{Int, Float64}()
    θ_diff_min_pu = OrderedDict{Int, Float64}()
    clamp_deg = limits.ang_diff_clamp_dc_deg   # mirrors dc_const_angle_differences! on the primal
    for (pair, pair_data) in pair_info
        circ_id = pair_circ_map_dict[pair]
        θ_diff_max_pu[circ_id] = deg2rad(min(pair_data.max_ang,  clamp_deg))
        θ_diff_min_pu[circ_id] = deg2rad(max(pair_data.min_ang, -clamp_deg))
    end

    dual_dict = OrderedDict{Symbol, Any}()
    dual_dict[:vars]       = OrderedDict{Symbol, Any}()
    dual_dict[:eq_const]   = OrderedDict{Symbol, OrderedDict{Int, JuMP.ConstraintRef}}()
    dual_dict[:ineq_const] = OrderedDict{Symbol, OrderedDict{Int, JuMP.ConstraintRef}}()
    dual_dict[:susceptance_model] = susc_model

    λ = var_dcopf_dual_lambda!(model, DBUS.bus)
    dual_dict[:vars][:λ] = λ

    μ_ref = var_dcopf_dual_mu_ref!(model)
    dual_dict[:vars][:μ_ref] = μ_ref

    α_lower = α_upper = nothing
    if bound_P_g
        α_lower, α_upper = var_dcopf_dual_alpha!(model, DGEN.id[active_gen])
        dual_dict[:vars][:α_lower] = α_lower
        dual_dict[:vars][:α_upper] = α_upper
    end

    β_lower = β_upper = nothing
    if bound_θ
        β_lower, β_upper = var_dcopf_dual_beta!(model, DBUS.bus)
        dual_dict[:vars][:β_lower] = β_lower
        dual_dict[:vars][:β_upper] = β_upper
    end

    μ_plus = μ_minus = nothing
    if include_branch_capacity && !isempty(rated_active_branches)
        μ_plus, μ_minus = var_dcopf_dual_mu_branch!(
            model,
            DCIR.id[rated_active_branches],
            DCIR.from_bus[rated_active_branches],
            DCIR.to_bus[rated_active_branches],
        )
        dual_dict[:vars][:μ_plus]  = μ_plus
        dual_dict[:vars][:μ_minus] = μ_minus
    end

    ν_plus = ν_minus = nothing
    if include_ang_diff
        ν_plus, ν_minus = var_dcopf_dual_nu_anglediff!(
            model,
            DCIR.id[active_branches],
            DCIR.from_bus[active_branches],
            DCIR.to_bus[active_branches],
        )
        dual_dict[:vars][:ν_plus]  = ν_plus
        dual_dict[:vars][:ν_minus] = ν_minus
    end

    begin
        d_obj = JuMP.AffExpr(0.0)
        for bus in keys(bus_gen_circ_dict_ON)
            p_d_pu = bus_gen_circ_dict_ON[bus][:pd_tot] / base_MVA
            JuMP.add_to_expression!(d_obj, -p_d_pu, λ[bus])
        end

        if bound_P_g
            for (idx, gen) in enumerate(DGEN.id[active_gen])
                pg_min_pu = DGEN.pg_min[active_gen[idx]] / base_MVA
                pg_max_pu = DGEN.pg_max[active_gen[idx]] / base_MVA
                JuMP.add_to_expression!(d_obj,  pg_min_pu, α_lower[gen])
                JuMP.add_to_expression!(d_obj, -pg_max_pu, α_upper[gen])
            end
        end

        if bound_θ
            # Primal box θ_min ≤ θ ≤ θ_max contributes +θ_min·β_lower − θ_max·β_upper.
            for bus in DBUS.bus
                JuMP.add_to_expression!(d_obj,  limits.θ_min_rad, β_lower[bus])
                JuMP.add_to_expression!(d_obj, -limits.θ_max_rad, β_upper[bus])
            end
        end

        if include_branch_capacity && μ_plus !== nothing
            for branch in keys(μ_plus)
                JuMP.add_to_expression!(d_obj, -P_max_pu[branch], μ_plus[branch])
                JuMP.add_to_expression!(d_obj, -P_max_pu[branch], μ_minus[branch])
            end
        end

        if include_ang_diff
            for branch in active_branches
                if haskey(θ_diff_max_pu, branch)
                    θ_max = θ_diff_max_pu[branch]
                    θ_min = θ_diff_min_pu[branch]
                    JuMP.add_to_expression!(d_obj, -θ_max, ν_plus[branch])
                    JuMP.add_to_expression!(d_obj,  θ_min, ν_minus[branch])
                end
            end
        end

        @objective(model, Max, d_obj)
    end

    eq_const_stat_Pg = eq_const_dcopf_dual_stationarity_Pg!(
        model, λ, α_lower, α_upper, DGEN, active_gen, bus_gen_circ_dict_ON, bound_P_g)
    dual_dict[:eq_const][:eq_const_stat_Pg] = eq_const_stat_Pg

    eq_const_stat_θ = eq_const_dcopf_dual_stationarity_theta!(
        model, λ, μ_ref, β_lower, β_upper, μ_plus, μ_minus, ν_plus, ν_minus,
        DBUS, DCIR, active_branches, slack_bus, susc_model,
        bound_θ, include_branch_capacity, include_ang_diff)
    dual_dict[:eq_const][:eq_const_stat_θ] = eq_const_stat_θ

    println("--------------------------------------------------------------------------------------------------------------------------------------")
    println("  DC-OPF DUAL MODEL SUMMARY")
    println("  Buses: $nBUS  |  Generators (ON): $(length(active_gen))  |  Branches (ON): $(length(active_branches))")
    println("  Dual variables : $(num_variables(model))")
    println("  Equality const : $(length(eq_const_stat_Pg) + length(eq_const_stat_θ))")
    println("--------------------------------------------------------------------------------------------------------------------------------------")

    Export_DCOPF_Dual_Model(model, path_names, d_obj, dual_dict)

    return model, dual_dict
end


# ==============================================================================
#               SECTION 4 — STRONG DUALITY VERIFICATION
# ==============================================================================

function Verify_Strong_Duality_DCOPF!(
    primal_model::Model,
    primal_dict::OrderedDict{Symbol, Any},
    dual_model::Model,
    dual_dict::OrderedDict{Symbol, Any},
    DBUS::DataFrame,
    DGEN::DataFrame,
    DCIR::DataFrame;
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

    λ_vals       = OrderedDict(bus => JuMP.value(v) for (bus, v) in dual_dict[:vars][:λ])
    α_lo_vals    = haskey(dual_dict[:vars], :α_lower) ?
        OrderedDict(g => JuMP.value(v) for (g, v) in dual_dict[:vars][:α_lower]) :
        OrderedDict{Int, Float64}()
    α_up_vals    = haskey(dual_dict[:vars], :α_upper) ?
        OrderedDict(g => JuMP.value(v) for (g, v) in dual_dict[:vars][:α_upper]) :
        OrderedDict{Int, Float64}()
    μ_plus_vals  = haskey(dual_dict[:vars], :μ_plus) ?
        OrderedDict(br => JuMP.value(v) for (br, v) in dual_dict[:vars][:μ_plus]) :
        OrderedDict{Int, Float64}()
    μ_minus_vals = haskey(dual_dict[:vars], :μ_minus) ?
        OrderedDict(br => JuMP.value(v) for (br, v) in dual_dict[:vars][:μ_minus]) :
        OrderedDict{Int, Float64}()

    lmp = OrderedDict(bus => -λ_vals[bus] for bus in keys(λ_vals))

    primal_λ = OrderedDict{Int, Float64}()
    if cross_check_primal_duals && haskey(primal_dict[:eq_const], :eq_const_p_balance)
        for (bus, cref) in primal_dict[:eq_const][:eq_const_p_balance]
            primal_λ[bus] = JuMP.dual(cref)
        end
    end

    println()
    println("══════════════════════════════════════════════════════════════════════")
    println("           DC-OPF STRONG DUALITY VERIFICATION REPORT")
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
    println("  ── Locational Marginal Prices π_k = -λ_k ──")
    for (bus, π) in sort(collect(lmp), by = x -> x[2], rev = true)
        primal_λ_val = get(primal_λ, bus, NaN)
        @printf("    Bus %3d : π = %10.4f  [λ_dual = %10.4f  λ_primal = %10.4f]\n",
            bus, π, λ_vals[bus], primal_λ_val)
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

    println()
    println("  ── Congested Branches  (μ⁺ or μ⁻ > 1e-4) ──")
    any_congested = false
    active_branches = findall(x -> x == 1, DCIR.l_status)
    for branch in active_branches
        μp = get(μ_plus_vals,  branch, 0.0)
        μm = get(μ_minus_vals, branch, 0.0)
        if μp > 1e-4 || μm > 1e-4
            i = DCIR.from_bus[branch]
            k = DCIR.to_bus[branch]
            @printf("    Branch %3d (%d→%d) : μ⁺ = %.4f  μ⁻ = %.4f  |π_i - π_k| = %.4f\n",
                branch, i, k, μp, μm, abs(lmp[i] - lmp[k]))
            any_congested = true
        end
    end
    !any_congested && println("    (none)")
    println("══════════════════════════════════════════════════════════════════════")
    println()

    return (;
        f_star, d_star, gap_abs, gap_rel, lmp, λ = λ_vals,
        α_lower = α_lo_vals, α_upper = α_up_vals,
        μ_plus = μ_plus_vals, μ_minus = μ_minus_vals,
    )
end


# ==============================================================================
#               SECTION 5 — EXPORT / SAVE DUAL RESULTS
# ==============================================================================

function Export_DCOPF_Dual_Model(
    model::Model,
    path_names::OrderedDict{Symbol, String},
    obj_function::JuMP.AbstractJuMPScalar,
    dual_dict::OrderedDict{Symbol, Any},
)
    open(joinpath(path_names[:pf_dispatch], "model_summary.txt"), "w") do io
        show(io, model)
    end

    var_groups = Symbol[]
    for sym in (:λ, :μ_ref, :α_lower, :α_upper, :β_lower, :β_upper, :μ_plus, :μ_minus, :ν_plus, :ν_minus)
        haskey(dual_dict[:vars], sym) && push!(var_groups, sym)
    end

    open(joinpath(path_names[:pf_dispatch], "model_details.txt"), "w") do io
        println(io, "=========")
        println(io, "Objective (max dual)")
        println(io, "=========")
        println(io, obj_function)
        println(io, "\n")
        println(io, "=========")
        println(io, "Dual variables")
        println(io, "=========")
        for sym in var_groups
            entry = dual_dict[:vars][sym]
            if entry isa JuMP.VariableRef
                println(io, "$sym: ", entry)
            else
                for (j, info) in entry
                    println(io, "$sym[$j]: ", info)
                end
            end
        end
        println(io, "\n")
        for (label, key) in (("D1 stationarity P_g", :eq_const_stat_Pg),
                             ("D2 stationarity θ", :eq_const_stat_θ))
            haskey(dual_dict[:eq_const], key) || continue
            println(io, "=========")
            println(io, label)
            println(io, "=========")
            for (j, cref) in dual_dict[:eq_const][key]
                println(io, "$j: ", cref)
            end
            println(io, "\n")
        end
    end
    println("DC-OPF Dual model saved as TXT in: ", path_names[:pf_dispatch])
end

function Save_Duals_DC_OPF_Dual_Model(
    path_names::OrderedDict{Symbol, String},
    dual_model::Model,
    dual_dict::OrderedDict{Symbol, Any},
    verification::NamedTuple,
)
    λ_vals = OrderedDict(bus => JuMP.value(v) for (bus, v) in dual_dict[:vars][:λ])
    μ_ref_val = JuMP.value(dual_dict[:vars][:μ_ref])
    α_lo_vals = haskey(dual_dict[:vars], :α_lower) ?
        OrderedDict(g => JuMP.value(v) for (g, v) in dual_dict[:vars][:α_lower]) :
        OrderedDict{Int, Float64}()
    α_up_vals = haskey(dual_dict[:vars], :α_upper) ?
        OrderedDict(g => JuMP.value(v) for (g, v) in dual_dict[:vars][:α_upper]) :
        OrderedDict{Int, Float64}()
    β_lo_vals = haskey(dual_dict[:vars], :β_lower) ?
        OrderedDict(b => JuMP.value(v) for (b, v) in dual_dict[:vars][:β_lower]) :
        OrderedDict{Int, Float64}()
    β_up_vals = haskey(dual_dict[:vars], :β_upper) ?
        OrderedDict(b => JuMP.value(v) for (b, v) in dual_dict[:vars][:β_upper]) :
        OrderedDict{Int, Float64}()
    μ_plus_vals = haskey(dual_dict[:vars], :μ_plus) ?
        OrderedDict(br => JuMP.value(v) for (br, v) in dual_dict[:vars][:μ_plus]) :
        OrderedDict{Int, Float64}()
    μ_minus_vals = haskey(dual_dict[:vars], :μ_minus) ?
        OrderedDict(br => JuMP.value(v) for (br, v) in dual_dict[:vars][:μ_minus]) :
        OrderedDict{Int, Float64}()
    ν_plus_vals = haskey(dual_dict[:vars], :ν_plus) ?
        OrderedDict(br => JuMP.value(v) for (br, v) in dual_dict[:vars][:ν_plus]) :
        OrderedDict{Int, Float64}()
    ν_minus_vals = haskey(dual_dict[:vars], :ν_minus) ?
        OrderedDict(br => JuMP.value(v) for (br, v) in dual_dict[:vars][:ν_minus]) :
        OrderedDict{Int, Float64}()

    lmp_vals = OrderedDict(bus => -λ_vals[bus] for bus in keys(λ_vals))

    open(joinpath(path_names[:pf_dispatch], "dual_DC_OPF_dual_model.txt"), "w") do io
        println(io, "DC-OPF EXPLICIT DUAL — optimal multipliers")
        println(io, "  d* = $(objective_value(dual_model))  |  f* = $(verification.f_star)  |  gap = $(verification.gap_rel)")
        println(io, "  μ_ref = $μ_ref_val")
        println(io)
        for bus in sort(collect(keys(λ_vals)))
            @printf(io, "  Bus %3d : λ = %12.6f   π = %12.6f\n", bus, λ_vals[bus], lmp_vals[bus])
        end
    end

    open(joinpath(path_names[:pf_dispatch], "duality_verification.txt"), "w") do io
        @printf(io, "f* (primal) = %.8f\n", verification.f_star)
        @printf(io, "d* (dual)   = %.8f\n", verification.d_star)
        @printf(io, "gap_abs     = %.3e\n", verification.gap_abs)
        @printf(io, "gap_rel     = %.3e\n", verification.gap_rel)
    end

    df_lmp = DataFrame(Bus_ID = collect(keys(lmp_vals)), LMP = collect(values(lmp_vals)), λ = collect(values(λ_vals)))
    CSV.write(joinpath(path_names[:pf_dispatch_CSV_duals], "dual_dcopf_LMP.csv"), df_lmp; delim = ';')

    df_alpha = DataFrame()
    df_beta = DataFrame()
    df_mu = DataFrame()
    df_nu = DataFrame()

    if !isempty(α_lo_vals)
        df_alpha = DataFrame(
            Gen_ID = collect(keys(α_lo_vals)),
            α_lower = collect(values(α_lo_vals)),
            α_upper = [α_up_vals[g] for g in keys(α_lo_vals)],
        )
        CSV.write(joinpath(path_names[:pf_dispatch_CSV_duals], "dual_dcopf_alpha.csv"), df_alpha; delim = ';')
    end

    if !isempty(β_lo_vals)
        df_beta = DataFrame(
            Bus_ID = collect(keys(β_lo_vals)),
            β_lower = collect(values(β_lo_vals)),
            β_upper = [β_up_vals[b] for b in keys(β_lo_vals)],
        )
        CSV.write(joinpath(path_names[:pf_dispatch_CSV_duals], "dual_dcopf_beta.csv"), df_beta; delim = ';')
    end

    if !isempty(μ_plus_vals)
        df_mu = DataFrame(
            Branch_ID = collect(keys(μ_plus_vals)),
            μ_plus = collect(values(μ_plus_vals)),
            μ_minus = [μ_minus_vals[br] for br in keys(μ_plus_vals)],
        )
        CSV.write(joinpath(path_names[:pf_dispatch_CSV_duals], "dual_dcopf_mu_branch.csv"), df_mu; delim = ';')
    end

    if !isempty(ν_plus_vals)
        df_nu = DataFrame(
            Branch_ID = collect(keys(ν_plus_vals)),
            ν_plus = collect(values(ν_plus_vals)),
            ν_minus = [ν_minus_vals[br] for br in keys(ν_plus_vals)],
        )
        CSV.write(joinpath(path_names[:pf_dispatch_CSV_duals], "dual_dcopf_nu_anglediff.csv"), df_nu; delim = ';')
    end

    df_opt = DataFrame(
        metric = ["primal_objective_pu", "dual_objective_pu", "gap_abs", "gap_rel"],
        value = [verification.f_star, verification.d_star, verification.gap_abs, verification.gap_rel],
    )
    CSV.write(joinpath(path_names[:pf_dispatch_CSV], "optimization_report.csv"), df_opt; delim = ';')
    open(joinpath(path_names[:pf_dispatch], "optimization_report.txt"), "w") do io
        println(io, "DC-OPF explicit dual optimisation report")
        @printf(io, "Primal f* (p.u.) : %.8f\n", verification.f_star)
        @printf(io, "Dual   d* (p.u.) : %.8f\n", verification.d_star)
        @printf(io, "Gap (relative)   : %.3e\n", verification.gap_rel)
    end

    sheets = Pair{String, DataFrame}[]
    push!(sheets, "LMP" => df_lmp)
    !isempty(α_lo_vals) && push!(sheets, "Alpha_Gen" => df_alpha)
    !isempty(β_lo_vals) && push!(sheets, "Beta_Bus" => df_beta)
    !isempty(μ_plus_vals) && push!(sheets, "Mu_Branch" => df_mu)
    !isempty(ν_plus_vals) && push!(sheets, "Nu_AngDiff" => df_nu)
    push!(sheets, "Duality" => df_opt)

    excel_file = joinpath(path_names[:pf_dispatch], "DC_OPF_Dual_Results.xlsx")
    XLSX.writetable(excel_file, sheets...; overwrite = true)
    println("DC-OPF explicit dual results saved in: ", path_names[:pf_dispatch])
end


# ==============================================================================
#               SECTION 6 — ORCHESTRATION (build → solve → verify → save)
# ==============================================================================

"""Pick an LP solver for the explicit DC-OPF dual (Gurobi preferred, HiGHS fallback)."""
function lp_solver_for_dcopf_dual(primary_solver::String)::String
  if primary_solver in ("Gurobi", "HiGHS")
    return primary_solver
  end
  return gurobi_available() ? "Gurobi" : "HiGHS"
end

function Run_Explicit_DC_OPF_Dual!(
    path_names::OrderedDict{Symbol, String},
    dispatch::DispatchConfig,
    primal_model::Model,
    primal_dict::OrderedDict{Symbol, Any},
    DBUS::DataFrame,
    DGEN::DataFrame,
    DCIR::DataFrame,
    bus_gen_circ_dict_ON::OrderedDict,
    base_MVA::Float64,
    nBUS::Int64,
    nGEN::Int64,
    nCIR::Int64;
    primal_f_star::Union{Nothing, Float64} = nothing,
    solver_name::String = "Gurobi",
    silent_solver::Bool = false,
    highs::HiGHSSolverConfig = HiGHSSolverConfig(),
    gurobi::GurobiSolverConfig = GurobiSolverConfig(),
    cross_check_primal_duals::Bool = true,
)
    dual_paths = dispatch_dual_path_names(path_names)
    opf_input_param = build_opf_input_param(dispatch)

    lp_solver = lp_solver_for_dcopf_dual(solver_name)
    dual_model = if lp_solver == "Gurobi"
        Setup_Optim_Model(lp_solver; gurobi=gurobi, silent=silent_solver)
    else
        Setup_Optim_Model(lp_solver; highs=highs, silent=silent_solver)
    end
    dual_model, dual_dict = Make_DC_OPF_Dual_Model!(
        dual_model, dual_paths, DBUS, DGEN, DCIR, bus_gen_circ_dict_ON,
        base_MVA, nBUS, nGEN, nCIR, opf_input_param)

    dual_log = joinpath(dual_paths[:pf_dispatch], "solver_log.txt")
    set_solver_log_path!(dual_model, lp_solver, dual_log)
    JuMP.optimize!(dual_model)

    verification = Verify_Strong_Duality_DCOPF!(
        primal_model, primal_dict, dual_model, dual_dict, DBUS, DGEN, DCIR;
        primal_f_star = primal_f_star,
        cross_check_primal_duals = cross_check_primal_duals,
    )

    Save_Duals_DC_OPF_Dual_Model(dual_paths, dual_model, dual_dict, verification)

    release_solver_backend!(dual_model)

    return dual_model, dual_dict, verification
end
