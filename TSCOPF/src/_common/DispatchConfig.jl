#=
================================================================================
 DispatchConfig.jl  —  steady-state OPF builder configuration
================================================================================
 Replaces `Define_Inputs_4OPF()` in `config_file.jl`.  User-facing knobs for
 which constraints, bounds, and objective form the dispatch layer builds.

 Mandatory equalities (slack angle, P/Q balance, branch flows where required by
 the formulation) are always built inside `Make_*OPF_Model*!` — not toggled here.
================================================================================
=#

# Branch susceptance convention for the DC-OPF (ignored by AC/ED).
#   SIMPLE      : b = 1/x                       (textbook DC, ignores series resistance)
#   POWERMODELS : b = imag(inv(r + x·im))       = -x/(r²+x²)  (PowerModels/MATPOWER variant)
# The two coincide when r = 0 and diverge otherwise; see scripts/generate_pm_dc_reference.jl.
"""
    SusceptanceModel

DC-OPF branch susceptance convention (ignored by AC-OPF and ED):

- `SIMPLE`: `b = 1/x` (textbook DC, ignores series resistance).
- `POWERMODELS`: `b = imag(inv(r + x*im))` (matches PowerModels `DCPPowerModel`).

The two agree when `r = 0` and diverge otherwise. See the user guide, section 2.1.
"""
@enum SusceptanceModel SIMPLE POWERMODELS

"""
    DispatchLimitsConfig

Steady-state limit values consumed by the OPF builders (`DispatchConfig.limits`).

Only limits that are NOT case data live here; `V`, `P_g`, `Q_g`, branch thermal
ratings, and per-branch angle-difference limits always come from the input CSVs.

- `θ_min_rad` / `θ_max_rad` — bus voltage-angle box when `bound_θ=true`.
- `ang_diff_clamp_ac_deg` / `ang_diff_clamp_dc_deg` — safety clamp applied to
  `DCIR.ang_min`/`ang_max` by the AC / DC angle-difference constraints (values
  beyond the clamp are tightened to it, with a printed warning).

Branch thermal / flow limits use `DCIR.l_cap_1` only when
`DispatchConfig.ineq_sbranch_upper` (or optional `bound_P_ik` / …) is on.
Unrated branches (`l_cap_1 == 0`) are skipped (no finite sentinel). The DC dual
LP mirrors the same rated-branch set for strong duality.
"""
Base.@kwdef struct DispatchLimitsConfig
    θ_min_rad::Float64 = -π
    θ_max_rad::Float64 = π
    ang_diff_clamp_ac_deg::Float64 = 60.0
    ang_diff_clamp_dc_deg::Float64 = 30.0
end

function validate_dispatch_limits!(limits::DispatchLimitsConfig)
    limits.θ_min_rad < limits.θ_max_rad ||
        throw(ArgumentError("DispatchLimitsConfig: require θ_min_rad < θ_max_rad."))
    limits.ang_diff_clamp_ac_deg > 0 && limits.ang_diff_clamp_dc_deg > 0 ||
        throw(ArgumentError("DispatchLimitsConfig: angle-diff clamps must be positive."))
    return nothing
end

"""
    branch_flow_limit_vectors(DCIR, base_MVA) -> (lower, upper)

Per-branch flow / thermal limit vectors in p.u. Rated branches (`l_cap_1 != 0`)
get `±l_cap_1 / base_MVA`. Unrated branches get `±Inf` so attach helpers skip
them (same inactive policy as `TsBoundLimitPair`).
"""
function branch_flow_limit_vectors(DCIR::DataFrame, base_MVA::Float64)
    nCIR = length(DCIR.id)
    lower = fill(-Inf, nCIR)
    upper = fill(Inf, nCIR)
    for i in DCIR.id
        if DCIR.l_cap_1[i] != 0.0
            cap = DCIR.l_cap_1[i] / base_MVA
            lower[i] = -cap
            upper[i] = cap
        end
    end
    return lower, upper
end

"""
    DispatchConfig

Steady-state OPF builder configuration (`RunConfig.dispatch`, avenue 1).

Selects the dispatch formulation (`type_model`), the fuel-cost form, which
variable bounds and inequality families are built as explicit `@constraint`s
(dual-friendly, so `JuMP.dual` extraction stays clean), and the DC branch
susceptance convention (`susceptance_model`).

Mandatory equalities (slack angle, P/Q balance, branch flow where required
by the formulation) are always built and are not part of this struct.

See the user guide, section 2, and the configuration map, section 2, for the
full field reference.
"""
Base.@kwdef struct DispatchConfig
    type_model::String = "ACOPF"       # "ACOPF" | "DCOPF" | "ED" | "UC"
    use_matrix::Bool = true           # Ybus/Bbus matrix form for power balance
    cost_type::String = "quadratic"   # "quadratic" (MATPOWER) | "linear"
    susceptance_model::SusceptanceModel = SIMPLE  # DC-OPF only: SIMPLE = 1/x ; POWERMODELS = imag(inv(r+jx))

    # Explicit ≤-form variable bounds (dual-friendly)
    bound_V::Bool = true
    bound_θ::Bool = true
    bound_P_g::Bool = true
    bound_Q_g::Bool = true
    bound_P_ik::Bool = false
    bound_Q_ik::Bool = false
    bound_P_ki::Bool = false
    bound_Q_ki::Bool = false

    # Optional inequality constraint families
    ineq_sg_upper::Bool = false       # generator capability curve
    ineq_sbranch_upper::Bool = true     # branch thermal limits
    ineq_ang_diff_branch::Bool = true   # branch angle-difference limits

    # Promote the bus active demand from an inlined Float64 to a JuMP variable pinned
    # by `P_d[k] - p_d⁰[k] == 0`. The feasible set and the optimum are unchanged, but
    # the pin's dual σ_k IS ∂L/∂P_d(k): stationarity in P_d[k] sums the parameter
    # derivative over EVERY row that happens to contain the demand, so the demand-side
    # price needs no census of which constraint families those are. ACOPF only, active
    # power only (q_d stays a constant). Exported as `dual_Pd_pin.csv`.
    pd_as_variable::Bool = false

    # Solve the in-house dual LP after the primal (results → Dispatch_Dual/).
    # Enforced by validate_dispatch_config!: ED or DCOPF + linear cost; DCOPF also needs use_matrix=true.
    # Quadratic primal → use save_duals (JuMP.dual on primal constraints) instead.
    solve_explicit_dual::Bool = false

    # Non-case-data limit values (θ box, angle-diff clamps, unbounded-flow sentinel).
    limits::DispatchLimitsConfig = DispatchLimitsConfig()

    # How simple box limits are attached: CONSTRAINT (explicit ≤ inequalities) or
    # VARIABLE (JuMP lower_bound/upper_bound at creation). UC ignores this field.
    bound_encoding::BoundEncoding = CONSTRAINT
end

"""Default steady-state builder settings (matches legacy `Define_Inputs_4OPF`)."""
default_dispatch_config() = DispatchConfig()

function copy_dispatch_config(dc::DispatchConfig)::DispatchConfig
    return DispatchConfig(; (f => getfield(dc, f) for f in fieldnames(DispatchConfig))...)
end

"""
    build_opf_input_param(dispatch::DispatchConfig) -> OrderedDict

Build the legacy `opf_input_param` container consumed by `Make_*OPF_Model*!`.
Equality constraints are formulation-mandatory and are not part of this dict.
"""
function build_opf_input_param(dispatch::DispatchConfig)::OrderedDict{Symbol, Any}
    opf_input_param = OrderedDict{Symbol, Any}()
    opf_input_param[:var_names] = OrderedDict(
        :V => "V", :θ => "θ", :P_g => "P_g", :Q_g => "Q_g",
        :P_ik => "P_ik", :Q_ik => "Q_ik", :P_ki => "P_ki", :Q_ki => "Q_ki",
    )
    opf_input_param[:var_bounds] = OrderedDict(
        :V => dispatch.bound_V, :θ => dispatch.bound_θ,
        :P_g => dispatch.bound_P_g, :Q_g => dispatch.bound_Q_g,
        :P_ik => dispatch.bound_P_ik, :Q_ik => dispatch.bound_Q_ik,
        :P_ki => dispatch.bound_P_ki, :Q_ki => dispatch.bound_Q_ki,
    )
    opf_input_param[:obj_function] = OrderedDict(:type => dispatch.cost_type)
    opf_input_param[:susceptance_model] = dispatch.susceptance_model
    opf_input_param[:pd_as_variable] = dispatch.pd_as_variable
    opf_input_param[:ineq_cons] = OrderedDict(
        :sg_upper => dispatch.ineq_sg_upper,
        :sbranch_upper => dispatch.ineq_sbranch_upper,
        :ang_diff_branch => dispatch.ineq_ang_diff_branch,
    )
    opf_input_param[:limits] = dispatch.limits
    opf_input_param[:bound_encoding] = dispatch.bound_encoding
    return opf_input_param
end

"""
    validate_dispatch_config!(dc::DispatchConfig)

Fail fast on an incoherent steady-state builder configuration: unknown
`type_model` or `cost_type`, `UC` without linear cost, or
`solve_explicit_dual=true` combined with a formulation/cost/matrix-form it
does not support (see the user guide, section 2.3).
"""
function validate_dispatch_config!(dc::DispatchConfig)
    dc.type_model in ("ACOPF", "DCOPF", "ED", "UC") ||
        throw(ArgumentError("Unknown type_model \"$(dc.type_model)\"."))
    dc.cost_type in ("linear", "quadratic") ||
        throw(ArgumentError("cost_type must be \"linear\" or \"quadratic\"."))
    if dc.type_model == "UC"
        dc.cost_type == "linear" ||
            throw(ArgumentError("UC requires cost_type=\"linear\" (MILP)."))
    end
    if dc.solve_explicit_dual
        dc.type_model in ("DCOPF", "ED") ||
            throw(ArgumentError(
                "solve_explicit_dual=true requires type_model=\"DCOPF\" or \"ED\"."))
        dc.cost_type == "linear" ||
            throw(ArgumentError("solve_explicit_dual=true requires cost_type=\"linear\" (LP dual)."))
        if dc.type_model == "DCOPF" && !dc.use_matrix
            throw(ArgumentError(
                "solve_explicit_dual=true for DCOPF requires use_matrix=true (Bbus primal; " *
                "branch-flow primal has a different dual structure)."))
        end
    end
    validate_dispatch_limits!(dc.limits)
    return nothing
end
