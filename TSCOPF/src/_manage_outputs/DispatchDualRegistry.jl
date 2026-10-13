#=
================================================================================
 DispatchDualRegistry.jl  —  registry-driven export of steady-state OPF duals
================================================================================
 Steady-state counterpart to _transient_stability/DynDualRegistry.jl.

 Every steady-state constraint container in `opf_dict[:eq_const]` / `[:ineq_const]`
 has the same shape — `OrderedDict{Int, ConstraintRef}` (even `eq_const_angle_sw`
 is `[1] => ref`) — so extraction is always
     ids  = collect(keys(container))
     vals = [JuMP.dual(c) for (_, c) in container]
 The only thing that varies per constraint family is metadata: which id column to
 label it with (bus / generator / branch), the price unit, and the CSV/XLSX target
 names. That metadata lives here, once, in `STEADY_STATE_DUAL_SPECS`.

 `Save_Duals_OPF_Model` (functions_2_save_dispatch_duals.jl) loops this table,
 keeps the specs whose constraint key is present in the solved model, and writes
 TXT / CSV / XLSX. ACOPF, DCOPF and ED all share that one function — they differ
 only in which keys exist, and `spec_present` filters accordingly.

 To export a NEW steady-state constraint family: add one row below.

 UC restricted-pricing LP reuses this table via Save_Duals_OPF_Model (same keys as
 ED). UC-only artifacts (fixed u*, future ramp/min-up/down duals) go through
 Save_Duals_UC_Model (duals_preamble / extra_csv / extra_sheets). When UC gains
 new lp_dict families that mirror steady-state keys, add registry rows here.
================================================================================
=#

# --- metadata enums -------------------------------------------------------------

@enum DualIdKind BUS GEN BRANCH SWING        # → Bus_ID / Gen_ID / Branch_ID / (no CSV)
@enum DualUnit   U_MW U_MVAR U_MVA U_PU U_RAD

# --- one export spec ------------------------------------------------------------

struct DispatchDualSpec
    source::Symbol                      # :eq_const or :ineq_const
    key::Symbol                         # constraint family key in opf_dict[source]
    txt_name::String                    # TXT section header (preserves legacy duals.txt labels)
    id_kind::DualIdKind                 # id column kind (and TXT-only flag via SWING)
    unit::DualUnit                      # price unit (TXT prefix + docs)
    csv_file::Union{Nothing, String}    # basename under pf_dispatch_CSV_duals (nothing → TXT only)
    sheet::Union{Nothing, String}       # XLSX sheet name (nothing → no sheet)
    value_col::Symbol                   # CSV / XLSX value column header
end

# ==================================================================================
# Canonical steady-state export table
# ==================================================================================
# Order matches the legacy duals.txt section order so existing TXT reports are
# unchanged. `haskey` filtering at save time means ACOPF / DCOPF / ED each emit
# only the rows whose constraint exists in their model.
#
# `eq_const_angle_sw` is a single swing-bus scalar with no natural id column, so it
# stays TXT-only (csv = sheet = nothing). All other families now emit CSV + XLSX —
# including the branch flow duals that were previously TXT-only.

const STEADY_STATE_DUAL_SPECS = DispatchDualSpec[
    # --- equality: swing reference + nodal balance ------------------------------
    DispatchDualSpec(:eq_const, :eq_const_angle_sw, "dual_eq_const_θ_SW",
        SWING, U_RAD, nothing, nothing, :Dual_Angle_SW),
    DispatchDualSpec(:eq_const, :eq_const_p_balance, "dual_eq_const_P_balance",
        BUS, U_MW, "dual_P_balance.csv", "P_Balance", :Dual_P_Balance),
    DispatchDualSpec(:eq_const, :eq_const_q_balance, "dual_eq_const_Q_balance",
        BUS, U_MVAR, "dual_Q_balance.csv", "Q_Balance", :Dual_Q_Balance),

    # Demand pin (`DispatchConfig.pd_as_variable`, ACOPF only). σ_k is ∂L/∂P_d(k) itself:
    # the pin's stationarity sums the demand's parameter derivative over every row that
    # contains it, so no family census is needed. On a dispatch-only ACOPF it must come
    # back equal to `dual_P_balance.csv`; when transient rows also carry the demand the
    # gap between the two IS the stability content of the price.
    DispatchDualSpec(:eq_const, :eq_const_pd_pin, "dual_eq_const_Pd_pin",
        BUS, U_MW, "dual_Pd_pin.csv", "Pd_Pin", :Dual_Pd_Pin),

    # --- equality: branch flow definitions (promoted to CSV/XLSX) ---------------
    DispatchDualSpec(:eq_const, :eq_const_p_ik, "dual_eq_const_Pik",
        BRANCH, U_MW, "dual_eq_Pik.csv", "Eq_Pik", :Dual_Pik),
    DispatchDualSpec(:eq_const, :eq_const_q_ik, "dual_eq_const_Qik",
        BRANCH, U_MVAR, "dual_eq_Qik.csv", "Eq_Qik", :Dual_Qik),
    DispatchDualSpec(:eq_const, :eq_const_p_ki, "dual_eq_const_Pki",
        BRANCH, U_MW, "dual_eq_Pki.csv", "Eq_Pki", :Dual_Pki),
    DispatchDualSpec(:eq_const, :eq_const_q_ki, "dual_eq_const_Qki",
        BRANCH, U_MVAR, "dual_eq_Qki.csv", "Eq_Qki", :Dual_Qki),

    # --- inequality: apparent-power / capability limits -------------------------
    DispatchDualSpec(:ineq_const, :ineq_const_sg_upper, "dual_ineq_const_Sg_upper",
        GEN, U_MVA, "dual_Sg_Upper.csv", "Sg_Upper", :Dual_Sg_Upper),
    DispatchDualSpec(:ineq_const, :ineq_const_s_ik, "dual_ineq_const_Sik_upper",
        BRANCH, U_MVA, "dual_Sik_Upper.csv", "Sik_Upper", :Dual_Sik_Upper),
    DispatchDualSpec(:ineq_const, :ineq_const_s_ki, "dual_ineq_const_Ski_upper",
        BRANCH, U_MVA, "dual_Ski_Upper.csv", "Ski_Upper", :Dual_Ski_Upper),

    # --- inequality: grid-forming converter capability (allow_gfm ACOPF only) ---
    # Attached by `attach_gfm_acopf_limits!` and keyed by generator id, so they read
    # like any other GEN family. Both are squared forms, exactly as `ineq_const_sg_upper`
    # is; the unit tags follow the quantity being limited (nameplate current → MVA,
    # internal EMF → p.u.) rather than the algebraic square.
    DispatchDualSpec(:ineq_const, :ineq_const_gfm_Imax, "dual_ineq_const_GFM_Imax",
        GEN, U_MVA, "dual_gfm_Imax.csv", "GFM_Imax", :Dual_GFM_Imax),
    DispatchDualSpec(:ineq_const, :ineq_const_gfm_Emax, "dual_ineq_const_GFM_Emax",
        GEN, U_PU, "dual_gfm_Emax.csv", "GFM_Emax", :Dual_GFM_Emax),

    # --- inequality: angle-difference limits ------------------------------------
    DispatchDualSpec(:ineq_const, :ineq_const_ang_diff_lower, "dual_ineq_const_ang_diff_lower",
        BRANCH, U_RAD, "dual_diff_ang_Lower.csv", "Ang_Diff_Lo", :Dual_diff_ang_Lower),
    DispatchDualSpec(:ineq_const, :ineq_const_ang_diff_upper, "dual_ineq_const_ang_diff_upper",
        BRANCH, U_RAD, "dual_diff_ang_Upper.csv", "Ang_Diff_Up", :Dual_diff_ang_Upper),

    # --- inequality: variable bounds (V, θ) -------------------------------------
    DispatchDualSpec(:ineq_const, :ineq_const_volt_mag_lower, "dual_ineq_const_V_lower",
        BUS, U_PU, "dual_V_lower.csv", "V_Mag_Lo", :Dual_V_Lower),
    DispatchDualSpec(:ineq_const, :ineq_const_volt_mag_upper, "dual_ineq_const_V_upper",
        BUS, U_PU, "dual_V_upper.csv", "V_Mag_Up", :Dual_V_Upper),
    DispatchDualSpec(:ineq_const, :ineq_const_volt_ang_lower, "dual_ineq_const_θ_lower",
        BUS, U_RAD, "dual_θ_lower.csv", "V_Ang_Lo", :Dual_θ_Lower),
    DispatchDualSpec(:ineq_const, :ineq_const_volt_ang_upper, "dual_ineq_const_θ_upper",
        BUS, U_RAD, "dual_θ_upper.csv", "V_Ang_Up", :Dual_θ_Upper),

    # --- inequality: generation bounds (P, Q) -----------------------------------
    DispatchDualSpec(:ineq_const, :ineq_const_pg_lower, "dual_ineq_const_Pg_lower",
        GEN, U_MW, "dual_Pg_lower.csv", "Pg_Lo", :Dual_Pg_Lower),
    DispatchDualSpec(:ineq_const, :ineq_const_pg_upper, "dual_ineq_const_Pg_upper",
        GEN, U_MW, "dual_Pg_upper.csv", "Pg_Up", :Dual_Pg_Upper),
    DispatchDualSpec(:ineq_const, :ineq_const_qg_lower, "dual_ineq_const_Qg_lower",
        GEN, U_MVAR, "dual_Qg_lower.csv", "Qg_Lo", :Dual_Qg_Lower),
    DispatchDualSpec(:ineq_const, :ineq_const_qg_upper, "dual_ineq_const_Qg_upper",
        GEN, U_MVAR, "dual_Qg_upper.csv", "Qg_Up", :Dual_Qg_Upper),

    # --- inequality: branch flow bounds (promoted to CSV/XLSX) ------------------
    DispatchDualSpec(:ineq_const, :ineq_const_pik_lower, "dual_ineq_const_Pik_lower",
        BRANCH, U_MW, "dual_Pik_lower.csv", "Pik_Lo", :Dual_Pik_Lower),
    DispatchDualSpec(:ineq_const, :ineq_const_pik_upper, "dual_ineq_const_Pik_upper",
        BRANCH, U_MW, "dual_Pik_upper.csv", "Pik_Up", :Dual_Pik_Upper),
    DispatchDualSpec(:ineq_const, :ineq_const_qik_lower, "dual_ineq_const_Qik_lower",
        BRANCH, U_MVAR, "dual_Qik_lower.csv", "Qik_Lo", :Dual_Qik_Lower),
    DispatchDualSpec(:ineq_const, :ineq_const_qik_upper, "dual_ineq_const_Qik_upper",
        BRANCH, U_MVAR, "dual_Qik_upper.csv", "Qik_Up", :Dual_Qik_Upper),
    DispatchDualSpec(:ineq_const, :ineq_const_pki_lower, "dual_ineq_const_Pki_lower",
        BRANCH, U_MW, "dual_Pki_lower.csv", "Pki_Lo", :Dual_Pki_Lower),
    DispatchDualSpec(:ineq_const, :ineq_const_pki_upper, "dual_ineq_const_Pki_upper",
        BRANCH, U_MW, "dual_Pki_upper.csv", "Pki_Up", :Dual_Pki_Upper),
    DispatchDualSpec(:ineq_const, :ineq_const_qki_lower, "dual_ineq_const_Qki_lower",
        BRANCH, U_MVAR, "dual_Qki_lower.csv", "Qki_Lo", :Dual_Qki_Lower),
    DispatchDualSpec(:ineq_const, :ineq_const_qki_upper, "dual_ineq_const_Qki_upper",
        BRANCH, U_MVAR, "dual_Qki_upper.csv", "Qki_Up", :Dual_Qki_Upper),
]

# ==================================================================================
# Helpers
# ==================================================================================

"""CSV/XLSX id-column header for a constraint family's index kind."""
function _id_column(kind::DualIdKind)::Symbol
    kind == BUS    && return :Bus_ID
    kind == GEN    && return :Gen_ID
    kind == BRANCH && return :Branch_ID
    return :Index   # SWING / fallback (TXT-only families never reach CSV/XLSX)
end

"""TXT price-unit prefix for a dual value."""
function _unit_prefix(unit::DualUnit)::String
    unit == U_MW   && return "€/MW"
    unit == U_MVAR && return "€/MVAr"
    unit == U_MVA  && return "€/MVA"
    unit == U_PU   && return "€/p.u."
    return "€/rad"  # U_RAD
end

"""`true` when the constraint family or bound-manifest entry for `spec` exists."""
function spec_present(opf_dict::OrderedDict{Symbol, Any}, spec::DispatchDualSpec)::Bool
    if haskey(opf_dict, spec.source) && haskey(opf_dict[spec.source], spec.key)
        return true
    end
    return bound_manifest_entry_present(opf_dict, spec.key)
end

"""
    extract_spec_duals(opf_dict, spec) -> (ids, vals)

Flatten the constraint container or bound-manifest entry for `spec`.
"""
function extract_spec_duals(opf_dict::OrderedDict{Symbol, Any}, spec::DispatchDualSpec)
    if haskey(opf_dict, spec.source) && haskey(opf_dict[spec.source], spec.key)
        container = opf_dict[spec.source][spec.key]
        ids  = collect(keys(container))
        vals = [JuMP.dual(c) for (_, c) in container]
        return ids, vals
    end
    return extract_bound_manifest_duals(opf_dict, spec.key)
end
