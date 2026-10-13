#=
================================================================================
 DynDualRegistry.jl  —  registry-driven extraction of TS constraint duals
================================================================================
 Phase 1.3: after the dynamic model is built, `build_dual_registry!` records
 which constraint families are present.  `Save_Duals_Dynamic_Model_tsred` then
 iterates that list instead of hard-coded `haskey` blocks.

 Workflow
 --------
 1. Builder finishes → `build_dual_registry!(dyn_model_dict)`
 2. Solver returns optimal → `extract_registered_duals(dyn_model_dict)`
 3. Save layer writes TXT / CSV / XLSX from the returned Dict

 To export duals for a NEW constraint family: add one row to
 `CLASSICAL_KRON_DUAL_SPECS` (and store the refs in `dyn_model_dict` as today).
================================================================================
=#

# --- layout tags: how JuMP constraint containers are flattened to dual vectors ---

@enum DualExtractLayout begin
    GEN_INDEXED           # OrderedDict{gen_id, ConstraintRef}           → Vector (pre-fault init)
    TIME_INDEXED          # OrderedDict{time_step, ConstraintRef}      → Vector (COI scalars)
    TIME_INDEXED_MERGE    # fault-on + post-fault COI keys             → vcat(tf, tpf)
    PER_GEN_TIME          # OrderedDict{gen, OrderedDict{t, CR}}       → OrderedDict{gen, Vector}
    PER_GEN_TIME_MERGE    # swing / Pe trajectories: tf + tpf per gen  → merged OrderedDict
    PER_BUS_TIME_MERGE    # nodal balance duals: tf + tpf per bus
end

# --- one export spec: links dyn_model_dict keys → save-layer symbol names --------

struct DualRegistryEntry
    export_name::Symbol              # local name in Save_Duals (e.g. :dual_Pe)
    source::Symbol                   # :eq_const or :ineq_const in dyn_model_dict
    key_tf::Symbol                   # constraint key during fault-on window
    key_tpf::Union{Nothing, Symbol}  # matching post-fault key (nothing if single-window)
    layout::DualExtractLayout
    csv_file::Union{Nothing, String} # basename under pf_TS_CSV_duals (nothing → TXT/XLSX only)
    xlsx_kw::Symbol                  # keyword for Save_Duals_2_Excel_tsred
end

# ==================================================================================
# Canonical export table — classical 2nd-order Kron (TSC-ACOPF + TSC-DCOPF)
# ==================================================================================
# Rows are filtered at build time: optional families (Pm_init, ΔωCOI, E bounds, …)
# only enter dyn_model_dict[:dual_registry] when the corresponding constraints exist.

const CLASSICAL_KRON_DUAL_SPECS = DualRegistryEntry[
    # --- pre-fault initial conditions (equality) --------------------------------
    DualRegistryEntry(:dual_Pe_init, :eq_const, :eq_const_P_init, nothing,
        GEN_INDEXED, "dual_Pe_init.csv", :dual_Pe_init_xlsx),
    DualRegistryEntry(:dual_Qe_init, :eq_const, :eq_const_Q_init, nothing,
        GEN_INDEXED, "dual_Qe_init.csv", :dual_Qe_init_xlsx),
    DualRegistryEntry(:dual_Pm_init, :eq_const, :eq_const_Pm_init, nothing,
        GEN_INDEXED, "dual_Pm_init.csv", :dual_Pm_init_xlsx),  # USE_PM only

    # --- COI reference trajectories (equality, fault + post-fault merged) -------
    DualRegistryEntry(:dual_δCOI, :eq_const, :eq_const_δCOI_tf, :eq_const_δCOI_tpf,
        TIME_INDEXED_MERGE, "dual_delta_COI.csv", :dual_δCOI_xlsx),
    DualRegistryEntry(:dual_ΔωCOI, :eq_const, :eq_const_ΔωCOI_tf, :eq_const_ΔωCOI_tpf,
        TIME_INDEXED_MERGE, "dual_Delta_Omega_COI.csv", :dual_ΔωCOI_xlsx),  # constrain_Δω

    # --- swing / electrical power ODE discretization (equality, per generator) --
    DualRegistryEntry(:dual_Pe, :eq_const, :eq_const_Pe_tf, :eq_const_Pe_tpf,
        PER_GEN_TIME_MERGE, "dual_Pe.csv", :dual_Pe_xlsx),
    DualRegistryEntry(:dual_δ, :eq_const, :eq_const_δ_tf, :eq_const_δ_tpf,
        PER_GEN_TIME_MERGE, "dual_delta.csv", :dual_δ_xlsx),
    DualRegistryEntry(:dual_Δω, :eq_const, :eq_const_Δω_tf, :eq_const_Δω_tpf,
        PER_GEN_TIME_MERGE, "dual_Delta_Omega.csv", :dual_Δω_xlsx),

    # --- transient-stability bounds w.r.t. COI (inequality, per generator) ------
    DualRegistryEntry(:dual_δ_COI_lower, :ineq_const, :ineq_const_δ_COI_tf_lower,
        :ineq_const_δ_COI_tpf_lower, PER_GEN_TIME_MERGE,
        "dual_delta_COI_lower.csv", :dual_δ_COI_lower_xlsx),
    DualRegistryEntry(:dual_δ_COI_upper, :ineq_const, :ineq_const_δ_COI_tf_upper,
        :ineq_const_δ_COI_tpf_upper, PER_GEN_TIME_MERGE,
        "dual_delta_COI_upper.csv", :dual_δ_COI_upper_xlsx),
    DualRegistryEntry(:dual_Δω_COI_lower, :ineq_const, :ineq_const_Δω_COI_tf_lower,
        :ineq_const_Δω_COI_tpf_lower, PER_GEN_TIME_MERGE,
        "dual_Delta_Omega_COI_lower.csv", :dual_Δω_COI_lower_xlsx),
    DualRegistryEntry(:dual_Δω_COI_upper, :ineq_const, :ineq_const_Δω_COI_tf_upper,
        :ineq_const_Δω_COI_tpf_upper, PER_GEN_TIME_MERGE,
        "dual_Delta_Omega_COI_upper.csv", :dual_Δω_COI_upper_xlsx),

    # --- transient-stability bounds w.r.t. a reference machine ------------------
    # Built instead of the COI rows above when bound_style_δ is :highest_H or :ref_gen.
    # Separate export names on purpose: the shadow price of "stay within δ_tol of the COI"
    # and of "stay within δ_tol of machine k" are different quantities, and a study that
    # compares the two conventions must be able to tell the columns apart. Only one of the
    # two families exists in any given run, so nothing is exported twice.
    DualRegistryEntry(:dual_δ_ref_lower, :ineq_const, :ineq_const_δ_ref_tf_lower,
        :ineq_const_δ_ref_tpf_lower, PER_GEN_TIME_MERGE,
        "dual_delta_ref_lower.csv", :dual_δ_ref_lower_xlsx),
    DualRegistryEntry(:dual_δ_ref_upper, :ineq_const, :ineq_const_δ_ref_tf_upper,
        :ineq_const_δ_ref_tpf_upper, PER_GEN_TIME_MERGE,
        "dual_delta_ref_upper.csv", :dual_δ_ref_upper_xlsx),

    # --- pairwise angle-separation corridor (bound_style_δ = :pairwise) --------
    # Same PER_GEN_TIME_MERGE shape, but the outer key is the pair `(i, j)` rather than a
    # generator id — the only two-index family in the model. The extractor is key-generic;
    # the CSV/XLSX writers turn the tuple into a `Pair_G<i>_G<j>` column. A non-zero column
    # names the two units the corridor is holding apart, which is the whole reason this
    # family is exported separately from δ_ref rather than merged into it.
    DualRegistryEntry(:dual_δ_pair_lower, :ineq_const, :ineq_const_δ_pair_tf_lower,
        :ineq_const_δ_pair_tpf_lower, PER_GEN_TIME_MERGE,
        "dual_delta_pair_lower.csv", :dual_δ_pair_lower_xlsx),
    DualRegistryEntry(:dual_δ_pair_upper, :ineq_const, :ineq_const_δ_pair_tf_upper,
        :ineq_const_δ_pair_tpf_upper, PER_GEN_TIME_MERGE,
        "dual_delta_pair_upper.csv", :dual_δ_pair_upper_xlsx),

    # --- absolute speed corridor (bound_style_Δω = :abs) ------------------------
    # Priced against a fixed band rather than against the COI, so this dual is the
    # shadow price of the fleet's common-mode frequency excursion as well as of the
    # spread between machines. Kept apart from the COI columns for the same reason.
    DualRegistryEntry(:dual_Δω_abs_lower, :ineq_const, :ineq_const_Δω_abs_tf_lower,
        :ineq_const_Δω_abs_tpf_lower, PER_GEN_TIME_MERGE,
        "dual_Delta_Omega_abs_lower.csv", :dual_Δω_abs_lower_xlsx),
    DualRegistryEntry(:dual_Δω_abs_upper, :ineq_const, :ineq_const_Δω_abs_tf_upper,
        :ineq_const_Δω_abs_tpf_upper, PER_GEN_TIME_MERGE,
        "dual_Delta_Omega_abs_upper.csv", :dual_Δω_abs_upper_xlsx),

    # --- explicit variable bound duals (inequality, pre-fault) ------------------
    DualRegistryEntry(:dual_LB_E, :ineq_const, :ineq_const_E_lower, nothing,
        GEN_INDEXED, "dual_LB_E.csv", :dual_LB_E_xlsx),   # TSC-ACOPF (E magnitude)
    DualRegistryEntry(:dual_UB_E, :ineq_const, :ineq_const_E_upper, nothing,
        GEN_INDEXED, "dual_UB_E.csv", :dual_UB_E_xlsx),
    DualRegistryEntry(:dual_LB_δ, :ineq_const, :ineq_const_δ_lower, nothing,
        GEN_INDEXED, "dual_LB_δ.csv", :dual_LB_δ_xlsx),   # TSC-DCOPF linear path
    DualRegistryEntry(:dual_UB_δ, :ineq_const, :ineq_const_δ_upper, nothing,
        GEN_INDEXED, "dual_UB_δ.csv", :dual_UB_δ_xlsx),
    DualRegistryEntry(:dual_LB_Pm, :ineq_const, :ineq_const_P_m_lower, nothing,
        GEN_INDEXED, "dual_LB_Pm.csv", :dual_LB_Pm_xlsx), # USE_PM only
    DualRegistryEntry(:dual_UB_Pm, :ineq_const, :ineq_const_P_m_upper, nothing,
        GEN_INDEXED, "dual_UB_Pm.csv", :dual_UB_Pm_xlsx),

    # --- optional explicit tf/tpf variable-bound duals (toggle-gated) -----------
    # Every family below is attached in both windows by `attach_fault_tf_var_bounds!`
    # and `attach_postfault_tpf_var_bounds!`, so both keys belong on the row: the
    # `nothing` that used to sit in the tpf slot silently exported the fault half only.
    DualRegistryEntry(:dual_LB_δ_tf, :ineq_const, :ineq_const_δ_tf_lower,
        :ineq_const_δ_tpf_lower, PER_GEN_TIME_MERGE,
        "dual_LB_δ_tf.csv", :dual_LB_δ_tf_xlsx),
    DualRegistryEntry(:dual_UB_δ_tf, :ineq_const, :ineq_const_δ_tf_upper,
        :ineq_const_δ_tpf_upper, PER_GEN_TIME_MERGE,
        "dual_UB_δ_tf.csv", :dual_UB_δ_tf_xlsx),
    DualRegistryEntry(:dual_LB_Pe_tf, :ineq_const, :ineq_const_Pe_tf_lower,
        :ineq_const_Pe_tpf_lower, PER_GEN_TIME_MERGE,
        "dual_LB_Pe_tf.csv", :dual_LB_Pe_tf_xlsx),
    DualRegistryEntry(:dual_UB_Pe_tf, :ineq_const, :ineq_const_Pe_tf_upper,
        :ineq_const_Pe_tpf_upper, PER_GEN_TIME_MERGE,
        "dual_UB_Pe_tf.csv", :dual_UB_Pe_tf_xlsx),
    DualRegistryEntry(:dual_LB_Δω_tf, :ineq_const, :ineq_const_Δω_tf_lower,
        :ineq_const_Δω_tpf_lower, PER_GEN_TIME_MERGE,
        "dual_LB_Δω_tf.csv", :dual_LB_Δω_tf_xlsx),
    DualRegistryEntry(:dual_UB_Δω_tf, :ineq_const, :ineq_const_Δω_tf_upper,
        :ineq_const_Δω_tpf_upper, PER_GEN_TIME_MERGE,
        "dual_UB_Δω_tf.csv", :dual_UB_Δω_tf_xlsx),
    # Qe boxes target a JuMP expression on the Kron path, so they are always
    # CONSTRAINT-encoded (`ts_bound_limits.jl` announces the downgrade with @info).
    DualRegistryEntry(:dual_LB_Qe_tf, :ineq_const, :ineq_const_Qe_tf_lower,
        :ineq_const_Qe_tpf_lower, PER_GEN_TIME_MERGE,
        "dual_LB_Qe_tf.csv", :dual_LB_Qe_tf_xlsx),
    DualRegistryEntry(:dual_UB_Qe_tf, :ineq_const, :ineq_const_Qe_tf_upper,
        :ineq_const_Qe_tpf_upper, PER_GEN_TIME_MERGE,
        "dual_UB_Qe_tf.csv", :dual_UB_Qe_tf_xlsx),
    # δ_COI reference box: a COI scalar per time step, hence TIME_INDEXED_MERGE.
    DualRegistryEntry(:dual_LB_δCOI_tf, :ineq_const, :ineq_const_δCOI_tf_lower,
        :ineq_const_δCOI_tpf_lower, TIME_INDEXED_MERGE,
        "dual_LB_δCOI_tf.csv", :dual_LB_δCOI_tf_xlsx),
    DualRegistryEntry(:dual_UB_δCOI_tf, :ineq_const, :ineq_const_δCOI_tf_upper,
        :ineq_const_δCOI_tpf_upper, TIME_INDEXED_MERGE,
        "dual_UB_δCOI_tf.csv", :dual_UB_δCOI_tf_xlsx),

    # --- turbine governor (equality; only present when include_governor) --------
    DualRegistryEntry(:dual_Pref_init, :eq_const, :eq_const_Pref_init, nothing,
        GEN_INDEXED, "dual_Pref_init.csv", :dual_Pref_init_xlsx),
    DualRegistryEntry(:dual_gov_valve, :eq_const, :eq_const_gov_valve_tf, :eq_const_gov_valve_tpf,
        PER_GEN_TIME_MERGE, "dual_gov_valve.csv", :dual_gov_valve_xlsx),
    DualRegistryEntry(:dual_gov_mech, :eq_const, :eq_const_gov_mech_tf, :eq_const_gov_mech_tpf,
        PER_GEN_TIME_MERGE, "dual_gov_mech.csv", :dual_gov_mech_xlsx),
    # Valve limiter families — present only for the matching `governor_limiter` mode.
    # GOV_SMOOTH: the clamp equality `Pv = smoothclip(Pv_raw)` (same role as the AVR
    # field softsat, `dual_avr_E_fd_sat`). GOV_HARD_BOUND: the explicit ≤-form pair on
    # the raw valve state. GOV_NO_LIMIT adds neither. Keys are built by
    # `_store_gov_limit!` in `functions_4_TS_governor.jl`.
    DualRegistryEntry(:dual_gov_valve_sat, :eq_const, :eq_const_gov_valve_limit_tf,
        :eq_const_gov_valve_limit_tpf, PER_GEN_TIME_MERGE,
        "dual_gov_valve_sat.csv", :dual_gov_valve_sat_xlsx),
    DualRegistryEntry(:dual_LB_gov_valve, :ineq_const, :ineq_const_gov_valve_tf_lower,
        :ineq_const_gov_valve_tpf_lower, PER_GEN_TIME_MERGE,
        "dual_LB_gov_valve.csv", :dual_LB_gov_valve_xlsx),
    DualRegistryEntry(:dual_UB_gov_valve, :ineq_const, :ineq_const_gov_valve_tf_upper,
        :ineq_const_gov_valve_tpf_upper, PER_GEN_TIME_MERGE,
        "dual_UB_gov_valve.csv", :dual_UB_gov_valve_xlsx),
]

# Full-bus classical (Phase 3) — extends Kron specs with Q and nodal-balance families.
# (`dual_Qe_init` is already in CLASSICAL_KRON_DUAL_SPECS when eq_const_Q_init is built.)
const CLASSICAL_FULL_BUS_DUAL_SPECS = DualRegistryEntry[
    CLASSICAL_KRON_DUAL_SPECS...,
    DualRegistryEntry(:dual_Qe, :eq_const, :eq_const_Qe_tf, :eq_const_Qe_tpf,
        PER_GEN_TIME_MERGE, "dual_Qe.csv", :dual_Qe_xlsx),
    DualRegistryEntry(:dual_Pbalance, :eq_const, :eq_const_Pbalance_tf, :eq_const_Pbalance_tpf,
        PER_BUS_TIME_MERGE, "dual_Pbalance.csv", :dual_Pbalance_xlsx),
    DualRegistryEntry(:dual_Qbalance, :eq_const, :eq_const_Qbalance_tf, :eq_const_Qbalance_tpf,
        PER_BUS_TIME_MERGE, "dual_Qbalance.csv", :dual_Qbalance_xlsx),
    DualRegistryEntry(:dual_LB_V, :ineq_const, :ineq_const_V_tf_lower, :ineq_const_V_tpf_lower,
        PER_BUS_TIME_MERGE, "dual_V_lower.csv", :dual_LB_V_xlsx),
    DualRegistryEntry(:dual_LB_P_ref, :ineq_const, :ineq_const_P_ref_lower, nothing,
        GEN_INDEXED, "dual_LB_P_ref.csv", :dual_LB_P_ref_xlsx),
    DualRegistryEntry(:dual_UB_P_ref, :ineq_const, :ineq_const_P_ref_upper, nothing,
        GEN_INDEXED, "dual_UB_P_ref.csv", :dual_UB_P_ref_xlsx),
    DualRegistryEntry(:dual_LB_Pv, :ineq_const, :ineq_const_Pv_tf_lower, :ineq_const_Pv_tpf_lower,
        PER_GEN_TIME_MERGE, "dual_LB_Pv.csv", :dual_LB_Pv_xlsx),
    DualRegistryEntry(:dual_UB_Pv, :ineq_const, :ineq_const_Pv_tf_upper, :ineq_const_Pv_tpf_upper,
        PER_GEN_TIME_MERGE, "dual_UB_Pv.csv", :dual_UB_Pv_xlsx),
    DualRegistryEntry(:dual_LB_Pm_tf, :ineq_const, :ineq_const_Pm_tf_lower, :ineq_const_Pm_tpf_lower,
        PER_GEN_TIME_MERGE, "dual_LB_Pm_tf.csv", :dual_LB_Pm_tf_xlsx),
    DualRegistryEntry(:dual_UB_Pm_tf, :ineq_const, :ineq_const_Pm_tf_upper, :ineq_const_Pm_tpf_upper,
        PER_GEN_TIME_MERGE, "dual_UB_Pm_tf.csv", :dual_UB_Pm_tf_xlsx),
]

# DQ_4TH FULL_BUS — extends classical full-bus specs with dq-machine init, algebra, and EMF duals.
# Classical `dual_LB_E` / `dual_UB_E` are omitted (DQ uses `E_fd` bound keys instead).
const DQ_4TH_FULL_BUS_DUAL_SPECS = DualRegistryEntry[
    [e for e in CLASSICAL_FULL_BUS_DUAL_SPECS
        if e.export_name ∉ (:dual_LB_E, :dual_UB_E)]...,
    DualRegistryEntry(:dual_Ed_init, :eq_const, :eq_const_Ed_init, nothing,
        GEN_INDEXED, "dual_Ed_init.csv", :dual_Ed_init_xlsx),
    DualRegistryEntry(:dual_Eq_init, :eq_const, :eq_const_Eq_init, nothing,
        GEN_INDEXED, "dual_Eq_init.csv", :dual_Eq_init_xlsx),
    DualRegistryEntry(:dual_Vd_init, :eq_const, :eq_const_Vd_init, nothing,
        GEN_INDEXED, "dual_Vd_init.csv", :dual_Vd_init_xlsx),
    DualRegistryEntry(:dual_Vq_init, :eq_const, :eq_const_Vq_init, nothing,
        GEN_INDEXED, "dual_Vq_init.csv", :dual_Vq_init_xlsx),
    DualRegistryEntry(:dual_LB_E_fd, :ineq_const, :ineq_const_E_fd_lower, nothing,
        GEN_INDEXED, "dual_LB_E_fd.csv", :dual_LB_E_fd_xlsx),
    DualRegistryEntry(:dual_UB_E_fd, :ineq_const, :ineq_const_E_fd_upper, nothing,
        GEN_INDEXED, "dual_UB_E_fd.csv", :dual_UB_E_fd_xlsx),
    DualRegistryEntry(:dual_Te, :eq_const, :eq_const_Te_tf, :eq_const_Te_tpf,
        PER_GEN_TIME_MERGE, "dual_Te.csv", :dual_Te_xlsx),
    DualRegistryEntry(:dual_Vd, :eq_const, :eq_const_Vd_tf, :eq_const_Vd_tpf,
        PER_GEN_TIME_MERGE, "dual_Vd.csv", :dual_Vd_xlsx),
    DualRegistryEntry(:dual_Vq, :eq_const, :eq_const_Vq_tf, :eq_const_Vq_tpf,
        PER_GEN_TIME_MERGE, "dual_Vq.csv", :dual_Vq_xlsx),
    DualRegistryEntry(:dual_Ed, :eq_const, :eq_const_Ed_tf, :eq_const_Ed_tpf,
        PER_GEN_TIME_MERGE, "dual_Ed.csv", :dual_Ed_xlsx),
    DualRegistryEntry(:dual_Eq, :eq_const, :eq_const_Eq_tf, :eq_const_Eq_tpf,
        PER_GEN_TIME_MERGE, "dual_Eq.csv", :dual_Eq_xlsx),
    DualRegistryEntry(:dual_Vref_init, :eq_const, :eq_const_Efd_init, nothing,
        GEN_INDEXED, "dual_Vref_init.csv", :dual_Vref_init_xlsx),
    # SEXS lead-lag stage — absent whenever every generator sets Ta_exc = Tb_exc = 0
    # (the pass-through bypass in `_avr_leadlag_gens`), so the presence filter drops it.
    # Rows are normalized to a unit E_LL[t] coefficient, which keeps these duals
    # comparable across machines with different Tb_exc.
    DualRegistryEntry(:dual_avr_leadlag, :eq_const,
        :eq_const_avr_leadlag_tf, :eq_const_avr_leadlag_tpf,
        PER_GEN_TIME_MERGE, "dual_avr_leadlag.csv", :dual_avr_leadlag_xlsx),
    DualRegistryEntry(:dual_avr_E_fd, :eq_const, :eq_const_E_fd_unlim_tf, :eq_const_E_fd_unlim_tpf,
        PER_GEN_TIME_MERGE, "dual_avr_E_fd.csv", :dual_avr_E_fd_xlsx),
    DualRegistryEntry(:dual_avr_E_fd_sat, :eq_const, :eq_const_E_fd_tf, :eq_const_E_fd_tpf,
        PER_GEN_TIME_MERGE, "dual_avr_E_fd_sat.csv", :dual_avr_E_fd_sat_xlsx),
    DualRegistryEntry(:dual_LB_Ed, :ineq_const, :ineq_const_Ed_lower, nothing,
        GEN_INDEXED, "dual_LB_Ed.csv", :dual_LB_Ed_xlsx),
    DualRegistryEntry(:dual_UB_Ed, :ineq_const, :ineq_const_Ed_upper, nothing,
        GEN_INDEXED, "dual_UB_Ed.csv", :dual_UB_Ed_xlsx),
    DualRegistryEntry(:dual_LB_Eq, :ineq_const, :ineq_const_Eq_lower, nothing,
        GEN_INDEXED, "dual_LB_Eq.csv", :dual_LB_Eq_xlsx),
    DualRegistryEntry(:dual_UB_Eq, :ineq_const, :ineq_const_Eq_upper, nothing,
        GEN_INDEXED, "dual_UB_Eq.csv", :dual_UB_Eq_xlsx),
    DualRegistryEntry(:dual_LB_Id, :ineq_const, :ineq_const_Id_lower, nothing,
        GEN_INDEXED, "dual_LB_Id.csv", :dual_LB_Id_xlsx),
    DualRegistryEntry(:dual_UB_Id, :ineq_const, :ineq_const_Id_upper, nothing,
        GEN_INDEXED, "dual_UB_Id.csv", :dual_UB_Id_xlsx),
    DualRegistryEntry(:dual_LB_Iq, :ineq_const, :ineq_const_Iq_lower, nothing,
        GEN_INDEXED, "dual_LB_Iq.csv", :dual_LB_Iq_xlsx),
    DualRegistryEntry(:dual_UB_Iq, :ineq_const, :ineq_const_Iq_upper, nothing,
        GEN_INDEXED, "dual_UB_Iq.csv", :dual_UB_Iq_xlsx),
    DualRegistryEntry(:dual_LB_V_ref, :ineq_const, :ineq_const_V_ref_lower, nothing,
        GEN_INDEXED, "dual_LB_V_ref.csv", :dual_LB_V_ref_xlsx),
    DualRegistryEntry(:dual_UB_V_ref, :ineq_const, :ineq_const_V_ref_upper, nothing,
        GEN_INDEXED, "dual_UB_V_ref.csv", :dual_UB_V_ref_xlsx),
    DualRegistryEntry(:dual_LB_E_fd_unlim, :ineq_const, :ineq_const_E_fd_unlim_tf_lower,
        :ineq_const_E_fd_unlim_tpf_lower, PER_GEN_TIME_MERGE,
        "dual_LB_E_fd_unlim.csv", :dual_LB_E_fd_unlim_xlsx),
    DualRegistryEntry(:dual_UB_E_fd_unlim, :ineq_const, :ineq_const_E_fd_unlim_tf_upper,
        :ineq_const_E_fd_unlim_tpf_upper, PER_GEN_TIME_MERGE,
        "dual_UB_E_fd_unlim.csv", :dual_UB_E_fd_unlim_xlsx),
    DualRegistryEntry(:dual_LB_Ed_tf, :ineq_const, :ineq_const_Ed_tf_lower, :ineq_const_Ed_tpf_lower,
        PER_GEN_TIME_MERGE, "dual_LB_Ed_tf.csv", :dual_LB_Ed_tf_xlsx),
    DualRegistryEntry(:dual_UB_Ed_tf, :ineq_const, :ineq_const_Ed_tf_upper, :ineq_const_Ed_tpf_upper,
        PER_GEN_TIME_MERGE, "dual_UB_Ed_tf.csv", :dual_UB_Ed_tf_xlsx),
    # `attach_dq_fault_tf_var_bounds!` wires Ed, Eq, Id, Iq and Te alike; only Ed and
    # Te had rows here, so an Eq / Id / Iq box could bind with nothing exported.
    DualRegistryEntry(:dual_LB_Eq_tf, :ineq_const, :ineq_const_Eq_tf_lower, :ineq_const_Eq_tpf_lower,
        PER_GEN_TIME_MERGE, "dual_LB_Eq_tf.csv", :dual_LB_Eq_tf_xlsx),
    DualRegistryEntry(:dual_UB_Eq_tf, :ineq_const, :ineq_const_Eq_tf_upper, :ineq_const_Eq_tpf_upper,
        PER_GEN_TIME_MERGE, "dual_UB_Eq_tf.csv", :dual_UB_Eq_tf_xlsx),
    DualRegistryEntry(:dual_LB_Id_tf, :ineq_const, :ineq_const_Id_tf_lower, :ineq_const_Id_tpf_lower,
        PER_GEN_TIME_MERGE, "dual_LB_Id_tf.csv", :dual_LB_Id_tf_xlsx),
    DualRegistryEntry(:dual_UB_Id_tf, :ineq_const, :ineq_const_Id_tf_upper, :ineq_const_Id_tpf_upper,
        PER_GEN_TIME_MERGE, "dual_UB_Id_tf.csv", :dual_UB_Id_tf_xlsx),
    DualRegistryEntry(:dual_LB_Iq_tf, :ineq_const, :ineq_const_Iq_tf_lower, :ineq_const_Iq_tpf_lower,
        PER_GEN_TIME_MERGE, "dual_LB_Iq_tf.csv", :dual_LB_Iq_tf_xlsx),
    DualRegistryEntry(:dual_UB_Iq_tf, :ineq_const, :ineq_const_Iq_tf_upper, :ineq_const_Iq_tpf_upper,
        PER_GEN_TIME_MERGE, "dual_UB_Iq_tf.csv", :dual_UB_Iq_tf_xlsx),
    DualRegistryEntry(:dual_LB_Te_tf, :ineq_const, :ineq_const_Te_tf_lower, :ineq_const_Te_tpf_lower,
        PER_GEN_TIME_MERGE, "dual_LB_Te_tf.csv", :dual_LB_Te_tf_xlsx),
    DualRegistryEntry(:dual_UB_Te_tf, :ineq_const, :ineq_const_Te_tf_upper, :ineq_const_Te_tpf_upper,
        PER_GEN_TIME_MERGE, "dual_UB_Te_tf.csv", :dual_UB_Te_tf_xlsx),
    # --- GFM pre-fault init (allow_gfm; GEN_INDEXED over GFM ids only) ----------
    DualRegistryEntry(:dual_gfm_Pmeas_init, :eq_const, :eq_const_gfm_Pmeas_init, nothing,
        GEN_INDEXED, "dual_gfm_Pmeas_init.csv", :dual_gfm_Pmeas_init_xlsx),
    DualRegistryEntry(:dual_gfm_Qmeas_init, :eq_const, :eq_const_gfm_Qmeas_init, nothing,
        GEN_INDEXED, "dual_gfm_Qmeas_init.csv", :dual_gfm_Qmeas_init_xlsx),
    DualRegistryEntry(:dual_gfm_Vmeas_init, :eq_const, :eq_const_gfm_Vmeas_init, nothing,
        GEN_INDEXED, "dual_gfm_Vmeas_init.csv", :dual_gfm_Vmeas_init_xlsx),
    DualRegistryEntry(:dual_gfm_Vset_init, :eq_const, :eq_const_gfm_Vset_init, nothing,
        GEN_INDEXED, "dual_gfm_Vset_init.csv", :dual_gfm_Vset_init_xlsx),
    # --- GFM pre-fault boxes (allow_gfm; `Attach_GFM_init!`) --------------------
    DualRegistryEntry(:dual_LB_gfm_δ, :ineq_const, :ineq_const_gfm_δ_lower, nothing,
        GEN_INDEXED, "dual_LB_gfm_δ.csv", :dual_LB_gfm_δ_xlsx),
    DualRegistryEntry(:dual_UB_gfm_δ, :ineq_const, :ineq_const_gfm_δ_upper, nothing,
        GEN_INDEXED, "dual_UB_gfm_δ.csv", :dual_UB_gfm_δ_xlsx),
    # --- GFM transient physics (allow_gfm; PER_GEN_TIME over GFM ids only) -------
    # Every family in `_GFM_EQ_FAMILIES` now has a row. The measurement filters, the
    # δ integration step and the raw-E definitions used to be excluded as "bookkeeping",
    # but that does not survive comparison with the SG path: `delta` is the converter's
    # counterpart of `eq_const_δ`, and `Eint_raw` / `Edroop_raw` are the counterparts of
    # the AVR's `eq_const_E_fd_unlim` — both of which are exported. The filter rows carry
    # the marginal value of the measured P/Q/V a converter droops on, which is exactly
    # the quantity a GFM economics study asks about.
    DualRegistryEntry(:dual_gfm_filter_P, :eq_const,
        :eq_const_gfm_filter_P_tf, :eq_const_gfm_filter_P_tpf,
        PER_GEN_TIME_MERGE, "dual_gfm_filter_P.csv", :dual_gfm_filter_P_xlsx),
    DualRegistryEntry(:dual_gfm_filter_Q, :eq_const,
        :eq_const_gfm_filter_Q_tf, :eq_const_gfm_filter_Q_tpf,
        PER_GEN_TIME_MERGE, "dual_gfm_filter_Q.csv", :dual_gfm_filter_Q_xlsx),
    DualRegistryEntry(:dual_gfm_filter_V, :eq_const,
        :eq_const_gfm_filter_V_tf, :eq_const_gfm_filter_V_tpf,
        PER_GEN_TIME_MERGE, "dual_gfm_filter_V.csv", :dual_gfm_filter_V_xlsx),
    DualRegistryEntry(:dual_gfm_delta, :eq_const,
        :eq_const_gfm_delta_tf, :eq_const_gfm_delta_tpf,
        PER_GEN_TIME_MERGE, "dual_gfm_delta.csv", :dual_gfm_delta_xlsx),
    DualRegistryEntry(:dual_gfm_Eint_raw, :eq_const,
        :eq_const_gfm_Eint_raw_tf, :eq_const_gfm_Eint_raw_tpf,
        PER_GEN_TIME_MERGE, "dual_gfm_Eint_raw.csv", :dual_gfm_Eint_raw_xlsx),
    DualRegistryEntry(:dual_gfm_Edroop_raw, :eq_const,
        :eq_const_gfm_Edroop_raw_tf, :eq_const_gfm_Edroop_raw_tpf,
        PER_GEN_TIME_MERGE, "dual_gfm_Edroop_raw.csv", :dual_gfm_Edroop_raw_xlsx),
    DualRegistryEntry(:dual_gfm_droop, :eq_const,
        :eq_const_gfm_droop_tf, :eq_const_gfm_droop_tpf,
        PER_GEN_TIME_MERGE, "dual_gfm_droop.csv", :dual_gfm_droop_xlsx),
    DualRegistryEntry(:dual_gfm_limiter_Id, :eq_const,
        :eq_const_gfm_limiter_Id_tf, :eq_const_gfm_limiter_Id_tpf,
        PER_GEN_TIME_MERGE, "dual_gfm_limiter_Id.csv", :dual_gfm_limiter_Id_xlsx),
    DualRegistryEntry(:dual_gfm_limiter_Iq, :eq_const,
        :eq_const_gfm_limiter_Iq_tf, :eq_const_gfm_limiter_Iq_tpf,
        PER_GEN_TIME_MERGE, "dual_gfm_limiter_Iq.csv", :dual_gfm_limiter_Iq_xlsx),
    DualRegistryEntry(:dual_gfm_Eint_clip, :eq_const,
        :eq_const_gfm_Eint_clip_tf, :eq_const_gfm_Eint_clip_tpf,
        PER_GEN_TIME_MERGE, "dual_gfm_Eint_clip.csv", :dual_gfm_Eint_clip_xlsx),
    DualRegistryEntry(:dual_gfm_Edroop_clip, :eq_const,
        :eq_const_gfm_Edroop_clip_tf, :eq_const_gfm_Edroop_clip_tpf,
        PER_GEN_TIME_MERGE, "dual_gfm_Edroop_clip.csv", :dual_gfm_Edroop_clip_xlsx),
    DualRegistryEntry(:dual_gfm_Pe, :eq_const,
        :eq_const_gfm_Pe_tf, :eq_const_gfm_Pe_tpf,
        PER_GEN_TIME_MERGE, "dual_gfm_Pe.csv", :dual_gfm_Pe_xlsx),
    DualRegistryEntry(:dual_gfm_Qe, :eq_const,
        :eq_const_gfm_Qe_tf, :eq_const_gfm_Qe_tpf,
        PER_GEN_TIME_MERGE, "dual_gfm_Qe.csv", :dual_gfm_Qe_xlsx),
    # P/Q overload integrators. Present only when `gfm_pq_overload` is on; a missing
    # key is skipped by `registry_entry_present`, so the knob-off export is unchanged.
    DualRegistryEntry(:dual_gfm_Pmax_raw, :eq_const,
        :eq_const_gfm_Pmax_raw_tf, :eq_const_gfm_Pmax_raw_tpf,
        PER_GEN_TIME_MERGE, "dual_gfm_Pmax_raw.csv", :dual_gfm_Pmax_raw_xlsx),
    DualRegistryEntry(:dual_gfm_Pmax_clip, :eq_const,
        :eq_const_gfm_Pmax_clip_tf, :eq_const_gfm_Pmax_clip_tpf,
        PER_GEN_TIME_MERGE, "dual_gfm_Pmax_clip.csv", :dual_gfm_Pmax_clip_xlsx),
    DualRegistryEntry(:dual_gfm_Pmin_raw, :eq_const,
        :eq_const_gfm_Pmin_raw_tf, :eq_const_gfm_Pmin_raw_tpf,
        PER_GEN_TIME_MERGE, "dual_gfm_Pmin_raw.csv", :dual_gfm_Pmin_raw_xlsx),
    DualRegistryEntry(:dual_gfm_Pmin_clip, :eq_const,
        :eq_const_gfm_Pmin_clip_tf, :eq_const_gfm_Pmin_clip_tpf,
        PER_GEN_TIME_MERGE, "dual_gfm_Pmin_clip.csv", :dual_gfm_Pmin_clip_xlsx),
    DualRegistryEntry(:dual_gfm_Qmax_raw, :eq_const,
        :eq_const_gfm_Qmax_raw_tf, :eq_const_gfm_Qmax_raw_tpf,
        PER_GEN_TIME_MERGE, "dual_gfm_Qmax_raw.csv", :dual_gfm_Qmax_raw_xlsx),
    DualRegistryEntry(:dual_gfm_Qmax_clip, :eq_const,
        :eq_const_gfm_Qmax_clip_tf, :eq_const_gfm_Qmax_clip_tpf,
        PER_GEN_TIME_MERGE, "dual_gfm_Qmax_clip.csv", :dual_gfm_Qmax_clip_xlsx),
    DualRegistryEntry(:dual_gfm_Qmin_raw, :eq_const,
        :eq_const_gfm_Qmin_raw_tf, :eq_const_gfm_Qmin_raw_tpf,
        PER_GEN_TIME_MERGE, "dual_gfm_Qmin_raw.csv", :dual_gfm_Qmin_raw_xlsx),
    DualRegistryEntry(:dual_gfm_Qmin_clip, :eq_const,
        :eq_const_gfm_Qmin_clip_tf, :eq_const_gfm_Qmin_clip_tpf,
        PER_GEN_TIME_MERGE, "dual_gfm_Qmin_clip.csv", :dual_gfm_Qmin_clip_xlsx),
    # --- GFM boxes. Under CONSTRAINT these are ≤-rows in :ineq_const; under VARIABLE
    # the same keys resolve through the bound manifest, so one row serves both.
    # The P/Q/V measurement boxes are wide guard-rails that should stay slack — which
    # is precisely why they are exported: a non-zero dual there says a guard-rail is
    # carrying load it was never meant to carry, and the run needs a second look.
    DualRegistryEntry(:dual_LB_gfm_P_meas_tf, :ineq_const,
        :ineq_const_gfm_P_meas_tf_lower, :ineq_const_gfm_P_meas_tpf_lower,
        PER_GEN_TIME_MERGE, "dual_LB_gfm_P_meas_tf.csv", :dual_LB_gfm_P_meas_tf_xlsx),
    DualRegistryEntry(:dual_UB_gfm_P_meas_tf, :ineq_const,
        :ineq_const_gfm_P_meas_tf_upper, :ineq_const_gfm_P_meas_tpf_upper,
        PER_GEN_TIME_MERGE, "dual_UB_gfm_P_meas_tf.csv", :dual_UB_gfm_P_meas_tf_xlsx),
    DualRegistryEntry(:dual_LB_gfm_Q_meas_tf, :ineq_const,
        :ineq_const_gfm_Q_meas_tf_lower, :ineq_const_gfm_Q_meas_tpf_lower,
        PER_GEN_TIME_MERGE, "dual_LB_gfm_Q_meas_tf.csv", :dual_LB_gfm_Q_meas_tf_xlsx),
    DualRegistryEntry(:dual_UB_gfm_Q_meas_tf, :ineq_const,
        :ineq_const_gfm_Q_meas_tf_upper, :ineq_const_gfm_Q_meas_tpf_upper,
        PER_GEN_TIME_MERGE, "dual_UB_gfm_Q_meas_tf.csv", :dual_UB_gfm_Q_meas_tf_xlsx),
    DualRegistryEntry(:dual_LB_gfm_V_meas_tf, :ineq_const,
        :ineq_const_gfm_V_meas_tf_lower, :ineq_const_gfm_V_meas_tpf_lower,
        PER_GEN_TIME_MERGE, "dual_LB_gfm_V_meas_tf.csv", :dual_LB_gfm_V_meas_tf_xlsx),
    DualRegistryEntry(:dual_UB_gfm_V_meas_tf, :ineq_const,
        :ineq_const_gfm_V_meas_tf_upper, :ineq_const_gfm_V_meas_tpf_upper,
        PER_GEN_TIME_MERGE, "dual_UB_gfm_V_meas_tf.csv", :dual_UB_gfm_V_meas_tf_xlsx),
    DualRegistryEntry(:dual_LB_gfm_E_int_raw_tf, :ineq_const,
        :ineq_const_gfm_E_int_raw_tf_lower, :ineq_const_gfm_E_int_raw_tpf_lower,
        PER_GEN_TIME_MERGE, "dual_LB_gfm_E_int_raw_tf.csv", :dual_LB_gfm_E_int_raw_tf_xlsx),
    DualRegistryEntry(:dual_UB_gfm_E_int_raw_tf, :ineq_const,
        :ineq_const_gfm_E_int_raw_tf_upper, :ineq_const_gfm_E_int_raw_tpf_upper,
        PER_GEN_TIME_MERGE, "dual_UB_gfm_E_int_raw_tf.csv", :dual_UB_gfm_E_int_raw_tf_xlsx),
    DualRegistryEntry(:dual_LB_gfm_E_droop_raw_tf, :ineq_const,
        :ineq_const_gfm_E_droop_raw_tf_lower, :ineq_const_gfm_E_droop_raw_tpf_lower,
        PER_GEN_TIME_MERGE, "dual_LB_gfm_E_droop_raw_tf.csv", :dual_LB_gfm_E_droop_raw_tf_xlsx),
    DualRegistryEntry(:dual_UB_gfm_E_droop_raw_tf, :ineq_const,
        :ineq_const_gfm_E_droop_raw_tf_upper, :ineq_const_gfm_E_droop_raw_tpf_upper,
        PER_GEN_TIME_MERGE, "dual_UB_gfm_E_droop_raw_tf.csv", :dual_UB_gfm_E_droop_raw_tf_xlsx),
    DualRegistryEntry(:dual_LB_gfm_Id_tf, :ineq_const,
        :ineq_const_gfm_Id_tf_lower, :ineq_const_gfm_Id_tpf_lower,
        PER_GEN_TIME_MERGE, "dual_LB_gfm_Id_tf.csv", :dual_LB_gfm_Id_tf_xlsx),
    DualRegistryEntry(:dual_UB_gfm_Id_tf, :ineq_const,
        :ineq_const_gfm_Id_tf_upper, :ineq_const_gfm_Id_tpf_upper,
        PER_GEN_TIME_MERGE, "dual_UB_gfm_Id_tf.csv", :dual_UB_gfm_Id_tf_xlsx),
    DualRegistryEntry(:dual_LB_gfm_Iq_tf, :ineq_const,
        :ineq_const_gfm_Iq_tf_lower, :ineq_const_gfm_Iq_tpf_lower,
        PER_GEN_TIME_MERGE, "dual_LB_gfm_Iq_tf.csv", :dual_LB_gfm_Iq_tf_xlsx),
    DualRegistryEntry(:dual_UB_gfm_Iq_tf, :ineq_const,
        :ineq_const_gfm_Iq_tf_upper, :ineq_const_gfm_Iq_tpf_upper,
        PER_GEN_TIME_MERGE, "dual_UB_gfm_Iq_tf.csv", :dual_UB_gfm_Iq_tf_xlsx),
    DualRegistryEntry(:dual_LB_gfm_E_int_tf, :ineq_const,
        :ineq_const_gfm_E_int_tf_lower, :ineq_const_gfm_E_int_tpf_lower,
        PER_GEN_TIME_MERGE, "dual_LB_gfm_E_int_tf.csv", :dual_LB_gfm_E_int_tf_xlsx),
    DualRegistryEntry(:dual_UB_gfm_E_int_tf, :ineq_const,
        :ineq_const_gfm_E_int_tf_upper, :ineq_const_gfm_E_int_tpf_upper,
        PER_GEN_TIME_MERGE, "dual_UB_gfm_E_int_tf.csv", :dual_UB_gfm_E_int_tf_xlsx),
    DualRegistryEntry(:dual_LB_gfm_E_droop_tf, :ineq_const,
        :ineq_const_gfm_E_droop_tf_lower, :ineq_const_gfm_E_droop_tpf_lower,
        PER_GEN_TIME_MERGE, "dual_LB_gfm_E_droop_tf.csv", :dual_LB_gfm_E_droop_tf_xlsx),
    DualRegistryEntry(:dual_UB_gfm_E_droop_tf, :ineq_const,
        :ineq_const_gfm_E_droop_tf_upper, :ineq_const_gfm_E_droop_tpf_upper,
        PER_GEN_TIME_MERGE, "dual_UB_gfm_E_droop_tf.csv", :dual_UB_gfm_E_droop_tf_xlsx),
    # Overload-PI guard-rails. Same presence rule as the equality rows above.
    DualRegistryEntry(:dual_LB_gfm_Pmax_raw_tf, :ineq_const,
        :ineq_const_gfm_Pmax_raw_tf_lower, :ineq_const_gfm_Pmax_raw_tpf_lower,
        PER_GEN_TIME_MERGE, "dual_LB_gfm_Pmax_raw_tf.csv", :dual_LB_gfm_Pmax_raw_tf_xlsx),
    DualRegistryEntry(:dual_UB_gfm_Pmax_raw_tf, :ineq_const,
        :ineq_const_gfm_Pmax_raw_tf_upper, :ineq_const_gfm_Pmax_raw_tpf_upper,
        PER_GEN_TIME_MERGE, "dual_UB_gfm_Pmax_raw_tf.csv", :dual_UB_gfm_Pmax_raw_tf_xlsx),
    DualRegistryEntry(:dual_LB_gfm_Pmax_clip_tf, :ineq_const,
        :ineq_const_gfm_Pmax_clip_tf_lower, :ineq_const_gfm_Pmax_clip_tpf_lower,
        PER_GEN_TIME_MERGE, "dual_LB_gfm_Pmax_clip_tf.csv", :dual_LB_gfm_Pmax_clip_tf_xlsx),
    DualRegistryEntry(:dual_UB_gfm_Pmax_clip_tf, :ineq_const,
        :ineq_const_gfm_Pmax_clip_tf_upper, :ineq_const_gfm_Pmax_clip_tpf_upper,
        PER_GEN_TIME_MERGE, "dual_UB_gfm_Pmax_clip_tf.csv", :dual_UB_gfm_Pmax_clip_tf_xlsx),
    DualRegistryEntry(:dual_LB_gfm_Pmin_raw_tf, :ineq_const,
        :ineq_const_gfm_Pmin_raw_tf_lower, :ineq_const_gfm_Pmin_raw_tpf_lower,
        PER_GEN_TIME_MERGE, "dual_LB_gfm_Pmin_raw_tf.csv", :dual_LB_gfm_Pmin_raw_tf_xlsx),
    DualRegistryEntry(:dual_UB_gfm_Pmin_raw_tf, :ineq_const,
        :ineq_const_gfm_Pmin_raw_tf_upper, :ineq_const_gfm_Pmin_raw_tpf_upper,
        PER_GEN_TIME_MERGE, "dual_UB_gfm_Pmin_raw_tf.csv", :dual_UB_gfm_Pmin_raw_tf_xlsx),
    DualRegistryEntry(:dual_LB_gfm_Pmin_clip_tf, :ineq_const,
        :ineq_const_gfm_Pmin_clip_tf_lower, :ineq_const_gfm_Pmin_clip_tpf_lower,
        PER_GEN_TIME_MERGE, "dual_LB_gfm_Pmin_clip_tf.csv", :dual_LB_gfm_Pmin_clip_tf_xlsx),
    DualRegistryEntry(:dual_UB_gfm_Pmin_clip_tf, :ineq_const,
        :ineq_const_gfm_Pmin_clip_tf_upper, :ineq_const_gfm_Pmin_clip_tpf_upper,
        PER_GEN_TIME_MERGE, "dual_UB_gfm_Pmin_clip_tf.csv", :dual_UB_gfm_Pmin_clip_tf_xlsx),
    DualRegistryEntry(:dual_LB_gfm_Qmax_raw_tf, :ineq_const,
        :ineq_const_gfm_Qmax_raw_tf_lower, :ineq_const_gfm_Qmax_raw_tpf_lower,
        PER_GEN_TIME_MERGE, "dual_LB_gfm_Qmax_raw_tf.csv", :dual_LB_gfm_Qmax_raw_tf_xlsx),
    DualRegistryEntry(:dual_UB_gfm_Qmax_raw_tf, :ineq_const,
        :ineq_const_gfm_Qmax_raw_tf_upper, :ineq_const_gfm_Qmax_raw_tpf_upper,
        PER_GEN_TIME_MERGE, "dual_UB_gfm_Qmax_raw_tf.csv", :dual_UB_gfm_Qmax_raw_tf_xlsx),
    DualRegistryEntry(:dual_LB_gfm_Qmax_clip_tf, :ineq_const,
        :ineq_const_gfm_Qmax_clip_tf_lower, :ineq_const_gfm_Qmax_clip_tpf_lower,
        PER_GEN_TIME_MERGE, "dual_LB_gfm_Qmax_clip_tf.csv", :dual_LB_gfm_Qmax_clip_tf_xlsx),
    DualRegistryEntry(:dual_UB_gfm_Qmax_clip_tf, :ineq_const,
        :ineq_const_gfm_Qmax_clip_tf_upper, :ineq_const_gfm_Qmax_clip_tpf_upper,
        PER_GEN_TIME_MERGE, "dual_UB_gfm_Qmax_clip_tf.csv", :dual_UB_gfm_Qmax_clip_tf_xlsx),
    DualRegistryEntry(:dual_LB_gfm_Qmin_raw_tf, :ineq_const,
        :ineq_const_gfm_Qmin_raw_tf_lower, :ineq_const_gfm_Qmin_raw_tpf_lower,
        PER_GEN_TIME_MERGE, "dual_LB_gfm_Qmin_raw_tf.csv", :dual_LB_gfm_Qmin_raw_tf_xlsx),
    DualRegistryEntry(:dual_UB_gfm_Qmin_raw_tf, :ineq_const,
        :ineq_const_gfm_Qmin_raw_tf_upper, :ineq_const_gfm_Qmin_raw_tpf_upper,
        PER_GEN_TIME_MERGE, "dual_UB_gfm_Qmin_raw_tf.csv", :dual_UB_gfm_Qmin_raw_tf_xlsx),
    DualRegistryEntry(:dual_LB_gfm_Qmin_clip_tf, :ineq_const,
        :ineq_const_gfm_Qmin_clip_tf_lower, :ineq_const_gfm_Qmin_clip_tpf_lower,
        PER_GEN_TIME_MERGE, "dual_LB_gfm_Qmin_clip_tf.csv", :dual_LB_gfm_Qmin_clip_tf_xlsx),
    DualRegistryEntry(:dual_UB_gfm_Qmin_clip_tf, :ineq_const,
        :ineq_const_gfm_Qmin_clip_tf_upper, :ineq_const_gfm_Qmin_clip_tpf_upper,
        PER_GEN_TIME_MERGE, "dual_UB_gfm_Qmin_clip_tf.csv", :dual_UB_gfm_Qmin_clip_tf_xlsx),
]

"""Return the dual-export catalog for the assembled dynamic model."""
function dual_spec_catalog(dyn_model_dict::OrderedDict{Symbol, Any})::Vector{DualRegistryEntry}
    meta = get(dyn_model_dict, :meta, OrderedDict{Symbol, Any}())
    if get(meta, :gen_order, "") == "DQ_4TH"
        return DQ_4TH_FULL_BUS_DUAL_SPECS
    elseif get(meta, :network_form, "KRON_REDUCED") == "FULL_BUS"
        return CLASSICAL_FULL_BUS_DUAL_SPECS
    else
        return CLASSICAL_KRON_DUAL_SPECS
    end
end

# ==================================================================================
# Internal helpers — map dyn_model_dict containers to dual value arrays
# ==================================================================================

"""Return the `:eq_const` or `:ineq_const` OrderedDict inside `dyn_model_dict`."""
function _constraint_store(dyn_model_dict::OrderedDict{Symbol, Any}, source::Symbol)
    source == :eq_const && return dyn_model_dict[:eq_const]
    source == :ineq_const && return dyn_model_dict[:ineq_const]
    throw(ArgumentError("Unknown dual registry source: $source"))
end

"""
    registry_entry_present(dyn_model_dict, entry) -> Bool

`true` when the fault-on constraint key required by `entry` exists, either as a
constraint container or through the bound manifest (VARIABLE encoding).

Presence is decided by `key_tf` **alone**, on purpose.  Single-window
disturbances (GL generator/load trip, OB open branch) build only the fault
window — `time_windows_gld` returns the same window twice and the builders call
`Define_Fault_*` without a post-fault stage — so no `*_tpf` key is ever created.
Requiring `key_tpf` here silently dropped every stability dual on those runs.
`extract_dual_entry` already emits tf-only data when the tpf key is absent, so
the gate must not be stricter than the extractor it guards.  The same applies to
the four independent δ-COI / Δω-COI toggles, which may enable the fault-on box
without its post-fault twin.
"""
function registry_entry_present(
    dyn_model_dict::OrderedDict{Symbol, Any},
    entry::DualRegistryEntry,
)::Bool
    store = _constraint_store(dyn_model_dict, entry.source)
    haskey(store, entry.key_tf) && return true
    return bound_manifest_entry_present(dyn_model_dict, entry.key_tf)
end

# One dual per generator (pre-fault init constraints, variable bound duals).
function _duals_gen_indexed(container)::Vector{Float64}
    return [JuMP.dual(c) for (_, c) in container]
end

# One dual per time step (COI equality constraints).
function _duals_time_indexed(container)::Vector{Float64}
    return [JuMP.dual(c) for (_, c) in container]
end

# Time series per generator (swing, Pe, stability-bound duals).
# The key type is left to the container: every family but one is keyed by a generator or bus
# id, while the pairwise angle corridor is keyed by a `(i, j)` tuple. Nothing here reads the
# key, so keeping the signature generic is what lets that family reuse PER_GEN_TIME_MERGE
# instead of needing a layout tag of its own.
function _duals_per_gen_time(container)
    return OrderedDict(
        i => [JuMP.dual(v) for (_, v) in inner] for (i, inner) in container)
end

"""Concatenate fault-on and post-fault per-gen trajectories along the time axis."""
function _merge_per_gen_time(tf::AbstractDict, tpf::AbstractDict)
    out = OrderedDict{keytype(tf), Vector{Float64}}()
    for k in keys(tf)
        out[k] = vcat(tf[k], tpf[k])
    end
    return out
end

"""
    extract_dual_entry(dyn_model_dict, entry)

Read `JuMP.dual` for every constraint in the family described by `entry`.
Returns `nothing` if the keys are absent (optional constraint family).
"""
function extract_dual_entry(
    dyn_model_dict::OrderedDict{Symbol, Any},
    entry::DualRegistryEntry,
)
    store = _constraint_store(dyn_model_dict, entry.source)
    use_manifest_tf = !haskey(store, entry.key_tf) &&
        bound_manifest_entry_present(dyn_model_dict, entry.key_tf)

    if use_manifest_tf
        tf_data = extract_bound_manifest_duals(dyn_model_dict, entry.key_tf)
        if entry.layout == GEN_INDEXED
            _, tf_vec = tf_data
            return tf_vec
        elseif entry.layout == TIME_INDEXED
            _, tf_vec = tf_data
            return tf_vec
        elseif entry.layout == TIME_INDEXED_MERGE
            _, tf_vec = tf_data
            if entry.key_tpf !== nothing && bound_manifest_entry_present(dyn_model_dict, entry.key_tpf)
                _, tpf_vec = extract_bound_manifest_duals(dyn_model_dict, entry.key_tpf)
                return vcat(tf_vec, tpf_vec)
            end
            return tf_vec
        elseif entry.layout == PER_GEN_TIME
            return tf_data
        elseif entry.layout == PER_GEN_TIME_MERGE
            tf_dict = tf_data
            if entry.key_tpf !== nothing && bound_manifest_entry_present(dyn_model_dict, entry.key_tpf)
                tpf_dict = extract_bound_manifest_duals(dyn_model_dict, entry.key_tpf)
                return _merge_per_gen_time(tf_dict, tpf_dict)
            end
            return tf_dict
        elseif entry.layout == PER_BUS_TIME_MERGE
            tf_dict = tf_data
            if entry.key_tpf !== nothing && bound_manifest_entry_present(dyn_model_dict, entry.key_tpf)
                tpf_dict = extract_bound_manifest_duals(dyn_model_dict, entry.key_tpf)
                return _merge_per_gen_time(tf_dict, tpf_dict)
            end
            return tf_dict
        end
    end

    registry_entry_present(dyn_model_dict, entry) || return nothing
    store = _constraint_store(dyn_model_dict, entry.source)
    tf_data = store[entry.key_tf]

    if entry.layout == GEN_INDEXED
        return _duals_gen_indexed(tf_data)
    elseif entry.layout == TIME_INDEXED
        return _duals_time_indexed(tf_data)
    elseif entry.layout == TIME_INDEXED_MERGE
        tf_vec = _duals_time_indexed(tf_data)
        if entry.key_tpf !== nothing && haskey(store, entry.key_tpf)
            return vcat(tf_vec, _duals_time_indexed(store[entry.key_tpf]))
        end
        return tf_vec
    elseif entry.layout == PER_GEN_TIME
        return _duals_per_gen_time(tf_data)
    elseif entry.layout == PER_GEN_TIME_MERGE
        tf_dict = _duals_per_gen_time(tf_data)
        if entry.key_tpf !== nothing && haskey(store, entry.key_tpf)
            return _merge_per_gen_time(tf_dict, _duals_per_gen_time(store[entry.key_tpf]))
        end
        return tf_dict
    elseif entry.layout == PER_BUS_TIME_MERGE
        tf_dict = _duals_per_gen_time(tf_data)  # same OrderedDict{bus, Vector} shape
        if entry.key_tpf !== nothing && haskey(store, entry.key_tpf)
            return _merge_per_gen_time(tf_dict, _duals_per_gen_time(store[entry.key_tpf]))
        end
        return tf_dict
    else
        throw(ArgumentError("Unhandled DualExtractLayout: $(entry.layout)"))
    end
end

"""
    extract_dual_entry_ids(dyn_model_dict, entry) -> Union{Nothing, Vector{Int}}

Generator (or bus) ids behind a `GEN_INDEXED` family, in the same order as the
values returned by [`extract_dual_entry`](@ref).

`extract_dual_entry` flattens those containers to a bare vector, which is what the
CSV writer used to dump — a file with no id column at all. That is safe to read only
if every dual family covers exactly the same units in the same order, and it does
not: on a mixed fleet `dual_Pe_init.csv` spans SG + GFM while `dual_Pe.csv` spans the
machines alone, so aligning them by row silently attributes one unit's dual to
another. Returns `nothing` when the family is absent.
"""
function extract_dual_entry_ids(
    dyn_model_dict::OrderedDict{Symbol, Any},
    entry::DualRegistryEntry,
)::Union{Nothing, Vector{Int}}
    entry.layout == GEN_INDEXED || return nothing
    store = _constraint_store(dyn_model_dict, entry.source)
    haskey(store, entry.key_tf) && return Int[Int(k) for k in keys(store[entry.key_tf])]
    if bound_manifest_entry_present(dyn_model_dict, entry.key_tf)
        ids, _ = extract_bound_manifest_duals(dyn_model_dict, entry.key_tf)
        return Int[Int(i) for i in ids]
    end
    return nothing
end

# ==================================================================================
# Public API — called from builders (register) and save layer (extract)
# ==================================================================================

"""
    build_dual_registry!(dyn_model_dict)

Scan `CLASSICAL_KRON_DUAL_SPECS` and store the applicable entries in
`dyn_model_dict[:dual_registry]`.  Call once at the end of
`Make_Dynamic_Model_tsred!` / `Make_Dynamic_Model_tsredlinear!`.
"""
function build_dual_registry!(dyn_model_dict::OrderedDict{Symbol, Any})
    catalog = dual_spec_catalog(dyn_model_dict)
    specs = DualRegistryEntry[e for e in catalog if registry_entry_present(dyn_model_dict, e)]
    dyn_model_dict[:dual_registry] = specs
    return dyn_model_dict
end

"""
    extract_registered_duals(dyn_model_dict) -> Dict{Symbol, Any}

Collect all dual trajectories keyed by `export_name` (e.g. `:dual_Pe`).
Rebuilds the registry on the fly if the builder did not call `build_dual_registry!`.
"""
function extract_registered_duals(dyn_model_dict::OrderedDict{Symbol, Any})
    registry = get(dyn_model_dict, :dual_registry, DualRegistryEntry[])
    if isempty(registry)
        build_dual_registry!(dyn_model_dict)
        registry = dyn_model_dict[:dual_registry]
    end
    duals = Dict{Symbol, Any}()
    for entry in registry
        val = extract_dual_entry(dyn_model_dict, entry)
        val !== nothing && (duals[entry.export_name] = val)
    end
    return duals
end

"""
    duals_to_xlsx_kwargs(duals, dyn_model_dict) -> Dict{Symbol, Any}

Map extracted duals to the keyword names expected by `Save_Duals_2_Excel_tsred`.
Uses the full-bus or Kron spec catalog according to `dyn_model_dict[:meta][:network_form]`.
"""
function duals_to_xlsx_kwargs(
    duals::Dict{Symbol, Any},
    dyn_model_dict::OrderedDict{Symbol, Any},
)
    catalog = dual_spec_catalog(dyn_model_dict)
    kwargs = Dict{Symbol, Any}()
    for entry in catalog
        haskey(duals, entry.export_name) || continue
        kwargs[entry.xlsx_kw] = duals[entry.export_name]
    end
    return kwargs
end
