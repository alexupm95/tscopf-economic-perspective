#=
================================================================================
 TSCOPF.jl  —  Transient-Stability-Constrained Optimal Power Flow (package root)
================================================================================
 Stage A (Phase 6): module wrapper around the former `setup_includes.jl` tree.
 Include order matches the pre-package layout; paths are relative to `src/`.
================================================================================
=#

module TSCOPF

using Logging
using PiecewiseLinearOpt
using Dualization

using LinearAlgebra, SparseArrays
using Dates, NumericIO, DataFrames, Printf, CSV, DataStructures, XLSX, TOML
using JuMP, Ipopt, HiGHS
import MathOptInterface as MOI
using Trapz

# --- source includes (dependency order) ---------------------------------------

include("_common/constants.jl")
include("_common/bound_encoding.jl")
include("_manage_inputs/SolverConfigCommon.jl")
include("_manage_inputs/IpoptSolverConfig.jl")
include("_common/ipopt_hsl.jl")
include("_manage_inputs/HiGHSSolverConfig.jl")
include("_manage_inputs/GurobiSolverConfig.jl")
include("_manage_inputs/MadNLPSolverConfig.jl")
include("_common/madnlp_linear_solvers.jl")
include("solver_registry.jl")
include("plots_dispatch.jl")
include("_common/DispatchConfig.jl")

include("_acopf/functions_make_acopf_model.jl")
include("_acopf/functions_4_AC_power_flow.jl")

include("_common/function_build_dispatch_model.jl")
include("_common/functions_4_sanity_checks.jl")
include("_common/auxiliar_functions.jl")
include("_common/functions_4_derivatives.jl")
include("_common/function_get_angle_limits.jl")
include("_common/functions_4_eqconst.jl")
include("_common/functions_4_ineqconst.jl")
include("_common/functions_4_obj_func.jl")
include("_common/bound_builders.jl")
include("_common/functions_4_variables.jl")
include("_common/functions_4_admittance_matrices.jl")
include("_common/functions_2_manipulate_optimal_matrices.jl")

include("_dcopf/functions_make_dcopf_model.jl")
include("_dcopf/functions_4_DC_power_flow.jl")
include("_dcopf/functions_build_DC_OPF_Dual.jl")

include("_ed/functions_make_ed_model.jl")
include("_ed/functions_build_ED_Dual.jl")

include("_uc/uc_toy_data.jl")
include("_uc/uc_variables.jl")
include("_uc/uc_eq_constraints.jl")
include("_uc/uc_ineq_constraints.jl")
include("_uc/functions_make_uc_model.jl")

include("_manage_outputs/functions_2_save_TS_results.jl")
include("_transient_stability/numerics.jl")   # numerical constants; before every TS builder
include("_transient_stability/FaultConfig.jl")
include("_transient_stability/DynModelConfig.jl")
include("_transient_stability/coupling_init.jl")
include("_transient_stability/ts_bound_limits.jl")
include("_transient_stability/TsConfig.jl")
include("_transient_stability/AbstractDynamicGenModel.jl")
include("_transient_stability/DynDualRegistry.jl")
include("_transient_stability/functions_2_build_TS_model_common.jl")
include("_transient_stability/functions_2_build_TS_model_w_Kron_shared.jl")
include("_transient_stability/functions_4_TS_kron_variables.jl")
include("_transient_stability/functions_4_TS_kron_eqconst.jl")
include("_transient_stability/functions_4_TS_kron_ineqconst.jl")
include("_transient_stability/functions_4_TS_fullbus_helpers.jl")
include("_transient_stability/functions_4_TS_fullbus_variables.jl")
include("_transient_stability/functions_4_TS_fullbus_eqconst.jl")
include("_transient_stability/functions_4_TS_fullbus_ineqconst.jl")
include("_transient_stability/functions_4_TS_governor.jl")
include("_transient_stability/functions_4_TS_avr.jl")
# DQ 4th-order machine (FULL_BUS only; milestone 1 = core, no AVR/governor)
include("_transient_stability/functions_4_TS_dq_helpers.jl")
include("_transient_stability/functions_4_TS_dq_variables.jl")
include("_transient_stability/functions_4_TS_dq_eqconst.jl")
include("_transient_stability/functions_4_TS_gfm.jl")
include("_transient_stability/functions_2_build_TS_model_w_FullBus.jl")
include("_transient_stability/functions_2_build_TS_model_w_DqFullBus.jl")
include("_transient_stability/functions_2_build_TS_model_w_Kron.jl")
include("_transient_stability/functions_2_build_TS_model_w_Kron_Linear.jl")

include("_manage_inputs/functions_2_read_input_files.jl")
include("_manage_inputs/functions_2_parse_matpower.jl")
include("_manage_inputs/functions_2_setup_optim.jl")

include("_manage_outputs/functions_2_save_input_parameters.jl")
include("_manage_outputs/function_2_save_dispatch_common.jl")
include("_manage_outputs/functions_2_save_dispatch_model.jl")
include("_manage_outputs/functions_2_save_dispatch_results.jl")
include("_manage_outputs/DispatchDualRegistry.jl")
include("_manage_outputs/functions_2_save_dispatch_duals.jl")
# Needs both dual registries in scope (the [exports] table lists what they wrote).
include("_manage_outputs/functions_2_save_run_manifest.jl")

include("engine.jl")

# --- public API ---------------------------------------------------------------

export RunConfig, reconfigure, reconfigure_transient, SystemData, load_system, run_case!
export validate_run_config!, validate_dyn_config!, validate_δ_reference!,
       validate_gfm_dispatch_band!
export validate_all_gfm_fleet!
export DispatchConfig, DispatchLimitsConfig, BoundEncoding, TransientConfig, TsSimulationConfig, TsBuilderConfig
export TsBoundLimitsConfig, TsBoundLimitPair
export validate_dispatch_config!
export DynModelConfig, FaultConfig, SteadyStateHints
export apply_steady_state_hints_to_opf!
export GenOrder, NetworkForm, MechPowerMode, FaultType, SusceptanceModel, GovernorLimiter
export SC, GL, OB, CLASSICAL_2ND, DQ_4TH, KRON_REDUCED, FULL_BUS, USE_PG, USE_PM, SIMPLE, POWERMODELS
export GOV_NO_LIMIT, GOV_SMOOTH, GOV_HARD_BOUND
export Parse_GFM_Dynamic_DataFrame, Read_GFM_Dynamic_Data, apply_gfm_base_conversion!
export Parse_Load_DataFrame, Read_Load_Data, remap_load_buses!, validate_load_table!, zip_by_bus
export validate_gen_dyn_partition!, validate_dgfm_ids!, validate_dgen_dyn_subset_ids!
export gfm_id_set, sg_active_gens, gfm_active_gens, dgfm_row, gfm_machine_warmstart
export Attach_GFM_init!, Attach_GFM_fault!, Attach_GFM_postf!, register_gfm_meta!
export attach_gfm_acopf_limits!, resolve_sg_gfm_gens!
export Save_Prefault_Coupling_Starts!
export Save_Run_Manifest!
export build_constraint_evaluator
export Export_Variable_Bounds!
export resolve_gfm_bound_limits, smooth_max_expr, smooth_min_expr, smooth_clip_expr
export validate_fault_config!, build_fault_details, copy_fault_config, apply_gl_load_scaling!
export Setup_Optim_Model, set_solver_log_path!, release_solver_backend!, Clean_Terminal
export IpoptSolverConfig, apply_ipopt_options!, validate_ipopt_solver_config!
export HiGHSSolverConfig, apply_highs_options!, validate_highs_solver_config!
export GurobiSolverConfig, apply_gurobi_options!, validate_gurobi_solver_config!
export MadNLPSolverConfig, apply_madnlp_options!, validate_madnlp_solver_config!
export MADNLP_LOG_INFO, MADNLP_LOG_ERROR, madnlp_hessian_type_name
export FAULT_BUS_SHUNT
export hsl_jll_available, is_ipopt_backend_solver, IPOPT_HSL_LINEAR_SOLVERS
export IPOPT_PARDISO_SOLVER, pardiso_available, configure_ipopt_linear_solver!
export is_madnlp_backend_solver, MADNLP_BACKEND_SOLVERS, MADNLP_HSL_LINEAR_SOLVERS
export madnlp_hsl_ext_available, configure_madnlp_linear_solver!
export Check_Coherence_Input_Data, resolve_save_optim_matrices, gurobi_available, plots_extension_loaded
export load_plots_extension!
export project_root, default_results_dir, input_files_dir
export build_results_path_names, build_results_paths, results_folder_keys
export Import_Matpower_Case, build_opf_input_param
export Make_UC_Model!, Solve_UC_Restricted_Pricing!, uc_toy_system_data, build_uc_toy_system
export uc_commitment_values
export Save_Solution_UC_Model

function __init__()
    _register_core_solver_builders!()
end

end # module TSCOPF
