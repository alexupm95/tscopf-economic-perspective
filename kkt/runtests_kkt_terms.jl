#=
 Residual assertion for the KKT stationarity export.

   julia --project=TSCOPF kkt/runtests_kkt_terms.jl

 Deliberately NOT registered in test/runtests.jl: that file ships in the public
 snapshot, where this path does not exist.

 The load-bearing check is `|RESIDUAL| < TOL` for every generator. It pins the sign
 convention of every family at once and fails the moment a builder adds a `P_g`
 coupling the family table does not know about — which is exactly the failure mode
 the analytic coefficients cannot detect on their own.
=#

using Pkg
Pkg.activate(joinpath(@__DIR__, "..", ".."))

using Test
using TSCOPF
using CSV, DataFrames
using TOML

include(joinpath(@__DIR__, "kkt_pg_terms.jl"))

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const SOLVED = (OPTIMAL, LOCALLY_SOLVED, ALMOST_LOCALLY_SOLVED)
# Loose next to the Ipopt tolerance (1e-8) but tight enough that a dropped family
# — the smallest of which is worth O(1) €/MW — cannot hide under it.
const RESIDUAL_TOL = 1e-6

"""Run one case with the KKT hook attached; return `(result, kkt_table, manifest)`."""
function run_with_kkt(cfg::RunConfig)
    hook = (; cfg, sys, path_names, opf_dict, dyn_model_dict, status_solved, kwargs...) -> begin
        (cfg.save_duals && status_solved) &&
            Save_KKT_Pg_Terms!(path_names, cfg, sys, opf_dict, dyn_model_dict)
        nothing
    end
    cfg = reconfigure(cfg; post_solve_hook = hook)

    sys = load_system(cfg, PROJECT_ROOT)
    result = run_case!(cfg, sys, PROJECT_ROOT, default_results_dir())

    run_dir = result.path_names[:pf_results_date]
    manifest = TOML.parsefile(joinpath(run_dir, "run_manifest.toml"))
    kkt_path = joinpath(run_dir, "Dispatch", "CSV_duals", "kkt_Pg_terms.csv")
    table = isfile(kkt_path) ? CSV.read(kkt_path, DataFrame; delim = ';') : nothing
    return result, table, manifest
end

"""Every assertion that does not depend on the model path."""
function check_common(result, table, manifest, cfg)
    @test result.status in SOLVED
    @test manifest["run"]["load_factor"] == cfg.load_factor
    @test manifest["run"]["case"] == cfg.case
    @test manifest["transient"]["dyn_model"]["gen_order"] ==
          string(cfg.transient.dyn_model.gen_order)
    @test manifest["transient"]["dyn_model"]["network_form"] ==
          string(cfg.transient.dyn_model.network_form)
    @test haskey(manifest["resolved"], "t_step")
    @test !isempty(manifest["exports"]["ts_dual_files"])

    @test table !== nothing
    @test names(table) == ["gen_id", "bus_id", "family", "t_index",
                           "dual", "coeff", "contribution"]

    residuals = table[table.family .== "RESIDUAL", :]
    @test nrow(residuals) == length(unique(table.gen_id))
    for row in eachrow(residuals)
        @test abs(row.contribution) < RESIDUAL_TOL
    end

    # The decomposition is only meaningful if the families it found are the ones the
    # run actually built: a silent `continue` everywhere would also give residual 0.
    families = Set(table.family)
    @test "p_balance" in families
    @test "cost_gradient" in families
    @test "Pe_init" in families

    return residuals
end

kron_cfg = RunConfig(
    trans_stab = true,
    case = "9bus",
    load_factor = 1.5,
    solver_name = "Ipopt",
    silent_solver = true,
    save_duals = true,
    save_matrices = false,
    save_ts_plots = false,
    dispatch = DispatchConfig(type_model = "ACOPF", use_matrix = true,
                              ineq_sbranch_upper = true),
    transient = TransientConfig(
        simulation = TsSimulationConfig(δ_tol_deg = 100.0, t_end_sim = 0.6,
                                        t_step = 0.02, t_start_fault = 0.01,
                                        clearing_time = 0.15),
        builder = TsBuilderConfig(
            ineq_δ_COI_tf_lower = true, ineq_δ_COI_tf_upper = true,
            ineq_δ_COI_tpf_lower = true, ineq_δ_COI_tpf_upper = true),
        dyn_model = DynModelConfig(
            gen_order = CLASSICAL_2ND,
            network_form = KRON_REDUCED,
            mech_power_mode = USE_PM,      # exercises the Pm_init family
            bound_style_δ = :coi_box,
            constrain_Δω = false,
            fault = FaultConfig(fault_type = SC, contingency_id = 2),
        ),
    ),
)

gfm_cfg = RunConfig(
    trans_stab = true,
    case = "9bus_gfm",
    load_factor = 1.5,
    solver_name = "Ipopt",
    silent_solver = true,
    save_duals = true,
    save_matrices = false,
    save_ts_plots = false,
    dispatch = DispatchConfig(type_model = "ACOPF", use_matrix = true),
    transient = TransientConfig(
        gen_dynamic_filename = "gen_dynamic_data_full.csv",
        gfm_dynamic_filename = "gfm_dynamic_data.csv",
        simulation = TsSimulationConfig(δ_tol_deg = 100.0, t_end_sim = 0.6,
                                        t_step = 0.02, t_start_fault = 0.01,
                                        clearing_time = 0.15),
        builder = TsBuilderConfig(),
        dyn_model = DynModelConfig(
            allow_gfm = true,
            gen_order = DQ_4TH,
            network_form = FULL_BUS,
            mech_power_mode = USE_PM,
            bound_style_δ = :coi_box,
            zip_load_p = (1.0, 0.0, 0.0),
            zip_load_q = (1.0, 0.0, 0.0),
            dq_speed_dev_in_algebra = true,
            constrain_Δω = false,
            fault = FaultConfig(fault_type = SC, contingency_id = 2),
        ),
    ),
)

@testset "KKT P_g stationarity export" begin

    @testset "Kron classical (USE_PM)" begin
        result, table, manifest = run_with_kkt(kron_cfg)
        residuals = check_common(result, table, manifest, kron_cfg)
        @test "Pm_init" in Set(table.family)
        println("Kron residuals: ", round.(residuals.contribution; sigdigits = 3))
    end

    @testset "DQ FULL_BUS + GFM" begin
        result, table, manifest = run_with_kkt(gfm_cfg)
        residuals = check_common(result, table, manifest, gfm_cfg)
        # The mixed fleet is the only path that exercises the converter families.
        @test !isempty(manifest["resolved"]["gfm_ids"])
        @test "gfm_Pmeas_init" in Set(table.family)
        # FULL_BUS + :trapezoidal ⇒ the first swing step carries P_g.
        @test "swing_first_step" in Set(table.family)
        println("GFM residuals: ", round.(residuals.contribution; sigdigits = 3))
    end

end
