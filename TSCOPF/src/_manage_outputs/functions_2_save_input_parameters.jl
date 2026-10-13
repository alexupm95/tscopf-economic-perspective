"""
    Copy_Input_CSVs_To_Results!(path_names; trans_stab, gen_dynamic_filename,
                                 matpower_file=nothing)

Copy the input files used for this run into `path_names[:pf_inputs]` (the
`Inputs/` folder under the timestamped results directory) so results stay
paired with the exact network / contingency data that was read.

Two modes, mirroring `load_system`:
  CSV mode     (matpower_file=nothing): archives bus_data.csv, generators_data.csv,
               line_data.csv from the case folder.
  MATPOWER mode (matpower_file set):   archives the .m file instead of those CSVs.
In both modes, gen_dynamic_data.csv and contingencies.csv are archived when
trans_stab=true, and load_data.csv whenever the case defines one (they always live in
the case folder regardless of input mode).
"""
function Copy_Input_CSVs_To_Results!(
    path_names::OrderedDict{Symbol, String};
    trans_stab::Bool,
    gen_dynamic_filename::String,
    matpower_file::Union{Nothing, String} = nothing,
    gfm_dynamic_filename::Union{Nothing, String} = nothing,
)
    # Source: case folder under INPUT_FILES; destination: Results/.../Inputs/
    src_dir = path_names[:pf_input_files]
    dst_dir = path_names[:pf_inputs]
    mkpath(dst_dir)

    copied = String[]

    if isnothing(matpower_file)
        # ── CSV mode: archive the three steady-state CSVs ─────────────────────
        for fname in String["bus_data.csv", "generators_data.csv", "line_data.csv"]
            src = joinpath(src_dir, fname)
            isfile(src) || throw(ArgumentError("Input CSV not found for archival copy: $src"))
            cp(src, joinpath(dst_dir, fname); force=true)
            push!(copied, fname)
        end
    else
        # ── MATPOWER mode: archive the .m file instead ────────────────────────
        src_m = isabspath(matpower_file) ? matpower_file : joinpath(src_dir, matpower_file)
        isfile(src_m) || throw(ArgumentError(
            "MATPOWER file not found for archival copy: $src_m"))
        dst_name = basename(matpower_file)
        cp(src_m, joinpath(dst_dir, dst_name); force=true)
        push!(copied, dst_name)
    end

    # Dynamic machine data and contingency table: always archived for TSC runs,
    # regardless of input mode (these always live in the case folder).
    if trans_stab
        for fname in String[gen_dynamic_filename, "contingencies.csv"]
            src = joinpath(src_dir, fname)
            isfile(src) || throw(ArgumentError("Input CSV not found for archival copy: $src"))
            cp(src, joinpath(dst_dir, fname); force=true)
            push!(copied, fname)
        end
        if gfm_dynamic_filename !== nothing
            src = joinpath(src_dir, gfm_dynamic_filename)
            if isfile(src)
                cp(src, joinpath(dst_dir, gfm_dynamic_filename); force=true)
                push!(copied, gfm_dynamic_filename)
            end
        end
    end

    # Optional per-load demand split. Archived in both input modes (it lives in the case
    # folder either way) and unscaled, like bus_data.csv — the KKT post-processing reads
    # it from here and re-applies the manifest's load_factor.
    let src = joinpath(src_dir, "load_data.csv")
        if isfile(src)
            cp(src, joinpath(dst_dir, "load_data.csv"); force=true)
            push!(copied, "load_data.csv")
        end
    end

    println("Input files archived to: $dst_dir  ($(join(copied, ", ")))")
    return copied
end

function Print_Input_Parameters(
    path_names::OrderedDict{Symbol, String},
    overwrite_results::Bool,
    trans_stab::Bool,
    type_model::String,
    solver_name::String,
    use_matrix::Bool,
    silent_solver::Bool,
    case::String,
    base_MVA::Float64,
    load_factor::Float64,
    save_optim_matrices::Bool;
    fault::Union{Nothing, FaultConfig}=nothing,
    dyn_parameters_dict::Union{Nothing, OrderedDict{Symbol, Any}}=nothing,
)
    ts_fault_details = nothing
    if trans_stab && !isnothing(dyn_parameters_dict) && fault !== nothing
        ts_fault_details = build_fault_details(fault, path_names[:pf_input_files])
    end

    has_dynamic_parameters = !isnothing(dyn_parameters_dict)

    filename = joinpath(path_names[:pf_results_date], "input_parameters.txt")
    open(filename, "w") do io
        println(io, "******** Simulation Input Parameters ********")
        println(io, "==============================================")
        println(io, "Optimize Transient Stability?: $trans_stab")
        println(io, "Case:                          $case")
        println(io, "Base MVA:                      $base_MVA")
        
        # Accessing common dynamic parameters safely
        if has_dynamic_parameters
            common = dyn_parameters_dict[:common]
            if haskey(common, :δ_tol)
                # Assuming δ_tol is a vector or tuple [lower, upper]
                println(io, "δ tolerance (deg):             [$(rad2deg(common[:δ_tol][1])); +$(rad2deg(common[:δ_tol][2]))]")
            end
            # Only meaningful under bound_style_δ = :pairwise, but printed whenever it is
            # present so the file states every corridor knob the run carried.
            if haskey(common, :δ_pair_tol)
                println(io, "δ pair separation (deg):       [$(rad2deg(common[:δ_pair_tol][1])); +$(rad2deg(common[:δ_pair_tol][2]))]")
            end
            if haskey(common, :f_syn)
                println(io, "Synchronous frequency (Hz):    $(common[:f_syn])")
            end
        end

        # --- Specific Fault Details ---
        if trans_stab && has_dynamic_parameters && ts_fault_details !== nothing
            println(io, "----------------------------------------------")
            println(io, "Reduced Model:                 $(ts_fault_details[:reduced_model])")
            println(io, "Fault Category:                $(ts_fault_details[:fault_type])")

            # 1. Short Circuit Logic (SC)
            if ts_fault_details[:fault_type] == "SC"
                sc = ts_fault_details[:sc]
                println(io, "Fault Location Type:           $(sc[:fault_location])")
                
                if sc[:fault_location] == "bus"
                    println(io, "Fault Bus ID:                  $(sc[:bus][:bus_id])")
                    println(io, "Disconnect Branch?:            $(sc[:bus][:disconnect_branch])")
                    println(io, "Branches to Disconnect:        $(sc[:bus][:branch_id_2_disconnect])")
                
                elseif sc[:fault_location] == "branch"
                    br = sc[:branch]
                    println(io, "Fault Branch ID:               $(br[:branch_id])")
                    println(io, "Between Buses:                 $(br[:adjacent_buses])")
                    println(io, "Location (% from 'from' bus):  $(br[:fault_location_percent])%")
                end

            # 2. Generator/Load Disconnection Logic (GL)
            elseif ts_fault_details[:fault_type] == "GL"
                gl = ts_fault_details[:gl]
                println(io, "Element to Disconnect:         $(gl[:element_2_disconnect])")

                if gl[:element_2_disconnect] == "gen"
                    println(io, "Generator ID(s):               $(gl[:gen][:gen_id])")
                    println(io, "Connected Bus(es):             $(gl[:gen][:bus_id])")
                
                elseif gl[:element_2_disconnect] == "load"
                    ld = gl[:load]
                    println(io, "Load Bus ID(s):                $(ld[:bus_id])")
                    println(io, "Demand scaling:                p_d_new = p_d * (1 + alpha),  q_d_new = q_d * (1 + alpha)")
                    for (bus_id, α) in zip(ld[:bus_id], ld[:percent_power])
                        retained_pct = round(gl_load_demand_retained(α) * 100; digits=4)
                        println(io, "  Bus $bus_id:  alpha = $α  ->  $(retained_pct)% of pre-fault demand retained")
                    end
                end

            # 3. Open Branch (OB) — no short-circuit stage
            elseif ts_fault_details[:fault_type] == "OB"
                ob = ts_fault_details[:ob]
                println(io, "Open Branch (no short-circuit): $(ob[:branch_id])")
                if haskey(ob, :from_bus) && haskey(ob, :to_bus)
                    for (br, fb, tb) in zip(ob[:branch_id], ob[:from_bus], ob[:to_bus])
                        println(io, "  Circuit $br:  bus $fb — bus $tb")
                    end
                end
            end

            # Add time-domain parameters from dyn_parameters_dict if they exist
            if haskey(dyn_parameters_dict, :time)
                t = dyn_parameters_dict[:time]
                println(io, "----------------------------------------------")
                println(io, "Simulation start time (s):     $(t[:t_start_sim])")
                println(io, "Simulation end time (s):       $(t[:t_end_sim])")
                println(io, "Time step (s):                 $(t[:t_step])")
                println(io, "Fault start time (s):          $(t[:t_start_fault])")
                if haskey(t, :t_clear_fault)
                    println(io, "Clearing Time (s):            $(t[:t_clear_fault])")
                end
            end
        end

        println(io, "==============================================")
    end

    println("Input parameters successfully saved in: ", path_names[:pf_results_date])
end