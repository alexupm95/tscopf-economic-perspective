# ===================================================================================
#         Registry-driven steady-state OPF dual export (AC-OPF / DC-OPF / ED)
# ===================================================================================
# A single generic writer over `STEADY_STATE_DUAL_SPECS` (DispatchDualRegistry.jl)
# replaces the former hand-written Save_Duals_{ACOPF,DCOPF,ED}_Model trio. Specs
# whose constraint key is absent in the solved model are skipped via `spec_present`,
# so the same function serves all three model types: ED exports only P-balance and
# Pg bounds, DC-OPF adds branch/angle families, AC-OPF adds the Q / V / Qg families.
#
# The spec table is the single source of truth for id columns, units, CSV filenames
# and XLSX sheet names — which removes the label-drift that produced the historical
# mislabels (Ski column, ang-diff-upper, Bus_ID-for-generator).

"""Write one TXT dual section: header banner, then one `[i] = <unit> value` per entry."""
function write_dual_section(io::IO, name::AbstractString, vals::Vector{Float64}, unit_prefix::AbstractString)
    println(io, "======================================")
    println(io, "          $name:")
    println(io, "======================================")
    for (i, val) in enumerate(vals)
        println(io, "[$i] =\t $unit_prefix $val")
    end
    println(io)  # empty line between sections
end

"""
    Save_Duals_OPF_Model(path_names, opf_dict, base_MVA; kwargs...)

Export steady-state duals present in `opf_dict` via `STEADY_STATE_DUAL_SPECS`.

Optional kwargs (UC / extensions):
  `model_label` — console messages only
  `model` — solved JuMP model; when given, its termination / primal / dual status
    is stamped at the top of `duals.txt` and a solve with no dual certificate is
    skipped entirely rather than exported as if it carried prices
  `duals_preamble` — text prepended to `duals.txt` (e.g. UC u* header)
  `excel_file` — XLSX path (default `Dispatch/Dispatch_Duals.xlsx`)
  `extra_csv` — vector of `(basename, DataFrame)` written under `CSV_duals/`
  `extra_sheets` — vector of `(sheet_name, DataFrame)` merged into the XLSX
"""
function Save_Duals_OPF_Model(
    path_names::OrderedDict{Symbol, String},
    opf_dict::OrderedDict{Symbol, Any},
    base_MVA::Float64;
    model_label::String = "OPF",
    model = nothing,
    duals_preamble::Union{Nothing, AbstractString} = nothing,
    excel_file::Union{Nothing, String} = nothing,
    extra_csv::Vector{Pair{String, DataFrame}} = Pair{String, DataFrame}[],
    extra_sheets::Vector{Pair{String, DataFrame}} = Pair{String, DataFrame}[],
    )

    if model !== nothing && !_has_usable_duals(model)
        @warn "Dual solution unavailable for the $model_label model " *
              "(termination: $(JuMP.termination_status(model)), " *
              "dual status: $(JuMP.dual_status(model))). Skipping the dispatch dual export."
        return nothing
    end

    present = [spec for spec in STEADY_STATE_DUAL_SPECS if spec_present(opf_dict, spec)]

    # Extract each present family once (JuMP.dual is not free); reuse for all formats.
    extracted = OrderedDict{Symbol, Tuple{Vector, Vector{Float64}}}()
    for spec in present
        extracted[spec.key] = extract_spec_duals(opf_dict, spec)
    end

    # ========== WRITE TO TXT FILE ==========
    open(joinpath(path_names[:pf_dispatch], "duals.txt"), "w") do io
        model !== nothing &&
            _write_solve_status_header!(io, model, "$model_label dual solution")
        if duals_preamble !== nothing
            print(io, duals_preamble)
            endswith(duals_preamble, "\n") || println(io)
        end
        for spec in present
            _, vals = extracted[spec.key]
            write_dual_section(io, spec.txt_name, vals, _unit_prefix(spec.unit))
        end
    end
    println("Duals of the $model_label model successfully saved as TXT file in: ", path_names[:pf_dispatch])

    # ========== WRITE TO CSV FILES ==========
    for spec in present
        spec.csv_file === nothing && continue
        ids, vals = extracted[spec.key]
        df = DataFrame(_id_column(spec.id_kind) => ids, spec.value_col => vals)
        CSV.write(joinpath(path_names[:pf_dispatch_CSV_duals], spec.csv_file), df; delim = ';')
    end
    for (fname, df) in extra_csv
        CSV.write(joinpath(path_names[:pf_dispatch_CSV_duals], fname), df; delim = ';')
    end
    println("Duals of the $model_label model successfully saved as CSV files in: ", path_names[:pf_dispatch_CSV_duals])

    # ========== WRITE TO XLSX FILE ==========
    excel_file = something(excel_file, joinpath(path_names[:pf_dispatch], "Dispatch_Duals.xlsx"))
    sheets_to_save = Pair{String, DataFrame}[]
    for spec in present
        spec.sheet === nothing && continue
        ids, vals = extracted[spec.key]
        df = DataFrame(_id_column(spec.id_kind) => ids, spec.value_col => vals)
        push!(sheets_to_save, spec.sheet => df)
    end
    append!(sheets_to_save, extra_sheets)
    if !isempty(sheets_to_save)
        XLSX.writetable(excel_file, sheets_to_save...; overwrite = true)
        println("All $model_label Duals successfully saved to: ", excel_file)
    else
        println("No duals were found to save.")
    end

end

# ==============================================================================
# UC duals: restricted-pricing LP (u* fixed from MILP).
# LP shadow prices reuse Save_Duals_OPF_Model + STEADY_STATE_DUAL_SPECS; commitment
# and any future UC-only dual artifacts are layered via preamble / extra_* kwargs.
# When new UC inequality families land in lp_dict, add rows to the registry (or a
# future UC_DUAL_SPECS table) — no duplicate hand-written export loops here.
# ==============================================================================
function _uc_duals_preamble(u_star::OrderedDict{Int, Float64})::String
    io = IOBuffer()
    println(io, "======================================")
    println(io, "  UC restricted pricing (u* fixed from MILP)")
    println(io, "  Shadow prices = JuMP.dual on restricted LP (registry export below)")
    println(io, "======================================")
    for (g, val) in sort(collect(u_star))
        println(io, "u*[$g] = $val")
    end
    println(io)
    return String(take!(io))
end

function Save_Duals_UC_Model(
    path_names::OrderedDict{Symbol, String},
    lp_dict::OrderedDict{Symbol, Any},
    base_MVA::Float64,
    u_star::OrderedDict{Int, Float64},
)
    mkpath(path_names[:pf_dispatch])
    mkpath(path_names[:pf_dispatch_CSV_duals])

    df_uc = DataFrame(Gen_ID = collect(keys(u_star)), u_star_fixed = collect(values(u_star)))

    Save_Duals_OPF_Model(path_names, lp_dict, base_MVA;
        model_label = "UC restricted pricing LP",
        duals_preamble = _uc_duals_preamble(u_star),
        excel_file = joinpath(path_names[:pf_dispatch], "UC_Dispatch_Duals.xlsx"),
        extra_csv = ["dual_uc_commitment.csv" => df_uc],
        extra_sheets = ["Commitment_u_star" => df_uc],
    )

    println("UC duals (restricted LP) saved in: ", path_names[:pf_dispatch])
end
